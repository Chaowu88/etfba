'''Define classes for variability analysis'''


import re
import platform
import numpy as np
from multiprocess import Pool
from pyomo.environ import Objective, maximize, minimize, Constraint
from .optim import FBAOptimizer, TFBAOptimizer, EFBAOptimizer
from ..io.results import FVAResults, TVAResults, EVAResults


class FVAOptimizer(FBAOptimizer):
    '''
    FVA calculates the variability of net fluxes under the constraints of mass 
    balance. 
    
    It's advisable to run FBA initially to obtain the optimal objective. Improper 
    objective values or gamma settings may result in the failure of estimating flux 
    ranges.
    '''

    def __init__(
            self,
            model,
            objective, 
            direction,
            obj_value,
            gamma, 
            flux_bound, 
            spec_flux_bound, 
            preset_flux, 
            irr_reactions, 
            ex_mass_bal_cons,
            **kwargs
    ):
        '''
        Parameters
        ----------
        model: Model
            The model that calls FVAOptimizer.
        objective : dict
            Mapping of reaction IDs to coefficients in the objective expression.
        direction : {'max', 'min'}
            Direction of objective optimization.
        obj_value : float
            Optimal objective obtained through flux analysis.
        gamma : float in (0, 1)
            The expression required to be no less than gamma*obj_value, or no 
            greater than (1+gamma)*obj_value, based on the direction.
        flux_bound : tuple
            Lower and upper bounds of metabolic flux.
        spec_flux_bound : dict
            Mapping of reaction IDs to specific flux bounds (lb, ub).
        preset_flux : dict
            Mapping of reaction IDs to fixed metabolic fluxes.
        irr_reactions : list
            List of irreversible reaction IDs.
        ex_conc : list
            List of metabolite concentrations excluded from optimization.
        ex_mass_bal_cons : list
            List of metabolites excluded from mass balance constraints.
        '''

        super().__init__(
            model=model, 
            objective=objective, 
            direction=direction, 
            flux_bound=flux_bound, 
            spec_flux_bound=spec_flux_bound, 
            preset_flux=preset_flux, 
            irr_reactions=irr_reactions, 
            ex_mass_bal_cons=ex_mass_bal_cons,
            parsimonious=False,
            slack=1e-3,
            **kwargs
        )

        self.obj_value = obj_value
        self.gamma = gamma


    def _build_objective(self, rxnid, direction):
        def obj_rule(model):
            if self.model.reactions[rxnid].rev:
                return model.fluxes[rxnid+'_f'] - model.fluxes[rxnid+'_b']
            else:
                return model.fluxes[rxnid]
            
        self.pyoModel.obj = Objective(rule=obj_rule, sense=direction)

    
    # def _remove_objective(self):

    #     self.pyoModel.del_component(self.pyoModel.obj)    


    def _build_objective_constraint(self):
        if self.gamma < 0 or self.gamma > 1:
            raise ValueError('gamma should be a float in [0, 1]')

        def obj_cstr_rule(model):
            obj_expr = sum(
                coe*model.fluxes[fluxid] 
                for fluxid, coe in self.objective.items()
            )
            if self.direction.lower() == 'max':
                return obj_expr >= self.gamma*self.obj_value
            elif self.direction.lower() == 'min':
                return obj_expr <= (1+self.gamma)*self.obj_value
        
        self.pyoModel.OBJcstr = Constraint(rule=obj_cstr_rule)


    def _individual_solve(self, solver, rxnids):
        if platform.system() == 'Linux':
            import os
            os.sched_setaffinity(os.getpid(), range(os.cpu_count()))
        
        self._build_flux_variables()
        self._build_mass_balance_contraints()
        self._build_objective_constraint()

        flux_range = {}
        for rxnid in rxnids:
            self._build_objective(rxnid, maximize)
            solver.solve(self.pyoModel, report_timing=False)
            flux_max = self._get_opt_obj()
            self._remove_objective()

            self._build_objective(rxnid, minimize)
            solver.solve(self.pyoModel, report_timing=False)
            flux_min = self._get_opt_obj()
            self._remove_objective()

            flux_range[rxnid] = [flux_min, flux_max]

        return flux_range
    

    def solve(self, solver='glpk', n_jobs=1):
        '''
        Parameters
        ----------
        solver: {"glpk", "gurobi"}
            "gurobi" is highly recommended for large models.
        n_jobs: int
            Number of jobs to run in parallel.
        '''

        sol = self._get_solver(solver)

        pool = Pool(processes=n_jobs)

        rxnid_chunks = np.array_split(self.rxnIDs, n_jobs)   
        
        async_res = []
        for rxnid_chunk in rxnid_chunks:
            res = pool.apply_async(
                func = self._individual_solve,
                args = (sol, rxnid_chunk)
            )
            async_res.append(res)

        pool.close()
        pool.join()

        async_res = [res.get() for res in async_res]
        flux_ranges = {rxnid: flux_range 
                       for res in async_res 
                       for rxnid, flux_range in res.items()}
        
        return FVAResults(self.obj_value, self.gamma, flux_ranges)
    

class TFVAOptimizer(FVAOptimizer, TFBAOptimizer):
    '''
    TFVA calculates the variability of net fluxes under the constraints of mass 
    balance and thermodynamic feasibility. It's important to note that not all 
    reactions are subject to thermodynamic constraints.

    It's advisable to run FBA initially to obtain the optimal objective. Improper 
    objective values or gamma settings may result in the failure of estimating flux 
    ranges.
    '''

    def __init__(
            self, 
            model, 
            objective, 
            direction,
            obj_value,
            gamma, 
            flux_bound, 
            conc_bound, 
            spec_flux_bound, 
            spec_conc_bound, 
            preset_flux, 
            preset_conc,
            preset_conc_ratio, 
            irr_reactions, 
            ex_conc, 
            ex_mass_bal_cons, 
            ex_thermo_cons, 
            dgpm_conf_level, 
            **kwargs
    ):
        '''
        Parameters
        ----------
        model: Model
            The model that calls TFVAOptimizer.
        objective : dict
            Mapping of reaction IDs to coefficients in the objective expression.
        direction : {'max', 'min'}
            Direction of objective optimization.
        obj_value : float
            Optimal objective obtained by flux analysis.
        gamma : float in (0, 1)
            The expression required to be no less than gamma*obj_value, or no 
            greater than (1+gamma)*obj_value, based on the direction.
        flux_bound : tuple
            Lower and upper bounds of metabolic flux.
        conc_bound : tuple
            Lower and upper bounds of metabolite concentration.
        spec_flux_bound : dict
            Mapping of reaction IDs to specific flux bounds (lb, ub).
        spec_conc_bound : dict
            Mapping of metabolite IDs to specific concentration bounds (lb, ub).
        preset_flux : dict
            Mapping of reaction IDs to fixed metabolic fluxes.
        preset_conc : dict
            Mapping of metabolite IDs to fixed metabolite concentrations.
        preset_conc_ratio : dict
            Mapping of pairs of metabolite IDs to fixed ratios of metabolites.
        irr_reactions : list
            List of irreversible reaction IDs.
        ex_conc : list
            List of excluded metabolite concentrations for optimization.
        ex_mass_bal_cons : list
            List of metabolites excluded from mass balance constraints.
        ex_thermo_cons : list
            List of reactions excluded from thermodynamics constraints.
        dgpm_conf_level : float
            Confidence level considered if uncertainty of standard reaction Gibbs 
            energy is taken into account.
        '''

        super().__init__(
            model=model, 
            objective=objective, 
            direction=direction,
            obj_value=obj_value,
            gamma=gamma, 
            flux_bound=flux_bound, 
            conc_bound=conc_bound, 
            spec_flux_bound=spec_flux_bound, 
            spec_conc_bound=spec_conc_bound, 
            preset_flux=preset_flux, 
            preset_conc=preset_conc,
            preset_conc_ratio=preset_conc_ratio, 
            irr_reactions=irr_reactions, 
            ex_conc=ex_conc, 
            ex_mass_bal_cons=ex_mass_bal_cons, 
            ex_thermo_cons=ex_thermo_cons, 
            dgpm_conf_level=dgpm_conf_level,
            **kwargs
        )


    def _individual_solve(self, solver, rxnids):
        if platform.system() == 'Linux':
            import os
            os.sched_setaffinity(os.getpid(), range(os.cpu_count()))
        
        self._build_flux_variables()
        self._build_conc_variables()
        self._build_binary_variables()
        self._build_mass_balance_contraints()
        self._build_flux_bound_constraints()
        self._build_ratio_constraint()
        self._build_thermodynamics_constraints()
        self._build_objective_constraint()

        flux_range = {}
        for rxnid in rxnids:
            self._build_objective(rxnid, maximize)
            solver.solve(self.pyoModel, report_timing=False)
            flux_max = self._get_opt_obj()
            self._remove_objective()

            self._build_objective(rxnid, minimize)
            solver.solve(self.pyoModel, report_timing=False)
            flux_min = self._get_opt_obj()
            self._remove_objective()

            flux_range[rxnid] = [flux_min, flux_max]

        return flux_range


class EFVAOptimizer(FVAOptimizer, EFBAOptimizer):
    '''
    EFVA calculates the variability of net fluxes under the constraints of mass 
    balance and enzyme protein allocation. It's important to note that not all 
    reactions are subject to enzyme protein constraints.
    
    It's advisable to run FBA initially to obtain the optimal objective. Improper 
    objective values or gamma settings may result in the failure of estimating flux 
    ranges.
    '''

    def __init__(
            self,
            model,
            objective, 
            direction,
            obj_value,
            gamma, 
            flux_bound, 
            spec_flux_bound, 
            preset_flux, 
            irr_reactions, 
            ex_mass_bal_cons,
            inc_enz_cons,
            enz_prot_lb,
            **kwargs
    ):
        '''
        Parameters
        ----------
        model: Model
            The model that calls FVAOptimizer.
        objective : dict
            Mapping of reaction IDs to coefficients in the objective expression.
        direction : {'max', 'min'}
            Direction of objective optimization.
        obj_value : float
            Optimal objective obtained by flux analysis.
        gamma : float in (0, 1)
            The expression required to be no less than gamma*obj_value, or no 
            greater than (1+gamma)*obj_value, based on the direction.
        flux_bound : tuple
            Lower and upper bounds of metabolic flux.
        spec_flux_bound : dict
            Mapping of reaction IDs to specific flux bounds (lb, ub).
        preset_flux : dict
            Mapping of reaction IDs to fixed metabolic fluxes.
        irr_reactions : list
            List of irreversible reaction IDs.
        ex_conc : list
            List of metabolite concentrations excluded from optimization.
        ex_mass_bal_cons : list
            List of metabolites excluded from mass balance constraints.
        inc_enz_cons : list
            List of reactions excluded from enzyme protein cost constraints.
        enz_prot_lb : float
            Upper bound of enzyme protein fraction.
        '''

        super().__init__(
            model=model,
            objective=objective, 
            direction=direction,
            obj_value=obj_value,
            gamma=gamma, 
            flux_bound=flux_bound, 
            spec_flux_bound=spec_flux_bound, 
            preset_flux=preset_flux, 
            irr_reactions=irr_reactions, 
            ex_mass_bal_cons=ex_mass_bal_cons,
            inc_enz_cons=inc_enz_cons,
            enz_prot_lb=enz_prot_lb,
            **kwargs
        )


    def _individual_solve(self, solver, rxnids):
        if platform.system() == 'Linux':
            import os
            os.sched_setaffinity(os.getpid(), range(os.cpu_count()))
        
        self._build_flux_variables()
        self._build_mass_balance_contraints()
        self._build_objective_constraint()
        self._build_enzyme_cost_constraint()

        flux_range = {}
        for rxnid in rxnids:
            self._build_objective(rxnid, maximize)
            solver.solve(self.pyoModel, report_timing=False)
            flux_max = self._get_opt_obj()
            self._remove_objective()

            self._build_objective(rxnid, minimize)
            solver.solve(self.pyoModel, report_timing=False)
            flux_min = self._get_opt_obj()
            self._remove_objective()

            flux_range[rxnid] = [flux_min, flux_max]

        return flux_range


class ETFVAOptimizer(TFVAOptimizer, EFVAOptimizer):
    '''
    ETFVA calculates the variability of net fluxes under the constraints of mass 
    balance, thermodynamic feasibility, and enzyme protein allocation. It's 
    important to note that not all reactions are subject to enzyme protein 
    constraints or thermodynamic constraints.
    
    It's advisable to run FBA initially to obtain the optimal objective. Improper 
    objective values or gamma settings may result in the failure of estimating flux 
    ranges.
    '''

    def __init__(
            self,
            model, 
            objective, 
            direction,
            obj_value,
            gamma, 
            flux_bound, 
            conc_bound, 
            spec_flux_bound, 
            spec_conc_bound, 
            preset_flux, 
            preset_conc,
            preset_conc_ratio, 
            irr_reactions, 
            ex_conc, 
            ex_mass_bal_cons, 
            ex_thermo_cons,
            inc_enz_cons, 
            enz_prot_lb, 
            dgpm_conf_level,
    ):
        '''
        Parameters
        ----------
        model : Model
            The model that calls ETVAOptimizer.
        objective : dict
            Mapping of reaction IDs to coefficients in the objective expression.
        direction : {'max', 'min'}
            Direction of objective optimization.
        obj_value : float
            Optimal objective obtained by flux analysis.
        gamma : float in (0, 1)
            The expression required to be no less than gamma*obj_value, or no 
            greater than (1+gamma)*obj_value, based on the direction.
        flux_bound : tuple
            Lower and upper bounds of metabolic flux.
        conc_bound : tuple
            Lower and upper bounds of metabolite concentration.
        spec_flux_bound : dict
            Mapping of reaction IDs to specific flux bounds (lb, ub).
        spec_conc_bound : dict
            Mapping of metabolite IDs to specific concentration bounds (lb, ub).
        preset_flux : dict
            Mapping of reaction IDs to fixed metabolic fluxes.
        preset_conc : dict
            Mapping of metabolite IDs to fixed metabolite concentrations.
        preset_conc_ratio : dict
            Mapping of pairs of metabolite IDs to fixed ratios of metabolites.
        irr_reactions : list
            List of irreversible reaction IDs.
        ex_conc : list
            List of excluded metabolite concentrations for optimization.
        ex_mass_bal_cons : list
            List of metabolites excluded from mass balance constraints.
        ex_thermo_cons : list
            List of reactions excluded from thermodynamics constraints.
        inc_enz_cons : list
            List of reactions excluded from enzyme protein cost constraints.
        enz_prot_lb : float
            Upper bound of enzyme protein fraction.
        dgpm_conf_level : float
            Confidence level considered if uncertainty of standard reaction Gibbs 
            energy is taken into account.
        '''

        super().__init__(
            model=model, 
            objective=objective, 
            direction=direction,
            obj_value=obj_value,
            gamma=gamma, 
            flux_bound=flux_bound, 
            conc_bound=conc_bound, 
            spec_flux_bound=spec_flux_bound, 
            spec_conc_bound=spec_conc_bound, 
            preset_flux=preset_flux, 
            preset_conc=preset_conc,
            preset_conc_ratio=preset_conc_ratio, 
            irr_reactions=irr_reactions, 
            ex_conc=ex_conc, 
            ex_mass_bal_cons=ex_mass_bal_cons, 
            ex_thermo_cons=ex_thermo_cons,
            inc_enz_cons=inc_enz_cons,
            enz_prot_lb=enz_prot_lb, 
            dgpm_conf_level=dgpm_conf_level,
        )


    def _individual_solve(self, solver, rxnids):
        if platform.system() == 'Linux':
            import os
            os.sched_setaffinity(os.getpid(), range(os.cpu_count()))
        
        self._build_flux_variables()
        self._build_conc_variables()
        self._build_binary_variables()
        self._build_mass_balance_contraints()
        self._build_flux_bound_constraints()
        self._build_ratio_constraint()
        self._build_thermodynamics_constraints()
        self._build_enzyme_cost_constraint()
        self._build_objective_constraint()

        flux_range = {}
        for rxnid in rxnids:
            self._build_objective(rxnid, maximize)
            solver.solve(self.pyoModel, report_timing=False)
            flux_max = self._get_opt_obj()
            self._remove_objective()

            self._build_objective(rxnid, minimize)
            solver.solve(self.pyoModel, report_timing=False)
            flux_min = self._get_opt_obj()
            self._remove_objective()

            flux_range[rxnid] = [flux_min, flux_max]

        return flux_range    


class TVAOptimizer(TFVAOptimizer):
    '''
    TVA calculates the variability of Gibbs energy for reactions under the 
    constraints of thermodynamic feasibility and mass balance.
    
    It's advisable to run FBA initially to obtain the optimal objective. Improper 
    objective values or gamma settings may result in the failure of estimating 
    Gibbs energy ranges.
    '''

    def __init__(
            self, 
            model, 
            objective, 
            direction,
            obj_value,
            gamma, 
            flux_bound, 
            conc_bound, 
            spec_flux_bound, 
            spec_conc_bound, 
            preset_flux, 
            preset_conc,
            preset_conc_ratio, 
            irr_reactions, 
            ex_conc, 
            ex_mass_bal_cons, 
            ex_thermo_cons, 
            dgpm_conf_level, 
            **kwargs
    ):
        '''
        Parameters
        ----------
        model: Model
            The model that calls TVAOptimizer.
        objective : dict
            Mapping of reaction IDs to coefficients in the objective expression.
        direction : {'max', 'min'}
            Direction of objective optimization.
        obj_value : float
            Optimal objective obtained by flux analysis.
        gamma : float in (0, 1)
            The expression required to be no less than gamma*obj_value, or no 
            greater than (1+gamma)*obj_value, based on the direction.
        flux_bound : tuple
            Lower and upper bounds of metabolic flux.
        conc_bound : tuple
            Lower and upper bounds of metabolite concentration.
        spec_flux_bound : dict
            Mapping of reaction IDs to specific flux bounds (lb, ub).
        spec_conc_bound : dict
            Mapping of metabolite IDs to specific concentration bounds (lb, ub).
        preset_flux : dict
            Mapping of reaction IDs to fixed metabolic fluxes.
        preset_conc : dict
            Mapping of metabolite IDs to fixed metabolite concentrations.
        preset_conc_ratio : dict
            Mapping of pairs of metabolite IDs to fixed ratios of metabolites.
        irr_reactions : list
            List of irreversible reaction IDs.
        ex_conc : list
            List of excluded metabolite concentrations for optimization.
        ex_mass_bal_cons : list
            List of metabolites excluded from mass balance constraints.
        ex_thermo_cons : list
            List of reactions excluded from thermodynamics constraints.
        dgpm_conf_level : float
            Confidence level considered if uncertainty of standard reaction Gibbs 
            energy is taken into account.
        '''

        super().__init__(
            model=model, 
            objective=objective, 
            direction=direction,
            obj_value=obj_value,
            gamma=gamma, 
            flux_bound=flux_bound, 
            conc_bound=conc_bound, 
            spec_flux_bound=spec_flux_bound, 
            spec_conc_bound=spec_conc_bound, 
            preset_flux=preset_flux, 
            preset_conc=preset_conc,
            preset_conc_ratio=preset_conc_ratio, 
            irr_reactions=irr_reactions, 
            ex_conc=ex_conc, 
            ex_mass_bal_cons=ex_mass_bal_cons, 
            ex_thermo_cons=ex_thermo_cons, 
            dgpm_conf_level=dgpm_conf_level,
            **kwargs
        )


    def _build_objective(self, fluxid, direction):
        def obj_rule(model):
            return self._calculate_gibbs_energy(model, fluxid)

        self.pyoModel.obj = Objective(rule=obj_rule, sense=direction)
    

    # can use super()._individual_solve directly, here _individual_solve is 
    # explicitly defined just for clarification.
    # def _individual_solve(self, solver, fluxids):

    #     if platform.system() == 'Linux':
    #         import os
    #         os.sched_setaffinity(os.getpid(), range(os.cpu_count()))
        
    #     self._build_flux_variables()
    #     self._build_conc_variables()
    #     self._build_binary_variables()
    #     self._build_mass_balance_contraints()
    #     self._build_flux_bound_constraints()
    #     self._build_ratio_constraint()
    #     self._build_thermodynamics_constraints()
    #     self._build_objective_constraint()

    #     dgp_range = {}
    #     for fluxid in fluxids:
    #         self._build_objective(fluxid, maximize)
    #         solver.solve(self.pyoModel, report_timing = False)
    #         dgp_max = self._get_opt_obj()
    #         self._remove_objective()

    #         self._build_objective(fluxid, minimize)
    #         solver.solve(self.pyoModel, report_timing = False)
    #         dgp_min = self._get_opt_obj()
    #         self._remove_objective()

    #         dgp_range[fluxid] = [dgp_min, dgp_max]

    #     return dgp_range


    def solve(self, solver='glpk', n_jobs=1):
        '''
        Parameters
        ----------
        solver: {"glpk", "gurobi"}
            "gurobi" is highly recommended for large models.
        n_jobs: int
            Number of jobs to run in parallel.
        '''

        sol = self._get_solver(solver)

        pool = Pool(processes=n_jobs)

        fluxids_filtered = list(
            filter(lambda fluxid: not re.match(r'.+_b$', fluxid), self.cstrFluxIDs)
        )   # no need to estimate backward reactions
        # only estimate the deltaGprime of reactions that subject to thermodynamic constraints
        fluxid_chunks = np.array_split(fluxids_filtered, n_jobs)
        
        async_res = []
        for fluxid_chunk in fluxid_chunks:
            res = pool.apply_async(
                func=self._individual_solve,
                args=(sol, fluxid_chunk)
            )
            async_res.append(res)

        pool.close()
        pool.join()    

        async_res = [res.get() for res in async_res]
        dgp_ranges = {re.sub(r'_[fb]$', '', fluxid): dgp_range 
                      for res in async_res 
                      for fluxid, dgp_range in res.items()}
        
        return TVAResults(self.obj_value, self.gamma, dgp_ranges)

        
class ETVAOptimizer(TVAOptimizer, EFBAOptimizer):
    '''
    ETVA calculates the variability of Gibbs energy for reactions under the 
    constraints of thermodynamic feasibility, mass balance, and enzyme protein 
    allocation.
    
    It's advisable to run FBA initially to obtain the optimal objective. Improper 
    objective values or gamma settings may result in the failure of estimating 
    Gibbs energy ranges.
    '''

    def __init__(
            self, 
            model, 
            objective, 
            direction,
            obj_value,
            gamma, 
            flux_bound, 
            conc_bound, 
            spec_flux_bound, 
            spec_conc_bound, 
            preset_flux, 
            preset_conc,
            preset_conc_ratio, 
            irr_reactions, 
            ex_conc, 
            ex_mass_bal_cons, 
            ex_thermo_cons,
            inc_enz_cons, 
            enz_prot_lb, 
            dgpm_conf_level, 
    ):
        '''
        Parameters
        ----------
        model: Model
            The model that calls ETVAOptimizer.
        objective : dict
            Mapping of reaction IDs to coefficients in the objective expression.
        direction : {'max', 'min'}
            Direction of objective optimization.
        obj_value : float
            Optimal objective obtained by flux analysis.
        gamma : float in (0, 1)
            The expression required to be no less than gamma*obj_value, or no 
            greater than (1+gamma)*obj_value, based on the direction.
        flux_bound : tuple
            Lower and upper bounds of metabolic flux.
        conc_bound : tuple
            Lower and upper bounds of metabolite concentration.
        spec_flux_bound : dict
            Mapping of reaction IDs to specific flux bounds (lb, ub).
        spec_conc_bound : dict
            Mapping of metabolite IDs to specific concentration bounds (lb, ub).
        preset_flux : dict
            Mapping of reaction IDs to fixed metabolic fluxes.
        preset_conc : dict
            Mapping of metabolite IDs to fixed metabolite concentrations.
        preset_conc_ratio : dict
            Mapping of pairs of metabolite IDs to fixed ratios of metabolites.
        irr_reactions : list
            List of irreversible reaction IDs.
        ex_conc : list
            List of excluded metabolite concentrations for optimization.
        ex_mass_bal_cons : list
            List of metabolites excluded from mass balance constraints.
        ex_thermo_cons : list
            List of reactions excluded from thermodynamics constraints.
        inc_enz_cons : list
            List of reactions excluded from enzyme protein cost constraints.
        enz_prot_lb : float
            Upper bound of enzyme protein fraction.
        dgpm_conf_level : float
            Confidence level considered if uncertainty of standard reaction Gibbs 
            energy is taken into account.
        '''

        super().__init__(
            model=model, 
            objective=objective, 
            direction=direction,
            obj_value=obj_value,
            gamma=gamma, 
            flux_bound=flux_bound, 
            conc_bound=conc_bound, 
            spec_flux_bound=spec_flux_bound, 
            spec_conc_bound=spec_conc_bound, 
            preset_flux=preset_flux, 
            preset_conc=preset_conc,
            preset_conc_ratio=preset_conc_ratio, 
            irr_reactions=irr_reactions, 
            ex_conc=ex_conc, 
            ex_mass_bal_cons=ex_mass_bal_cons, 
            ex_thermo_cons=ex_thermo_cons,
            inc_enz_cons=inc_enz_cons,
            enz_prot_lb=enz_prot_lb, 
            dgpm_conf_level=dgpm_conf_level,
        )


    def _individual_solve(self, solver, fluxids):
        if platform.system() == 'Linux':
            import os
            os.sched_setaffinity(os.getpid(), range(os.cpu_count()))
        
        self._build_flux_variables()
        self._build_conc_variables()
        self._build_binary_variables()
        self._build_mass_balance_contraints()
        self._build_flux_bound_constraints()
        self._build_ratio_constraint()
        self._build_thermodynamics_constraints()
        self._build_enzyme_cost_constraint()
        self._build_objective_constraint()

        dgp_range = {}
        for fluxid in fluxids:
            self._build_objective(fluxid, maximize)
            solver.solve(self.pyoModel, report_timing=False)
            dgp_max = self._get_opt_obj()
            self._remove_objective()

            self._build_objective(fluxid, minimize)
            solver.solve(self.pyoModel, report_timing=False)
            dgp_min = self._get_opt_obj()
            self._remove_objective()

            dgp_range[fluxid] = [dgp_min, dgp_max]

        return dgp_range
    

class EVAOptimizer(EFVAOptimizer):
    '''
    EVA estimates the variability of enzyme protein costs under the constraints of 
    mass balance and total enzyme protein allocation.
    
    It's advisable to run FBA initially to obtain the optimal objective. Improper 
    objective values or gamma settings may result in the failure of estimating 
    enzyme protein cost ranges.
    '''

    def __init__(
            self,
            model,
            objective, 
            direction,
            obj_value,
            gamma, 
            flux_bound, 
            spec_flux_bound, 
            preset_flux,
            irr_reactions, 
            ex_mass_bal_cons, 
            inc_enz_cons, 
            enz_prot_lb,
            **kwargs
    ):
        '''
        Parameters
        ----------
        model: Model
            Model that calls EVAOptimizer.
        objective : dict
            Mapping of reaction IDs to coefficients in the objective expression.
        direction : {'max', 'min'}
            Direction of objective optimization.
        obj_value : float
            Optimal objective obtained by flux analysis.
        gamma : float in (0, 1)
            The expression required to be no less than gamma*obj_value, or no 
            greater than (1+gamma)*obj_value, based on the direction.
        flux_bound : tuple
            Lower and upper bounds of metabolic flux.
        spec_flux_bound : dict
            Mapping of reaction IDs to specific flux bounds (lb, ub).
        preset_flux : dict
            Mapping of reaction IDs to fixed metabolic fluxes.
        irr_reactions : list
            List of irreversible reaction IDs.
        ex_mass_bal_cons : list
            List of metabolites excluded from mass balance constraints.
        inc_enz_cons : list
            List of reactions excluded from enzyme protein cost constraints.
        enz_prot_lb : float
            Upper bound of enzyme protein fraction.
        '''

        super().__init__(
            model=model, 
            objective=objective, 
            direction=direction,
            obj_value=obj_value,
            gamma=gamma, 
            flux_bound=flux_bound, 
            spec_flux_bound=spec_flux_bound, 
            preset_flux=preset_flux, 
            irr_reactions=irr_reactions, 
            ex_mass_bal_cons=ex_mass_bal_cons, 
            inc_enz_cons=inc_enz_cons,
            enz_prot_lb=enz_prot_lb,
            **kwargs
        )

        
    def _build_objective(self, rxnid, direction):
        def obj_rule(model):
            return self._calculate_enzyme_cost(model, rxnid)
        
        self.pyoModel.obj = Objective(rule=obj_rule, sense=direction)


    def solve(self, solver='glpk', n_jobs=1):
        '''
        Parameters
        ----------
        solver: {"glpk", "gurobi"}
            "gurobi" is highly recommended for large models.
        n_jobs: int
            Number of jobs to run in parallel.
        '''

        sol = self._get_solver(solver)

        pool = Pool(processes=n_jobs)

        rxnid_chunks = np.array_split(self.inc_enz_cons, n_jobs)
        # only estimate the enzyme protein cost of reactions that subject to protein allocation 
        # constraint
        
        async_res = []
        for rxnid_chunk in rxnid_chunks:
            res = pool.apply_async(
                func = self._individual_solve,
                args = (sol, rxnid_chunk)
            )
            async_res.append(res)

        pool.close()
        pool.join()

        async_res = [res.get() for res in async_res]
        epc_ranges = {rxnid: epc_range 
                      for res in async_res 
                      for rxnid, epc_range in res.items()}
        
        return EVAResults(self.obj_value, self.gamma, epc_ranges)


class TEVAOptimizer(EVAOptimizer, TFBAOptimizer):
    '''
    TEVA estimates the variability of enzyme protein costs under the constraints of 
    thermodynamic feasibility, mass balance, and enzyme protein allocation.
    
    It's advisable to run FBA initially to obtain the optimal objective. Improper 
    objective values or gamma settings may result in the failure of estimating 
    enzyme protein cost ranges.
    '''
            
    def __init__(
            self, 
            model, 
            objective, 
            direction,
            obj_value,
            gamma, 
            flux_bound, 
            conc_bound, 
            spec_flux_bound, 
            spec_conc_bound, 
            preset_flux, 
            preset_conc,
            preset_conc_ratio, 
            irr_reactions, 
            ex_conc, 
            ex_mass_bal_cons, 
            ex_thermo_cons,
            inc_enz_cons, 
            enz_prot_lb, 
            dgpm_conf_level, 
    ):
        '''
        Parameters
        ----------
        model: Model
            Model that calls TEVAOptimizer.
        objective : dict
            Mapping of reaction IDs to coefficients in the objective expression.
        direction : {'max', 'min'}
            Direction of objective optimization.
        obj_value : float
            Optimal objective obtained by flux analysis.
        gamma : float in (0, 1)
            The expression required to be no less than gamma*obj_value, or no 
            greater than (1+gamma)*obj_value, based on the direction.
        flux_bound : tuple
            Lower and upper bounds of metabolic flux.
        conc_bound : tuple
            Lower and upper bounds of metabolite concentration.
        spec_flux_bound : dict
            Mapping of reaction IDs to specific flux bounds (lb, ub).
        spec_conc_bound : dict
            Mapping of metabolite IDs to specific concentration bounds (lb, ub).
        preset_flux : dict
            Mapping of reaction IDs to fixed metabolic fluxes.
        preset_conc : dict
            Mapping of metabolite IDs to fixed metabolite concentrations.
        preset_conc_ratio : dict
            Mapping of metabolite ratio IDs to fixed ratios of metabolites (metabid:metabid => float).
        irr_reactions : list
            List of irreversible reaction IDs.
        ex_conc : list
            List of metabolite concentrations excluded from optimization.
        ex_mass_bal_cons : list
            List of metabolites excluded from mass balance constraints.
        ex_thermo_cons : list
            List of reactions excluded from thermodynamics constraints.
        inc_enz_cons : list
            List of reactions included in enzyme protein cost constraints.
        enz_prot_lb : float
            Upper bound of enzyme protein fraction.
        dgpm_conf_level : float
            Confidence level if uncertainty of standard reaction Gibbs energy is 
            considered.
        '''

        super().__init__(
            model=model, 
            objective=objective, 
            direction=direction,
            obj_value=obj_value,
            gamma=gamma, 
            flux_bound=flux_bound, 
            conc_bound=conc_bound, 
            spec_flux_bound=spec_flux_bound, 
            spec_conc_bound=spec_conc_bound, 
            preset_flux=preset_flux, 
            preset_conc=preset_conc,
            preset_conc_ratio=preset_conc_ratio, 
            irr_reactions=irr_reactions, 
            ex_conc=ex_conc, 
            ex_mass_bal_cons=ex_mass_bal_cons, 
            ex_thermo_cons=ex_thermo_cons,
            inc_enz_cons=inc_enz_cons,
            enz_prot_lb=enz_prot_lb, 
            dgpm_conf_level=dgpm_conf_level,
        )

    def _individual_solve(self, solver, rxnids):
        if platform.system() == 'Linux':
            import os
            os.sched_setaffinity(os.getpid(), range(os.cpu_count()))
        
        self._build_flux_variables()
        self._build_conc_variables()
        self._build_binary_variables()
        self._build_mass_balance_contraints()
        self._build_flux_bound_constraints()
        self._build_ratio_constraint()
        self._build_thermodynamics_constraints()
        self._build_enzyme_cost_constraint()
        self._build_objective_constraint()

        epc_range = {}
        for rxnid in rxnids:
            self._build_objective(rxnid, maximize)
            solver.solve(self.pyoModel, report_timing=False)
            epc_max = self._get_opt_obj()
            self._remove_objective()

            self._build_objective(rxnid, minimize)
            solver.solve(self.pyoModel, report_timing=False)
            epc_min = self._get_opt_obj()
            self._remove_objective()

            epc_range[rxnid] = [epc_min, epc_max]

        return epc_range

                        