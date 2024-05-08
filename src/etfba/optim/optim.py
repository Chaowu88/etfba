'''Difine classes for constrained flux analysis.'''


import re
import numpy as np
from scipy.stats import norm
from pyomo.environ import (ConcreteModel, Set, Var, Objective, Constraint, 
                           SolverFactory, NonNegativeReals, Reals, Binary, 
                           value, maximize, minimize, log)
from pyomo.opt import SolverStatus, TerminationCondition
from pyomo.core.expr.numeric_expr import LinearExpression
import logging
logging.basicConfig(level = logging.INFO, format = '%(levelname)s: %(message)s')
from ..io.results import FBAResults, TFBAResults, EFBAResults, ETFBAResults


R = 8.315e-3         # Gas constant in kJ/mol/K
T = 298.15           # Absolute temperature in K, equivalent to 25 C
K = 1e6              # A sufficiently large constant
EPSILON = 1e-3       # Tolerance ensuring reactions proceed with Gibbs energy 
                     # dissipation. Use 1 for stricter constraint
# THDM_THLD = 0.5      # threshold of Gibbs energy change
# FLUX_THLD = 1e-7     # threshold of a zero-flux
# MAXITER = 1000       # maximum iterations for NLP problem


class FBAOptimizer():
    '''
    FBA computes the fluxes by optimizing the objective (e.g., biomass formation) 
    under the constraints of mass balance.
    '''
    
    def __init__(
            self, 
            model, 
            objective, 
            direction, 
            flux_bound, 
            spec_flux_bound, 
            preset_flux, 
            irr_reactions, 
            ex_mass_bal_cons,
            parsimonious,
            slack
    ):
        '''
        Parameters
        ----------
        model: Model
            The model that calls FBAOptimizer.
        objective : dict
            Mapping of reaction IDs to coefficients in the objective expression.
        direction : {'max', 'min'}
            Direction of optimization.
        flux_bound : tuple
            Lower and upper bounds of metabolic flux.
        spec_flux_bound : dict
            Mapping of reaction IDs to their specific flux bounds (lb, ub).
        preset_flux : dict
            Mapping of reaction IDs to fixed metabolic fluxes.
        irr_reactions : list
            List of irreversible reaction IDs.
        ex_mass_bal_cons : list
            List of metabolites excluded from mass balance constraints.
        '''
        
        self.model = model
        
        self.rxnIDs = self.model.stoichiometric_matrix.columns.tolist()
        self.metabIDs = self.model.total_stoichiometric_matrix.index.tolist()
        
        self.objective = objective
        self.direction = direction
        
        self.flux_bound = flux_bound
        
        if spec_flux_bound is None:
            self.spec_flux_bound = {}  
        else:
            self.spec_flux_bound = spec_flux_bound
        
        if preset_flux is None:
            self.preset_flux = {}  
        else:
            self.preset_flux = preset_flux
        
        self.irr_reactions = irr_reactions
        if self.irr_reactions is not None:   # irr_reactions is prioritized
            for rxnid in self.rxnIDs:
                if rxnid in self.irr_reactions:
                    self.model.reactions[rxnid].rev = False
                else:
                    self.model.reactions[rxnid].rev = True
        else:
            self.irr_reactions = [rxnid for rxnid in self.rxnIDs 
                                  if not self.model.reactions[rxnid].rev]
        
        self.varFluxIDs = self.model.total_stoichiometric_matrix.columns.tolist()
        # varFluxIDs denote total fluxes as variables
        
        if ex_mass_bal_cons is None:
            self.ex_mass_bal_cons = []  
        else:
            self.ex_mass_bal_cons = list(set(ex_mass_bal_cons))
        
        self.cstrMetabIDs = []
        for metabid in self.metabIDs:
            #if metabid not in self.model.end_metabolites and metabid not in self.ex_mass_bal_cons:
            if metabid not in self.ex_mass_bal_cons:
                self.cstrMetabIDs.append(metabid)
                self.model.metabolites[metabid].is_constrained_by_mass_balance = True
        # cstrMetabIDs denote metabolites involved in mass balance
        # metabolites in ex_mass_bal_cons are excluded
        
        self.parsimonious = parsimonious
        self.slack = slack

        self.pyoModel = ConcreteModel()
        self.pyoModel.varFluxIDs = Set(initialize = self.varFluxIDs)
        self.pyoModel.cstrMetabIDs = Set(initialize = self.cstrMetabIDs)
        
    
    def _build_flux_variables(self, initial = None):
        '''
        Parameters
        ----------
        initial: dict
            Mapping of reaction IDs to their initial values.
        '''
        
        if not set(self.preset_flux.keys()).issubset(self.varFluxIDs):
            logging.warning(
                'some preset fluxes are not used, note "_f" and "_b" should '
                'be added as suffix for reversible reactions'
            )

        def flux_bounds_rule(model, fluxid):
            if fluxid in self.preset_flux:
                return (self.preset_flux[fluxid],)*2
            elif fluxid in self.spec_flux_bound:
                return self.spec_flux_bound[fluxid]
            else:
                return self.flux_bound
          
        if initial is None:
            self.pyoModel.fluxes = Var(
                self.pyoModel.varFluxIDs, 
                within = NonNegativeReals, 
                bounds = flux_bounds_rule
            )
        else:
            logging.info('load initial flux values')

            initial = {fluxid: flux for fluxid, flux in initial.items() 
                       if fluxid in self.pyoModel.varFluxIDs}
            self.pyoModel.fluxes = Var(
                self.pyoModel.varFluxIDs, 
                within = NonNegativeReals,
                bounds = flux_bounds_rule, 
                initialize = initial
            )
        
    
    def _build_objective(self):
        for fluxid in self.objective:
            if fluxid not in self.varFluxIDs:
                raise KeyError(f'{fluxid} in objective not exist in the model')
        
        if self.direction.lower() == 'max':
            direction = maximize
        elif self.direction.lower() == 'min':
            direction = minimize
        else:
            raise ValueError("only 'max' or 'min' is acceptable")
        
        def obj_rule(model):
            return sum(coe*model.fluxes[fluxid] 
                       for fluxid, coe in self.objective.items())
            
        self.pyoModel.obj = Objective(rule = obj_rule, sense = direction)    
        
    
    def _build_parsimonious_objective(self):
        def obj_rule(model):
            return sum(model.fluxes[fluxid] for fluxid in self.varFluxIDs)
        
        self.pyoModel.obj = Objective(rule = obj_rule, sense = minimize)
        

    def _remove_objective(self):
        self.pyoModel.del_component(self.pyoModel.obj)


    def _build_mass_balance_contraints(self):
        def mb_rule(model, metabid):
            mb_cstr = LinearExpression(
                constant = 0, 
                linear_coefs = self.model.total_stoichiometric_matrix.loc[metabid, :].tolist(),
                linear_vars = [model.fluxes[fluxid] for fluxid in self.varFluxIDs]
            )
            # LinearExpression will greatly speed up building the constraints
            return mb_cstr == 0
            
        self.pyoModel.MBcstrs = Constraint(
            self.pyoModel.cstrMetabIDs, 
            rule = mb_rule
        )


    def _build_objective_constraint(self, opt_obj):
        def obj_cstr_rule(model):
            obj_expr = sum(
                coe*model.fluxes[fluxid] for fluxid, coe in self.objective.items()
            )
            if self.direction.lower() == 'max':
                return obj_expr >= (1-self.slack)*opt_obj
            elif self.direction.lower() == 'min':
                return obj_expr <= (1+self.slack)*opt_obj
        
        self.pyoModel.OBJcstr = Constraint(rule = obj_cstr_rule)
        

    def _get_solver(self, solver):
        if solver == 'glpk':
            sol = SolverFactory('glpk')
        elif solver == 'gurobi':
            sol = SolverFactory('gurobi_direct')
            # sol = SolverFactory('gurobi', solver_io = 'python')
        else:
            raise ValueError('solver should be "glpk" or "gurobi"')
        
        return sol


    def _optimization_successful(self):
        return (
            self.res.solver.status == SolverStatus.ok and 
            self.res.solver.termination_condition in [
                TerminationCondition.optimal, 
                TerminationCondition.feasible
            ]
        )
    

    def _get_opt_obj(self):
        return value(self.pyoModel.obj)
        
        
    def _get_opt_fluxes(self):
        '''
        Return net fluxes.
        '''
        
        optTotalFluxes = {}
        for fluxid in self.pyoModel.varFluxIDs:
            optTotalFluxes[fluxid] = value(self.pyoModel.fluxes[fluxid])
            
        optNetFluxes = (
            self.model.transformation_matrix@list(optTotalFluxes.values())
        ).to_dict()
        
        return optNetFluxes

    
    def solve(self, solver = 'glpk'):
        '''
        Parameters
        ----------
        solver: {"glpk", "gurobi"}
            "gurobi" is highly recommended for large models.
        '''    
        
        self._build_flux_variables()
        self._build_objective()
        self._build_mass_balance_contraints()
        
        sol = self._get_solver(solver)
        self.res = sol.solve(self.pyoModel, report_timing = False)
        #self.pyoModel.pprint()#!!!
        #self.pyoModel.display()#!!!
        optObj = self._get_opt_obj()
        optFluxes = self._get_opt_fluxes()
        optSuccess = self._optimization_successful()

        if self.parsimonious and optSuccess:
            logging.info('estimating parsimonious fluxes')
            
            self._remove_objective()
            self._build_parsimonious_objective()
            self._build_objective_constraint(optObj)

            self.res = sol.solve(self.pyoModel, report_timing = False)

            optFluxes = self._get_opt_fluxes()
            optSuccess = self._optimization_successful()
        
        return FBAResults(
            optObj, 
            optFluxes, 
            optSuccess, 
            self.model.stoichiometric_matrix
        )    


class TFBAOptimizer(FBAOptimizer):
    '''
    TFBA computes the fluxes by optimizing the objective (e.g., biomass formation) 
    under the constraints of mass balance and thermodynamics.
    '''
    
    def __init__(
            self, 
            model, 
            objective, 
            direction, 
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
            parsimonious,
            slack, 
            dgpm_conf_level,
            **kwargs
    ):
        '''
        Parameters
        ----------
        model: Model
            The model that calls TFBAOptimizer.
        objective : dict
            Mapping of reaction IDs to coefficients in the objective expression.
        direction : {'max', 'min'}
            Direction of optimization.
        flux_bound : tuple
            Lower and upper bounds of metabolic flux.
        conc_bound : tuple
            Lower and upper bounds of metabolite concentration.
        spec_flux_bound : dict
            Mapping of reaction IDs to their specific flux bounds (lb, ub).
        spec_conc_bound : dict
            Mapping of metabolite IDs to their specific concentration bounds 
            (lb, ub).
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
            Confidence level for considering uncertainty in standard reaction Gibbs 
            energy.
        '''
        
        super().__init__(
            model = model, 
            objective = objective, 
            direction = direction, 
            flux_bound = flux_bound, 
            spec_flux_bound = spec_flux_bound, 
            preset_flux = preset_flux, 
            irr_reactions = irr_reactions, 
            ex_mass_bal_cons = ex_mass_bal_cons,
            parsimonious = parsimonious,
            slack = slack, 
            **kwargs
        )

        self.conf_level = dgpm_conf_level
        
        self.conc_bound = conc_bound
        self.lnconc_bounds = tuple(np.log(self.conc_bound))
        
        if spec_conc_bound is None:
            self.spec_conc_bound = {}  
        else:
            self.spec_conc_bound = spec_conc_bound
        self.spec_lnconc_bounds = {metabid: tuple(np.log(bounds)) 
                                   for metabid, bounds 
                                   in self.spec_conc_bound.items()}
        
        if preset_conc is None:
            self.preset_conc = {}  
        else:
            self.preset_conc = preset_conc

        if preset_conc_ratio is None:
            self.preset_conc_ratio = {} 
        else:
            self.preset_conc_ratio = preset_conc_ratio
        
        if ex_thermo_cons is None:
            self.ex_thermo_cons = []
        else:
            self.ex_thermo_cons = list(set(ex_thermo_cons))
        
        self.cstrFluxIDs = []
        for rxnid in self.rxnIDs:
            if all([not self.model.reactions[rxnid].is_h2o_transport,
                    not self.model.reactions[rxnid].is_biomass_formation,
                    not self.model.reactions[rxnid].is_exch_reaction,
                    rxnid not in self.ex_thermo_cons]):

                if self.model.reactions[rxnid].rev:
                    self.cstrFluxIDs.append(rxnid+'_f')
                    self.cstrFluxIDs.append(rxnid+'_b')
                else:
                    self.cstrFluxIDs.append(rxnid)
                self.model.reactions[rxnid].is_constrained_by_thermodynamics = True
        # cstrFluxIDs denotes total fluxes subject to thermodynamic constraints
        # H2O transport, biomass formation, exchange reactions, and reactions 
        # in ex_thermo_cons are excluded
        # it is a subset of varFluxIDs

        if ex_conc is None:
            self.ex_conc = []
        else: 
            self.ex_conc = list(set(ex_conc))

        totalStoyMat_reduced = self.model.total_stoichiometric_matrix[
            self.cstrFluxIDs
        ]
        self.totalStoyMat_reduced = totalStoyMat_reduced[
            ~(totalStoyMat_reduced == 0).all(axis = 1)
        ]
        # totalStoyMat_reduced only includes reactions subject to thermodynamic 
        # constraints as well as related metabolites
        self.varMetabIDs = []
        for metabid in self.totalStoyMat_reduced.index:
            if all([not self.model.metabolites[metabid].is_h2o,   # is_h
                    metabid not in self.ex_conc]):
                self.varMetabIDs.append(metabid)
        # varMetabIDs represent variable metabolites found in reactions that are 
        # subjected to thermodynamic constraints.
        # H2O and metabolites in ex_conc are excluded

        self.pyoModel.varMetabIDs = Set(initialize = self.varMetabIDs)
        self.pyoModel.cstrFluxIDs = Set(initialize = self.cstrFluxIDs)
        
    
    def _build_conc_variables(self, initial = None):
        '''
        Parameters
        ----------
        initial: dict
            Mapping of reaction IDs to their corresponding initial values.
        '''
        
        def conc_bounds_rule(model, metabid):
            if metabid in self.preset_conc:   
                # set the bounds of fixed concentrations
                return (np.log(self.preset_conc[metabid]),)*2
            elif metabid in self.spec_lnconc_bounds:
                return self.spec_lnconc_bounds[metabid]
            else:
                return self.lnconc_bounds
            
        if initial is None:
            self.pyoModel.lnconcs = Var(
                self.pyoModel.varMetabIDs, 
                bounds = conc_bounds_rule
            )
        else:
            logging.info('load initial conc. values')
            
            initial = {metabid: conc for metabid, conc in initial.items() 
                       if metabid in self.pyoModel.varMetabIDs}
            self.pyoModel.lnconcs = Var(
                self.pyoModel.varMetabIDs, 
                bounds = conc_bounds_rule, 
                initialize = initial
            )
                
    
    def _build_binary_variables(self):
        self.pyoModel.xs = Var(self.pyoModel.cstrFluxIDs, within = Binary)


    def _build_error_variables(self):
        z = norm.ppf((1+self.conf_level) / 2)
        
        def error_bounds_rule(model, fluxid):
            rxnid = re.sub(r'_[fb]$', '', fluxid)
            dgpm_err = self.model.reactions[rxnid].dgpm_error
            return [-z*dgpm_err, z*dgpm_err]

        self.pyoModel.errors = Var(
            self.pyoModel.cstrFluxIDs, 
            within = Reals, 
            bounds = error_bounds_rule
        )
        
        
    def _calculate_gibbs_energy(self, model, fluxid):
        '''
        Parameters
        ----------
        model: pyomo model
            The Pyomo model object.
        fluxid: str
            Flux ID.
        '''
        
        rxnid = re.sub(r'_[fb]$', '', fluxid)

        subs = self.model.reactions[rxnid].substrates
        pros = self.model.reactions[rxnid].products
        dgpm = self.model.reactions[rxnid].dgpm

        subsSum = sum([subs[subid].coe*model.lnconcs[subid] for subid in subs 
                       if subid in self.varMetabIDs])
        prosSum = sum([pros[proid].coe*model.lnconcs[proid] for proid in pros 
                       if proid in self.varMetabIDs])

        dgp = dgpm + (prosSum - subsSum)*R*T

        if re.match(r'.+_b$', fluxid):   # backward reaction
            return -dgp
        else:   # forward reaction and irreversible reaction
            return dgp
    

    def _build_flux_bound_constraints(self):
        def bound_rule(model, fluxid):
            return (
                model.fluxes[fluxid] 
                <= model.xs[fluxid]*model.fluxes[fluxid].bounds[1]
            )

        self.pyoModel.FLUXBNDcstr = Constraint(
            self.pyoModel.cstrFluxIDs, 
            rule = bound_rule
        )


    def _build_thermodynamics_constraints(self):
        def thmd_rule(model, fluxid):
            return (
                self._calculate_gibbs_energy(model, fluxid) 
                <= K*(1-model.xs[fluxid]) - EPSILON
            )

        self.pyoModel.THMDcstr = Constraint(
            self.pyoModel.cstrFluxIDs, 
            rule = thmd_rule
        )


    def _build_thermodynamics_constraints_with_uncertainty(self):
        def thmd_rule(model, fluxid):
            return (
                self._calculate_gibbs_energy(model, fluxid) + model.errors[fluxid] 
                <= K*(1-model.xs[fluxid]) - EPSILON
            )

        self.pyoModel.THMDcstr = Constraint(
            self.pyoModel.cstrFluxIDs, 
            rule = thmd_rule
        )    

    
    def _build_ratio_constraint(self):
        if self.preset_conc_ratio:
            def ratio_rule(model, ratioid):
                num, den = ratioid.split(':')
                return (
                    model.lnconcs[num] - model.lnconcs[den] 
                    == log(self.preset_conc_ratio[ratioid])
                )
            
            self.pyoModel.RATIOcstr = Constraint(
                self.preset_conc_ratio.keys(), 
                rule = ratio_rule
            )


    def _get_opt_lnconcs(self):
        optLnconcs = {metabid: value(self.pyoModel.lnconcs[metabid]) 
                      for metabid in self.pyoModel.varMetabIDs}
                
        return optLnconcs
    

    def _get_opt_gibbis_energies(self, error = False):
        fluxids_filtered = list(
            filter(lambda fluxid: not re.match(r'.+_b$', fluxid), self.cstrFluxIDs)
        )   # no need to estimate backward reactions
        optDgps = {}
        for fluxid in fluxids_filtered:
            opt_dgp = value(self._calculate_gibbs_energy(self.pyoModel, fluxid))
            if error:
                opt_dgp_err = value(self.pyoModel.errors[fluxid])
            else:
                opt_dgp_err = 0
                    
            rxnid = re.sub(r'_[fb]$', '', fluxid)
            optDgps[rxnid] = opt_dgp + opt_dgp_err
                        

            
        # optDgps = {}
        # for fluxid in self.pyoModel.cstrFluxIDs:
            
        #     rxnid = re.sub(r'_[fb]$', '', fluxid)
        #     if rxnid in optDgps:
        #         continue
        #     else:
        #         if self.model.reactions[rxnid].rev:
        #             direction = '_f' if optFluxes[rxnid] >= 0 else '_b'
        #             opt_dgp = value(self._calculate_gibbs_energy(self.pyoModel, rxnid+direction))
        #             if error:
        #                 opt_dgp_err = value(self.pyoModel.errors[rxnid+direction])
        #             else:
        #                 opt_dgp_err = 0
        #         else:
        #             opt_dgp = value(self._calculate_gibbs_energy(self.pyoModel, rxnid))
        #             if error:
        #                 opt_dgp_err = value(self.pyoModel.errors[rxnid])
        #             else:
        #                 opt_dgp_err = 0
                
        #         optDgps[rxnid] = opt_dgp + opt_dgp_err

        return optDgps
            
                
    def solve(self, solver = 'glpk'):    
        '''
        Parameters
        ----------
        solver: {"glpk", "gurobi"}
            "gurobi" is highly recommended for large models.
        '''

        if_error = not bool(self.conf_level is None)

        self._build_flux_variables()
        self._build_conc_variables()
        self._build_binary_variables()
        self._build_objective()
        self._build_mass_balance_contraints()
        self._build_flux_bound_constraints()
        self._build_ratio_constraint()
        if self.conf_level is None:
            self._build_thermodynamics_constraints()
        else:
            self._build_error_variables()
            self._build_thermodynamics_constraints_with_uncertainty()
        
        sol = self._get_solver(solver)
        self.res = sol.solve(self.pyoModel, report_timing = False)
        #self.pyoModel.pprint()#!!!
        #self.pyoModel.display()#!!!
        optObj = self._get_opt_obj()
        optFluxes = self._get_opt_fluxes()
        optLnconcs = self._get_opt_lnconcs()    
        optDgps = self._get_opt_gibbis_energies(if_error)
        optSuccess = self._optimization_successful()
            
        if self.parsimonious and optSuccess:
            logging.info('estimating parsimonious fluxes')
            
            self._remove_objective()
            self._build_parsimonious_objective()
            self._build_objective_constraint(optObj)

            self.res = sol.solve(self.pyoModel, report_timing = False)

            optFluxes = self._get_opt_fluxes()
            optLnconcs = self._get_opt_lnconcs()    
            optDgps = self._get_opt_gibbis_energies(if_error)
            optSuccess = self._optimization_successful()
        
        return TFBAResults(
            optObj, 
            optFluxes, 
            optLnconcs, 
            optDgps, optSuccess, 
            self.model.stoichiometric_matrix
        )


class EFBAOptimizer(FBAOptimizer):
    '''
    EFBA computes the fluxes by optimizing the objective (e.g., biomass formation) 
    under the constraints of mass balance and the enzyme protein allocation.
    '''

    def __init__(
            self, 
            model, 
            objective, 
            direction, 
            flux_bound, 
            spec_flux_bound, 
            preset_flux,
            irr_reactions, 
            ex_mass_bal_cons, 
            inc_enz_cons, 
            enz_prot_lb,
            parsimonious,
            slack, 
            **kwargs
    ):
        '''
        Parameters
        ----------
        model: Model
            The model that calls EFBAOptimizer.
        objective : dict
            Mapping of reaction IDs to coefficients in the objective expression.
        direction : {'max', 'min'}
            Direction of optimization.
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
            List of reactions included in the enzyme protein cost constraint.
        enz_prot_lb : float
            Upper bound of enzyme protein fraction.
        '''
        
        super().__init__(
            model, 
            objective = objective, 
            direction = direction, 
            flux_bound = flux_bound, 
            spec_flux_bound = spec_flux_bound, 
            preset_flux = preset_flux, 
            irr_reactions = irr_reactions, 
            ex_mass_bal_cons = ex_mass_bal_cons,
            parsimonious = parsimonious,
            slack = slack, 
            **kwargs
        )
        
        if inc_enz_cons is None:
            self.inc_enz_cons = []
        else:
            self.inc_enz_cons = list(set(inc_enz_cons))
        self.q = enz_prot_lb


    def _calculate_enzyme_cost(self, model, rxnid):
        '''
        Parameters
        ----------
        model: pyomo model
            The Pyomo model object.
        rxnid: str
            Reaction ID.
        '''
        
        if self.model.reactions[rxnid].rev:
            fflux = model.fluxes[rxnid+'_f']
            bflux = model.fluxes[rxnid+'_b']
            fkcat = self.model.reactions[rxnid].fkcat
            bkcat = self.model.reactions[rxnid].bkcat
            e = fflux/fkcat + bflux/bkcat
        else:
            flux = model.fluxes[rxnid]
            kcat = self.model.reactions[rxnid].fkcat
            e = flux/kcat
        
        mw = self.model.reactions[rxnid].mw
        cost = 1/3600*mw*e

        return cost

    
    def _build_enzyme_cost_constraint(self):
        for rxnid in self.inc_enz_cons:
            if self.model.reactions[rxnid].is_biomass_formation:
                raise ValueError(
                    "biomass formation can't be included in enzyme protein cost"
                )
                
            if self.model.reactions[rxnid].is_exch_reaction:
                raise ValueError(
                    f"exchange reaction {rxnid} can't be included in "
                    "enzyme protein cost"
                )
                
        def epc_rule(model):
            costs = [self._calculate_enzyme_cost(model, rxnid) 
                     for rxnid in self.inc_enz_cons]
            return (0, sum(costs), self.q)
        
        self.pyoModel.EPCcstr = Constraint(rule = epc_rule)


    def _get_opt_enzyme_protein_cost(self):
        optEcosts = {
            rxnid: value(self._calculate_enzyme_cost(self.pyoModel, rxnid)) 
            for rxnid in self.inc_enz_cons
        }
        optTotalEcost = sum(optEcosts.values())

        return optTotalEcost, optEcosts
    

    def solve(self, solver = 'glpk'):
        '''
        Parameters
        ----------
        solver: {"glpk", "gurobi"}
            "gurobi" is highly recommended for large models.
        '''

        self._build_flux_variables()
        self._build_objective()
        self._build_mass_balance_contraints()
        self._build_enzyme_cost_constraint()
        
        sol = self._get_solver(solver)
        self.res = sol.solve(self.pyoModel, report_timing = False, tee = False)
        
        #self.pyoModel.display()#!!!
        #self.pyoModel.pprint()#!!!
        optObj = self._get_opt_obj()
        optFluxes = self._get_opt_fluxes()
        optTotalEcost, optEcosts = self._get_opt_enzyme_protein_cost()
        optSuccess = self._optimization_successful()

        if self.parsimonious and optSuccess:
            logging.info('estimating parsimonious fluxes')
            
            self._remove_objective()
            self._build_parsimonious_objective()
            self._build_objective_constraint(optObj)

            self.res = sol.solve(self.pyoModel, report_timing = False)

            optFluxes = self._get_opt_fluxes()    
            optTotalEcost, optEcosts = self._get_opt_enzyme_protein_cost()
        
        return EFBAResults(
            optObj, 
            optFluxes, 
            optTotalEcost, 
            optEcosts, 
            optSuccess,
            self.model.stoichiometric_matrix
        )


class ETFBAOptimizer(TFBAOptimizer, EFBAOptimizer):
    '''
    ETFBA computes the fluxes by optimizing the objective (e.g., biomass formation) 
    under the constraints of mass balance, thermodynamics, and enzyme protein 
    allocation.
    '''

    def __init__(
            self, 
            model, 
            objective, 
            direction, 
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
            parsimonious,
            slack,
            dgpm_conf_level
    ):
        '''
        Parameters
        ----------
        model: Model
            The model that calls ETFBAOptimizer.
        objective : dict
            Mapping of reaction IDs to coefficients in the objective expression.
        direction : {'max', 'min'}
            Direction of optimization.
        flux_bound : tuple
            Lower and upper bounds of metabolic flux.
        conc_bound : tuple
            Lower and upper bounds of metabolite concentration.
        spec_flux_bound : dict
            Mapping of reaction IDs to their specific flux bounds (lb, ub).
        spec_conc_bound : dict
            Mapping of metabolite IDs to their specific concentration bounds 
            (lb, ub).
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
            List of reactions included in the enzyme protein cost constraint.
        enz_prot_lb : float
            Upper bound of enzyme protein fraction.
        dgpm_conf_level : float
            Confidence level for considering uncertainty in standard reaction Gibbs 
            energy.
        '''

        super().__init__(
            model = model, 
            objective = objective, 
            direction = direction,
            flux_bound = flux_bound, 
            conc_bound = conc_bound, 
            spec_flux_bound = spec_flux_bound, 
            spec_conc_bound = spec_conc_bound, 
            preset_flux = preset_flux, 
            preset_conc = preset_conc,
            preset_conc_ratio = preset_conc_ratio, 
            irr_reactions = irr_reactions, 
            ex_conc = ex_conc, 
            ex_mass_bal_cons = ex_mass_bal_cons, 
            ex_thermo_cons = ex_thermo_cons, 
            inc_enz_cons = inc_enz_cons, 
            enz_prot_lb = enz_prot_lb,
            parsimonious = parsimonious,
            slack = slack,
            dgpm_conf_level = dgpm_conf_level
        )


    def solve(self, solver = 'glpk'):
        '''
        Parameters
        ----------
        solver: {"glpk", "gurobi"}
            "gurobi" is highly recommended for large models.
        '''
        
        self._build_flux_variables()
        self._build_conc_variables()
        self._build_binary_variables()
        self._build_objective()
        self._build_mass_balance_contraints()
        self._build_flux_bound_constraints()
        self._build_ratio_constraint()
        if self.conf_level is None:
            self._build_thermodynamics_constraints()
        else:
            self._build_error_variables()
            self._build_thermodynamics_constraints_with_uncertainty()
        self._build_enzyme_cost_constraint()
        
        sol = self._get_solver(solver)
        self.res = sol.solve(self.pyoModel, report_timing = False, tee = False)
        
        #self.pyoModel.display()#!!!
        #self.pyoModel.pprint()#!!!
        optObj = self._get_opt_obj()
        optFluxes = self._get_opt_fluxes()
        optLnconcs = self._get_opt_lnconcs()    
        optDgps = self._get_opt_gibbis_energies(not self.conf_level is None)
        optTotalEcost, optEcosts = self._get_opt_enzyme_protein_cost()
        optSuccess = self._optimization_successful()

        if self.parsimonious and optSuccess:
            logging.info('estimating parsimonious fluxes')
            
            self._remove_objective()
            self._build_parsimonious_objective()
            self._build_objective_constraint(optObj)

            self.res = sol.solve(self.pyoModel, report_timing = False)

            optFluxes = self._get_opt_fluxes()
            optLnconcs = self._get_opt_lnconcs()    
            optDgps = self._get_opt_gibbis_energies(not self.conf_level is None)    
            optTotalEcost, optEcosts = self._get_opt_enzyme_protein_cost()
        
        return ETFBAResults(
            optObj, 
            optFluxes, 
            optLnconcs, 
            optDgps, 
            optTotalEcost, 
            optEcosts, 
            optSuccess,
            self.model.stoichiometric_matrix
        )
