'''Difine classes of analysis results.'''


from math import exp
from .io import save_values


class PrettyDict(dict):
    
    def __init__(self, *args, ndigits=3, **kwargs):
        super().__init__(*args, **kwargs)
        self.ndigits = ndigits
        
    
    def save(self, file, exp_transform=False):
        '''
        Parameters
        ----------
        file: str
            File ending with .xlsx, .tsv, or .bin format.
        log_transform: bool
            Determines whether to apply natural exponential transformation.
        '''
        
        save_values(file, self, exp_transform)
        
    
    def __repr__(self):
        itemsStr = []
        for key, value in self.items():
            if isinstance(value, list):
                itemsStr.append(f'{key}: [{round(value[0], self.ndigits)}, '
                                f'{round(value[1], self.ndigits)}]')
            elif isinstance(value, self.__class__):
                itemsStr.append(f'{key}:\n{self.__class__(value)}')
            else:
                if isinstance(value, float):
                    itemsStr.append(f'{key}: {round(value, self.ndigits)}')
                else:
                    itemsStr.append(f'{key}: {value}')
        
        return '\n'.join(itemsStr)
        
        
class FBAResults:
    '''
    Attributes
    ----------
    opt_objective: float
        Optimal objective value achieved during optimization.
    opt_fluxes: dict
        Dictionary mapping reaction ID to its optimal flux value.
    optimization_successful: bool
        Indicates whether the optimization process was successful.
    '''
    
    def __init__(self, opt_obj, opt_fluxes, opt_success, stoy_mat):
        '''
        Parameters
        ----------
        opt_obj: float
            Optimal objective value achieved during optimization.
        opt_fluxes: dict
            Dictionary mapping reaction ID to its optimal flux value.
        opt_success: bool
            Indicates whether the optimization process was successful (True if 
            successful).
        stoy_mat: DataFrame
            Stoichiometric matrix.
        '''
        
        self._opt_obj = opt_obj
        self._opt_fluxes = opt_fluxes
        self._opt_success = opt_success
        self._stoy_mat = stoy_mat
    
    
    @property
    def opt_objective(self):
        return round(self._opt_obj, 3)
        
    
    @property
    def opt_fluxes(self):
        return PrettyDict(self._opt_fluxes)   
        
        
    @property
    def optimization_successful(self):
        return self._opt_success

    
    def statement(self, metabid):
        '''
        Parameters
        ----------
        metabid: str
            Metabolite ID.
        '''

        productions = PrettyDict()
        consumptions = PrettyDict()
        for rxnid in self._stoy_mat.columns[self._stoy_mat.loc[metabid, :] != 0]:
            
            if self._stoy_mat.loc[metabid, rxnid]*self._opt_fluxes[rxnid] > 0:
                productions[rxnid] = self._opt_fluxes[rxnid]
            elif self._stoy_mat.loc[metabid, rxnid]*self._opt_fluxes[rxnid] < 0:
                consumptions[rxnid] = self._opt_fluxes[rxnid]

        return PrettyDict({
            'productions': productions, '\nconsumptions': consumptions
        })

    
    # def __repr__(self):
    
    #     return (
    #         f'optimal objective: {self.opt_objective}\n'
    #         f'\noptimal fluxes\n{str(self.opt_fluxes)}\n'
    #     )
        

class TFBAResults(FBAResults):
    '''
    Attributes
    ----------
    ...
    opt_concentrations: dict
        Dictionary mapping metabolite ID to its optimal concentration.
    opt_directions: dict
        Dictionary mapping reaction ID to its direction ("f" for forward, "b" for 
        backward).
    opt_gibbs_energy: dict
        Dictionary mapping reaction ID to its optimal deltaGprime value.
    '''
    
    def __init__(
            self, 
            opt_obj, 
            opt_fluxes, 
            opt_lnconcs, 
            opt_dgps, 
            opt_success,
            stoy_mat, 
            **kwargs
    ):
        '''
        Parameters
        ----------
        opt_obj: float
            Optimal objective value achieved by the optimization process.
        opt_fluxes: dict
            Dictionary mapping reaction ID to its optimal flux value.
        opt_lnconcs: dict
            Dictionary mapping metabolite ID to its optimal natural logarithm 
            concentration.
        opt_dgps: dict
            Dictionary mapping reaction ID to its optimal Gibbs energy change.
        opt_success: bool
            Boolean indicating whether the optimization was successful (True) or not 
            (False).
        stoy_mat: DataFrame
            Stoichiometric matrix.
        '''
        
        super().__init__(
            opt_obj=opt_obj, 
            opt_fluxes=opt_fluxes, 
            opt_success=opt_success,
            stoy_mat=stoy_mat,
            **kwargs
        )
        
        self._opt_lnconcs = opt_lnconcs
        self._opt_concs = {metabid: exp(lnconc) 
                           for metabid, lnconc in self._opt_lnconcs.items()}
        self._opt_dgps = opt_dgps
        
        
    @property
    def opt_concentrations(self):
        return PrettyDict(self._opt_concs)
    
        
    @property
    def opt_directions(self):
        directions = {}
        for rxnid, flux in self.opt_fluxes.items():
            if flux > 0:
                directions[rxnid] = 'forward'
            elif flux < 0:
                directions[rxnid] = 'reverse'
            else:
                directions[rxnid] = 'zero flux'
        
        return PrettyDict(directions)
    
    
    @property
    def opt_gibbs_energy(self):
        return PrettyDict(self._opt_dgps)
        
    
    # def __repr__(self):

    #     return (
    #         f'optimal objective: {self.opt_objective}\n'
    #         f'\noptimal fluxes\n{str(self.opt_fluxes)}\n'
    #         f'\noptimal concentrations\n{str(self.opt_concentrations)}\n'
    #     )       
        
        
class EFBAResults(FBAResults):
    '''
    Attributes
    ----------
    ...
    opt_total_enzyme_cost: float
        Optimal total enzyme protein cost achieved.
    opt_enzyme_costs: dict
        Dictionary mapping reaction ID to its optimal enzyme protein abundance.
    '''

    def __init__(
            self, 
            opt_obj, 
            opt_fluxes, 
            opt_total_epc, 
            opt_epcs, 
            opt_success, 
            stoy_mat,
            **kwargs
    ):
        '''
        Parameters
        ----------
        opt_obj: float
            Optimal objective value obtained from the optimization.
        opt_fluxes: dict
            Dictionary containing reaction ID mapped to its optimal flux value.
        opt_total_epc: float
            Optimal total enzyme protein abundance achieved.
        opt_epcs: dict
            Dictionary mapping reaction ID to its optimal enzyme protein abundance.
        opt_success: bool
            Boolean indicating whether the optimization was successful (True) or not 
            (False).
        stoy_mat: DataFrame
            Stoichiometric matrix.
        '''

        super().__init__(
            opt_obj=opt_obj, 
            opt_fluxes=opt_fluxes, 
            opt_success=opt_success,
            stoy_mat=stoy_mat, 
            **kwargs
        )

        self._opt_total_epc = opt_total_epc
        self._opt_epcs = opt_epcs


    @property
    def opt_total_enzyme_cost(self):
        return round(self._opt_total_epc, 5)
    
    
    @property
    def opt_enzyme_costs(self):
        return PrettyDict(self._opt_epcs, ndigits=5)   
    

class ETFBAResults(TFBAResults, EFBAResults):
    '''
    Attributes
    ----------
    ...    
    '''
    
    def __init__(
            self, 
            opt_obj, 
            opt_fluxes, 
            opt_lnconcs, 
            opt_dgps, 
            opt_total_epc, 
            opt_epcs, 
            opt_success,
            stoy_mat
    ):
        '''
        Parameters
        ----------
        opt_obj: float
            Optimal objective value achieved.
        opt_fluxes: dict
            Dictionary mapping reaction ID to its optimal flux value.
        opt_lnconcs: dict
            Dictionary mapping metabolite ID to its optimal natural logarithm 
            concentration.
        opt_dgps: dict
            Dictionary mapping reaction ID to its optimal reaction Gibbs energy 
            change.
        opt_total_epc: float
            Optimal total enzyme protein abundance achieved.
        opt_epcs: dict
            Dictionary mapping reaction ID to its optimal enzyme protein abundance.
        opt_success: bool
            Indicates whether the optimization process was successful.
        stoy_mat: df
            Stoichiometric matrix.    
        '''
        
        super().__init__(
            opt_obj=opt_obj, 
            opt_fluxes=opt_fluxes, 
            opt_lnconcs=opt_lnconcs, 
            opt_dgps=opt_dgps, 
            opt_total_epc=opt_total_epc, 
            opt_epcs=opt_epcs,
            opt_success=opt_success,
            stoy_mat=stoy_mat
        )
    
"""
class EPCMResults(TFBAResults):
    '''
    Attributes
    ----------
    ...
    opt_total_enzyme_cost: float
        optimal total enzyme protein abundance
    opt_enzyme_costs: dict
        reaction ID => optimal enzyme protein abundance
    opt_enzyme_base_costs
        reaction ID => optimal enzyme protein abundance (base term)
    opt_enzyme_thermodynamic_costs
        reaction ID => optimal enzyme protein abundance (thermodynamic term)
    opt_enzyme_kinetic_costs
        reaction ID => optimal enzyme protein abundance (kinetic term)    
    '''
    
    def __init__(
            self, 
            opt_obj, 
            opt_fluxes, 
            opt_lnconcs, 
            opt_dgps, 
            opt_total_epc, 
            opt_epcs, 
            opt_epcs_base, 
            opt_epcs_thmd, 
            opt_epcs_kin, 
            opt_success
    ):
        '''
        Parameters
        ----------
        opt_obj: float
            optimal objective
        opt_fluxes: dict
            reaction ID => optimal flux value
        opt_lnconcs: dict
            metabolite ID => optimal concentration
        opt_dgps: dict
            reaction ID => optimal reaction Gibbs energy change
        opt_total_epc: float
            optimal total enzyme protein abundance
        opt_epcs: dict
            reaction ID => optimal enzyme protein abundance
        opt_epcs_base: dict
            reaction ID => optimal enzyme protein abundance (base term)
        opt_epcs_thmd: dict
            reaction ID => optimal enzyme protein abundance (thermodynamic term)
        opt_epcs_kin: dict
            reaction ID => optimal enzyme protein abundance (kinetic term)
        opt_success: bool
            True if optimization is successful    
        '''
        
        super().__init__(
            opt_obj, 
            opt_fluxes, 
            opt_lnconcs, 
            opt_dgps, 
            opt_success
        )
        
        self._opt_total_epc = opt_total_epc
        self._opt_epcs = opt_epcs

        self._opt_epcs_base = opt_epcs_base
        self._opt_epcs_thmd = opt_epcs_thmd
        self._opt_epcs_kin = opt_epcs_kin
        
        
    @property
    def opt_total_enzyme_cost(self):
        
        return round(self._opt_total_epc, 5)
    
    
    @property
    def opt_enzyme_costs(self):
        
        return PrettyDict(self._opt_epcs, ndigits = 5)


    @property
    def opt_enzyme_base_costs(self):
        
        return PrettyDict(self._opt_epcs_base, ndigits = 5)


    @property
    def opt_enzyme_thermodynamic_costs(self):
        
        return PrettyDict(self._opt_epcs_thmd, ndigits = 5)


    @property
    def opt_enzyme_kinetic_costs(self):
        
        return PrettyDict(self._opt_epcs_kin, ndigits = 5)
"""

class VariabilityResults():
    '''
    Attributes
    ----------
    objective_value: float
        Objective value obtained from TFBA, EFBA, or ETFBA analysis.
    gamma: float
        Control parameter for the objective value.
    '''

    def __init__(self, obj_value, gamma, ranges):
        '''
        Parameters
        ----------
        obj_value: float
            Objective value obtained from the analysis.
        gamma: float
            The value indicating the constraint on the objective. The objective 
            should be no less than gamma*obj_value, or no greater than 
            (1+gamma)*obj_value in the variability analysis, based on the direction.
        ranges: dict
            Dictionary containing reaction ID mapped to its corresponding flux 
            ranges [lb, ub].
        '''
        
        self._obj_value = obj_value
        self._gamma = gamma
        self._ranges = ranges


    @property
    def objective_value(self):
        return round(self._obj_value)


    @property
    def gamma(self):
        return self._gamma
    

class FVAResults(VariabilityResults):
    '''
    Attributes
    ----------
    ...
    flux_ranges: dict
        Dictionary containing reaction IDs mapped to their corresponding net flux 
        ranges [lb, ub].
    '''

    @property
    def flux_ranges(self):
        return PrettyDict(self._ranges)
    

class TVAResults(VariabilityResults):
    '''
    Attributes
    ----------
    ...
    gibbs_energy_ranges: dict
        Dictionary containing reaction IDs mapped to their corresponding ranges of 
        deltaGprime [lb, ub].
    '''

    @property
    def gibbs_energy_ranges(self):
        return PrettyDict(self._ranges)
    

class EVAResults(VariabilityResults):
    '''
    Attributes
    ----------
    ...
    protein_cost_ranges: dict
        Dictionary mapping reaction IDs to their respective ranges of protein cost 
        [lb, ub].
    '''

    @property
    def protein_cost_ranges(self):
        return PrettyDict(self._ranges)