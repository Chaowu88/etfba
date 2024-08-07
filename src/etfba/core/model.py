'''Define the Model class.'''


import re
from functools import lru_cache
from collections.abc import Iterable
import pandas as pd
from .reaction import Reaction, DEFAULT_MW, DEFAULT_KCAT, DEFAULT_KM, DEFAULT_DGPM
from .metabolite import Metabolite
from ..optim.optim import FBAOptimizer, TFBAOptimizer, EFBAOptimizer, ETFBAOptimizer
from ..optim.variability import (FVAOptimizer, TFVAOptimizer, EFVAOptimizer, 
                                 ETFVAOptimizer, TVAOptimizer, ETVAOptimizer, 
                                 EVAOptimizer, TEVAOptimizer)
from ..io.results import PrettyDict
from ..io.io import load_model, save_model


class Model():
    '''
    Attributes
    ----------
    metabolites : PrettyDict
        A dictionary mapping metabolite IDs to corresponding Metabolite objects.
    reactions : PrettyDict
        A dictionary mapping reaction IDs to corresponding Reaction objects.
    end_metabolites : PrettyDict
        A dictionary mapping metabolite IDs to Metabolite objects representing 
        initial substrates or final products within the model.
    stoichiometric_matrix : DataFrame
        Represents the stoichiometric matrix where rows correspond to metabolites 
        and columns correspond to net reactions. Negative values indicate 
        substrates, while positive values denote products.
    total_stoichiometric_matrix : DataFrame
        Represents the stoichiometric matrix where rows correspond to metabolites 
        and columns correspond to total reactions (reversible reactions are split 
        into forward and backward reactions). Negative values denote substrates, 
        while positive values denote products.
    transformation_matrix : DataFrame
        Matrix facilitating the conversion of total fluxes into net fluxes.
    '''
    
    def __init__(self, name=None):
        '''
        Parameters
        ----------
        name: str
            Model name.
        '''
        
        self.name = name
        
        self._metabolites = PrettyDict()
        self._reactions = PrettyDict()
        
    
    @classmethod
    def load(cls, filename):
        '''
        Parameters
        ----------
        filename: str
            Filename of the model file. The filename should end with the extension 
            '.bin'.
        '''
        
        return load_model(filename)


    def save(self, filename):
        '''
        Parameters
        ----------
        filename: str
            Filename of the model to save. The filename should end with the 
            extension '.bin'.
        '''

        save_model(filename, self)


    @staticmethod
    def _set_value(value, default):
        return default if value == '' else float(value)
        

    def _build_reactant(self, rxn, reacsStr, reackms, label):
        reacStrLst = reacsStr.split(';')
        reackmLst = reackms.split(';')
        if len(reacStrLst) != len(reackmLst) and not rxn.is_biomass_formation:
            if label == 'substrate':
                raise ValueError(
                    f'the number of subtrates in {rxn.rxnid} does not match '
                    f'the number of subtrate Km values'
                )
            elif label == 'product':
                raise ValueError(
                    f'the number of products in {rxn.rxnid} does not match '
                    f'the number of product Km values'
                )

        for idx, reacStr in enumerate(reacStrLst):
                
            coe_reac = reacStr.split()
            if len(coe_reac) == 0:   # no substrate or product in exchange reaction
                continue
            elif len(coe_reac) == 1:
                coe, reacid = 1.0, coe_reac[0]
            else:
                coe, reacid = coe_reac
            
            reac = self._metabolites.setdefault(reacid, Metabolite(reacid))
            if label == 'substrate':
                reac.coes[rxn.rxnid] = -float(coe)
            elif label == 'product':
                reac.coes[rxn.rxnid] = float(coe)

            if rxn.is_biomass_formation or rxn.is_exch_reaction:   
            # reactants in biomass formation and exchange reactions don't have Km 
            # values
                reac.kms[rxn.rxnid] = None
            else:
                reac.kms[rxn.rxnid] = self._set_value(reackmLst[idx], DEFAULT_KM)

            if re.match(r'^h($|\.[\w\._]+$)', reacid, flags=re.I):
                reac.is_h = True

            if re.match(r'^h2o($|\.[\w\._]+$)', reacid, flags=re.I):
                reac.is_h2o = True 

            if label == 'substrate':
                rxn.substrates[reacid] = reac
            elif label == 'product':
                rxn.products[reacid] = reac


    def read_from_excel(self, filename):
        '''
        Parameters
        ----------
        filename: str
            Filename of an Excel file containing the following fields: 
            Enzyme, Substrates, Products, Sub Kms (mM), Pro Kms (mM), 
            Fwd kcat (1/s), Bwd kcat (1/s), MW (kDa), and ΔrG'm (kJ/mol).
        '''
        
        data = pd.read_excel(
            filename, 
            header=0, 
            index_col=0, 
            comment='#'
        ).fillna('').astype(str)
        
        for rxnid, rowInfos in data.iterrows():
            subsStr, prosStr, rev, subkms, prokms, fkcat, bkcat, mw, dgpm = rowInfos
            
            fkcat = self._set_value(fkcat, DEFAULT_KCAT)
            bkcat = self._set_value(bkcat, DEFAULT_KCAT)
            rev = float(rev)
            mw = self._set_value(mw, DEFAULT_MW)
            dgpm = self._set_value(dgpm, DEFAULT_DGPM)
            rxn = Reaction(
                rxnid, 
                forward_kcat=fkcat, 
                backward_kcat=bkcat, 
                reversible=bool(rev),
                molecular_weight=mw, 
                standard_gibbs_energy=dgpm
            )
            
            if re.search(r'biom', rxnid, flags=re.I):
                rxn.is_biomass_formation = True    
            
            if subsStr == '' or prosStr == '':
                rxn.is_exch_reaction = True
            
            if re.match(r'^h\.e$', subsStr) or re.match(r'^h\.e$', prosStr):
                rxn.is_h_transport = True

            if re.match(r'^h2o\.e$', subsStr) or re.match(r'^h2o\.e$', prosStr):
                rxn.is_h2o_transport = True

            self._build_reactant(rxn, subsStr, subkms, 'substrate')
            self._build_reactant(rxn, prosStr, prokms, 'product')
            
            self._reactions[rxnid] = rxn
    
    
    def add_reactions(self, reactions):
        '''
        Parameters
        ----------
        reactions: Reaction or list of Reactions
            Reaction objects to be added into the model.
        '''

        if not isinstance(reactions, Iterable):
            reactions = [reactions]

        for rxn in reactions:
            self._reactions[rxn.rxnid] = rxn

            for subid in rxn.substrates:
                self._metabolites[subid] = rxn.substrates[subid]

            for proid in rxn.products:
                self._metabolites[proid] = rxn.products[proid]


    def remove_reactions(self, reactions):
        '''
        Parameters
        ----------
        reactions: Reaction or list of Reactions
            Reaction objects to be removed from the model.
        '''
        
        if not isinstance(reactions, Iterable):
            reactions = [reactions]
        
        for rxn in reactions:
            del self._reactions[rxn.rxnid]

            for subid in rxn.substrates:
                del self._metabolites[subid]

            for proid in rxn.products:
                del self._metabolites[proid]


    @property
    def metabolites(self):
        if len(self._metabolites) == 0:
            raise AttributeError('no metabolite found, model empty')
        else:
            return self._metabolites
            
            
    @property
    def reactions(self):
        if len(self._reactions) == 0:
            raise AttributeError('no reaction found, model empty')
        else:
            return self._reactions
            
    
    @lru_cache()   
    # Re-estimates the stoichiometric matrix if changes occur in 
    # self._metabolites or self._reactions, otherwise uses cached results.
    def _get_stoichiometric_matrix(self, metabolites, reactions):
        '''
        Parameters
        ----------
        metabolites: tuple
            A tuple of metabolite IDs.
        reactions: tuple
            A tuple of reaction IDs.
        '''
        
        stoyMat_net = pd.DataFrame(0.0, index=metabolites, columns=reactions)
        for rxnid in self.reactions:
            
            for metabid in self.reactions[rxnid].substrates:
                sub_coe = -self.reactions[rxnid].substrates[metabid].coe
                stoyMat_net.loc[metabid, rxnid] = sub_coe
                
            for metabid in self.reactions[rxnid].products:
                pro_coe = self.reactions[rxnid].products[metabid].coe
                stoyMat_net.loc[metabid, rxnid] = pro_coe
                
        return stoyMat_net    
    
    
    @property
    def stoichiometric_matrix(self):
        if len(self._metabolites) == 0 and len(self._reactions) == 0:
            raise AttributeError(
                "can't compute stoichiometric matrix, "
                "no metabolite or reaction found, model empty"
            )
        
        stoyMat_net = self._get_stoichiometric_matrix(
            tuple(sorted(self._metabolites)), 
            tuple(self._reactions.keys())
        )
            
        return stoyMat_net
            
    
    @lru_cache()
    # Re-estimates the total stoichiometric matrix if changes occur in 
    # self._metabolites or self._reactions, otherwise uses cached results.
    def _get_total_stoichiometric_matrix(self, metabolites, reactions):
        '''
        Parameters
        ----------
        metabolites: tuple
            A tuple of metabolite IDs.
        reactions: tuple
            A tuple of reaction IDs.
        '''

        stoyMat_net = self._get_stoichiometric_matrix(metabolites, reactions)
        
        stoyMat_total = []
        total_rxnids = []
        for rxnid, col in stoyMat_net.items():
            if self.reactions[rxnid].rev:
                stoyMat_total.append(col)
                stoyMat_total.append(-col+0.0)   # avoid -0.0
                total_rxnids.append(rxnid+'_f')
                total_rxnids.append(rxnid+'_b')
            else:
                stoyMat_total.append(col)
                total_rxnids.append(rxnid)
        
        stoyMat_total = pd.DataFrame(stoyMat_total, index=total_rxnids).T
        
        return stoyMat_total    
    
    
    @property
    def total_stoichiometric_matrix(self):
        if len(self._metabolites) == 0 and len(self._reactions) == 0:
            raise AttributeError(
                "can't compute total stoichiometric matrix, "
                "no metabolite or reaction found, model empty"
            )
        
        stoyMat_total = self._get_total_stoichiometric_matrix(
            tuple(sorted(self._metabolites)), 
            tuple(self._reactions.keys())
        )
            
        return stoyMat_total    
    

    @lru_cache()
    # Re-estimates the transformation matrix if changes occur in 
    # self._metabolites or self._reactions, otherwise uses cached results.
    def _get_transformation_matrix(self, metabolites, reactions):
        '''
        Parameters
        ----------
        metabolites: tuple
            A tuple of metabolite IDs.
        reactions: tuple
            A tuple of reaction IDs.
        '''

        transMat = pd.DataFrame(
            0, 
            index=self.stoichiometric_matrix.columns,
            columns=self.total_stoichiometric_matrix.columns
        )
        for rxnid in transMat.index:
            if self.reactions[rxnid].rev:
                transMat.loc[rxnid, rxnid+'_f'] = 1
                transMat.loc[rxnid, rxnid+'_b'] = -1
            else:
                transMat.loc[rxnid, rxnid] = 1

        return transMat


    @property
    def transformation_matrix(self):
        if len(self._metabolites) == 0 and len(self._reactions) == 0:
            raise AttributeError(
                "can't transformation matrix, "
                "no metabolite or reaction found, model empty"
            )

        transMat = self._get_transformation_matrix(
            tuple(sorted(self._metabolites)), 
            tuple(self._reactions.keys())
        )

        return transMat

    
    @property
    def end_metabolites(self):
        endsDict = PrettyDict()
        for metabid, row in self.stoichiometric_matrix.iterrows():
            if row[row!=0].size == 1:
                endsDict[metabid] = self.metabolites[metabid]
                
        return endsDict
    
        
    def optimize(
            self, 
            kind, 
            *, 
            objective=None, 
            flux_bound=(0, 1000), 
            conc_bound=(0.001, 100), 
            spec_flux_bound=None, 
            spec_conc_bound=None, 
            preset_flux=None, 
            preset_conc=None,
            preset_conc_ratio=None, 
            irr_reactions=None, 
            ex_conc=None, 
            ex_mass_bal_cons=None, 
            ex_thermo_cons=None, 
            inc_enz_cons=None, 
            enz_prot_lb=None,
            parsimonious=False,
            slack=1e-3, 
            dgpm_conf_level=None
    ):
        '''
        Perform constraint-based optimization considering various constraints such 
        as mass balance, thermodynamic, and enzyme protein allocation constraints.

        Parameters
        ----------
        
        kind: {'fba', 'tfba', 'efba', 'etfba'}
            Type of optimization to perform:
            - 'fba': Flux Balance Analysis with mass balance constraints only.
            - 'tfba': Flux Balance Analysis with both mass balance and 
            thermodynamic constraints.
            - 'efba': Flux Balance Analysis with both mass balance and enzyme 
            protein cost constraints.
            - 'etfba': Flux Balance Analysis with mass balance, thermodynamics, and 
            enzyme protein cost constraints.
        objective: dict
            Mapping of flux IDs to coefficients in the objective expression, e.g., 
            {'r1': 2, 'r2': -1} defines the expression "2*r1 - 1*r2". Suffix of 
            '_f' or '_b' is required to indicate forward or backward flux for 
            reversible reactions. Valid in 'fba', 'tfba', 'efba', and 'etfba'.
        flux_bound: 2-tuple
            Lower and upper bounds of metabolic fluxes in mmol/gCDW/h. 
            Valid in 'fba', 'tfba', 'efba' and 'etfba'.
        conc_bound: 2-tuple
            Lower and upper bounds of metabolite concentrations in mM. 
            Valid in 'tfba' and 'etfba'.
        spec_flux_bound: dict
            Mapping of flux IDs to their bounds (lb, ub), where spec_flux_bound 
            takes priority over flux_bound. Suffix of '_f' or '_b' is required to 
            indicate forward or backward flux for reversible reactions. Valid in 
            'fba', 'tfba', 'efba', and 'etfba'.
        spec_conc_bound: dict
            Mapping of metabolite IDs to their bounds (lb, ub), where 
            spec_conc_bound takes priority over conc_bound. Valid in 'tfba' and 
            'etfba'.
        preset_flux: dict
            Mapping of flux IDs to fixed metabolic fluxes, e.g., substrate update 
            rates. preset_flux takes priority over flux_bound. Suffix of '_f' or 
            '_b' is required to indicate forward or backward flux for reversible 
            reactions. Valid in 'fba', 'tfba', 'efba', and 'etfba'.
        preset_conc: dict
            Mapping of metabolite IDs to fixed metabolite concentrations, e.g., 
            substrate concentrations in the media. preset_conc takes priority over 
            conc_bound. Valid in 'tfba' and 'etfba'.
        preset_conc_ratio: dict
            Mapping of ratio IDs to fixed ratios of metabolites. Ratio ID has the 
            format "metabid:metabid". Valid in 'tfba' and 'etfba'.
        irr_reactions: list
            List of irreversible reaction IDs. irr_reactions is prioritized in 
            defining reversibilities. Valid in 'fba', 'tfba', 'efba', and 'etfba'.
        ex_conc: list
            List of metabolite IDs excluded from optimization. Valid in 'tfba' and 
            'etfba'.
        ex_mass_bal_cons: list
            List of metabolite IDs excluded from mass balance constraints. Valid in 
            'fba', 'tfba', 'efba', and 'etfba'.
        ex_thermo_cons: list
            List of reaction IDs excluded from thermodynamics constraints. For 
            reversible reactions, both forward and backward reactions are excluded. 
            Valid in 'tfba' and 'etfba'.
        inc_enz_cons: list
            List of reaction IDs included in enzyme protein cost constraints. 
            Default values are used for missing kinetic parameters of Km, kcat, and 
            MW. Valid in 'efba' and 'etfba'. 
        enz_prot_lb: float
            Upper bound of enzyme protein fraction in g/gCDW. Valid in 'efba' and 
            'etfba'.
        parsimonious: bool
            Whether to further calculate parsimonious fluxes while maintaining the 
            objective no worse than that obtained by conventional flux balance 
            analysis. It may take longer to compute but could be useful for 
            eliminating loops. Valid in 'fba', 'tfba', 'efba', and 'etfba'.
        slack: float
            Small nonnegative constant used to relax the objective constraint in 
            parsimonious FBA. The objective is required to be no less than 
            (1-slack)*opt_obj. Considering adjusting slack if parsimonious 
            FBA encounters difficulties in finding feasible solutions. Valid in 
            'fba', 'tfba', 'efba' and 'etfba'.
        dgpm_conf_level: float
            If dgpm_conf_level is set, the standard reaction Gibbs energy is 
            allowed to vary between [dgpm - z*error, dgpm + z*error], where dgpm 
            and error are the standard reaction Gibbs energy and its error, and z 
            is determined by dgpm_conf_level. Valid in 'tfba' and 'etfba'.
        use_initial_fluxes (deprecated): dict or file (.xlsx, .tsv, .bin) or None
            Initial values of fluxes; optimization will be performed without 
            initials if None. Suffix of '_f' or '_b' is required to indicate 
            forward or backward flux for reversible reactions. Valid only in 
            'etfba'.
        use_initial_concs (deprecated): dict or file (.xlsx, .tsv, .bin) or None
            Initial values of concentrations; optimization will be performed 
            without initials if None. Valid only in 'etfba'.
        '''
        
        direction = 'max'

        if kind.lower() == 'fba':
            return FBAOptimizer(
                self, 
                objective, 
                direction, 
                flux_bound, 
                spec_flux_bound, 
                preset_flux, 
                irr_reactions, 
                ex_mass_bal_cons,
                parsimonious,
                slack
            )
            
        elif kind.lower() == 'tfba':
            return TFBAOptimizer(
                self, 
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
                dgpm_conf_level
            )
        
        elif kind.lower() == 'efba':
            if inc_enz_cons is None:
                raise TypeError('inc_enz_cons argument should be set for ETFBA')
            
            if enz_prot_lb is None:
                raise TypeError(
                    'enz_prot_lb argument should be set for ETFBA'
                )
            
            return EFBAOptimizer(
                self, 
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
                slack
            ) 

        elif kind.lower() == 'etfba':
            if inc_enz_cons is None:
                raise TypeError('inc_enz_cons argument should be set for ETFBA')
            
            if enz_prot_lb is None:
                raise TypeError(
                    'enz_prot_lb argument should be set for ETFBA'
                )    

            return ETFBAOptimizer(
                self, 
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
            )    
        
        else:
            raise ValueError(
                'kind should be one of {"fba", "tfba", "efba", "etfba"}'
            )
        

    def evaluate_variability(
            self, 
            kind, 
            *, 
            objective=None,
            obj_value=None,
            gamma=1,
            flux_bound=(0, 100), 
            conc_bound=(0.001, 10), 
            spec_flux_bound=None, 
            spec_conc_bound=None, 
            preset_flux=None, 
            preset_conc=None,
            preset_conc_ratio=None, 
            irr_reactions=None, 
            ex_conc=None, 
            ex_mass_bal_cons=None, 
            ex_thermo_cons=None,
            inc_enz_cons=None, 
            enz_prot_lb=None,
            dgpm_conf_level=None,
    ):
        '''
        Perform variability analysis to assess the feasible range of derived 
        fluxes, reaction Gibbs energies and enzyme protein costs.

        Parameters
        ----------
        kind: {'fva', 'tfva', 'efva', 'etfva', 'tva', 'etva', 'eva', 'teva'}
            Type of variability analysis to perform:
            - 'fva': Flux Variability Analysis.
            - 'tfva': Thermodynamically constrained Flux Variability Analysis.
            - 'efva': Enzyme protein allocation constrained Flux Variability 
            Analysis.
            - 'etfva': Enzyme protein allocation and thermodynamically constrained 
            Flux 
            Variability Analysis.
            - 'tva': Thermodynamic Variability Analysis.
            - 'etva': Enzyme protein allocation constrained Thermodynamic 
            Variability Analysis.
            - 'eva': Enzyme Protein Variability Analysis.
            - 'teva': Thermodynamically constrained Enzyme Protein Variability 
            Analysis.
        objective: dict
            Objective function used in corresponding 'fba', 'tfba', 'efba', or 
            'etfba' analysis. Mapping of flux IDs to coefficients in the objective 
            expression, e.g., {'r1': 2, 'r2': -1} defines the expression "2*r1 - 
            1*r2". Suffix of '_f' or '_b' is required to indicate forward or 
            backward flux for reversible reactions. 
            Valid in 'fva', 'tfva', 'efva', 'etfva', 'tva', 'etva', 'eva', 'teva'.
        obj_value: non-negative float
            Optimal objective of corresponding 'fba', 'tfba', 'efba', or 'etfba' 
            analysis.
        gamma: float in [0, 1]
            Objective expression required to be no less than gamma*obj_value. 
            Consider adjusting gamma if etfba encounters difficulties in finding 
            feasible solutions.
        flux_bound: 2-tuple
            Lower and upper bound of metabolic fluxes in mmol/gCDW/h. Valid in 
            'fva', 'tfva', 'efva', 'etfva', 'tva', 'etva', 'eva', 'teva'.
        conc_bound: 2-tuple
            Lower and upper bound of metabolite concentrations in mM. Valid in 
            'tfva', 'etfva', 'tva', 'etva', and 'teva'.
        spec_flux_bound: dict
            Mapping of flux IDs to their bounds (lb, ub), where spec_flux_bound 
            takes priority over flux_bound. Suffix of '_f' or '_b' is required to 
            indicate forward or backward flux for reversible reactions. Valid in 
            'fva', 'tfva', 'efva', 'etfva', 'tva', 'etva', 'eva', 'teva'.
        spec_conc_bound: dict
            Mapping of metabolite IDs to their bounds (lb, ub), where 
            spec_conc_bound takes priority over conc_bound. Valid in 'tfva', 
            'etfva', 'tva', 'etva' and 'teva'.
        preset_flux: dict
            Mapping of flux IDs to fixed metabolic fluxes, e.g., substrate update 
            rates. preset_flux takes priority over flux_bound. Suffix of '_f' or 
            '_b' is required to indicate forward or backward flux for reversible 
            reactions. Valid in 'fva', 'tfva', 'efva', 'etfva', 'tva', 'etva', 
            'eva', 'teva'.
        preset_conc: dict
            Mapping of metabolite IDs to fixed metabolite concentrations, e.g., 
            substrate concentrations in the media. preset_conc takes priority over 
            conc_bound. Valid in 'tfva', 'etfva', 'tva', 'etva', and 'teva'.
        preset_conc_ratio: dict
            Mapping of ratio IDs to fixed ratios of metabolites. Ratio ID has the 
            format "metabid:metabid". Valid in 'tfva', 'etfva', 'tva', 'etva' and 
            'teva'.
        irr_reactions: list of reaction ID
            List of irreversible reaction IDs. irr_reactions is prioritized in 
            defining reversibilities. Valid in 'fva', 'tfva', 'efva', 'etfva', 
            'tva', 'etva', 'eva', 'teva'.
        ex_conc: list of metabolite ID
            List of metabolite IDs excluded from optimization. Valid in 'tfva', 
            'etfva', 'tva', 'etva' and 'teva'.
        ex_mass_bal_cons: list of metabolite ID
            List of metabolite IDs excluded from mass balance constraints. Valid in 
            'fva', 'tfva', 'efva', 'etfva', 'tva', 'etva', 'eva', 'teva'.
        ex_thermo_cons: list of reaction ID
            List of reaction IDs excluded from thermodynamics constraints. For 
            reversible reactions, both forward and backward reactions are excluded. 
            Valid in 'tfva', 'etfva', 'tva', 'etva' and 'teva'.
        inc_enz_cons: list of reaction ID
            List of reaction IDs included in enzyme protein cost constraints. 
            Default values are used for missing kinetic parameters of Km, kcat, and 
            MW. Valid in 'etva', 'teva', 'efva' and 'etfva'. 
        enz_prot_lb: float
            Upper bound of enzyme protein fraction in g/gCDW. Valid in 'etva', 
            'teva', 'efva' and 'etfva'.
        dgpm_conf_level: float
            If dgpm_conf_level is set, the standard reaction Gibbs energy is 
            allowed to vary between [dgpm - z*error, dgpm + z*error], where dgpm 
            and error are the standard reaction Gibbs energy and its error, and z 
            is determined by dgpm_conf_level. Valid in 'tfva', 'etfva', 'tva', 
            'etva' and 'teva'.
        '''
        
        direction = 'max'

        if kind.lower() == 'tva':
            if obj_value is None:
                raise ValueError(
                    'obj_value should be provided for variability analysis, '
                    'call optimize(kind = "tfba") first to get this value'
                ) 

            return TVAOptimizer( 
                self, 
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
                dgpm_conf_level
            )

        elif kind.lower() == 'etva':
            if obj_value is None:
                raise ValueError(
                    'obj_value should be provided for variability analysis, '
                    'call optimize(kind = "etfba") first to get this value'
                )
            
            if inc_enz_cons is None:
                raise TypeError('inc_enz_cons argument should be set for ETVA')
            
            if enz_prot_lb is None:
                raise TypeError('enz_prot_lb argument should be set for ETVA')

            return ETVAOptimizer(
                self, 
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
                dgpm_conf_level
            )   

        elif kind.lower() == 'eva':
            if obj_value is None:
                raise ValueError(
                    'obj_value should be provided for variability analysis, '
                    'call optimize(kind = "efba") first to get this value'
                )
            
            if inc_enz_cons is None:
                raise TypeError('inc_enz_cons argument should be set for EVA')
            
            if enz_prot_lb is None:
                raise TypeError('enz_prot_lb argument should be set for EVA')
            
            return EVAOptimizer(
                self,
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
                enz_prot_lb
            )
        
        elif kind.lower() == 'teva':
            if obj_value is None:
                raise ValueError(
                    'obj_value should be provided for variability analysis, '
                    'call optimize(kind = "etfba") first to get this value'
                )
            
            if inc_enz_cons is None:
                raise TypeError('inc_enz_cons argument should be set for TEVA')
            
            if enz_prot_lb is None:
                raise TypeError('enz_prot_lb argument should be set for TEVA')
            
            return TEVAOptimizer(
                self,
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
                dgpm_conf_level
            )

        elif kind.lower() == 'fva':
            if obj_value is None:
                raise ValueError(
                    'obj_value should be provided for variability analysis, '
                    'call optimize(kind = "fba") first to get this value'
                )
            
            return FVAOptimizer(
                self,
                objective, 
                direction,
                obj_value,
                gamma, 
                flux_bound, 
                spec_flux_bound, 
                preset_flux, 
                irr_reactions, 
                ex_mass_bal_cons,
            )
        
        elif kind.lower() == 'tfva':
            if obj_value is None:
                raise ValueError(
                    'obj_value should be provided for variability analysis, '
                    'call optimize(kind = "tfba") first to get this value'
                )
            
            return TFVAOptimizer(
                self, 
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
                dgpm_conf_level
            )
        
        elif kind.lower() == 'efva':
            if obj_value is None:
                raise ValueError(
                    'obj_value should be provided for variability analysis, '
                    'call optimize(kind = "efba") first to get this value'
                )

            if inc_enz_cons is None:
                raise TypeError('inc_enz_cons argument should be set for EFVA')
            
            if enz_prot_lb is None:
                raise TypeError(
                    'enz_prot_lb argument should be set for EFVA'
                ) 
            
            return EFVAOptimizer(
                self,
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
                enz_prot_lb
            )
        
        elif kind.lower() == 'etfva':
            if obj_value is None:
                raise ValueError(
                    'obj_value should be provided for variability analysis, '
                    'call optimize(kind = "etfba") first to get this value'
                )

            if inc_enz_cons is None:
                raise TypeError('inc_enz_cons argument should be set for ETFVA')
            
            if enz_prot_lb is None:
                raise TypeError(
                    'enz_prot_lb argument should be set for ETFVA'
                )

            return ETFVAOptimizer(
                self, 
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
                dgpm_conf_level
            )
    

    def __repr__(self):
        if len(self._metabolites) != 0 and len(self._reactions) != 0:
            rxn_plural = 's' if len(self._reactions) > 1 else ''
            metab_plural = 's' if len(self._metabolites) > 1 else ''

            return (
                f'model {self.name if self.name else "unknown"} with '
                f'{len(self._reactions)} reaction{rxn_plural} and '
                f'{len(self._metabolites)} metabolite{metab_plural}'
            )
        else:
            return f'model {self.name if self.name else "unknown"} empty'    
    