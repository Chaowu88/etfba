'''Difine the Reaction class.'''


from collections.abc import Iterable
from ..io.results import PrettyDict


DEFAULT_MW = 40       # Default enzyme molecular weight in kDa
DEFAULT_KCAT = 200    # Default reaction catalytic rate constant in 1/s
DEFAULT_KM = 0.2      # Default reactant Michaelis constant in mM
DEFAULT_DGPM = 0      # Default reaction standard Gibbs energy in KJ/mol


class ReactantDict(PrettyDict):
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        self.caller = None   # caller is Reaction
        self.role = None   # {'substrate', 'product'}
        
    
    def __getitem__(self, key):
        item = super().__getitem__(key)   # item is Metabolite

        if self.caller is not None:
            item.host = self.caller
        
        if self.role is not None:
            item.role = self.role
            
        return item   
        
        
class Reaction():
    '''
    Attributes
    ----------
    substrates: PrettyDict
        Key is the metabolite ID, value is the Metabolite instances as substrates.
    
    products: PrettyDict
        Key is the metabolite ID, value is the Metabolite instances as products.
    
    forward_kcat (fkcat): float
        Catalytic rate constant in the forward direction.
    
    backward_kcat (bkcat): float
        Catalytic rate constant in the backward direction .   
    
    molecular_weight (mw): float
        Molecular weight of the catalytic enzyme.
    
    standard_gibbs_energy (dgpm): float
        Standard Gibbs energy of the reaction.
    
    standard_gibbs_energy_error (dgpm_error): float
        Standard error of the standard reaction Gibbs energy.    
    
    reversible (rev): bool
        Reaction reversibility.
    
    is_biomass_formation: bool
        Indicates if it's the biomass formation reaction.
    
    is_exch_reaction: bool
        Indicates if it's an exchange reaction.
    
    is_h_transport: bool
        Indicates if it's a proton transport reaction.
    
    is_h2o_transport: bool
        Indicates if it's a water transport reaction.
    
    is_constrained_by_thermodynamics: bool
        Indicates if it's constrained by thermodynamics.
    '''
    
    def __init__(
            self, 
            rxnid, 
            enzyme_name = None, 
            category = None, 
            *, 
            forward_kcat = None, 
            backward_kcat = None, 
            molecular_weight = None, 
            standard_gibbs_energy = None, 
            standard_gibbs_energy_error = None, 
            reversible = True,
            is_biomass_formation = False, 
            is_exch_reaction = False, 
            is_h_transport = False, 
            is_h2o_transport = False
    ):
        '''
        Parameters
        ----------
        rxnid: str
            Reaction ID.
        enzyme_name: str
            Enzyme name.
        category: str
            Reaction category.    
        forward_kcat: float
            Kcat value (1/s) in the forward direction.
        backward_kcat: float
            Kcat value (1/s) in the backward direction.
        molecular_weight: float
            Enzyme molecular weight in kDa.
        standard_gibbs_energy: float
            Standard reaction Gibbs energy in kJ/mol with metabolite 
            concentrations in mM.
        standard_gibbs_energy_error: float
            Standard error of the standard Gibbs energy.
        reversible: bool
            Reaction reversibility.
        is_biomass_formation: bool
            Indicates if it's the biomass formation reaction.
        is_exch_reaction: bool
            Indicates if it's an exchange reaction.
        is_h_transport: bool
            Indicates if it's a proton transport reaction.
        is_h2o_transport: bool
            Indicates if it's a water transport reaction.    
        '''
    
        self.rxnid = rxnid
        self.enzyme = enzyme_name
        self.category = category

        self.fkcat = forward_kcat
        self.bkcat = backward_kcat
        self.mw = molecular_weight
        self.dgpm = standard_gibbs_energy
        self.dgpm_error = standard_gibbs_energy_error
        self.rev = reversible
        
        self.is_biomass_formation = is_biomass_formation
        self.is_exch_reaction = is_exch_reaction
        self.is_h_transport = is_h_transport
        self.is_h2o_transport = is_h2o_transport

        self._substrates = ReactantDict()
        self._products = ReactantDict()
        
        self.is_constrained_by_thermodynamics = False


    def _add_reactants(self, coes, kms, label):
        if label == 'substrate':
            for reac, coe in coes.items():
                reac.coes[self.rxnid] = -coe
                self.substrates[reac.metabid] = reac
        elif label == 'products':
            for reac, coe in coes.items():
                reac.coes[self.rxnid] = coe
                self.products[reac.metabid] = reac
        
        if not self.is_biomass_formation and not self.is_exch_reaction:
            if kms is None:
                for reac in coes.keys():
                    reac.kms[self.rxnid] = DEFAULT_KM
            else:
                if set(kms).issubset(coes):    
                    for reac in coes.keys():
                        if reac in kms:
                            reac.kms[self.rxnid] = kms[reac]
                        else:
                            reac.kms[self.rxnid] = DEFAULT_KM
                else:
                    raise ValueError('keys in kms should be subset of those in coes')
        else:   
            # reactants in biomass formation and exchange reactions don't have Km 
            # values
            for reac in coes.keys():
                reac.kms[self.rxnid] = None
                

    def add_substrates(self, coes, kms = None):
        '''
        Parameters
        ----------
        coes: dict
            Dictionary mapping substrates to their corresponding 
            stoichiometric coefficients (positive values).
            All substrates should have a provided stoichiometric 
            coefficient.
        kms: dict
            Dictionary mapping substrates to their corresponding Km 
            values.
            Km values should only be provided as a subset of 
            substrates present in the 'coes' dictionary.
            An empty dict is acceptable if no Km values are available.
        '''

        self._add_reactants(coes, kms, 'substrate')


    def add_products(self, coes, kms = None):
        '''
        Parameters
        ----------
        coes: dict
            Dictionary mapping products to their corresponding 
            stoichiometric coefficients (positive values).
            All products should have a provided stoichiometric 
            coefficient.
        kms: dict
            Dictionary mapping products to their corresponding Km 
            values.
            Km values should only be provided as a subset of 
            products present in the 'coes' dictionary.
            An empty dict is acceptable if no Km values are available.
        '''

        self._add_reactants(coes, kms, 'products')    


    def _remove_reactants(self, reactants, label):
        if not isinstance(reactants, Iterable):
            reactants = [reactants]
        
        for reac in reactants:
            del reac.coes[self.rxnid]
            del reac.kms[self.rxnid]

        if label == 'substrate':
            for reac in reactants:
                del self.substrates[reac.metabid]
        elif label == 'product':
            for reac in reactants:
                del self.products[reac.metabid]


    def remove_substrates(self, substrates):
        '''
        Parameters
        ----------
        substrates: list of Metabolite
            List of substrates to be removed from the current 
            reaction.
        '''

        self._remove_reactants(substrates, 'substrate')


    def remove_products(self, products):
        '''
        Parameters
        ----------
        products: list of Metabolites
            List of products to be removed from current reaction.
        '''

        self._remove_reactants(products, 'product')


    @property
    def substrates(self):
        self._substrates.caller = self
        self._substrates.role = 'substrate'
        
        return self._substrates
        
        
    @property
    def products(self):
        self._products.caller = self
        self._products.role = 'product'
        
        return self._products
    

    @property
    def forward_kcat(self):
        return self.fkcat
    

    @forward_kcat.setter
    def forward_kcat(self, value):
        self.fkcat = value


    @property
    def backward_kcat(self):
        return self.bkcat
    

    @backward_kcat.setter
    def backward_kcat(self, value):
        self.bkcat = value    


    @property
    def molecular_weight(self):
        return self.mw
    

    @molecular_weight.setter
    def molecular_weight(self, value):
        self.mw = value


    @property
    def standard_gibbs_energy(self):
        return self.dgpm
    

    @standard_gibbs_energy.setter
    def standard_gibbs_energy(self, value):
        self.dgpm = value
        

    @property
    def standard_gibbs_energy_error(self):
        return self.dgpm_error
    

    @standard_gibbs_energy_error.setter
    def standard_gibbs_energy_error(self, value):
        self.dgpm_error = value
    

    @property
    def reversible(self):
        return self.rev
    

    @reversible.setter
    def reversible(self, value):
        self.rev = value
    
        
    def __repr__(self):
        if self.substrates or self.products:
            coe_subs = []
            for subid, sub in sorted(self.substrates.items()):
                coe = -sub.coes[self.rxnid]
                coe_subs.append(f'{coe} {subid}')
            subsStr = ' + '.join(coe_subs)

            coe_pros = []
            for proid, pro in sorted(self.products.items()):
                coe = pro.coes[self.rxnid]
                coe_pros.append(f'{coe} {proid}')
            prosStr = ' + '.join(sorted(coe_pros))
            
            if self.rev:
                return f'{subsStr} <=> {prosStr}'
            else:
                return f'{subsStr} => {prosStr}'
        else:
            return 'reaction not constructed'