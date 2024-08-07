'''Difine the Metabolite class.'''


from ..io.results import PrettyDict


class Singleton(type):
    '''
    Make a singleton metaclass. 
    '''
    
    def __init__(cls, name, bases, attrs):
        cls._instances = {}
        
    def __call__(cls, *args, **kwargs):
        key = (sum(hash(argv) for argv in args) 
               + sum(hash(argv) for argv in kwargs.values()))
        instance = cls._instances.setdefault(
            key, super().__call__(*args, **kwargs)
        )
        
        return instance
    

# class Singleton(type):
#     '''
#     make a singleton metaclass
#     '''
    
#     _instances = {}
#     _init = {}
    
#     def __init__(cls, name, bases, attrs):
        
#         cls._init[cls] = attrs.get('__init__', None)
    
    
#     def __call__(cls, *args, **kwargs):
        
#         init = cls._init[cls]
#         if init is not None:
#             key = sum(hash(argv) for argv in args) + sum(hash(argv) for argv in kwargs.values())
#         else:
#             key = cls
        
#         if key not in cls._instances:
#             cls._instances[key] = super().__call__(*args, **kwargs)
        
#         return cls._instances[key]


class Metabolite(metaclass=Singleton):
    '''
    Metabolite instances created with identical IDs, names, and all other attributes 
    are treated as a singular instance.

    Attributes
    ----------
    kms: PrettyDict
        Key is the reaction ID where the metabolite is involved, value is the 
        corresponding Km value.
    km: positive float
        Km value of the metabolite in the host reaction.
    coes: PrettyDict
        Key is the reaction ID where the metabolite is involved, value is the 
        corresponding stoichiometric coefficient in the reaction. Negative for 
        substrates and positive for products.
    coe: positive float
        Stoichiometric coefficient of the metabolite in the host reaction.
    is_h: bool
        Indicates whether it is a proton.
    is_h2o: bool
        Indicates whether it is water.
    is_constrained_by_mass_balance: bool
        Indicates whether it appears in the mass balance constraints.
    host: Reaction
        The host reaction of the metabolite.
    role: {'substrate', 'product'}
        Role played by the metabolite in the host reaction.    
    '''
        
    def __init__(
            self, 
            metabid, 
            name=None, 
            compartment=None, 
            *, 
            is_h=False, 
            is_h2o=False
    ):
        '''
        Parameters
        ----------
        metabid: str
            Metabolite ID.
        name: str
            Metabolite name.
        compartment: str
            Compartment identifier: "c" for cytoplasm; "p" for 
            periplasm; "e" for extracellular.
        is_h: bool
            Indicates whether the metabolite is a proton.
        is_h2o: bool
            Indicates whether the metabolite is water.    
        '''
    
        self.metabid = metabid
        self.name = name
        self.compartment = compartment

        self.is_h = is_h
        self.is_h2o = is_h2o

        self.kms = PrettyDict()
        self.coes = PrettyDict()

        self.is_constrained_by_mass_balance = False
        
        self.host = None
        self.role = None
        
    
    @property
    def km(self):
        if self.host is not None:
            host = self.host
            self.host = None
            
            return self.kms[host.rxnid]
        else:
            raise AttributeError('host reaction not found, use kms instead')


    @km.setter
    def km(self, value):
        if self.host is not None:
            if value <= 0:
                raise ValueError('Km value should be positive')
            else:
                host = self.host
                self.host = None
                self.kms[host.rxnid] = value
        else:
            raise AttributeError('host reaction not found, can not set km')

    
    @property
    def coe(self):
        if self.host is not None:
            host = self.host
            self.host = None
            
            return abs(self.coes[host.rxnid])
        else:
            raise AttributeError('host reaction not found, use coes instead')


    @coe.setter
    def coe(self, value):
        if self.host is not None:
            if value <= 0:
                raise ValueError('coefficient should be postive')
            else:
                host = self.host
                self.host = None
                if self.role == 'substrate':
                    self.coes[host.rxnid] = -value
                elif self.role == 'product':
                    self.coes[host.rxnid] = value
        else:
            raise AttributeError('host reaction not found, can not set coe')
    

    # @property
    # def is_exch_metabolite(self):
        
    #     return bool(re.match(r'^[\w\._]+\.o$', self.metabid))


    # @property
    # def is_biomass(self):

    #     return bool(re.match(r'^biomass$', self.metabid, flags = re.I))    
        
    
    def __repr__(self):
        return self.name if self.name else self.metabid
        