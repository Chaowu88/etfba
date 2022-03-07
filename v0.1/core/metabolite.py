#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'



from .pdict import PrettyDict




class Singleton(type):
	'''
	make a singleton metaclass
	'''
	
	_instances = {}
	_init = {}
	
	def __init__(cls, name, bases, dct):
		
		cls._init[cls] = dct.get('__init__', None)
	
	
	def __call__(cls, *args, **kwargs):
		
		init = cls._init[cls]
		if init is not None:
			key = str(args) + str(kwargs)
		else:
			key = cls
		
		if key not in cls._instances:
			cls._instances[key] = super().__call__(*args, **kwargs)
		
		return cls._instances[key]




class Metabolite(metaclass = Singleton):
	'''
	Attributes
	kms: PrettyDict
		Km values in reactions, reaction ID => km value
		
	km: positive float
		km value of metabolite in host reaction
		
	coes: PrettyDict
		stoichiometric coefficients in reactions, reaction ID => coefficient
		
	coe: positive float
		stoichiometric coefficient of metabolite in host reaction
		
	is_constrained_by_mass_balance: bool
		if appears in the mass balance constraints
	'''
		
	def __init__(self, id, name = None):
	
		self.id = id
		self.name = name
		
		self.kms = PrettyDict()
		self.coes = PrettyDict()
		self.is_constrained_by_mass_balance = True
		
		self.host = None
		
	
	@property
	def km(self):
		
		if self.host != None:
			host = self.host
			self.host = None
			
			return self.kms[host.id]
		
		else:
			raise AttributeError('host reaction not found, use kms instead')
		
	
	@property
	def coe(self):
		
		if self.host != None:
			host = self.host
			self.host = None
			
			return self.coes[host.id]
		
		else:
			raise AttributeError('host reaction not found, use coes instead')
	
	
	def __repr__(self):
	
		return self.name if self.name else self.id
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	