#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'




from .pdict import SmartDict




class Reaction():
	'''
	Attributes
	substrates: PrettyDict
		metabolite ID => Metabolite 
	
	products: PrettyDict
		metabolite ID => Metabolite
		
	fkcat, bkcat: float
		kcat value in forward and backward direction
		
	mw: float
		molecular weight of the catalytic enzyme
		
	dgpm: float
		stand reaction Gibbs energy change
		
	rev: bool
		reaction revisibility
		
	is_biomass_formation: bool
		whether the biomass formation reaction
		
	is_constrained_by_thermodynamics: bool
		whether in the thermodynamics constraints
	'''
	
	def __init__(self, rxnid, enzyme_name = None, *, forward_kcat = None, backward_kcat = None, 
				 molecular_weight = None, standard_deltaG = None, reversible = True):
		'''
		Parameters
		rxnid: str
			reaction ID
		
		enzyme_name: str
			enzyme name
		
		forward_kcat: float
			kcat value (1/s) in forward direction
		
		backward kcat: float
			kcat value (1/s) in backward direction
		
		molecular weight: float
			enzyme molecular weight in kDa
		
		standard_deltaG: float
			standard reaction Gibbs energy change in kJ/mol with metabolite concentraions in mM
		
		reversible: bool
			reversibility
		'''
	
		self.rxnid = rxnid
		self.enzyme = enzyme_name
		self.fkcat = forward_kcat
		self.bkcat = backward_kcat
		self.mw = molecular_weight
		self.dgpm = standard_deltaG
		self.rev = reversible
		
		self._substrates = SmartDict()
		self._products = SmartDict()
		
		self.is_biomass_formation = False
		self.is_exch_reaction = False
		
		self.is_constrained_by_thermodynamics = False
		
	
	@property
	def substrates(self):
		
		self._substrates.caller = self
		
		return self._substrates
		
		
	@property
	def products(self):
		
		self._products.caller = self
		
		return self._products
	
		
	def __repr__(self):
		
		if self.substrates and self.products:
			subsStr = ' + '.join(sorted(self.substrates.keys()))
			prosStr = ' + '.join(sorted(self.products.keys()))
			
			return '%s => %s' % (subsStr, prosStr)
		
		else:
			return 'reaction not constructed'
			
			
			
			
			
			
			
		
		
		
		
		
	
	
	
