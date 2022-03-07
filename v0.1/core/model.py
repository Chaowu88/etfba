#!/usr/bin/env python
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'



import re
from functools import lru_cache
import numpy as np
import pandas as pd
from .reaction import Reaction
from .metabolite import Metabolite
from ..optim.optim import FBAOptimizer, TFBAOptimizer, ETFBAOptimizer
from .pdict import PrettyDict




class Model():
	'''
	Attributes
	metabolites: PrettyDict
		metabolite ID => Metabolite 
	
	reactions: PrettyDict
		reaction ID => Reaction
	
	end_metabolites: PrettyDict
		metabolite ID => Metabolite, initial substrates or final products of the model
	
	stoichiometric_matrix： DataFrame
		stoichiometric matrix, rows are metabolites, columns are forward reactions
		
	total_stoichiometric_matrix: DataFrame
		stoichiometric matrix, rows are metabolites, columns are total (forward and backward) reactions
	'''
	
	def __init__(self, name = None):
		'''
		Parameters
		name: str
			model name
		'''
		
		self.name = name
		
		self._metabolites = PrettyDict()
		self._reactions = PrettyDict()
		self._end_metabolites = PrettyDict()
		self._stoichiometric_matrix = pd.DataFrame()
		self._total_stoichiometric_matrix = pd.DataFrame()
		
	
	@staticmethod
	def _float(string):
		
		return np.nan if string == '' else float(string)
		
		
	def read_from_excel(self, filename):
		'''
		Parameters
		filename: str
			path of excel file with fields: 
			Enzyme, Substrates, Products, Sub Kms (mM), Pro Kms (mM), 
			Fwd kcat (1/s), Bwd kcat (1/s), MW (kDa) and ΔrG'm (kJ/mol)
		'''
		
		modelData = pd.read_excel(filename, header = 0, index_col = 0, comment = '#').fillna('').astype(str)
		
		for rxnid, rowInfos in modelData.iterrows():
			subsStr, prosStr, rev, subkms, prokms, fkcat, bkcat, mw, dgpm = rowInfos
			
			fkcat = self._float(fkcat)
			bkcat = self._float(bkcat)
			rev = float(rev)
			mw = self._float(mw)
			dgpm = self._float(dgpm)
			rxn = Reaction(rxnid, forward_kcat = fkcat, backward_kcat = bkcat, reversible = bool(rev),
						   molecular_weight = mw, standard_deltaG = dgpm)
			
			if re.search(r'biomass', prosStr, flags = re.I):
				rxn.is_biomass_formation = True
				
			if re.match(r'^[\w\._]+\.o$', subsStr) or re.match(r'^[\w\._]+\.o$', prosStr):
				rxn.is_exchange = True
			
			subStrLst = subsStr.split(';')
			subkmLst = subkms.split(';')
			if len(subStrLst) != len(subkmLst) and not rxn.is_biomass_formation:
				raise ValueError('the number of subtrates in %s does not match the number of subtrate Km values' % rxnid)
			
			proStrLst = prosStr.split(';')
			prokmLst = prokms.split(';')
			if len(proStrLst) != len(prokmLst) and not rxn.is_biomass_formation:
				raise ValueError('the number of products in %s does not match the number of product Km values' % rxnid)
			
			for idx, subStr in enumerate(subStrLst):
				
				coe_sub = subStr.split()
				if len(coe_sub) == 1:
					coe, subid = 1.0, coe_sub[0]
				else:
					coe, subid = coe_sub
				
				sub = self._metabolites.setdefault(subid, Metabolite(subid))
				sub.coes[rxnid] = float(coe)
				
				if rxn.is_biomass_formation or rxn.is_exchange:
					sub.kms[rxnid] = None
				else:
					sub.kms[rxnid] = self._float(subkmLst[idx])
				
				rxn.substrates[subid] = sub	
			
			for idx, proStr in enumerate(proStrLst):
			
				coe_pro = proStr.split()
				if len(coe_pro) == 1:
					coe, proid = 1.0, coe_pro[0]
				else:
					coe, proid = coe_pro
				
				pro = self._metabolites.setdefault(proid, Metabolite(proid))
				pro.coes[rxnid] = float(coe)
				
				if rxn.is_biomass_formation or rxn.is_exchange:
					pro.kms[rxnid] = None
				else:
					pro.kms[rxnid] = self._float(prokmLst[idx])
				
				rxn.products[proid] = pro
			
			self._reactions[rxnid] = rxn
	
	
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
	def _get_stoichiometric_matrix(self):
		
		stoyMat = pd.DataFrame(0, index = sorted(self._metabolites), columns = sorted(self._reactions))
		for rxnid in self.reactions:
			
			for metabid in self.reactions[rxnid].substrates:
				stoyMat.loc[metabid, rxnid] = -self.reactions[rxnid].substrates[metabid].coe
				
			for metabid in self.reactions[rxnid].products:
				stoyMat.loc[metabid, rxnid] = self.reactions[rxnid].products[metabid].coe
				
		return stoyMat	
	
	
	@property
	def stoichiometric_matrix(self):
		
		if len(self._metabolites) == 0 and len(self._reactions) == 0:
			raise AttributeError('no metabolite or reaction found, model empty')
		
		if self._stoichiometric_matrix.empty:
			self._stoichiometric_matrix = self._get_stoichiometric_matrix()
			
		return self._stoichiometric_matrix
			
	
	@lru_cache()
	def _get_total_stoichiometric_matrix(self):
		
		cols = ['']*len(self._reactions)*2
		cols[::2] = [r+'_f' for r in sorted(self._reactions)]
		cols[1::2] = [r+'_b' for r in sorted(self._reactions)]
		
		stoyExtMat = pd.DataFrame(0, index = sorted(self._metabolites), columns = cols)
		for rxnid in self.reactions:
			
			for metabid in self.reactions[rxnid].substrates:
				coe = self.reactions[rxnid].substrates[metabid].coe
				stoyExtMat.loc[metabid, rxnid+'_f'] = -coe
				stoyExtMat.loc[metabid, rxnid+'_b'] = coe
				
			for metabid in self.reactions[rxnid].products:
				coe = self.reactions[rxnid].products[metabid].coe
				stoyExtMat.loc[metabid, rxnid+'_f'] = coe
				stoyExtMat.loc[metabid, rxnid+'_b'] = -coe
				
		return stoyExtMat	
	
	
	@property
	def total_stoichiometric_matrix(self):
		
		if len(self._metabolites) == 0 and len(self._reactions) == 0:
			raise AttributeError('no metabolite or reaction found, model empty')
		
		if self._total_stoichiometric_matrix.empty:
			self._total_stoichiometric_matrix = self._get_total_stoichiometric_matrix()
			
		return self._total_stoichiometric_matrix
	
	
	@lru_cache()
	def _get_end_metabolites(self):
		
		endsDict = PrettyDict()
		for metabid, row in self.stoichiometric_matrix.iterrows():
			if row[row!=0].size == 1:
				endsDict[metabid] = self.metabolites[metabid]
				
		return endsDict
		
		
	@property
	def end_metabolites(self):
		
		if len(self._end_metabolites) == 0:   # there should be at least one end metabolite in the model
			self._end_metabolites = self._get_end_metabolites()
			
		return self._end_metabolites
	
		
	def optimize(self, kind, *, objective, direction = 'max', flux_bounds = (-100, 100), conc_bounds = None, 
				 preset_fluxes = None, preset_concs = None, irr_reactions = None, excluded_concs = None, 
				 excluded_mb = None, excluded_thmd = None, included_epc = None, use_fba_results = None):
		'''
		Parameters
		kind: str
			kind of optimization to perform:
			'fba', flux balance analysis with mass balance constraints only;
			'tfba', flux balance analysis with both mass balance and thermodynamics constraints;
			'etfba', flux balance analysis with mass balance, thermodynamics and enzyme protein cost constraints
		
		objective: dict,
			reaction ID => coefficient in the objective expression, e.g., {'r1': 2, 'r2': -1} defines the expression
			"2*r1 - 1*r2", objective forms the complete objective expression in 'fba' and 'tfba', and the denominator
			of objective expression in 'etfba'
		
		direction: str
			direction of optimization:
			'max', maximize;
			'min', minimize,
			tunable in 'fba' and 'tfba', the argument is ignored in 'etfba' because the objective is 
			enzyme_protein_cost/objective_flux or enzyme_protein_cost which will always be optimized in minimizing direction
		
		flux_bounds: tuple
			lower and upper bounds of metabolic flux in mmol/gCDW/h, valid in 'fba', 'tfba' and 'etfba'
			
		conc_bounds: tuple
			lower and upper bounds of metabolite concentration in mM, valid in 'tfba' and 'etfba'
			
		preset_fluxes: dict
			rxnid => float, fixed metabolic fluxes, valid in 'fba', 'tfba' and 'etfba'
		
		preset_concs: dict
			metabid => float, fixed metabolite concentrations, e.g., substrate concentrations in media, 
			valid in 'tfba' and 'etfba'
		
		irr_reactions: list of reaction ID
			irreversible reactions, backward flux will be set to 0. reversibilities set by irr_reactions will 
			overwrite those set in file. valid in 'fba', 'tfba' and 'etfba'
		
		excluded_concs: list of metabolite ID
			metabolite concentrations excluded from optimization, valid in 'tfba' and 'etfba'
		
		excluded_mb: list of metabolite ID
			metabolites excluded from mass balance constraints, valid in 'fba', 'tfba' and 'etfba'
			initial substrates and final products are excludeded in any case
			
		excluded_thmd: list of reaction	ID
			reactions excluded from thermodynamics constraints, valid in 'tfba' and 'etfba'
			reactions without available standard deltaG' are excluded in any case
			
		included_epc: list of reaction ID
			reactions included in enzyme protein cost constraints, valid only in 'etfba'
			if kinetic parameters of Km, kcat and M.W. are not available for included reactions, default values will be used
			
		use_fba_results: bool
			if True, enzyme_protein_cost will be minimized with precomputed fluxes by FBA;
			if False, enzyme_protein_cost/objective_flux will be minimized with variable fluxes, valid only in 'etfba'
		'''
		
		if kind.lower() == 'fba':
			if conc_bounds is not None:
				raise TypeError('FBA model does not accept conc_bounds argument')
				
			if preset_concs is not None:
				raise TypeError('FBA model does not accept preset_concs argument')
				
			if excluded_concs is not None:
				raise TypeError('FBA model does not accept excluded_concs argument')
			
			if excluded_thmd is not None:
				raise TypeError('FBA model does not accept excluded_thmd argument')
			
			if included_epc is not None:
				raise TypeError('FBA model does not accept included_epc argument')
				
			if use_fba_results is not None:
				raise TypeError('FBA model does not accept use_fba_results argument')
			
			return FBAOptimizer(self, objective, direction, flux_bounds, preset_fluxes, irr_reactions, excluded_mb)
			
		if kind.lower() == 'tfba':
			if conc_bounds is None:
				conc_bounds = (0.001, 10)
				
			if included_epc is not None:
				raise TypeError('TFBA model does not accept included_epc argument')
				
			if use_fba_results is not None:
				raise TypeError('FBA model does not accept use_fba_results argument')	
			
			return TFBAOptimizer(self, objective, direction, flux_bounds, conc_bounds, preset_fluxes, preset_concs,
								 irr_reactions, excluded_concs, excluded_mb, excluded_thmd)
			
		if kind.lower() == 'etfba':
			if conc_bounds is None:
				conc_bounds = (0.001, 10)
				
			if included_epc is None:
				raise TypeError('included_epc argument should be set for ETFBA')
				
			if use_fba_results is None:
				use_fba_results = False
				
			return ETFBAOptimizer(self, objective, direction, flux_bounds, conc_bounds, preset_fluxes, preset_concs,
								  irr_reactions, excluded_concs, excluded_mb, excluded_thmd, included_epc, use_fba_results)
			
		
	def __repr__(self):
		
		if len(self._metabolites) != 0 and len(self._reactions) != 0:
			return 'model %s with %s reactions and %s metabolites' % (self.name if self.name else 'unknown',
																	  len(self._reactions), len(self._metabolites))
		
		else:
			return 'model %s not constructed' % self.name if self.name else 'unknown'
		
		

		
		
		
		
		
		
		
		
		




