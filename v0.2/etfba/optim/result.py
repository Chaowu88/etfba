#!/usr/bin/env python
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'




import re
import numpy as np
from ..core.pdict import PrettyDict




class FBAResults():
	'''
	Attributes
	opt_objective: float
		optimal objective
		
	opt_fluxes: dict
		reaction ID => optimal flux value
	'''
	
	def __init__(self, opt_obj, opt_fluxes):
		'''
		Parameters
		opt_obj: float
			optimal objective
			
		opt_fluxes: dict
			reaction ID => optimal flux value
		'''
		
		self._opt_obj = opt_obj
		self._opt_fluxes = opt_fluxes
	
	
	@property
	def opt_objective(self):
		
		return round(self._opt_obj, 2)
		
	
	@property
	def opt_fluxes(self):
		
		return PrettyDict({rxnid: round(flux, 2) for rxnid, flux in self._opt_fluxes.items()})
	
	
	def __repr__(self):
		
		optObjStr = 'optimal objective: %s' % self.opt_objective
		optFluxesStr = str(self.opt_fluxes)
		
		return optObjStr + '\n\noptimal fluxes\n' + optFluxesStr + '\n'
		



class TFBAResults(FBAResults):
	'''
	Attributes
	...
		
	opt_concentrations: dict
		metabolite ID => optimal metabolite concentration
		
	opt_directions: dict
		reaction ID => "f"/"b"
		
	opt_gibbs_energy: dict
		reaction ID => optimal deltaGprime
	'''
	
	def __init__(self, opt_obj, opt_fluxes, opt_lnconcs, opt_dgps):
		'''
		Parameters
		opt_obj: float
			optimal objective
			
		opt_fluxes: dict
			reaction ID => optimal flux value
			
		opt_lnconcs: dict
			metabolite ID => optimal concentration
			
		opt_dgps: dict
			reaction ID => optimal reaction Gibbs energy change
		'''
		
		super().__init__(opt_obj, opt_fluxes)
		
		self._opt_lnconcs = opt_lnconcs
		self._opt_concs = dict(zip(self._opt_lnconcs.keys(), np.exp(list(self._opt_lnconcs.values()))))
		self._opt_dgps = opt_dgps
		
		
	@property
	def opt_concentrations(self):
		
		return PrettyDict({metabid: round(conc, 2) for metabid, conc in self._opt_concs.items()})
	
		
	@property
	def opt_directions(self):
	
		return PrettyDict({rxnid: 'f' if flux >= 0 else 'b' for rxnid, flux in self.opt_fluxes.items()})
	
	
	@property
	def opt_gibbs_energy(self):
		
		return PrettyDict({rxnid: round(dgp, 2) for rxnid, dgp in self._opt_dgps.items()})
		
	
	def __repr__(self):
		
		optObjStr = 'optimal objective: %s' % self.opt_objective
		optFluxesStr = str(self.opt_fluxes)
		optConcsStr = str(self.opt_concentrations)
		
		return optObjStr + '\n\noptimal fluxes\n' + optFluxesStr + '\n\noptimal concentrations\n' + optConcsStr + '\n'
		
		
		
		
class ETFBAResults(TFBAResults):
	'''
	Attributes
	...
		
	opt_total_enzyme_cost: float
		optimal total enzyme protein abundance
			
	opt_enzyme_costs: dict
		reaction IDs => optimal enzyme protein abundance	
	'''
	
	def __init__(self, opt_obj, opt_total_fluxes, opt_lnconcs, opt_dgps, opt_epc, opt_epcs):
		'''
		Parameters
		opt_obj: float
			optimal objective
			
		opt_fluxes: dict
			forward/backward reaction ID => optimal flux value
			
		opt_lnconcs: dict
			metabolite ID => optimal concentration
			
		opt_dgps: dict
			reaction ID => optimal reaction Gibbs energy change
			
		opt_epc: float
			optimal total enzyme protein abundance
			
		opt_epcs: dict
			reaction IDs => optimal enzyme protein abundance
		'''
		
		super().__init__(opt_obj, opt_total_fluxes, opt_lnconcs, opt_dgps)
		
		self._opt_epc = opt_epc
		self._opt_epcs = opt_epcs
		
		
	@property
	def opt_total_enzyme_cost(self):
		
		return round(self._opt_epc, 5)
	
	
	@property
	def opt_enzyme_costs(self):
		
		return PrettyDict({rxnid: round(epc, 5) for rxnid, epc in self._opt_epcs.items()})
	
			
			
	

