#!/usr/bin/env python
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'



R = 8.315e-3		# gas constant in kJ/mol/K
T = 298.15			# absolute temperature in K, or 25 C
K = 10000			# big enough constant
defaultMW = 40		# default enzyme molecular weight in kDa
defaultKcat = 200	# default reaction catalytic rate constant in 1/s
defaultKm = 0.2		# default reactant Michaelis constant in mM
maxIter = 100000




import re
import numpy as np
from pyomo.environ import (ConcreteModel, Set, Param, Var, Objective, Constraint, SolverFactory,
						   NonNegativeReals, Binary, value, maximize, minimize, log, exp)
from .result import FBAResults, TFBAResults, ETFBAResults




class FBAOptimizer():
	
	def __init__(self, model, objective, direction, flux_bounds, preset_fluxes, irr_reactions, excluded_mb):
		'''
		Parameters
		model: Model
			model that calls FBAOptimizer
			
		objective: dict
			reaction ID => coefficient in the objective expression
		
		direction: str
			direction of optimization
			
		flux_bounds: tuple
			lower and upper bounds of metabolic flux
		
		preset_fluxes: dict
			rxnid => float, fixed metabolic fluxes
		
		irr_reactions: list
			irreversible reactions
			
		excluded_mb: list
			metabolites excluded from mass balance constraints
		'''
		
		self.model = model
		self.objective = objective
		self.direction = direction
		self.flux_bounds = flux_bounds
		self.irr_reactions = irr_reactions
		if self.irr_reactions is not None:
			for rxnid, rxn in self.model.reactions.items():
				if rxnid in self.irr_reactions:
					rxn.rev = False
				else:
					rxn.rev = True
		else:
			self.irr_reactions = [rxnid for rxnid, rxn in self.model.reactions.items() if rxn.rev == False]
		self.preset_fluxes = preset_fluxes
		self.excluded_mb = [] if excluded_mb is None else excluded_mb
		
		self.pyoModel = ConcreteModel()
		
		self.pyoModel.fluxIDs = Set(initialize = self.model.stoichiometric_matrix.columns)
		self.pyoModel.metabIDs = Set(initialize = self.model.stoichiometric_matrix.index)
		
	
	def _build_flux_variables(self):
		
		if self.preset_fluxes is None:
			
			def flux_bounds_rule(model, rxnid):
				if rxnid in self.irr_reactions:
					return (max(0, self.flux_bounds[0]), self.flux_bounds[1])
				else:
					return self.flux_bounds
			
			self.pyoModel.fluxes = Var(self.pyoModel.fluxIDs, bounds = flux_bounds_rule)
		
		else:
			def flux_bounds_rule(model, rxnid):
				if rxnid in self.preset_fluxes:   # set the bounds of fixed fluxes
					return (self.preset_fluxes[rxnid],)*2
				elif rxnid in self.irr_reactions:
					return (max(0, self.flux_bounds[0]), self.flux_bounds[1])
				else:
					return self.flux_bounds
				
			self.pyoModel.fluxes = Var(self.pyoModel.fluxIDs, bounds = flux_bounds_rule)
	
	
	def _build_objective(self):
		
		for k in self.objective:
			if k not in self.model.reactions:
				raise KeyError('use valid reaction IDs in objective')
		
		if self.direction.lower() == 'max':
			direction = maximize
		elif self.direction.lower() == 'min':
			direction = minimize
		else:
			raise ValueError("only 'max' or 'min' is acceptable")
		
		def obj_rule(model):
			return sum(coe*model.fluxes[rxnid] for rxnid, coe in self.objective.items())
			
		self.pyoModel.obj = Objective(rule = obj_rule, sense = direction)	
		
	
	def _build_mass_balance_contraints(self):
		
		def mb_rule(model, metabid):
			stoyMat = self.model.stoichiometric_matrix
			
			return sum(stoyMat.loc[metabid, rxnid]*model.fluxes[rxnid] for rxnid in stoyMat.columns) == 0
			
		self.pyoModel.MBcstrs = Constraint(self.pyoModel.metabIDs, rule = mb_rule)
		
		# end metabolites and metabolites in excluded_mb are not balanced
		for metabid in self.pyoModel.metabIDs:
			if metabid in self.model.end_metabolites or metabid in self.excluded_mb:
				
				self.model.metabolites[metabid].is_constrained_by_mass_balance = False
				self.pyoModel.MBcstrs[metabid].deactivate()   
		
		
	def solve(self):	
		
		self._build_flux_variables()
		self._build_objective()
		self._build_mass_balance_contraints()
		
		solver = SolverFactory('glpk')
		solver.solve(self.pyoModel)
		
		optObj = value(self.pyoModel.obj)
		
		optFluxes = {}
		for rxnid in self.pyoModel.fluxes:
			optFluxes[rxnid] = value(self.pyoModel.fluxes[rxnid])
		
		return FBAResults(optObj, optFluxes)
	
	


class TFBAOptimizer(FBAOptimizer):
	
	def __init__(self, model, objective, direction, flux_bounds, conc_bounds, preset_fluxes, preset_concs, 
				 irr_reactions, excluded_concs, excluded_mb, excluded_thmd):
		'''
		Parameters
		model: Model
			model that calls FBAOptimizer
			
		objective: dict
			reaction ID => coefficient in the objective expression
		
		direction: str
			direction of optimization
			
		flux_bounds: tuple
			lower and upper bounds of metabolic flux
		
		conc_bounds: tuple
			lower and upper bounds of metabolite concentration
			
		preset_fluxes: dict
			rxnid => float, fixed metabolic fluxes
		
		preset_concs: dict
			metabid => float, fixed metabolite concentration
		
		irr_reactions: list
			irreversible reactions
			
		excluded_concs: list
			metabolite concentrations excluded from optimization
			
		excluded_mb: list
			metabolites excluded from mass balance constraints
			
		excluded_thmd: list
			metabolites excluded from thermodynamics constraints
		'''
		
		super().__init__(model, objective, direction, flux_bounds, preset_fluxes, irr_reactions, excluded_mb)
		
		self.conc_bounds = conc_bounds
		self.lnconc_bounds = tuple(np.log(self.conc_bounds))
		self.preset_concs = preset_concs
		self.excluded_thmd = [] if excluded_thmd is None else excluded_thmd
		self.excluded_concs = [] if excluded_concs is None else excluded_concs
		
		self.pyoModel.inMetabIDs = Set(initialize = [metabid for metabid in self.model.metabolites
															 if metabid.lower() != 'biomass' and 
															 not self._is_exch_metab(metabid) and
															 not self._is_excluded_conc(metabid)])   # exclude biomass, exchange metabolites, and metabolites excluded from optimization
		self.pyoModel.influxIDs = Set(initialize = [rxnid for rxnid, rxn in self.model.reactions.items()
														  if not rxn.is_biomass_formation and 
														  not rxn.is_exchange])   # exclude biomass formation and exchange reactions
		
	
	def _is_exch_metab(self, metabid):
		
		return True if re.match(r'^[\w\._]+\.o$', metabid) else False
	
	
	def _is_excluded_conc(self, metabid):
		
		return re.sub(r'(.+)\..+$', '\g<1>', metabid).lower() in self.excluded_concs
		
	
	def _build_conc_variables(self):
		
		if self.preset_concs is None:
			
			self.pyoModel.lnconcs = Var(self.pyoModel.inMetabIDs, bounds = self.lnconc_bounds)
		
		else:
			def conc_bounds_rule(model, metabid):
				if metabid in self.preset_concs:   # set the bounds of fixed concentrations
					return (np.log(self.preset_concs[metabid]),)*2
				else:
					return self.lnconc_bounds
				
			self.pyoModel.lnconcs = Var(self.pyoModel.inMetabIDs, bounds = conc_bounds_rule)
		
	
	def _calculate_gibbs_energy(self, model, rxnid):
		'''
		Parameters
		model: pyomo model
			pyomo model
		
		rxnid: str
			reaction ID
		'''
		
		subs = self.model.reactions[rxnid].substrates
		pros = self.model.reactions[rxnid].products
		
		subsSum = sum([subs[subid].coe*model.lnconcs[subid] for subid in subs if not self._is_excluded_conc(subid)])
		prosSum = sum([pros[proid].coe*model.lnconcs[proid] for proid in pros if not self._is_excluded_conc(proid)])
		
		return self.model.reactions[rxnid].dgpm + (prosSum - subsSum)*R*T

		
	def _build_thermodynamics_constraints(self):
		
		def thmd_rule(model, rxnid):
			return model.fluxes[rxnid]*self._calculate_gibbs_energy(model, rxnid) <= 0
		
		self.pyoModel.THMDcstr = Constraint(self.pyoModel.influxIDs, rule = thmd_rule)
		
		# reactions in excluded_thmd are excluded from thermodynamics constraints
		for rxnid in self.pyoModel.influxIDs:
			if rxnid in self.excluded_thmd:
				
				self.model.reactions[rxnid].is_constrained_by_thermodynamics = False
				self.pyoModel.THMDcstr[rxnid].deactivate()
		
		
	def solve(self):	
		
		self._build_flux_variables()
		self._build_conc_variables()
		self._build_objective()
		self._build_mass_balance_contraints()
		self._build_thermodynamics_constraints()
		
		solver = SolverFactory('ipopt')
		solver.options['max_iter']= maxIter
		solver.solve(self.pyoModel)
		
		optObj = value(self.pyoModel.obj)
		
		optFluxes = {}
		for rxnid in self.pyoModel.fluxes:
			optFluxes[rxnid] = value(self.pyoModel.fluxes[rxnid])
		
		optLnconcs = {}
		for metabid in self.pyoModel.lnconcs:
			optLnconcs[metabid] = value(self.pyoModel.lnconcs[metabid])
		
		optDgps = {}
		for rxnid in self.pyoModel.THMDcstr:
			optDgps[rxnid] = value(self._calculate_gibbs_energy(self.pyoModel, rxnid))
		
		return TFBAResults(optObj, optFluxes, optLnconcs, optDgps)
		



class ETFBAOptimizer(TFBAOptimizer):
	
	def __init__(self, model, objective, direction, flux_bounds, conc_bounds, preset_fluxes, preset_concs, 
				 irr_reactions, excluded_concs, excluded_mb, excluded_thmd, included_epc, use_fba_results):
		'''
		Parameters
		model: Model
			model that calls FBAOptimizer
			
		objective: dict
			reaction ID => coefficient in the objective expression
		
		direction: str
			direction of optimization
			
		flux_bounds: tuple
			lower and upper bounds of metabolic flux
		
		conc_bounds: tuple
			lower and upper bounds of metabolite concentration
			
		preset_fluxes: dict
			rxnid => float, fixed metabolic fluxes
			
		preset_concs: dict
			metabid => float, fixed metabolite concentration
			
		irr_reactions: list
			irreversible reactions
		
		excluded_concs: list
			metabolite concentrations excluded from optimization
		
		excluded_mb: list
			metabolites excluded from mass balance constraints
			
		excluded_thmd: list
			reactions excluded from thermodynamics constraints
			
		included_epc: list
			reactions excluded from enzyme protein cost constraints
			
		use_fba_results: bool
			whether to use precomputed fluxes by FBA
		'''
		
		super().__init__(model, objective, direction, flux_bounds, conc_bounds, preset_fluxes, preset_concs, 
						 irr_reactions, excluded_concs, excluded_mb, excluded_thmd)
		
		self.included_epc = included_epc
		self.use_fba_results = use_fba_results
		
		
	@staticmethod
	def _get_real(value, default):
		
		return value if value is not np.nan else default
		
	
	def _calculate_enzyme_cost(self, model, rxnids):
		'''
		Parameters
		model: pyomo model
			pyomo model
		
		rxnids: list
			reaction IDs
		'''
		
		costs = []
		for rxnid in rxnids:
			
			v = model.fluxes[rxnid]
			fkcat = self._get_real(self.model.reactions[rxnid].fkcat, defaultKcat)
			mw = self._get_real(self.model.reactions[rxnid].mw, defaultMW)
			dgpm = self.model.reactions[rxnid].dgpm	
				
			subs = self.model.reactions[rxnid].substrates
			pros = self.model.reactions[rxnid].products
			
			subsKmSum = sum([subs[subid].coe*log(self._get_real(self.model.metabolites[subid].kms[rxnid], defaultKm)) 
							 for subid in subs if not self._is_excluded_conc(subid)])
			prosKmSum = sum([pros[proid].coe*log(self._get_real(self.model.metabolites[proid].kms[rxnid], defaultKm)) 
							 for proid in pros if not self._is_excluded_conc(proid)])
			
			subsConcSum = sum([subs[subid].coe*model.lnconcs[subid] for subid in subs if not self._is_excluded_conc(subid)])
			prosConcSum = sum([pros[proid].coe*model.lnconcs[proid] for proid in pros if not self._is_excluded_conc(proid)])
		
			#e = v/fkcat*(1 + exp(subsKmSum - subsConcSum))/(1 - exp(prosConcSum - subsConcSum + dgpm/R/T))
			e = v/fkcat*(exp(subsConcSum - subsKmSum) + exp(prosConcSum - prosKmSum) + 1)/(exp(subsConcSum - subsKmSum)*(1 - exp(prosConcSum - subsConcSum + dgpm/R/T)))
			
			costs.append(1/3600*mw*e)
			
		return costs
	
	
	def _build_flux_parameters(self):
		
		fbaFluxes = FBAOptimizer(self.model, self.objective, self.direction, self.flux_bounds, self.preset_fluxes, 
								   self.irr_reactions, self.excluded_mb).solve().opt_fluxes
		
		self.pyoModel.fluxes = Param(self.pyoModel.fluxIDs, initialize = fbaFluxes)
		
	
	def _build_objective_without_fluxes(self):
		
		for rxnid in self.included_epc:
			if self.model.reactions[rxnid].is_biomass_formation:
				raise ValueError("biomass formation can't be included in enzyme protein cost")
				
			if self.model.reactions[rxnid].is_exchange:
				raise ValueError("exchange reaction %s can't be included in enzyme protein cost" % rxnid)
				
		def obj_rule(model):
			return sum(self._calculate_enzyme_cost(model, self.included_epc))
			
		self.pyoModel.obj = Objective(rule = obj_rule, sense = minimize)   # always minimize the objective
		
	
	def _build_objective_with_fluxes(self):
		
		for k in self.objective:
			if k not in self.model.reactions:
				raise KeyError('use valid reaction IDs in objective')
		
		for rxnid in self.included_epc:
			if self.model.reactions[rxnid].is_biomass_formation:
				raise ValueError("biomass formation can't be included in enzyme protein cost")
				
			if self.model.reactions[rxnid].is_exchange:
				raise ValueError("exchange reaction %s can't be included in enzyme protein cost" % rxnid)
		
		def obj_rule(model):
			flux_obj = sum(coe*model.fluxes[rxnid] for rxnid, coe in self.objective.items())
			epc_obj = sum(self._calculate_enzyme_cost(model, self.included_epc))
			
			return epc_obj/flux_obj
			
		self.pyoModel.obj = Objective(rule = obj_rule, sense = minimize)   # always minimize the objective
	
	
	def _build_thermodynamics_constraints_without_fluxes(self):
		
		def thmd_rule(model, rxnid):
			if model.fluxes[rxnid] >= 0:
				return self._calculate_gibbs_energy(model, rxnid) <= 0
			else:
				return self._calculate_gibbs_energy(model, rxnid) >= 0
		
		self.pyoModel.THMDcstr = Constraint(self.pyoModel.influxIDs, rule = thmd_rule)
		
		# reactions in excluded_thmd are excluded from thermodynamics constraints
		for rxnid in self.pyoModel.influxIDs:
			if rxnid in self.excluded_thmd:
				
				self.model.reactions[rxnid].is_constrained_by_thermodynamics = False
				self.pyoModel.THMDcstr[rxnid].deactivate()
				
	
	def solve(self):
		
		if self.use_fba_results:
			self._build_flux_parameters()
			self._build_conc_variables()
			self._build_objective_without_fluxes()
			self._build_thermodynamics_constraints_without_fluxes()
		
		else:
			self._build_flux_variables()
			self._build_conc_variables()
			self._build_objective_with_fluxes()
			self._build_mass_balance_contraints()
			self._build_thermodynamics_constraints()
		
		solver = SolverFactory('ipopt')
		solver.options['max_iter']= maxIter
		res = solver.solve(self.pyoModel)
		
		optObj = value(sum(coe*self.pyoModel.fluxes[rxnid] for rxnid, coe in self.objective.items()))
		
		optFluxes = {}
		for rxnid in self.pyoModel.fluxes:
			optFluxes[rxnid] = value(self.pyoModel.fluxes[rxnid])
		
		optLnconcs = {}
		for metabid in self.pyoModel.lnconcs:
			optLnconcs[metabid] = value(self.pyoModel.lnconcs[metabid])
			
		optDgps = {}
		for rxnid in self.pyoModel.THMDcstr:
			optDgps[rxnid] = value(self._calculate_gibbs_energy(self.pyoModel, rxnid))
		
		ecosts = self._calculate_enzyme_cost(self.pyoModel, self.included_epc)
		
		optTotalEcost = value(sum(ecosts))
		
		optEcosts = {}
		for rxnid, cost in zip(self.included_epc, ecosts):
			optEcosts[rxnid] = value(cost)
		
		return ETFBAResults(optObj, optFluxes, optLnconcs, optDgps, optTotalEcost, optEcosts)	
		

