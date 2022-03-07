#!/usr/bin/env python
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'



import os
import re
import numpy as np
import pandas as pd
from equilibrator_api import ComponentContribution, Q_, Reaction


MODEL_NAME = 'ecoli'
REACTION_FILE = 'path\to\reactions.xlsx'
OUT_FILE = 'path\to\dgpms.xlsx'
CYTOPLASMIC_PH_RANGE = (6.5, 7.5)
PERIPLASMIC_PH_RANGE = (5, 6)
NPOINTS = 3

CYTOPLASMIC_PMG = 3
PERIPLASMIC_PMG = 3
CYTOPLASMIC_IS = '250 mM'
PERIPLASMIC_IS = '200 mM'
V_DIFFERENCE = '0.15 V'
R = 8.315e-3
T = 298.15


class GibbsEnergyEstimator():
	'''
	Attributes
	name: str
	metabolites: dict
		metabolite ID => Compound
	transmembrane_reactions: dict
		reaction ID => Reaction
	nontransmembrane_reactions: dict
		reaction ID => Reaction
	nreactants: dict
		reaction ID => [# of substrates, # of products], excluding H+ and H2O
	'''
	
	def __init__(self, name):
		
		self.name = name
		self.metabolites = {}
		self.transmembrane_reactions = {}
		self.nontransmembrane_reactions = {}
		self.nreactants = {}
		
	
	def read_reactions(self, filename):
		'''
		Parameters
		filename: str
			path of excel file, required fields in "metabolites" sheet: "Metabolite ID" and "KEGG ID"
			required fields in "reactions" sheet: "Reaction ID", "Substrates", "Products" and "Transmembrane"
		'''
		
		modelData = pd.read_excel(filename, sheet_name = None, header = 0, index_col = 0, comment = '#', squeeze = True)
		metabInfo = modelData['metabolites']
		rxnInfo = modelData['reactions']
		
		self.model = ComponentContribution()
		
		for metabid, keggid in metabInfo.iteritems():
			self.metabolites[metabid] = self.model.get_compound(keggid)
		
		for rxnid, (subsStr, prosStr, trans) in rxnInfo.iterrows():
			
			subStrLst = subsStr.split(';')
			proStrLst = prosStr.split(';')
			
			rxnDict = {}
			inRxnDict = {}
			outRxnDict = {}
			nsubs = 0
			npros = 0
			for subStr in subStrLst:
				
				coe_sub = subStr.split()
				if len(coe_sub) == 1:
					coe, subid = 1, coe_sub[0]
				else:
					coe, subid = coe_sub
				
				if not trans:
					rxnDict[self.metabolites[subid]] = -float(coe)
				else:
					if re.search(r'.e$', subid):
						outRxnDict[self.metabolites[re.sub(r'(.*).e$', '\g<1>', subid)]] = -float(coe)
					else:
						inRxnDict[self.metabolites[subid]] = -float(coe)
						
				if self.metabolites[re.sub(r'(.*).e$', '\g<1>', subid)].id not in [4, 5]:   # H+ and H2O
					nsubs += 1
				
			for proStr in proStrLst:
				
				coe_pro = proStr.split()
				if len(coe_pro) == 1:
					coe, proid = 1.0, coe_pro[0]
				else:
					coe, proid = coe_pro
					
				if not trans:	
					rxnDict[self.metabolites[proid]] = float(coe)	
				else:
					if re.search(r'.e$', proid):
						outRxnDict[self.metabolites[re.sub(r'(.*).e$', '\g<1>', proid)]] = float(coe)
					else:
						inRxnDict[self.metabolites[proid]] = float(coe)
						
				if self.metabolites[re.sub(r'(.*).e$', '\g<1>', proid)].id not in [4, 5]:   # H+ and H2O
					npros += 1		
			
			if rxnDict:
				self.nontransmembrane_reactions[rxnid] = Reaction(rxnDict)
			if inRxnDict and outRxnDict:
				self.transmembrane_reactions[rxnid] = [Reaction(inRxnDict), Reaction(outRxnDict)]
				
			self.nreactants[rxnid] = [nsubs, npros]	
		
	
	def check_reaction_balance(self):
		'''
		check balance of non-transmembrane reactions
		'''
		
		unbalanced_found = False
		for rxnid, rxn in self.nontransmembrane_reactions.items():
			if not rxn.is_balanced():
				unbalanced_found = True
				print('%s is not balanced' % rxnid)
		
		if not unbalanced_found:
			print('all reactions balanced')
	
	
	def estimate_deltaGprime0(self, kind, *, cytoplasmic_ph = None, cytoplasmic_pmg = None, cytoplasmic_is = None,
							  periplasmic_ph = None, periplasmic_pmg = None, periplasmic_is = None,
							  e_potential_difference = None):
		'''
		Parameters
		kind: str
			't' or 'nt', representing transmembrane and non-transmembrane, respectively
		cytoplasmic_ph: float
			cytoplasmic pH
		cytoplasmic_pmg: float
			cytoplasmic pMg	
		cytoplasmic_is: str with unit
			cytoplasmic ionic strength
		periplasmic_ph: float
			periplasmic pH, required if kind == 't'
		periplasmic_pmg: float
			periplasmic pMg, required if kind == 't'	
		periplasmic_is: str with unit
			periplasmic ionic strength, required if kind == 't'
		e_potential_difference: str with unit
			electrical potential difference, required if kind == 't'
			
		Returns
		deltaGprime0s: dict
			rxnid => deltaGprime0
		'''
		
		deltaGprime0s = {}
		
		self.model.p_h = Q_(cytoplasmic_ph)
		self.model.p_mg = Q_(cytoplasmic_pmg)
		self.model.ionic_strength = Q_(cytoplasmic_is)
		
		if kind == 'nt':
			if periplasmic_ph is not None:
				raise TypeError('periplasmic_ph not required for non-transmembrane reactions')
				
			if periplasmic_is is not None:
				raise TypeError('periplasmic_is not required for non-transmembrane reactions')
			
			if periplasmic_pmg is not None:
				raise TypeError('periplasmic_pmg not required for non-transmembrane reactions')
				
			if e_potential_difference is not None:
				raise TypeError('e_potential_difference not required for non-transmembrane reactions')
			
			for rxnid, rxn in self.nontransmembrane_reactions.items():
				deltaGprime0s[rxnid] = self.model.standard_dg_prime(rxn).value.m_as('kJ/mol')
				
		if kind == 't':
			for rxnid, rxns in self.transmembrane_reactions.items():
				deltaGprime0s[rxnid] = self.model.multicompartmental_standard_dg_prime(
											reaction_inner = rxns[0], reaction_outer = rxns[1], 
											e_potential_difference = Q_(e_potential_difference),
											p_h_outer = Q_(periplasmic_ph), p_mg_outer = Q_(periplasmic_pmg),
											ionic_strength_outer = Q_(periplasmic_is)).value.m_as('kJ/mol')
			
		#print('\n'.join(['%s: %s' % (rxnid, dgp0) for rxnid, dgp0 in deltaGprime0s.items()]))
		
		return deltaGprime0s
		
		
	def estimate_deltaGprimem(self, kind, *, cytoplasmic_ph = None, cytoplasmic_pmg = None, cytoplasmic_is = None,
							  periplasmic_ph = None, periplasmic_pmg = None, periplasmic_is = None,
							  e_potential_difference = None):
		
		deltaGprime0s = self.estimate_deltaGprime0(kind, cytoplasmic_ph = cytoplasmic_ph, cytoplasmic_pmg = cytoplasmic_pmg,
												   cytoplasmic_is = cytoplasmic_is, periplasmic_ph = periplasmic_ph, 			
												   periplasmic_pmg = periplasmic_pmg, periplasmic_is = periplasmic_is, 
												   e_potential_difference = e_potential_difference)
		
		deltaGprimems = {}
		for rxnid, dgp0 in deltaGprime0s.items():
			nsubs, npros = self.nreactants[rxnid]
			deltaGprimems[rxnid] = dgp0 + 3*np.log(10)*R*T*(nsubs - npros)
		
		#print('\n'.join(['%s: %s' % (rxnid, dgpm) for rxnid, dgpm in deltaGprimems.items()]))
		
		return deltaGprimems




if __name__ == '__main__':
	
	os.makedirs(os.path.dirname(OUT_FILE), exist_ok = True)
	
	ecoli_gibbs = GibbsEnergyEstimator(MODEL_NAME)
	ecoli_gibbs.read_reactions(REACTION_FILE)
	#ecoli_gibbs.check_reaction_balance()
	
	pH_dgpms = []
	for i, phi in enumerate(np.linspace(*CYTOPLASMIC_PH_RANGE, NPOINTS)):
		dgpms_nt = ecoli_gibbs.estimate_deltaGprimem('nt', cytoplasmic_ph = phi, cytoplasmic_pmg = CYTOPLASMIC_PMG, 
													 cytoplasmic_is = CYTOPLASMIC_IS)
		
		for j, pho in enumerate(np.linspace(*PERIPLASMIC_PH_RANGE, NPOINTS)):
			dgpms_t = ecoli_gibbs.estimate_deltaGprimem('t', cytoplasmic_ph = phi, cytoplasmic_pmg = CYTOPLASMIC_PMG, 
													    cytoplasmic_is = CYTOPLASMIC_IS, periplasmic_ph = pho, 
													    periplasmic_pmg = PERIPLASMIC_PMG, periplasmic_is = PERIPLASMIC_IS, 
													    e_potential_difference = V_DIFFERENCE)
														
			pH_dgpms.append([phi, pho] + list(dgpms_nt.values()) + list(dgpms_t.values()))	  
	
	pH_dgpms = pd.DataFrame(pH_dgpms, columns = ['ph_in', 'ph_out'] + list(dgpms_nt.keys()) + list(dgpms_t.keys()))
	pH_dgpms.to_excel(OUT_FILE, header = True, index = True)
	
	
		
		
		
	
	
