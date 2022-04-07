#!/usr/bin/env python
# -*- coding: UTF-8 -*-


r'''
python C:\Users\cwu\Desktop\Software\TFBA-ext\etfba\test.py
'''


ETFBAPATH = 'path\to\etfba'
MODELFILE = 'path\to\E_coli_core.xlsx'


import sys
sys.path.append(ETFBAPATH)
from etfba import Model


def test_Ecoli_core_etfba():
	
	ecoli = Model('ecoli')
	ecoli.read_from_excel(MODELFILE)
	#print(ecoli.end_metabolites)
	#print(ecoli.reactions['r73'].substrates['3pg'].km)
	#ecoli.stoichiometric_matrix.to_excel(r'C:\Users\cwu\Desktop\S.xlsx')
	#ecoli.total_stoichiometric_matrix.to_excel(r'C:\Users\cwu\Desktop\tS.xlsx')
	'''
	res = ecoli.optimize('fba', objective = {'r94': 1}, preset_fluxes = {'r53': 10, 'r54': 0, 'r50': 8.39}, 
					     flux_bounds = (-100, 100),
					     irr_reactions = None,
					     excluded_mb = None).solve()
	print(res)
					 
	res = ecoli.optimize('fba', objective = {'biom': 1}, preset_fluxes = {'glc_e': 10, 'fru_e': 0, 'atpm': 8.39}, 
					     flux_bounds = (-100, 100),
					     irr_reactions = None,
					     excluded_mb = None).solve()
	print(res)
	
	res = ecoli.optimize('fba', objective = {'biom': 1}, preset_fluxes = {'glc_e': 10, 'fru_e': 0, 'o2_e': 0, 'atpm': 8.39}, 
					     flux_bounds = (-100, 100),
					     irr_reactions = None,
					     excluded_mb = None).solve()
	print(res)
	#'r50', 'r53', 'r54', 'r58'
	#'atpm', 'glc_e', 'fru_e', 'o2_e'
	'''
	'''
	res = ecoli.optimize('tfba', objective = {'r94': 1}, 
						 preset_fluxes = {'r53': 10, 'r54': 0, 'r50': 8.39}, preset_concs = {'glc.e': 5.6},   # 1 g/L glucose
					     flux_bounds = (-100, 100), conc_bounds = (0.001, 10),
					     irr_reactions = None,
					     excluded_mb = None).solve()
	print(res)
	#print(res.opt_total_fluxes)
	
	res = ecoli.optimize('tfba', objective = {'biom': 1}, 
						 preset_fluxes = {'glc_e': 10, 'fru_e': 0, 'atpm': 8.39}, preset_concs = {'glc.e': 5.6},   # 1 g/L glucose
					     flux_bounds = (-100, 100), conc_bounds = (0.001, 10),
					     irr_reactions = None,
					     excluded_mb = None).solve()
	print(res)
	#print(res.opt_directions)
	#print(res.opt_gibbs_energy)
	'''
	
	res = ecoli.optimize('etfba', objective = {'r94': 1}, 
						 preset_fluxes = {'r53': 10, 'r54': 0, 'r50': 8.39}, preset_concs = {'glc.e': 5.6},
						 flux_bounds = (-100, 100), conc_bounds = (0.001, 10),
					     excluded_mb = None,
					     included_epc = ['r01', 'r02', 'r03', 'r04', 'r05', 'r06', 'r07', 'r08', 'r09', 'r10', 'r11', 'r12', 'r13', 'r14', 'r15', 'r16', 'r17', 'r18', 'r19', 'r20', 'r21', 'r22', 'r23', 'r24', 'r25', 'r26', 'r27', 'r28', 'r29', 'r30', 'r31', 'r32', 'r33', 'r34', 'r35', 'r36', 'r37', 'r38', 'r39', 'r40', 'r41', 'r42', 'r43', 'r44', 'r45'],
						 use_fba_results = False, use_tfba_results = True).solve()
	print(res)
	print(res.opt_total_enzyme_cost)
	print(res.opt_enzyme_costs)
	'''
	res = ecoli.optimize('etfba', objective = {'biom': 1}, 
						 preset_fluxes = {'glc_e': 10, 'fru_e': 0, 'o2_e': 0, 'atpm': 8.39}, preset_concs = {'glc.e': 5.6},
						 flux_bounds = (-100, 100), conc_bounds = (0.001, 10),
					     excluded_mb = None,
					     included_epc = ['pgi', 'pfk', 'fbp', 'fba', 'tpi', 'gap', 'pgk', 'gpm', 'eno', 'pyk', 'pps', 'pdh', 'zwf', 'pgl', 'gnd', 'rpi', 'rpe', 'tkt1', 'tal', 'tkt2', 'cs', 'acn1', 'acn2', 'icd', 'kgd', 'suc', 'sdh', 'fum', 'mdh', 'gs', 'gdh', 'gls', 'gogat', 'aldh', 'adh', 'pta', 'ak', 'ldh', 'pfl', 'icl', 'mals', 'me1', 'me2', 'ppc', 'ppck'],
						  use_fba_results = False, use_tfba_results = True).solve()
	#Glycolysis/Gluconeogenesis
	#'pgi', 'pfk', 'fbp', 'fba', 'tpi', 'gap', 'pgk', 'gpm', 'eno', 'pyk', 'pps', 'pdh'
	#'r01', 'r02', 'r03', 'r04', 'r05', 'r06', 'r07', 'r08', 'r09', 'r10', 'r11', 'r12'
	#Pentose Phosphate Pathway
	#'zwf', 'pgl', 'gnd', 'rpi', 'rpe', 'tkt1', 'tal', 'tkt2'
	#'r13', 'r14', 'r15', 'r16', 'r17', 'r18', 'r19', 'r20'
	#Citric Acid Cycle
	#'cs', 'acn1', 'acn2', 'icd', 'kgd', 'suc', 'sdh', 'fum', 'mdh'
	#'r21', 'r22', 'r23', 'r24', 'r25', 'r26', 'r27', 'r28', 'r29'
	#Glutamate Metabolism
	#'gs', 'gdh', 'gls', 'gogat'
	#'r30', 'r31', 'r32', 'r33'
	#Pyruvate Metabolism
	#'aldh', 'adh', 'pta', 'ak', 'ldh', 'pfl'
	#'r34', 'r35', 'r36', 'r37', 'r38', 'r39'
	#Glyoxylate Shunt
	#'icl', 'mals'
	#'r40', 'r41'
	#Anaplerotic Reactions
	#'me1', 'me2', 'ppc', 'ppck'
	#'r42', 'r43', 'r44', 'r45'
	print(res)
	#print(res.opt_objective)
	print(res.opt_total_enzyme_cost)
	print(res.opt_enzyme_costs)
	'''
	
	
	
	
if __name__ == '__main__':

	test_Ecoli_core_etfba()
	
	