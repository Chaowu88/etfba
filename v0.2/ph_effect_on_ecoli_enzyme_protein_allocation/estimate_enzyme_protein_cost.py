#!/usr/bin/env python
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'


r'''
python C:\Users\cwu\Desktop\Software\Papers\pH_effect\plot_ph_effect_contour\estimate_enzyme_protein_cost.py
'''


import os
import pandas as pd
import sys
sys.path.append('path\to\etfba')
from etfba import Model


MODEL_NAME = 'ecoli'
MODEL_FILE = 'path\to\E_coli_core.xlsx'
GIBBSENERGY_FILE = 'path\to\dgpms.xlsx'
OUT_DIR = 'path\to\output'
OBJECTIVE = {'biom': 1}
PRESET_FLUXES = {'glc_e': 10, 'fru_e': 0, 'o2_e': 0, 'atpm': 8.39}
PRESET_CONCS = {'glc.e': 5.6}
INCLUDED_EPC = ['pgi', 'pfk', 'fbp', 'fba', 'tpi', 'gap', 'pgk', 'gpm', 'eno', 'pyk', 'pps', 'pdh', 'zwf', 'pgl', 'gnd', 
				'rpi', 'rpe', 'tkt1', 'tal', 'tkt2', 'cs', 'acn1', 'acn2', 'icd', 'kgd', 'suc', 'sdh', 'fum', 'mdh', 'gs', 'gdh', 'gls', 'gogat', 'aldh', 'adh', 'pta', 'ak', 'ldh', 'pfl', 'icl', 'mals', 'me1', 'me2', 'ppc', 'ppck']
USE_FBA_RESULTS = False
USE_TFBA_RESULTS = True
TRIAL_LIMIT = 20	
	
	
def main():
	
	model = Model(MODEL_NAME)
	model.read_from_excel(MODEL_FILE)
	
	dgpmsInfo = pd.read_excel(GIBBSENERGY_FILE, header = 0, index_col = 0)
	
	epcs = []
	fluxes = []
	concs = []
	dgps = []
	for i, row in dgpmsInfo.iterrows():
		
		phs, dgpms = row[:2], row[2:]
		
		for rxnid, dgpm in dgpms.items():
			model.reactions[rxnid].dgpm = dgpm
		
		count = 0
		while count < TRIAL_LIMIT:
			try:
				print('run %s' % i)
				
				res = model.optimize('etfba', objective = OBJECTIVE, 
									 preset_fluxes = PRESET_FLUXES, preset_concs = PRESET_CONCS,
									 flux_bounds = (-100, 100), conc_bounds = (0.001, 10),
									 excluded_mb = None, included_epc = INCLUDED_EPC,
									 use_fba_results = USE_FBA_RESULTS, use_tfba_results = USE_TFBA_RESULTS).solve()
				optcosts = res.opt_enzyme_costs
				optfluxes = res.opt_fluxes
				optconcs = res.opt_concentrations
				optdgps = res.opt_gibbs_energy
				break
			
			except:
				count += 1
				continue
			
		if count == TRIAL_LIMIT:	
			optcosts = {}
			optfluxes = {}
			optconcs = {}
			optdgps = {}
			
		epcs.append({**phs.to_dict(), **optcosts})
		fluxes.append({**phs.to_dict(), **optfluxes})
		concs.append({**phs.to_dict(), **optconcs})
		dgps.append({**phs.to_dict(), **optdgps})
	
	epcs = pd.DataFrame(epcs)
	fluxes = pd.DataFrame(fluxes)
	concs = pd.DataFrame(concs)
	dgps = pd.DataFrame(dgps)
	
	os.makedirs(OUT_DIR, exist_ok = True)
	epcs.to_excel(r'%s/epcs.xlsx' % OUT_DIR, header = True, index = True)
	fluxes.to_excel(r'%s/fluxes.xlsx' % OUT_DIR, header = True, index = True)
	concs.to_excel(r'%s/concs.xlsx' % OUT_DIR, header = True, index = True)
	dgps.to_excel(r'%s/dgps.xlsx' % OUT_DIR, header = True, index = True)



if __name__ == '__main__':
	
	main()
	
	
	
	
	