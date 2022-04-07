#!/usr/bin/env python
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'


r'''
python C:\Users\cwu\Desktop\Software\Papers\pH_effect\plot_ph_effect_contour\plot_contour.py
'''


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


EPCS_FILE = 'path\to\epcs.xlsx'
DGPS_FILE = 'path\to\dgps.xlsx'
OUT_DIR = 'path\to\output'
ENZYME_CAT = {'glycolysis': ['pgi', 'pfk', 'fbp', 'fba', 'tpi', 'gap', 'pgk', 'gpm', 'eno', 'pyk', 'pps', 'pdh'],
			  'pentose phosphate pathway': ['zwf', 'pgl', 'gnd', 'rpi', 'rpe', 'tkt1', 'tal', 'tkt2'],
			  'TCA cycle': ['cs', 'acn1', 'acn2', 'icd', 'kgd', 'suc', 'sdh', 'fum', 'mdh', 'icl', 'mals'],
			  'glutamate metabolism': ['gs', 'gdh', 'gls', 'gogat'],
			  'pyruvate metabolism': ['aldh', 'adh', 'pta', 'ak', 'ldh', 'pfl'],
			  'anaplerotic reactions': ['me1', 'me2', 'ppc', 'ppck'],
			  'ATP metabolism': ['atps4r', 'atpm']}
PATHWAYS = ['glycolysis', 'pentose phosphate pathway', 'TCA cycle', 'glutamate metabolism', 'pyruvate metabolism', 'anaplerotic reactions']			  
NPOINTS = 20


def plot_contour(out_dir, filename, data_file, cmap, constant_color):
	
	dataInfo = pd.read_excel(data_file, header = 0, index_col = 0)
	
	X_phin = dataInfo['ph_in'].values.reshape(NPOINTS, NPOINTS).T
	Y_phout = dataInfo['ph_out'].values.reshape(NPOINTS, NPOINTS).T
	
	data = dataInfo.iloc[:, 2:].copy()
	data = (data.fillna(method = 'bfill') + data.fillna(method = 'ffill'))/2   # impute by mean
	if filename == 'epcs':
		for i, row in data.iterrows():
			if (row < 0).any() or (row >= 1).any() or row.sum() >= 1:
				data.loc[i, :] = 0
	#data.to_excel(r'C:\Users\cwu\Desktop\all.xlsx')#!!!
	if filename == 'epcs':
		ndigits = 5
	elif filename == 'dgps':
		ndigits = 2
	
	# plot per enzyme
	for pathway in PATHWAYS:
		
		enzymes = ENZYME_CAT[pathway]
		
		ncols = 3
		nrows = int(np.ceil((len(enzymes)+1)/ncols))

		fig, axes = plt.subplots(nrows = nrows, ncols = ncols, figsize = (12, nrows*3), sharex = 'all', sharey = 'all')
		
		for i, enz in enumerate(enzymes+['sum']):
			
			if enz == 'sum':
				Z = data[enzymes].sum(axis = 1).values.reshape(NPOINTS, NPOINTS).T
			else:
				Z = data[enz].values.reshape(NPOINTS, NPOINTS).T
				
			if axes.ndim == 2:
				indexer = (i//ncols, i%ncols)
			elif axes.ndim == 1:
				indexer = i
			
			vmin = Z.min().min()
			
			if vmin == 0:
				vmin = 0.00001
				
			vmax = Z.max().max()
			levels = np.linspace(vmin, vmax, NPOINTS)
			
			if vmax - vmin > 0.0001:
				ctf = axes[indexer].contourf(X_phin, Y_phout, Z, vmin = vmin, vmax = vmax, levels = levels, 
											 cmap = plt.cm.get_cmap(cmap).reversed())
				
				cbar = fig.colorbar(mappable = ctf, ax = axes[indexer])
				cbarTicks = cbar.get_ticks()
				cbarTicksNew = np.linspace(cbarTicks.min(), cbarTicks.max(), 4)
				cbar.set_ticks(cbarTicksNew)
				cbar.ax.set_yticklabels(cbarTicksNew.round(ndigits))
				cbar.ax.tick_params(labelsize = 13)
			
			else:
				Z = np.full_like(Z, (vmax + vmin)/2)
				
				ctf = axes[indexer].contourf(X_phin, Y_phout, Z, NPOINTS, colors = constant_color)
				
				cbar = fig.colorbar(mappable = ctf, ax = axes[indexer])
				cbar.set_ticks([])
				cbar.ax.set_yticklabels([])
				cbar.set_label(round((vmax + vmin)/2, ndigits), horizontalalignment = 'left', rotation = 360, 
							   labelpad = 5, fontsize = 13)
			
			axes[indexer].locator_params(axis = 'x', nbins = 3)
			axes[indexer].locator_params(axis = 'y', nbins = 4)
			axes[indexer].tick_params(labelsize = 15)
			axes[indexer].set_xlabel(enz, fontsize = 25)
		
		ax_label = fig.add_subplot(111, frameon = False)
		ax_label.tick_params(labelcolor = 'none', top = False, bottom = False, left = False, right = False)
		ax_label.set_xlabel('Cytoplasmic pH', labelpad = 50, fontsize = 35)
		ax_label.set_ylabel('Periplasmic pH', labelpad = 30, fontsize = 35)
		
		for i in range(len(enzymes)+1, ncols*nrows):
			if axes.ndim == 2:
				indexer = (i//ncols, i%ncols)
			elif axes.ndim == 1:
				indexer = i
			
			fig.delaxes(ax = axes[indexer])
		
		os.makedirs(out_dir, exist_ok = True)
		#plt.tight_layout()
		fig.subplots_adjust(wspace = 0.4, hspace = 0.3)
		plt.savefig('%s/%s_%s.jpg' % (out_dir, pathway, filename), dpi = 300, bbox_inches = 'tight')
			
	
def main():
	
	plot_contour(OUT_DIR, 'epcs', EPCS_FILE, 'viridis', '#3C528B')
	plot_contour(OUT_DIR, 'dgps', DGPS_FILE, 'plasma', '#D8556C')
	
	
	
	
if __name__ == '__main__':
	
	main()




