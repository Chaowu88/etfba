'''
Estimate the sensitivity of growth rate to kcat downregulation and upregulation.

Usage:
python /path/to/etfba-main/scripts/estimate_growth_sensitivity.py
'''


import os
import numpy as np
import pandas as pd
from multiprocess import Pool
from etfba import Model
from helper import (get_internal_ph, add_H_leak_rxn, set_kcat_and_MW, 
                    get_default_KaC_and_KdC, set_deltaGprimem, get_H_leak_flux, 
                    set_adjusted_kcat, reset_kcat, set_cpu_affinity)
from config import (DIR, MODEL_FILE, MW_FILE, GIBBS_ENERGY_FILE, KCAT_FILE, 
                    CORR_KAPP, KD_FILE, SEP_RXNS, OBJECTIVE, PRESET_FLUXES, 
                    FLUX_BOUNDS, SPEC_FLUX_BOUNDS, PRESET_CONCS, CONC_BOUNDS, 
                    SPEC_CONC_BOUNDS, Q, GROWTH_TYPE, PERTURB_DIRECTION, 
                    PERTURB_FOLD, LOW_PH, NEUTRAL_PH, N_JOBS)


def _change_kcat(model, enzyme, fold):
    model.reactions[enzyme].forward_kcat *= fold
    if model.reactions[enzyme].rev:
        model.reactions[enzyme].backward_kcat *= fold


def perturb_kcat(model, enzyme):
    if PERTURB_DIRECTION == 'downregulation':
        _change_kcat(model, enzyme, 1/PERTURB_FOLD)
    elif PERTURB_DIRECTION == 'upregulation':
        _change_kcat(model, enzyme, PERTURB_FOLD)


def restore_kcat(model, enzyme):
    if PERTURB_DIRECTION == 'downregulation':
        _change_kcat(model, enzyme, PERTURB_FOLD)
    elif PERTURB_DIRECTION == 'upregulation':
        _change_kcat(model, enzyme, 1/PERTURB_FOLD)


def etfba_worker(
        out_dir, 
        model, 
        pert_enzymes, 
        dgpms_info, 
        kd_data, 
        KaC_default, 
        KdC_default, 
        inc_enz_cons, 
        ex_thermo_cons
    ):
    set_cpu_affinity()

    for enzyme in pert_enzymes:
        enz_out_dir = f'{out_dir}/{enzyme}'
        os.makedirs(enz_out_dir, exist_ok=True)

        perturb_kcat(model, enzyme)

        fluxes = []
        dgps = []
        epcs = []
        sel_dgpms_info = dgpms_info[[LOW_PH, NEUTRAL_PH]]
        for ph_out, dgpms in sel_dgpms_info.items():
            print(f'{enzyme} {PERTURB_DIRECTION} at external pH {ph_out}')

            set_deltaGprimem(model, dgpms)
            ph_in = get_internal_ph(ph_out)
            preset_flux = PRESET_FLUXES.copy()
            preset_flux.update(get_H_leak_flux(ph_in, ph_out))
            preset_conc = PRESET_CONCS.copy()
            preset_conc.update({'h_c': 10**(-ph_in+3)})   
            
            old_kcats = set_adjusted_kcat(
                model, kd_data, inc_enz_cons, SEP_RXNS, ph_in, 
                KaC_default, KdC_default
            )

            res = model.optimize(
                'etfba', 
                objective=OBJECTIVE, 
                flux_bound=FLUX_BOUNDS, 
                conc_bound=CONC_BOUNDS,
                spec_flux_bound=SPEC_FLUX_BOUNDS,
                spec_conc_bound=SPEC_CONC_BOUNDS,
                preset_flux=preset_flux,
                preset_conc=preset_conc,   
                ex_thermo_cons=ex_thermo_cons,
                inc_enz_cons=inc_enz_cons,
                enz_prot_lb=Q,
                parsimonious=True,
            ).solve(solver='gurobi')

            if res.optimization_successful:
                fluxes.append(pd.Series(res.opt_fluxes))
                dgps.append(pd.Series(res.opt_gibbs_energy))
                epcs.append(pd.Series(res.opt_enzyme_costs))

            reset_kcat(model, old_kcats, inc_enz_cons)

        restore_kcat(model, enzyme)

        pd.DataFrame(
            fluxes, index=sel_dgpms_info.columns
        ).T.to_excel(f'{enz_out_dir}/metabolic_fluxes.xlsx')
        pd.DataFrame(
            dgps, index=sel_dgpms_info.columns
        ).T.to_excel(f'{enz_out_dir}/gibbs_energies.xlsx')
        pd.DataFrame(
            epcs, index=sel_dgpms_info.columns
        ).T.to_excel(f'{enz_out_dir}/enzyme_protein_costs.xlsx')
    

def etfba_sensitivity(out_dir):
    model = Model.load(MODEL_FILE)
    add_H_leak_rxn(model)

    dgpms_info = pd.read_excel(GIBBS_ENERGY_FILE, header=0, index_col=0)
    ex_thermo_cons = [rxnid for rxnid in model.reactions 
                      if rxnid not in dgpms_info.index]
    kcat_data = pd.read_excel(KCAT_FILE, header=None, index_col=0).squeeze()
    for rxnid in CORR_KAPP:
        kcat_data.loc[rxnid] = CORR_KAPP[rxnid]
    mw_data = pd.read_excel(MW_FILE, header=None, index_col=0).squeeze()
    inc_enz_cons = set_kcat_and_MW(model, kcat_data, mw_data)
    kd_data = pd.read_excel(KD_FILE, header=0, index_col=0)
    KaC_default, KdC_default = get_default_KaC_and_KdC(kd_data)
    
    pool = Pool(processes=N_JOBS)
    for pert_enzymes in np.array_split(inc_enz_cons, N_JOBS):
        res = pool.apply_async(
            func=etfba_worker,
            args=(
                out_dir, 
                model, 
                pert_enzymes, 
                dgpms_info, 
                kd_data, 
                KaC_default, 
                KdC_default, 
                inc_enz_cons,
                ex_thermo_cons
            )
        )
    
    pool.close()
    pool.join()




if __name__ == '__main__':
    etfba_sensitivity(f'{DIR}/sensitivity/{GROWTH_TYPE}/{PERTURB_DIRECTION}')