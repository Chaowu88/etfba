'''
Estimate feasible range of reaction Gibbs energies at low and neutral pHs.

Usage:
python /path/to/etfba-main/scripts/estimate_gibbs_energy_ranges.py
'''


import os
import pandas as pd
from etfba import Model
from helper import (get_internal_ph, add_H_leak_rxn, set_kcat_and_MW, 
                    get_default_KaC_and_KdC, set_deltaGprimem, get_H_leak_flux, 
                    set_adjusted_kcat, reset_kcat)
from config import (DIR, MODEL_FILE, MW_FILE, GIBBS_ENERGY_FILE, KCAT_FILE, 
                    CORR_KAPP, KD_FILE, SEP_RXNS, ETFBA_FLUX_FILE, OBJECTIVE, 
                    GROWTH_TYPE, PRESET_FLUXES, FLUX_BOUNDS, SPEC_FLUX_BOUNDS, 
                    PRESET_CONCS, CONC_BOUNDS, SPEC_CONC_BOUNDS, Q, LOW_PH, 
                    NEUTRAL_PH, N_JOBS)

    
def etva(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    
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
    
    etfba_bioms = pd.read_excel(
        ETFBA_FLUX_FILE, header=0, index_col=0
    ).loc['BIOMASS_Ec_iML1515_core_75p37M',:]

    dgp_ranges = []
    for ph_out, dgpms in dgpms_info[[LOW_PH, NEUTRAL_PH]].items():
        print('\nEstimate Gibbs energy variability at external pH', ph_out)

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

        res = model.evaluate_variability(
            'etva',
            objective=OBJECTIVE,
            obj_value=etfba_bioms[ph_out],
            gamma=0.99,   # release gamma if a feasible solution cannot be found
            flux_bound=FLUX_BOUNDS, 
            conc_bound=CONC_BOUNDS,
            spec_flux_bound=SPEC_FLUX_BOUNDS,
            spec_conc_bound=SPEC_CONC_BOUNDS,    
            preset_flux=preset_flux, 
            preset_conc=preset_conc,   
            ex_thermo_cons=ex_thermo_cons,
            inc_enz_cons=inc_enz_cons,
            enz_prot_lb=Q,
        ).solve(solver='gurobi', n_jobs=N_JOBS)
        
        dgp_ranges.append(
            pd.DataFrame(
                res.gibbs_energy_ranges, 
                index=[[ph_out, ph_out], ['lb', 'ub']]
            )
        )

        reset_kcat(model, old_kcats, inc_enz_cons)

    pd.concat(dgp_ranges).T.to_excel(f'{out_dir}/gibbs_energy_ranges.xlsx')




if __name__ == '__main__':
    etva(f'{DIR}/variability/{GROWTH_TYPE}')
    
