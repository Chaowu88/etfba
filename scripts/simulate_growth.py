'''Simulate anaerobic and aerobic growth at various external pHs.

Usage:
python /path/to/etfba-main/scripts/simulate_growth.py
'''


import os
import pandas as pd
from etfba import Model
from helper import (get_internal_ph, add_H_leak_rxn, set_kcat_and_MW, 
                    get_default_KaC_and_KdC, set_deltaGprimem, get_H_leak_flux, 
                    set_adjusted_kcat, reset_kcat)
from config import (DIR, MODEL_FILE, MW_FILE, GIBBS_ENERGY_FILE, KCAT_FILE, 
                    CORR_KAPP, KD_FILE, SEP_RXNS, OBJECTIVE, 
                    GROWTH_TYPE, PRESET_FLUXES, FLUX_BOUNDS, SPEC_FLUX_BOUNDS, 
                    PRESET_CONCS, CONC_BOUNDS, SPEC_CONC_BOUNDS, Q)
    

def etfba(out_dir):
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
    
    fluxes = []
    dgps = []
    epcs = []
    for ph_out, dgpms in dgpms_info.items():
        print('\nRun simulation at external pH', ph_out)
        
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
        
        print(f"Growth rate: {res.opt_fluxes['BIOMASS_Ec_iML1515_core_75p37M']:.3f}")
        print(f"Glucose uptake: {res.opt_fluxes['EX_glc__D_e']:.3f}")
        print(f"Acetate formation: {res.opt_fluxes['EX_ac_e']:.3f}")
        print(f"formate formation: {res.opt_fluxes['EX_for_e']:.3f}")
        print(f"Ethanol formation: {res.opt_fluxes['EX_etoh_e']:.3f}")
        print(f"Lactate formation: {res.opt_fluxes['EX_lac__D_e']:.3f}")
        print(f"Succinate formation: {res.opt_fluxes['EX_succ_e']:.3f}")

        if res.optimization_successful:
            fluxes.append(pd.Series(res.opt_fluxes))
            dgps.append(pd.Series(res.opt_gibbs_energy))
            epcs.append(pd.Series(res.opt_enzyme_costs))

        reset_kcat(model, old_kcats, inc_enz_cons)
    
    pd.DataFrame(
        fluxes, index=dgpms_info.columns
    ).T.to_excel(f'{out_dir}/metabolic_fluxes.xlsx')
    pd.DataFrame(
        dgps, index=dgpms_info.columns
    ).T.to_excel(f'{out_dir}/gibbs_energies.xlsx')
    pd.DataFrame(
        epcs, index=dgpms_info.columns
    ).T.to_excel(f'{out_dir}/enzyme_protein_costs.xlsx')
    



if __name__ == '__main__':
    etfba(f'{DIR}/simulation/{GROWTH_TYPE}')