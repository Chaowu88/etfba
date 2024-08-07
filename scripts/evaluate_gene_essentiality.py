r'''
Estimate the growth by deleting each of the essential, substitutable and blocked reactions 
at single (low) pH value.

ssh -m hmac-sha2-512 cwu@kestrel.hpc.nrel.gov
Woyaolvkai123!
srun --time=60 --account=runflux --partition=debug --ntasks=100 --pty $SHELL
srun --time=240 --account=runflux --partition=short --ntasks=100 --pty $SHELL
srun --time=12:00:00 --account=runflux --partition=standard --ntasks=36 --pty $SHELL
conda activate etfba-py38

python /home/cwu/Projects/Papers/pH_effect/two_step_optimization_GSM/5_simulate_growth/6_estimate_growth_on_deletion.py

~ 20 min (etfba, 100 jobs)
'''


import os
import numpy as np
import pandas as pd
from multiprocess import Pool
import sys
# ETFBA_PATH = '/home/cwu/Projects/ETFBA'
# sys.path.append(ETFBA_PATH)
from etfba import Model
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from utils import (add_H_leak_rxn, set_kcat_and_MW, get_default_KaC_and_KdC, 
                   set_deltaGprimem, get_H_leak_flux, set_adjusted_kcat, reset_kcat,
                   set_cpu_affinity)
from constants import (WORKING_DIR, MODEL_FILE, MW_FILE, 
                       GIBBS_ENERGY_FILE,
                       KCAT_FILE, CORR_KAPP, ADJUST_KCAT, KD_FILE, SEP_RXNS,
                       OBJECTIVE, 
                       PRESET_FLUXES, FLUX_BOUNDS, SPEC_FLUX_BOUNDS,
                       PRESET_CONCS, CONC_BOUNDS, SPEC_CONC_BOUNDS,
                       Q,
                       AEROBIC,
                       INTER_ESSEN_GENES_INFO_FILE, TRADIS_ESSEN_GENES_INFO_FILE,
                       KEIO_ESSEN_GENES_INFO_FILE, PEC_ESSEN_GENES_INFO_FILE,
                       ALL_GENES_INFO_FILE2,
                       ALL_GENES_INFO_FILE,
                       ETFBA_FLUX_FILE,
                       ETFBA_ESSEN_RXNS_LOW_PH_FILE, ETFBA_SUBST_RXNS_LOW_PH_FILE, 
                       ETFBA_BLOCKED_RXNS_LOW_PH_FILE,
                       ETFBA_ESSEN_RXNS_MED_PH_FILE, ETFBA_SUBST_RXNS_MED_PH_FILE, 
                       ETFBA_BLOCKED_RXNS_MED_PH_FILE,
                       FBA_ESSEN_RXNS_MED_PH_FILE, FBA_SUBST_RXNS_MED_PH_FILE, 
                       FBA_BLOCKED_RXNS_MED_PH_FILE, 
                       N_JOBS)


PH_LEVEL = 'MED'   # 'LOW' 'MED'
TYPE = 'all_genes2'   
# 'from_genes' 'all_genes2' 'all_enzymes' 'classified_enzymes' 'all_genes'

if PH_LEVEL == 'MED':
    DGMP_COL_IDX = -1
    PH_LABEL = '7.5'
    if TYPE == 'classified_enzymes':
        FBA_RXN_FILES = [
            FBA_ESSEN_RXNS_MED_PH_FILE, 
            FBA_SUBST_RXNS_MED_PH_FILE, 
            FBA_BLOCKED_RXNS_MED_PH_FILE
        ]
        ETFBA_RXN_FILES = [
            ETFBA_ESSEN_RXNS_MED_PH_FILE, 
            ETFBA_SUBST_RXNS_MED_PH_FILE, 
            ETFBA_BLOCKED_RXNS_MED_PH_FILE
        ]
    elif TYPE == 'from_genes':
        ESSEN_GENE_INFO_FILES = [
            INTER_ESSEN_GENES_INFO_FILE, 
            TRADIS_ESSEN_GENES_INFO_FILE,
            KEIO_ESSEN_GENES_INFO_FILE, 
            PEC_ESSEN_GENES_INFO_FILE
        ]
elif PH_LEVEL == 'LOW':
    DGMP_COL_IDX = 0
    PH_LABEL = '4'
    if TYPE == 'classified_enzymes':
        ETFBA_RXN_FILES = [
            ETFBA_ESSEN_RXNS_LOW_PH_FILE, 
            ETFBA_SUBST_RXNS_LOW_PH_FILE, 
            ETFBA_BLOCKED_RXNS_LOW_PH_FILE
        ]


def set_deleted_rxn_flux(model, rxnids):
    if not isinstance(rxnids, list):
        rxnids = rxnids
    
    preset_flux = {}
    for rxnid in rxnids:
        if model.reactions[rxnid].rev:
            preset_flux[rxnid+'_f'] = 0.0
            preset_flux[rxnid+'_b'] = 0.0
        else:
            preset_flux[rxnid] = 0.0

    return preset_flux


def fba_worker(mut_rxnid_chunk, model, dgpms_info):
    set_cpu_affinity()

    # run fba
    fluxes_all = {}
    for mut_rxnid in mut_rxnid_chunk:
        
        cyto_ph, peri_ph = dgpms_info.iloc[:2, DGMP_COL_IDX]

        # set fluxes
        preset_flux = {}
        preset_flux.update(get_H_leak_flux(cyto_ph, peri_ph))
        preset_flux.update(set_deleted_rxn_flux(model, mut_rxnid))
        preset_flux.update(PRESET_FLUXES)
        
        try:
            res = model.optimize(
                'fba', 
                objective=OBJECTIVE, 
                flux_bound=FLUX_BOUNDS, 
                spec_flux_bound=SPEC_FLUX_BOUNDS,
                preset_flux=preset_flux,
                parsimonious=False,
            ).solve(solver='gurobi')
            
            fluxes = res.opt_fluxes
        
        except:   # exception indicates infeasible solution
            fluxes = dict.fromkeys(model.reactions.keys(), 0.0)

        finally:
            fluxes_all[mut_rxnid] = pd.Series(fluxes)
    
    return fluxes_all


def etfba_worker(
        mut_rxnid_chunk, 
        model, 
        dgpms_info, 
        kd_data, 
        KaC_default, 
        KdC_default, 
        inc_enz_cons,
        ex_thermo_cons
    ):

    set_cpu_affinity()
    
    # run etfba
    if isinstance(mut_rxnid_chunk, (list, pd.Index)):
        mut_rxnid_chunk = dict(zip(mut_rxnid_chunk, mut_rxnid_chunk))
    elif isinstance(mut_rxnid_chunk, pd.Series):
        pass 

    fluxes_all = {}
    for mut_id, rxnids in mut_rxnid_chunk.items():
        
        cyto_ph, peri_ph = dgpms_info.iloc[:2, DGMP_COL_IDX]
        dgpms = dgpms_info.iloc[2:, DGMP_COL_IDX]

        # set standard Gibbs energy
        set_deltaGprimem(model, dgpms)

        # set fluxes
        preset_flux = {}
        preset_flux.update(get_H_leak_flux(cyto_ph, peri_ph))
        preset_flux.update(set_deleted_rxn_flux(model, rxnids))
        preset_flux.update(PRESET_FLUXES)

        # set concentrations
        preset_conc = {}
        preset_conc = {'h_c': 10**(-cyto_ph+3)}   
        #! do not set h_p because periplasmic [H+] may not be equal to external [H+]
        preset_conc.update(PRESET_CONCS)

        # set adjusted Kcat
        if ADJUST_KCAT:
            old_kcats = set_adjusted_kcat(
                model, 
                kd_data, 
                inc_enz_cons, 
                SEP_RXNS,  
                cyto_ph, 
                KaC_default, 
                KdC_default
            )
        
        try:
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
            
            fluxes = res.opt_fluxes
        
        except:   # exception indicates infeasible solution
            fluxes = dict.fromkeys(model.reactions.keys(), 0.0)

        finally:
            fluxes_all[mut_id] = pd.Series(fluxes)

            # reset kcat
            if ADJUST_KCAT:
                reset_kcat(model, old_kcats, inc_enz_cons)
    
    return fluxes_all


def fba_mutant(out_dir, rxn_file):
    model = Model.load(MODEL_FILE)

    # get reaction list to knock out
    mut_rxnids = pd.read_excel(rxn_file, header=0, index_col=0).index

    # add H+ leak
    add_H_leak_rxn(model)

    # get pH info
    dgpms_info = pd.read_excel(GIBBS_ENERGY_FILE, header=0, index_col=0)

    pool = Pool(processes=N_JOBS)

    async_res = []
    mut_rxnid_chunks = np.array_split(mut_rxnids, N_JOBS)
    for mut_rxnid_chunk in mut_rxnid_chunks:
        res = pool.apply_async(
            func=fba_worker,
            args=(mut_rxnid_chunk, model, dgpms_info)
        )
        async_res.append(res)
    
    pool.close()
    pool.join()
    
    results = {}
    for res in async_res:
        results.update(res.get())
    
    os.makedirs(out_dir, exist_ok = True)
    pd.DataFrame(results).to_csv(   
        f'{out_dir}/mutant_fluxes.tsv', 
        header = True, index = True, sep = '\t'
    )   # rows are all reactions, columns are deleted reactions


def etfba_mutant(out_dir, rxn_file):
    model = Model.load(MODEL_FILE)

    # get reaction list to knock out
    if TYPE in ['from_genes', 'all_genes2', 'all_genes']:
        info = pd.read_excel(rxn_file, header=0, index_col=0)
        info_exist = info.dropna(subset=['essential_reactions'])
        mut_rxnids = info_exist['essential_reactions'].str.split(',')
    elif TYPE == 'classified_enzymes':
        mut_rxnids = pd.read_excel(rxn_file, header=0, index_col=0).index
    elif TYPE == 'all_enzymes':
        mut_rxnids = pd.read_excel(
            rxn_file, header=0, index_col=0
        ).iloc[2:,:].index

    # add H+ leak
    add_H_leak_rxn(model)

    # get reactions excluded from thermodynamic constraints
    dgpms_info = pd.read_excel(GIBBS_ENERGY_FILE, header=0, index_col=0)
    
    ex_thermo_cons = [rxnid for rxnid in model.reactions 
                      if rxnid not in dgpms_info.index]

    # set kcat and mw, get enzymes subject to enzyme protein cost constraint
    mw_data = pd.read_excel(MW_FILE, header=0, index_col=0).squeeze()
    
    kcat_data = pd.read_excel(KCAT_FILE, header=0, index_col=0)
    for rxnid in CORR_KAPP:
        kcat_data.loc[rxnid] = CORR_KAPP[rxnid]

    inc_enz_cons = set_kcat_and_MW(model, kcat_data, mw_data)

    # set default KaC and KdC
    if ADJUST_KCAT:
        kd_data = pd.read_excel(KD_FILE, header=0, index_col=0)
        KaC_default, KdC_default = get_default_KaC_and_KdC(kd_data)

    pool = Pool(processes=N_JOBS)

    async_res = []
    mut_rxnid_chunks = np.array_split(mut_rxnids, N_JOBS)
    for mut_rxnid_chunk in mut_rxnid_chunks:
        res = pool.apply_async(
            func=etfba_worker,
            args=(
                mut_rxnid_chunk, 
                model, 
                dgpms_info, 
                kd_data, 
                KaC_default, 
                KdC_default, 
                inc_enz_cons,
                ex_thermo_cons
            )
        )
        async_res.append(res)
    
    pool.close()
    pool.join()
    
    results = {}
    for res in async_res:
        results.update(res.get())
    
    os.makedirs(out_dir, exist_ok = True)
    pd.DataFrame(results).to_csv(   
        f'{out_dir}/mutant_fluxes.tsv', 
        header = True, 
        index = True, 
        sep = '\t'
    )   # rows are all reactions, columns are deleted reactions




if __name__ == '__main__':
    
    # flag = '_aerobic' if AEROBIC else ''
    
    # # for rxn_file, rxn_type in zip(
    # #     FBA_RXN_FILES,
    # #     ['essential', 'substitutable', 'blocked']
    # # ):
    # #     fba_mutant(
    # #         f'{WORKING_DIR}/5_simulate_growth/new_ph_curve'
    # #         f'/fba{flag}_ph=4-7.5_8/deletion/ph={PH_LABEL}/{rxn_type}',
    # #         rxn_file
    # #     )
    
    # for rxn_file, rxn_type in zip(
    #     ETFBA_RXN_FILES,
    #     ['essential', 'substitutable', 'blocked']
    # ):
    #     etfba_mutant(
    #         f'{WORKING_DIR}/5_simulate_growth/new_ph_curve'
    #         f'/complete_etfba{flag}_ph=4-7.5_8_Q=0.19_adjust_kcat_new_ph_curve'
    #         f'/deletion/ph={PH_LABEL}/{rxn_type}/test',
    #         rxn_file
    #     )

    flag = '_aer' if AEROBIC else '_ana'
    # compare with essential gene databases
    if TYPE == 'from_genes':
        for essen_genes_info_file, database in zip(
            ESSEN_GENE_INFO_FILES,
            ['intersect', 'TraDIS', 'Keio', 'PEC']
        ):
            etfba_mutant(
                f'{WORKING_DIR}/5_simulate_growth/new_ph_curve2'
                f'/complete_etfba{flag}_ph=4-7.5_8_Q=0.19_adjust_kcat'
                f'/deletion/ph={PH_LABEL}/{TYPE}/{database}',
                essen_genes_info_file
            )

    # compare with all genes from Monk paper
    if TYPE == 'all_genes2':
        etfba_mutant(
            f'{WORKING_DIR}/5_simulate_growth/new_ph_curve2'
            f'/complete_etfba{flag}_ph=4-7.5_8_Q=0.19_adjust_kcat'
            f'/deletion/ph={PH_LABEL}/all_genes_Monk_paper',
            ALL_GENES_INFO_FILE2
        )

    elif TYPE == 'all_genes':
        etfba_mutant(
            f'{WORKING_DIR}/5_simulate_growth/new_ph_curve2'
            f'/complete_etfba{flag}_ph=4-7.5_8_Q=0.19_adjust_kcat'
            f'/deletion/ph={PH_LABEL}/{TYPE}',
            ALL_GENES_INFO_FILE
        )

    elif TYPE == 'all_enzymes':
        etfba_mutant(
            f'{WORKING_DIR}/5_simulate_growth/new_ph_curve2'
            f'/complete_etfba{flag}_ph=4-7.5_8_Q=0.19_adjust_kcat'
            f'/deletion/ph={PH_LABEL}/{TYPE}',
            ETFBA_FLUX_FILE
        )
    
    elif TYPE == 'classified_enzymes':
        for rxn_file, rxn_type in zip(
            ETFBA_RXN_FILES,
            ['essential', 'substitutable', 'blocked']
        ):
            etfba_mutant(
                f'{WORKING_DIR}/5_simulate_growth/new_ph_curve2'
                f'/complete_etfba{flag}_ph=4-7.5_8_Q=0.19_adjust_kcat'
                f'/deletion/ph={PH_LABEL}/{TYPE}/{rxn_type}',
                rxn_file
            )
    