DIR = '/path/to/etfba-main/scripts'

GROWTH_TYPE = 'anaerobic'   # {'anaerobic', 'aerobic'}

MODEL_FILE = f'{DIR}/model.bin'
GIBBS_ENERGY_FILE = f'{DIR}/dgpms.xlsx'
KCAT_FILE = f'{DIR}/kcats.xlsx'
CORR_KAPP = {
    'GLCptspp': 240.3936579,   
}   # estimated based on the specific activity data from Brenda
KD_FILE = f'{DIR}/params.xlsx'
SEP_RXNS = ['DHORTS', 'M1PD', 'TRSARr']   
# these reactions have seperate forward and backward kcats
MW_FILE = f'{DIR}/mws.xlsx'

OBJECTIVE = {'BIOMASS_Ec_iML1515_core_75p37M': 1}
FLUX_BOUNDS = (0, 1000)
SPEC_FLUX_BOUNDS = {
    'EX_glc__D_e_b': (0, 1000)
}
CONC_BOUNDS = (0.0001, 200)
PRESET_CONCS = {
    'glc__D_p': 20,   
    'pi_p': 56,   
    'so4_p': 3,   
    'nh4_p': 19,   
    'na1_p': 160,   
    'k_p': 22,   
    'fe2_p': 62,   
}
if GROWTH_TYPE == 'aerobic':
    PRESET_FLUXES = {
        'EX_glc__D_e_f': 0, 
        'EX_fru_e': 0,
        'FHL': 0,   
        'F6PA_f': 0,   
        'F6PA_b': 0,   
        'DHAPT': 0, 
    }
    PRESET_CONCS.update({
        'o2_p': 0.0082,  
        'co2_p': 0.1,
    })
    SPEC_CONC_BOUNDS = {
    'o2_c': (0.0001, 0.0082),   
    'co2_c': (0.1, 200)   
    }
else:
    PRESET_FLUXES = {
        'EX_glc__D_e_f': 0, 
        'EX_fru_e': 0, 
        'EX_o2_e_f': 0,   
        'EX_o2_e_b': 0,   
        'F6PA_f': 0,   
        'F6PA_b': 0,   
        'DHAPT': 0, 
    }   
    PRESET_CONCS.update({
        'o2_p': 1e-9,   
        'o2_c': 1e-9,   
        'co2_p': 0.1,    
    })    
    SPEC_CONC_BOUNDS = {
        'co2_c': (0.1, 200)
    }   
Q = 0.19

LOW_PH = 4.0
NEUTRAL_PH = 7.5

ETFBA_FLUX_FILE = f'{DIR}/simulation/{GROWTH_TYPE}/metabolic_fluxes.xlsx'
# used for variability analysis

PERTURB_DIRECTION = 'downregulation'   # {'upregulation', 'downregulation'}
PERTURB_FOLD = 5

N_JOBS = 100   
