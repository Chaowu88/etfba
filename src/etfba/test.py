#!/usr/bin/env python
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '12/12/2021'
__version__ = '1.0'


r'''
python C:\Users\cwu\Desktop\Software\ETFBA\etfba\test.py
python /home/cwu/Projects/ETFBA/etfba/test.py
'''


import pandas as pd
import pyomo.environ as pyo
import platform  


system = platform.system()
if system == 'Linux':
    ETFBA_PATH = '/home/cwu/Projects/ETFBA'
    OUTPUT_PATH = '/home/cwu/Projects/ETFBA/results/test'
    ECOLI_MODEL_FILE1 = '/home/cwu/Projects/ETFBA/Data/ecoli_core/E_coli_core.xlsx'
    FIXED_FLUXES_FILE = OUTPUT_PATH+'/etfba_fluxes.xlsx'

elif system == 'Windows':
    ETFBA_PATH = r'C:\Users\cwu\Desktop\Software\ETFBA'
    SCIP_PATH = r'C:\Users\cwu\Desktop\Software\solvers\scipampl'
    OUTPUT_PATH = r'C:\Users\cwu\Desktop'

    TCA_MODEL_FILE = r'C:\Users\cwu\Desktop\Software\ETFBA\Data\tca\demo_TCA.xlsx'

    ECOLI_MODEL_FILE1 = r'C:\Users\cwu\Desktop\Software\ETFBA\Data\ecoli_core\E_coli_core.xlsx'
    ECOLI_MODEL_FILE_ORI = r'C:\Users\cwu\Desktop\Software\ETFBA\Data\ecoli_core\E_coli_core_ori.xlsx'
    FIXED_FLUXES_FILE = r'C:\Users\cwu\Desktop\etfba_fluxes.xlsx'
    
    ECOLI_MODEL_FILE2 = r'C:\Users\cwu\Desktop\Software\ETFBA\Data\ecoli_core\E_coli_core2.xlsx'
    ECOLI_MODEL_SPLIT_ATPS4R_FILE1 = r'C:\Users\cwu\Desktop\Software\ETFBA\Data\ecoli_core\E_coli_core_split_atps4r.xlsx'
    ECOLI_MODEL_SPLIT_ATPS4R_FILE2 = r'C:\Users\cwu\Desktop\Software\ETFBA\Data\ecoli_core\E_coli_core2_split_atps4r.xlsx'
    ECOLI_MODEL_NEW_FILE1 = r'C:\Users\cwu\Desktop\Software\ETFBA\Data\ecoli_core\E_coli_core_new1.xlsx'
    ECOLI_MODEL_P_TYPE = r'C:\Users\cwu\Desktop\Software\ETFBA\Data\ecoli_core\E_coli_core2_p_type.xlsx'
    ECOLI_AERO_INI_FLUXES_FILE1 = r'C:\Users\cwu\Desktop\Software\ETFBA\Data\ecoli_core\ini\aero_iniFluxes.bin'
    ECOLI_AERO_INI_FLUXES_FILE2 = r'C:\Users\cwu\Desktop\Software\ETFBA\Data\ecoli_core\ini\aero_iniFluxes2.bin'
    ECOLI_AERO_INI_CONCS_FILE = r'C:\Users\cwu\Desktop\Software\ETFBA\Data\ecoli_core\ini\aero_iniConcs.bin'
    ECOLI_AERO_INI_LNCONCS_FILE = r'C:\Users\cwu\Desktop\Software\ETFBA\Data\ecoli_core\ini\aero_iniLnConcs.bin'
    ECOLI_ANA_INI_FLUXES_FILE1 = r'C:\Users\cwu\Desktop\Software\ETFBA\Data\ecoli_core\ini\ana_iniFluxes.xlsx'
    ECOLI_ANA_INI_FLUXES_FILE2 = r'C:\Users\cwu\Desktop\Software\ETFBA\Data\ecoli_core\ini\ana_iniFluxes2.bin'
    ECOLI_ANA_INI_CONCS_FILE = r'C:\Users\cwu\Desktop\Software\ETFBA\Data\ecoli_core\ini\ana_iniConcs.xlsx'
    ECOLI_ANA_INI_LNCONCS_FILE = r'C:\Users\cwu\Desktop\Software\ETFBA\Data\ecoli_core\ini\ana_iniLnConcs.bin'

    ECOLI_MODEL_FILE_NEW = r'C:\Users\cwu\Desktop\Software\ETFBA\Data\ecoli_core\E_coli_core_new3.xlsx'




def test_glpk():
    
    model = pyo.ConcreteModel()

    model.x = pyo.Var([1,2,3,4], domain = pyo.NonNegativeReals)
    model.y = pyo.Var(domain = pyo.Binary)

    model.obj = pyo.Objective(expr = 2*model.x[1] + 3*model.x[2] + 3*model.x[3] + 4*model.x[4])

    model.c1 = pyo.Constraint(expr = 3*model.x[1] + 4*model.x[2] >= 1)
    model.c2 = pyo.Constraint(expr = 3*model.x[3] + 4*model.x[4] >= 1)
    model.c3 = pyo.Constraint(expr = model.x[2] <= 10000*(1 - model.y))
    model.c4 = pyo.Constraint(expr = model.x[2] <= 10000*(model.y))
    
    pyo.SolverFactory('glpk').solve(model)
    
    model.obj.display()
    model.pprint()
    

def test_mindtpy():
    
    model = pyo.ConcreteModel()

    model.x = pyo.Var(bounds = (1.0, 10.0), initialize = 5.0)
    model.y = pyo.Var(domain = pyo.Binary)

    model.c1 = pyo.Constraint(expr = (model.x - 4.0)**2 - model.x <= 50.0*(1 - model.y))
    model.c2 = pyo.Constraint(expr = model.x*pyo.log(model.x) + 5.0 <= 50.0*(model.y))

    model.objective = pyo.Objective(expr = model.x, sense = pyo.minimize)

    pyo.SolverFactory('mindtpy').solve(model, mip_solver='glpk', nlp_solver='ipopt')
    
    model.objective.display()
    #model.display()
    model.pprint()
    
    
def test_scip():
    
    model = pyo.ConcreteModel()

    model.x = pyo.Var(bounds = (1.0, 10.0), initialize = 5.0)
    model.y = pyo.Var(domain = pyo.Binary)

    model.c1 = pyo.Constraint(expr = (model.x - 4.0)**2 - model.x <= 50.0*(1 - model.y))
    model.c2 = pyo.Constraint(expr = model.x*pyo.log(model.x) + 5.0 <= 50.0*(model.y))

    model.objective = pyo.Objective(expr = model.x, sense = pyo.minimize)

    pyo.SolverFactory('scipampl', executable = SCIP_PATH).solve(model)
    
    model.objective.display()
    #model.display()
    model.pprint()

    
def test_apopt1():
    
    model = pyo.ConcreteModel()

    model.x = pyo.Var(bounds = (1.0, 10.0), initialize = 5.0)
    model.y = pyo.Var(domain = pyo.Binary)

    model.c1 = pyo.Constraint(expr = (model.x - 4.0)**2 - model.x <= 50.0*(1 - model.y))
    model.c2 = pyo.Constraint(expr = model.x*pyo.log(model.x) + 5.0 <= 50.0*(model.y))

    model.objective = pyo.Objective(expr = model.x, sense = pyo.minimize)

    #pyo.SolverFactory('mindtpy').solve(model, mip_solver='glpk', nlp_solver='ipopt')
    pyo.SolverFactory('apopt.py').solve(model)
    
    #model.objective.display()
    #model.display()
    model.pprint()


def test_apopt2():
    
    model = pyo.ConcreteModel()

    model.x = pyo.Var([1,2,3,4], domain = pyo.NonNegativeReals)
    model.y = pyo.Var(domain = pyo.Binary)

    model.obj = pyo.Objective(expr = 2*model.x[1] + 3*model.x[2] + 3*model.x[3] + 4*model.x[4])

    model.c1 = pyo.Constraint(expr = 3*model.x[1] + 4*model.x[2] >= 1)
    model.c2 = pyo.Constraint(expr = 3*model.x[3] + 4*model.x[4] >= 1)
    model.c3 = pyo.Constraint(expr = 1000*pyo.cos(model.x[3]) <= 1000)
    model.c4 = pyo.Constraint(expr = 1000*pyo.sin(model.x[4]) <= 1000)
    model.c5 = pyo.Constraint(expr = model.x[2] <= 10000*(1 - model.y))
    model.c6 = pyo.Constraint(expr = model.x[2] <= 10000*(model.y))
    
    pyo.SolverFactory('apopt.py').solve(model)
    
    model.obj.display()
    model.pprint()
    
    
def test_apopt3():
        
    model = pyo.ConcreteModel()
    
    model.x = pyo.Var([1,2,3,4], domain = pyo.Binary)
    
    model.obj = pyo.Objective(expr = model.x[1] + model.x[3], sense = pyo.maximize)
    
    model.c1 = pyo.Constraint(expr = model.x[1]**2 + model.x[2]**2 == 1)
    model.c2 = pyo.Constraint(expr = model.x[3]**2 + model.x[4]**2 == 1)
    
    pyo.SolverFactory('apopt.py').solve(model)
    
    model.obj.display()
    model.pprint()    
    
    
def test_apopt4():
        
    model = pyo.ConcreteModel()
    
    model.x = pyo.Var([1,2,3,4], domain = pyo.Binary)
    model.y = pyo.Var(['a','b'])
    
    model.obj = pyo.Objective(expr = model.x[1]*model.y['a'] - model.x[2]*model.y['b'], sense = pyo.maximize)
    
    model.c1 = pyo.Constraint(expr = model.x[1]**2 + model.x[2]**2 == 1)
    model.c2 = pyo.Constraint(expr = model.x[3]**2 + model.x[4]**2 == 1)
    model.c3 = pyo.Constraint(expr = model.y['a'] + model.y['b']**2 <= 10)
    
    pyo.SolverFactory('apopt.py').solve(model)
    
    model.obj.display()
    model.pprint()    
    

def test_TCA():
    # use binary variables
    
    Q = 0.01
    K = 10000
    R = 8.315e-3
    T = 298.15
    
    MW = 40
    kcat = 200
    Km = 0.2
    lnKm = pyo.log(Km)
    
    dGp0_2 = -19.3
    dGp0_3 = 2.5
    dGp0_4 = -18.8
    dGp0_5 = 4.6
    dGp0_6 = 5.5
    dGp0_7 = -19.1
    dGp0_8 = 4.1
    
    model = pyo.ConcreteModel()
    
    model.v = pyo.Var(['v1f', 'v1b', 'v2f', 'v2b', 'v3f', 'v3b', 'v4f', 'v4b', 'v5f', 'v5b', 'v6f', 'v6b', 
                       'v7f', 'v7b', 'v8f', 'v8b', 'v9f', 'v9b'], bounds = (0, 100))   # 0 - 100 mmol/gCDW/h
    model.lnc = pyo.Var(['accoa', 'icit', 'akg', 'suc', 'mal', 'oaa', 'gox'], bounds = (-6.9, 4.6))   # 0.001 mM - 100 mM
    model.x = pyo.Var(['v2f', 'v3f', 'v4f', 'v5f', 'v6f', 'v7f', 'v8f'], domain = pyo.Binary)   # no v1 and v9
    
    #model.obj = pyo.Objective(expr = model.v['v9f'] - model.v['v9b'], sense = pyo.maximize)
    model.obj = pyo.Objective(expr = (model.v['v9f'] - model.v['v9b'])/(1/3600*(
    model.x['v2f']*model.v['v2f']/kcat*(1 + pyo.exp(2*lnKm - model.lnc['accoa'] - model.lnc['oaa']))/(1 - pyo.exp(model.lnc['icit'] - model.lnc['accoa'] - model.lnc['oaa'] + dGp0_2/R/T)) + (1 - model.x['v2f'])*model.v['v2b']/kcat*(1 + pyo.exp(lnKm - model.lnc['icit']))/(1 - pyo.exp(model.lnc['accoa'] + model.lnc['oaa'] - model.lnc['icit'] - dGp0_2/R/T)) +
    model.x['v3f']*model.v['v3f']/kcat*(1 + pyo.exp(lnKm - model.lnc['icit']))/(1 - pyo.exp(model.lnc['akg'] - model.lnc['icit'] + dGp0_3/R/T)) + (1 - model.x['v3f'])*model.v['v3b']/kcat*(1 + pyo.exp(lnKm - model.lnc['akg']))/(1 - pyo.exp(model.lnc['icit'] - model.lnc['akg'] - dGp0_3/R/T)) +
    model.x['v4f']*model.v['v4f']/kcat*(1 + pyo.exp(lnKm - model.lnc['akg']))/(1 - pyo.exp(model.lnc['suc'] - model.lnc['akg'] + dGp0_4/R/T)) + (1 - model.x['v4f'])*model.v['v4b']/kcat*(1 + pyo.exp(lnKm - model.lnc['suc']))/(1 - pyo.exp(model.lnc['akg'] - model.lnc['suc'] - dGp0_4/R/T)) +
    model.x['v5f']*model.v['v5f']/kcat*(1 + pyo.exp(lnKm - model.lnc['suc']))/(1 - pyo.exp(model.lnc['mal'] - model.lnc['suc'] + dGp0_5/R/T)) + (1 - model.x['v5f'])*model.v['v5b']/kcat*(1 + pyo.exp(lnKm - model.lnc['mal']))/(1 - pyo.exp(model.lnc['suc'] - model.lnc['mal'] - dGp0_5/R/T)) +
    model.x['v6f']*model.v['v6f']/kcat*(1 + pyo.exp(lnKm - model.lnc['mal']))/(1 - pyo.exp(model.lnc['oaa'] - model.lnc['mal'] + dGp0_6/R/T)) + (1 - model.x['v6f'])*model.v['v6b']/kcat*(1 + pyo.exp(lnKm - model.lnc['oaa']))/(1 - pyo.exp(model.lnc['mal'] - model.lnc['oaa'] - dGp0_6/R/T)) +
    model.x['v7f']*model.v['v7f']/kcat*(1 + pyo.exp(lnKm - model.lnc['icit']))/(1 - pyo.exp(model.lnc['suc'] + model.lnc['gox'] - model.lnc['icit'] + dGp0_7/R/T)) + (1 - model.x['v7f'])*model.v['v7b']/kcat*(1 + pyo.exp(2*lnKm - model.lnc['suc'] - model.lnc['gox']))/(1 - pyo.exp(model.lnc['icit'] - model.lnc['suc'] - model.lnc['gox'] - dGp0_7/R/T)) +
    model.x['v8f']*model.v['v8f']/kcat*(1 + pyo.exp(2*lnKm - model.lnc['accoa'] - model.lnc['gox']))/(1 - pyo.exp(model.lnc['mal'] - model.lnc['accoa'] - model.lnc['gox'] + dGp0_8/R/T)) + (1 - model.x['v8f'])*model.v['v8b']/kcat*(1 + pyo.exp(lnKm - model.lnc['mal']))/(1 - pyo.exp(model.lnc['accoa'] + model.lnc['gox'] - model.lnc['mal'] - dGp0_8/R/T)))*MW), sense = pyo.maximize)
    
    # mass balance
    model.mb1 = pyo.Constraint(expr = model.v['v1f'] + model.v['v2b'] + model.v['v8b'] + 6.85*model.v['v9b'] - 
                                      model.v['v1b'] - model.v['v2f'] - model.v['v8f'] - 6.85*model.v['v9f'] == 0)   # accoa
    model.mb2 = pyo.Constraint(expr = model.v['v2f'] + model.v['v3b'] + model.v['v7b'] - 
                                      model.v['v2b'] - model.v['v3f'] - model.v['v7f'] == 0)   # icit
    model.mb3 = pyo.Constraint(expr = model.v['v3f'] + model.v['v4b'] + 0.75*model.v['v9b'] - 
                                      model.v['v3b'] - model.v['v4f'] - 0.75*model.v['v9f'] == 0)   # akg
    model.mb4 = pyo.Constraint(expr = model.v['v4f'] + model.v['v5b'] + model.v['v7f'] - 
                                      model.v['v4b'] - model.v['v5f'] - model.v['v7b'] == 0)   # suc
    model.mb5 = pyo.Constraint(expr = model.v['v5f'] + model.v['v6b'] + model.v['v8f'] - 
                                      model.v['v5b'] - model.v['v6f'] - model.v['v8b'] == 0)   # mal
    model.mb6 = pyo.Constraint(expr = model.v['v2b'] + model.v['v6f'] + 0.62*model.v['v9b'] - 
                                      model.v['v2f'] - model.v['v6b'] - 0.62*model.v['v9f']== 0)   # oaa
    model.mb7 = pyo.Constraint(expr = model.v['v7f'] + model.v['v8b'] - 
                                      model.v['v7b'] - model.v['v8f'] == 0)   # gox
    
    model.mb8 = pyo.Constraint(expr = model.v['v1f'] == 50)   # substrate uptake
    model.mb9 = pyo.Constraint(expr = model.v['v1b'] == 0)
    model.mb10 = pyo.Constraint(expr = model.v['v9b'] == 0)
    '''
    model.mb11 = pyo.Constraint(expr = model.v['v2f'] <= model.x['v2f']*100)   # v2f
    model.mb12 = pyo.Constraint(expr = model.v['v2b'] <= (1 - model.x['v2f'])*100)   # v2b
    model.mb13 = pyo.Constraint(expr = model.v['v3f'] <= model.x['v3f']*100)   # v3f
    model.mb14 = pyo.Constraint(expr = model.v['v3b'] <= (1 - model.x['v3f'])*100)   # v3b
    model.mb15 = pyo.Constraint(expr = model.v['v4f'] <= model.x['v4f']*100)   # v4f
    model.mb16 = pyo.Constraint(expr = model.v['v4b'] <= (1 - model.x['v4f'])*100)   # v4b
    model.mb17 = pyo.Constraint(expr = model.v['v5f'] <= model.x['v5f']*100)   # v5f
    model.mb18 = pyo.Constraint(expr = model.v['v5b'] <= (1 - model.x['v5f'])*100)   # v5b
    model.mb19 = pyo.Constraint(expr = model.v['v6f'] <= model.x['v6f']*100)   # v6f
    model.mb20 = pyo.Constraint(expr = model.v['v6b'] <= (1 - model.x['v6f'])*100)   # v6b
    model.mb21 = pyo.Constraint(expr = model.v['v7f'] <= model.x['v7f']*100)   # v7f
    model.mb22 = pyo.Constraint(expr = model.v['v7b'] <= (1 - model.x['v7f'])*100)   # v7b
    model.mb23 = pyo.Constraint(expr = model.v['v8f'] <= model.x['v8f']*100)   # v8f
    model.mb24 = pyo.Constraint(expr = model.v['v8b'] <= (1 - model.x['v8f'])*100)   # v8b
    '''
    '''
    # enzyme protein constraints
    model.epc = pyo.Constraint(expr = 1/3600*(
    model.x['v2f']*model.v['v2f']/kcat*(1 + pyo.exp(2*lnKm - model.lnc['accoa'] - model.lnc['oaa']))/(1 - pyo.exp(model.lnc['icit'] - model.lnc['accoa'] - model.lnc['oaa'] + dGp0_2/R/T)) + (1 - model.x['v2f'])*model.v['v2b']/kcat*(1 + pyo.exp(lnKm - model.lnc['icit']))/(1 - pyo.exp(model.lnc['accoa'] + model.lnc['oaa'] - model.lnc['icit'] - dGp0_2/R/T)) +
    model.x['v3f']*model.v['v3f']/kcat*(1 + pyo.exp(lnKm - model.lnc['icit']))/(1 - pyo.exp(model.lnc['akg'] - model.lnc['icit'] + dGp0_3/R/T)) + (1 - model.x['v3f'])*model.v['v3b']/kcat*(1 + pyo.exp(lnKm - model.lnc['akg']))/(1 - pyo.exp(model.lnc['icit'] - model.lnc['akg'] - dGp0_3/R/T)) +
    model.x['v4f']*model.v['v4f']/kcat*(1 + pyo.exp(lnKm - model.lnc['akg']))/(1 - pyo.exp(model.lnc['suc'] - model.lnc['akg'] + dGp0_4/R/T)) + (1 - model.x['v4f'])*model.v['v4b']/kcat*(1 + pyo.exp(lnKm - model.lnc['suc']))/(1 - pyo.exp(model.lnc['akg'] - model.lnc['suc'] - dGp0_4/R/T)) +
    model.x['v5f']*model.v['v5f']/kcat*(1 + pyo.exp(lnKm - model.lnc['suc']))/(1 - pyo.exp(model.lnc['mal'] - model.lnc['suc'] + dGp0_5/R/T)) + (1 - model.x['v5f'])*model.v['v5b']/kcat*(1 + pyo.exp(lnKm - model.lnc['mal']))/(1 - pyo.exp(model.lnc['suc'] - model.lnc['mal'] - dGp0_5/R/T)) +
    model.x['v6f']*model.v['v6f']/kcat*(1 + pyo.exp(lnKm - model.lnc['mal']))/(1 - pyo.exp(model.lnc['oaa'] - model.lnc['mal'] + dGp0_6/R/T)) + (1 - model.x['v6f'])*model.v['v6b']/kcat*(1 + pyo.exp(lnKm - model.lnc['oaa']))/(1 - pyo.exp(model.lnc['mal'] - model.lnc['oaa'] - dGp0_6/R/T)) +
    model.x['v7f']*model.v['v7f']/kcat*(1 + pyo.exp(lnKm - model.lnc['icit']))/(1 - pyo.exp(model.lnc['suc'] + model.lnc['gox'] - model.lnc['icit'] + dGp0_7/R/T)) + (1 - model.x['v7f'])*model.v['v7b']/kcat*(1 + pyo.exp(2*lnKm - model.lnc['suc'] - model.lnc['gox']))/(1 - pyo.exp(model.lnc['icit'] - model.lnc['suc'] - model.lnc['gox'] - dGp0_7/R/T)) +
    model.x['v8f']*model.v['v8f']/kcat*(1 + pyo.exp(2*lnKm - model.lnc['accoa'] - model.lnc['gox']))/(1 - pyo.exp(model.lnc['mal'] - model.lnc['accoa'] - model.lnc['gox'] + dGp0_8/R/T)) + (1 - model.x['v8f'])*model.v['v8b']/kcat*(1 + pyo.exp(lnKm - model.lnc['mal']))/(1 - pyo.exp(model.lnc['accoa'] + model.lnc['gox'] - model.lnc['mal'] - dGp0_8/R/T)))*MW <= Q)    
    '''
    
    # thermodynamics
    model.thmd1 = pyo.Constraint(expr = dGp0_2 + R*T*(model.lnc['icit'] - model.lnc['accoa'] - model.lnc['oaa']) 
                                        <= (1 - model.x['v2f'])*K-0.001)   # v2f
    model.thmd2 = pyo.Constraint(expr = -dGp0_2 + R*T*(model.lnc['accoa'] + model.lnc['oaa'] - model.lnc['icit']) 
                                        <= model.x['v2f']*K-0.001)   # v2b
    model.thmd3 = pyo.Constraint(expr = dGp0_3 + R*T*(model.lnc['akg'] - model.lnc['icit']) 
                                        <= (1 - model.x['v3f'])*K-0.001)   # v3f
    model.thmd4 = pyo.Constraint(expr = -dGp0_3 + R*T*(model.lnc['icit'] - model.lnc['akg']) 
                                        <= model.x['v3f']*K-0.001)   # v3b
    model.thmd5 = pyo.Constraint(expr = dGp0_4 + R*T*(model.lnc['suc'] - model.lnc['akg']) 
                                        <= (1 - model.x['v4f'])*K-0.001)   # v4f
    model.thmd6 = pyo.Constraint(expr = -dGp0_4 + R*T*(model.lnc['akg'] - model.lnc['suc']) 
                                        <= model.x['v4f']*K-0.001)   # v4b
    model.thmd7 = pyo.Constraint(expr = dGp0_5 + R*T*(model.lnc['mal'] - model.lnc['suc']) 
                                        <= (1 - model.x['v5f'])*K-0.001)   # v5f
    model.thmd8 = pyo.Constraint(expr = -dGp0_5 + R*T*(model.lnc['suc'] - model.lnc['mal']) 
                                        <= model.x['v5f']*K-0.001)   # v5b
    model.thmd9 = pyo.Constraint(expr = dGp0_6 + R*T*(model.lnc['oaa'] - model.lnc['mal']) 
                                        <= (1 - model.x['v6f'])*K-0.001)   # v6f
    model.thmd10 = pyo.Constraint(expr = -dGp0_6 + R*T*(model.lnc['mal'] - model.lnc['oaa']) 
                                         <= model.x['v6f']*K-0.001)   # v6b
    model.thmd11 = pyo.Constraint(expr = dGp0_7 + R*T*(model.lnc['suc'] + model.lnc['gox'] - model.lnc['icit']) 
                                         <= (1 - model.x['v7f'])*K-0.001)   # v7f
    model.thmd12 = pyo.Constraint(expr = -dGp0_7 + R*T*(model.lnc['icit'] - model.lnc['suc'] - model.lnc['gox']) 
                                         <= model.x['v7f']*K-0.001)   # v7b
    model.thmd13 = pyo.Constraint(expr = dGp0_8 + R*T*(model.lnc['mal'] + model.lnc['gox'] - model.lnc['accoa']) 
                                         <= (1 - model.x['v8f'])*K-0.001)   # v8f
    model.thmd14 = pyo.Constraint(expr = -dGp0_8 + R*T*(model.lnc['gox'] + model.lnc['accoa'] - model.lnc['mal']) 
                                         <= model.x['v8f']*K-0.001)   # v8b
    
    # 0-1 binary
    model.bin1 = pyo.Constraint(expr = model.v['v2f'] >= model.v['v2b'] - (1 - model.x['v2f'])*K)   # v2
    model.bin2 = pyo.Constraint(expr = model.v['v2f'] <= model.v['v2b'] + model.x['v2f']*K)   # v2
    model.bin3 = pyo.Constraint(expr = model.v['v3f'] >= model.v['v3b'] - (1 - model.x['v3f'])*K)   # v3
    model.bin4 = pyo.Constraint(expr = model.v['v3f'] <= model.v['v3b'] + model.x['v3f']*K)   # v3
    model.bin5 = pyo.Constraint(expr = model.v['v4f'] >= model.v['v4b'] - (1 - model.x['v4f'])*K)   # v4
    model.bin6 = pyo.Constraint(expr = model.v['v4f'] <= model.v['v4b'] + model.x['v4f']*K)   # v4
    model.bin7 = pyo.Constraint(expr = model.v['v5f'] >= model.v['v5b'] - (1 - model.x['v5f'])*K)   # v5
    model.bin8 = pyo.Constraint(expr = model.v['v5f'] <= model.v['v5b'] + model.x['v5f']*K)   # v5
    model.bin9 = pyo.Constraint(expr = model.v['v6f'] >= model.v['v6b'] - (1 - model.x['v6f'])*K)   # v6
    model.bin10 = pyo.Constraint(expr = model.v['v6f'] <= model.v['v6b'] + model.x['v6f']*K)   # v6
    model.bin11 = pyo.Constraint(expr = model.v['v7f'] >= model.v['v7b'] - (1 - model.x['v7f'])*K)   # v7
    model.bin12 = pyo.Constraint(expr = model.v['v7f'] <= model.v['v7b'] + model.x['v7f']*K)   # v7
    model.bin13 = pyo.Constraint(expr = model.v['v8f'] >= model.v['v8b'] - (1 - model.x['v8f'])*K)   # v8
    model.bin14 = pyo.Constraint(expr = model.v['v8f'] <= model.v['v8b'] + model.x['v8f']*K)   # v8
    
    #pyo.SolverFactory('scipampl', executable = SCIP_PATH).solve(model)
    pyo.SolverFactory('mindtpy').solve(model, mip_solver='glpk', nlp_solver='ipopt',
                                       strategy='OA', init_strategy = 'max_binary', tee = True)
    
    model.obj.display()
    model.display()     
    
    
def test_TCA_etfba():
    
    import sys
    sys.path.append(ETFBA_PATH)
    from scipy.linalg import null_space
    from etfba import Model
    
    tca = Model('TCA')
    #print(tca.metabolites)
    tca.read_from_excel(TCA_MODEL_FILE)
    #print(tca.metabolites)
    #print(tca)
    
    #print(tca.reactions['r1'])
    #print(tca.metabolites['AcCoA'])
    #print(tca.metabolites['accoa'].coes)
    #print(tca.metabolites['AcCoA'].kms)
    
    #print(tca.reactions['r2'].substrates['AcCoA'].km)
    #print(tca.metabolites['AcCoA'].kms['r2'])
    #print(tca.metabolites['AcCoA'].km)   # raise error
    #print(tca.reactions['r02'].substrates['accoa'].coe)
    
    #print(tca.stoichiometric_matrix)
    #print(null_space(tca.stoichiometric_matrix))
    #print(tca.total_stoichiometric_matrix)
    #print(tca.substrates)
    #print(tca.products)
    
    # res = tca.optimize('fba', 
    #                    objective = {'r11': 1}, 
    #                    preset_fluxes = {'r09': 10}, 
    #                    flux_bounds = (-50, 50),
    #                    excluded_mb = ['atp', 'amp', 'gtp', 'gdp', 'ppi', 'pi', 'nadh', 'nad', 'qh2', 'q', 'coa']
    #                    ).solve()
    # #print(res)
    # print('\nopt fluxes')
    # print(res.opt_fluxes)
    
    
    # res = tca.optimize('tfba', 
    #                    objective = {'r11': 1}, 
    #                    preset_fluxes = {'r09': 10}, preset_concs = {'ac.e': 5.6},
    #                    excluded_mb = ['atp', 'amp', 'gtp', 'gdp', 'ppi', 'pi', 'nadh', 'nad', 'qh2', 'q', 'coa'],
    #                    flux_bounds = (-50, 50)
    #                    ).solve()
    # #print(res)
    # print('\nopt fluxes')
    # print(res.opt_fluxes)
    # print('\nopt concs')
    # print(res.opt_concentrations)
    # print('\nopt gibbs energy')
    # print(res.opt_gibbs_energy)
    
    
    # res = tca.optimize('etfba', 
    #                    objective = {'r11': 1}, 
    #                    preset_fluxes = {'r09': 10}, preset_concs = None,   #{'ac.e': 5.6}
    #                    flux_bounds = (-50, 50), conc_bounds = (0.001, 10),
    #                    excluded_mb = ['atp', 'amp', 'gtp', 'gdp', 'ppi', 'pi', 'nadh', 'nad', 'qh2', 'q', 'coa'],
    #                    included_epcs = ['r01', 'r02', 'r03', 'r04', 'r05', 'r06', 'r07', 'r08', 'r09', 'r10'],
    #                    enzyme_protein_ub = 0.1
    #                    ).solve()
    # #print(res)
    # #print(res.opt_directions)
    # print('\nopt fluxes')
    # print(res.opt_fluxes)
    # '''
    # print('\nopt concs')
    # print(res.opt_concentrations)
    # print('\nopt gibbs energy')
    # print(res.opt_gibbs_energy)
    # '''
    # print('\nopt enzyme cost')
    # print(res.opt_total_enzyme_cost)
    # print(res.opt_enzyme_costs)
    
    

def test_Ecoli_core_etfba():
    
    import sys
    sys.path.append(ETFBA_PATH)
    from etfba import Model
    
    ecoli = Model('ecoli')
    ecoli.read_from_excel(ECOLI_MODEL_FILE1)
    #print(ecoli.end_metabolites)
    #print(ecoli.reactions['r73'].substrates['3pg'].km)
    #ecoli.stoichiometric_matrix.to_excel(OUTPUT_PATH+'\S.xlsx')
    #ecoli.total_stoichiometric_matrix.to_excel(OUTPUT_PATH+'\tS.xlsx')
    
    # res = ecoli.optimize('fba', 
    #                      objective = {'r94': 1}, 
    #                      preset_fluxes = {'r53': 10, 'r54': 0, 'r58': 0, 'r50': 8.39}, 
    #                      flux_bounds = (-100, 100),
    #                      irr_reactions = None,
    #                      excluded_mb = None).solve()
    # print(res)
    
    
    # res = ecoli.optimize('fba', 
    #                      objective = {'biom': 1}, 
    #                      preset_fluxes = {'glc_e': 10, 'fru_e': 0, 'atpm': 8.39}, 
    #                      flux_bounds = (-100, 100),
    #                      irr_reactions = None,
    #                      excluded_mb = None).solve()
    # print(res)
    
    res = ecoli.optimize('fba', 
                         objective = {'biom': 1}, 
                         preset_fluxes = {'glc_e': 10, 'fru_e': 0, 'o2_e': 0, 'atpm': 8.39}, 
                         flux_bounds = (-100, 100),
                         irr_reactions = None,
                         excluded_mb = None).solve()
    print(res)
    #'r50', 'r53', 'r54', 'r58'
    #'atpm', 'glc_e', 'fru_e', 'o2_e'
    
    
    # res = ecoli.optimize('tfba', 
    #                      objective = {'r94': 1}, 
    #                      preset_fluxes = {'r53': 10, 'r54': 0, 'r58': 0, 'r50': 8.39}, 
    #                      preset_concs = None,   # {'glc.e': 5.6} 1 g/L glucose
    #                      flux_bounds = (-100, 100), conc_bounds = (0.001, 10),
    #                      irr_reactions = None,
    #                      excluded_mb = None).solve()
    # print(res)
    # #print(res.opt_total_fluxes)
    # print('\nopt fluxes')
    # print(res.opt_fluxes)
    # #dump(res.opt_fluxes, open(OUTPUT_PATH+'\iniFluxes.bin', 'wb'))
    # print('\nopt concs')
    # print(res.opt_concentrations)
    # #dump(res.opt_concentrations, open(OUTPUT_PATH+'\iniConcs.bin', 'wb'))
    
    
    # res = ecoli.optimize('tfba', 
    #                      objective = {'biom': 1}, 
    #                      preset_fluxes = {'glcpts': 10, 'frupts': 0, 'o2_t': 0, 'atpm': 8.39}, 
    #                      preset_concs = None,   # {'glc.e': 5.6} 1 g/L glucose
    #                      flux_bounds = (-100, 100), conc_bounds = (0.001, 10),
    #                      irr_reactions = None,
    #                      excluded_mb = None).solve()
    # print(res)
    # #print(res.opt_directions)
    # #print(res.opt_gibbs_energy)
    
    
    # aerobic
    # res = ecoli.optimize('etfba', 
    #                      objective = {'r94': 1}, 
    #                      preset_fluxes = {'r53': 10, 'r54': 0},   #{'r53': 10, 'r54': 0, 'r50': 8.39}
    #                      preset_concs = None,   #{'glc.e': 5.6}
    #                      flux_bounds = (-100, 100), conc_bounds = (0.001, 10),
    #                      spec_flux_bounds = None,
    #                      excluded_mb = None, #excluded_thmd = ['r55', 'r57', 'r58'],
    #                      included_epcs = ['r01', 'r02', 'r03', 'r04', 'r05', 'r06', 'r07', 'r08', 'r09', 'r10', 'r11', 'r12', 'r13', 'r14', 'r15', 'r16', 'r17', 'r18', 'r19', 'r20', 'r21', 'r22', 'r23', 'r24', 'r25', 'r26', 'r27', 'r28', 'r29', 'r30', 'r31', 'r32', 'r33', 'r34', 'r35', 'r36', 'r37', 'r38', 'r39', 'r40', 'r41', 'r42', 'r43', 'r44', 'r45'],
    #                      enzyme_protein_ub = 0.1, 
    #                      use_initial_fluxes = None, #ECOLI_AERO_INI_FLUXES_FILE2, # OUTPUT_PATH+r'\optFluxes.bin'
    #                      use_initial_concs = None, #ECOLI_AERO_INI_CONCS_FILE, # OUTPUT_PATH+r'\optConcs.bin',
    #                      ).solve()
    
    
    # anaerobic
    # res = ecoli.optimize('etfba', 
    #                      objective = {'r94': 1}, 
    #                      preset_fluxes = {'r53': 10, 'r54': 0, 'r58': 0, 'r50': 8.39},
    #                      preset_concs = None,   #{'glc.e': 5.6}
    #                      flux_bounds = (-100, 100), conc_bounds = (0.001, 10),
    #                      excluded_mb = None, #excluded_thmd = ['r55', 'r57', 'r58'],
    #                      included_epcs = ['r01', 'r02', 'r03', 'r04', 'r05', 'r06', 'r07', 'r08', 'r09', 'r10', 'r11', 'r12', 'r13', 'r14', 'r15', 'r16', 'r17', 'r18', 'r19', 'r20', 'r21', 'r22', 'r23', 'r24', 'r25', 'r26', 'r27', 'r28', 'r29', 'r30', 'r31', 'r32', 'r33', 'r34', 'r35', 'r36', 'r37', 'r38', 'r39', 'r40', 'r41', 'r42', 'r43', 'r44', 'r45'],
    #                      enzyme_protein_ub = 0.1, 
    #                      use_initial_fluxes = OUTPUT_PATH+r'\fluxes.xlsx', #ECOLI_ANA_INI_FLUXES_FILE2, # OUTPUT_PATH+r'\fluxes.xlsx'
    #                      use_initial_concs = OUTPUT_PATH+r'\concs.xlsx', #ECOLI_ANA_INI_CONCS_FILE, # OUTPUT_PATH+r'\concs.xlsx'
    #                      ).solve()
    

    # res = ecoli.optimize('etfba', 
    #                      objective = {'biom': 1}, 
    #                      preset_fluxes = {'glcpts': 10, 'frupts': 0, 'o2_t': 0, 'atpm': 8.39},   #{'glcpts': 10, 'frupts': 0, 'o2_t': 0, 'atpm': 8.39} 
    #                      preset_concs = None, #{'glc.e': 5.6}
    #                      flux_bounds = (-100, 100), conc_bounds = (0.001, 10),
    #                      excluded_mb = None,
    #                      included_epcs = ['pgi', 'pfk', 'fbp', 'fba', 'tpi', 'gap', 'pgk', 'gpm', 'eno', 'pyk', 'pps', 'pdh', 'zwf', 'pgl', 'gnd', 'rpi', 'rpe', 'tkt1', 'tal', 'tkt2', 'cs', 'acn1', 'acn2', 'icd', 'kgd', 'suc', 'sdh', 'fum', 'mdh', 'gs', 'gdh', 'gls', 'gogat', 'aldh', 'adh', 'pta', 'ak', 'ldh', 'pfl', 'icl', 'mals', 'me1', 'me2', 'ppc', 'ppck'],
    #                      enzyme_protein_ub = 0.1, 
    #                      use_initial_fluxes = ECOLI_ANA_INI_FLUXES_FILE1, #OUTPUT_PATH+r'\fluxes.xlsx', #ECOLI_ANA_INI_FLUXES_FILE1,
    #                      use_initial_concs = ECOLI_ANA_INI_CONCS_FILE, #OUTPUT_PATH+r'\concs.xlsx' #ECOLI_ANA_INI_CONCS_FILE,
    #                      ).solve()

    
    # res = ecoli.optimize('epcm', 
    #                      objective = {'r94': 1}, 
    #                      preset_fluxes = {'r53': 10, 'r54': 0, 'r50': 8.39}, 
    #                      preset_concs = None,   #{'glc.e': 5.6}
    #                      flux_bounds = (-100, 100), conc_bounds = (0.001, 10),
    #                      excluded_mb = None, #excluded_thmd = ['r55', 'r57', 'r58'],
    #                      included_epcs = ['r01', 'r02', 'r03', 'r04', 'r05', 'r06', 'r07', 'r08', 'r09', 'r10', 'r11', 'r12', 'r13', 'r14', 'r15', 'r16', 'r17', 'r18', 'r19', 'r20', 'r21', 'r22', 'r23', 'r24', 'r25', 'r26', 'r27', 'r28', 'r29', 'r30', 'r31', 'r32', 'r33', 'r34', 'r35', 'r36', 'r37', 'r38', 'r39', 'r40', 'r41', 'r42', 'r43', 'r44', 'r45'],
    #                      use_fixed_fluxes = 'fba').solve()


    # res = ecoli.optimize('epcm', 
    #                      objective = {'biom': 1}, 
    #                      preset_fluxes = {'glcpts': 10, 'frupts': 0, 'o2_t': 0, 'atpm': 8.39}, 
    #                      preset_concs = None,   #{'glc.e': 5.6}
    #                      flux_bounds = (-100, 100), conc_bounds = (0.001, 10),
    #                      excluded_mb = None, #excluded_thmd = ['r55', 'r57', 'r58'],
    #                      included_epcs = ['pgi', 'pfk', 'fbp', 'fba', 'tpi', 'gap', 'pgk', 'gpm', 'eno', 'pyk', 'pps', 'pdh', 'zwf', 'pgl', 'gnd', 'rpi', 'rpe', 'tkt1', 'tal', 'tkt2', 'cs', 'acn1', 'acn2', 'icd', 'kgd', 'suc', 'sdh', 'fum', 'mdh', 'gs', 'gdh', 'gls', 'gogat', 'aldh', 'adh', 'pta', 'ak', 'ldh', 'pfl', 'icl', 'mals', 'me1', 'me2', 'ppc', 'ppck'],
    #                      use_fixed_fluxes = 'tfba').solve()
    
    
    print('optimization successful: ' + str(res.optimization_successful))
    #print(res)
    
    #print(res.opt_objective)
    print('\nopt fluxes')
    print(res.opt_fluxes)
    # if res.optimization_successful:
    #     res.opt_fluxes.save(OUTPUT_PATH+r'\fluxes.xlsx')
    
    print('\nopt concs')
    print(res.opt_concentrations)
    # if res.optimization_successful:
    #     res.opt_concentrations.save(OUTPUT_PATH+r'\concs.xlsx')
    
    #print('\nopt gibbs energy')
    #print(res.opt_gibbs_energy)
    
    print('\nopt enzyme cost')
    print(res.opt_total_enzyme_cost)
    print(res.opt_enzyme_costs)              
    

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
    

def test_Ecoli_core_etfba_ac():

    import sys
    sys.path.append(ETFBA_PATH)
    from etfba import Model
    
    ecoli = Model('ecoli')
    ecoli.read_from_excel(ECOLI_MODEL_SPLIT_ATPS4R_FILE1)

    #GIBBSENERGY_FILE = r'C:\Users\cwu\Desktop\Software\Papers\pH_effect\ph_effect_ATP_consumption_new\phin7.5_phout5.4-5.57_10\dgpms.xlsx'
    #GIBBSENERGY_FILE = r'C:\Users\cwu\Desktop\Software\Papers\pH_effect\ph_effect_ATP_consumption_new\phin7.5_phout5-9_50\dgpms.xlsx'
    GIBBSENERGY_FILE = r'C:\Users\cwu\Desktop\Software\Papers\pH_effect\ph_effect_ATP_consumption\phin7.5_phout3-8_10\dgpms.xlsx'
    
    dgpmsInfo = pd.read_excel(GIBBSENERGY_FILE, header = 0, index_col = 0)
    for i, row in dgpmsInfo.iterrows():
        if i == 8:
            (phin, phout), dgpms = row[:2], row[2:]
            for rxnid, dgpm in dgpms.items():
                ecoli.reactions[rxnid].dgpm = dgpm

            # res = ecoli.optimize('fba', 
            #                      objective = {'biom': 1}, 
            #                      preset_fluxes = {'frupts': 0, 'o2_t': 0, 'atpm': 8.39},   #'o2_t': 0
            #                      spec_flux_bounds = {'glcpts': (0, 10)},
            #                      flux_bounds = (-100, 100),
            #                      irr_reactions = None,
            #                      excluded_mb = None).solve()

            # res = ecoli.optimize('tfba', 
            #             objective = {'biom': 1}, 
            #             preset_fluxes = {'frupts': 0, 'o2_t': 0, 'atpm': 8.39},   #'o2_t': 0
            #             preset_concs = {'h': 10**(-phin), 'h.e': 10**(-phout)},   #{'ac.e': 60}
            #             spec_flux_bounds = {'glcpts': (0, 10)},
            #             flux_bounds = (-100, 100), 
            #             conc_bounds = (0.001, 10),   #(0.001, 10)
            #             ).solve()

            res = ecoli.optimize('etfba', 
                        objective = {'biom': 1}, 
                        preset_fluxes = {'frupts': 0, 'o2_t': 0},   #'o2_t': 0, 'atpm': 8.39
                        preset_concs = None, #{'h': 10**(-phin+3), 'h.e': 10**(-phout+3)},   # conc unit in mM
                        spec_flux_bounds = None, #{'glcpts': (0, 10)},
                        flux_bounds = (-100, 100), 
                        conc_bounds = (0.001, 10),   #(0.001, 10)
                        included_epcs = ['pgi', 'pfk', 'fbp', 'fba', 'tpi', 'gap', 'pgk', 'gpm', 'eno', 'pyk', 'pps', 'pdh', 'zwf', 'pgl', 'gnd', 'rpi', 'rpe', 'tkt1', 'tal', 'tkt2', 'cs', 'acn1', 'acn2', 'icd', 'kgd', 'suc', 'sdh', 'fum', 'mdh', 'gs', 'gdh', 'gls', 'gogat', 'aldh', 'adh', 'pta', 'ak', 'ldh', 'pfl', 'icl', 'mals', 'me1', 'me2', 'ppc', 'ppck', 'nadh16', 'cytbd', 'atps4r', 'atpase', 'adk', 'atpm', 'nnt1', 'nnt2', 'glcpts', 'frupts', 'pi_t', 'co2_t', 'o2_t', 'nh4_t', 'acald_t', 'ac_t', 'akg_t', 'lac_t', 'etoh_t', 'for_t1', 'for_t2', 'fum_t', 'mal_t', 'pyr_t', 'suc_t1', 'suc_t2', 'gln_t', 'glu_t'],
                        enzyme_protein_ub = 0.03,
                        use_initial_fluxes = None, #OUTPUT_PATH+r'\fluxes.xlsx', #ECOLI_ANA_INI_FLUXES_FILE1,
                        use_initial_concs = None #OUTPUT_PATH+r'\concs.xlsx', #ECOLI_ANA_INI_CONCS_FILE
                        ).solve()

            #'nadh16', 'cytbd', 'atps4r', 'atpase', 'adk', 'atpm', 'nnt1', 'nnt2', 'glcpts', 'frupts', 'pi_t', 'co2_t', 'o2_t', 'nh4_t', 'acald_t', 'ac_t', 'akg_t', 'lac_t', 'etoh_t', 'for_t1', 'for_t2', 'fum_t', 'mal_t', 'pyr_t', 'suc_t1', 'suc_t2', 'gln_t', 'glu_t'
            # exclude 'h2o_t'    

            # res = ecoli.optimize('epcm', 
            #             objective = {'biom': 1}, 
            #             preset_fluxes = {'frupts': 0, 'o2_t': 0, 'atpm': 8.39}, 
            #             preset_concs = {'h': 10**(-phin), 'h.e': 10**(-phout)},   #{'ac.e': 60}
            #             spec_flux_bounds = {'glcpts': (0, 10)},
            #             flux_bounds = (-100, 100), 
            #             conc_bounds = (0.001, 10),   #(0.001, 10)
            #             included_epcs = ['pgi', 'pfk', 'fbp', 'fba', 'tpi', 'gap', 'pgk', 'gpm', 'eno', 'pyk', 'pps', 'pdh', 'zwf', 'pgl', 'gnd', 'rpi', 'rpe', 'tkt1', 'tal', 'tkt2', 'cs', 'acn1', 'acn2', 'icd', 'kgd', 'suc', 'sdh', 'fum', 'mdh', 'gs', 'gdh', 'gls', 'gogat', 'aldh', 'adh', 'pta', 'ak', 'ldh', 'pfl', 'icl', 'mals', 'me1', 'me2', 'ppc', 'ppck'],
            #             use_fixed_fluxes = 'etfba',
            #             enzyme_protein_ub = 0.015,
            #             #use_initial_fluxes = None, #OUTPUT_PATH+r'\fluxes.xlsx', #ECOLI_ANA_INI_FLUXES_FILE1,
            #             #use_initial_concs = None #OUTPUT_PATH+r'\concs.xlsx', #ECOLI_ANA_INI_CONCS_FILE
            #             ).solve()

            print('optimization successful: ' + str(res.optimization_successful))
            print('\nopt fluxes')
            print(res.opt_fluxes)
            # print('\nopt concs')
            # print(res.opt_concentrations)
            # print('\nopt gibbs energy')
            # print(res.opt_gibbs_energy)
            print('\nopt enzyme cost')
            print(res.opt_total_enzyme_cost)
            print(res.opt_enzyme_costs)
            
            # if res.optimization_successful:
            #     res.opt_fluxes.save(OUTPUT_PATH+r'\fluxes.xlsx')
    
            # if res.optimization_successful:
            #     res.opt_concentrations.save(OUTPUT_PATH+r'\concs.xlsx')


def test_Ecoli_core_new():

    import sys
    sys.path.append(ETFBA_PATH)
    from etfba import Model
    
    ecoli = Model('ecoli')
    ecoli.read_from_excel(ECOLI_MODEL_FILE_ORI)

    print(ecoli.end_metabolites)
    
    #ecoli.transformation_matrix.to_excel(OUTPUT_PATH+'/T.xlsx')


    res = ecoli.optimize('fba', 
                         objective = {'biom': 1}, 
                         preset_fluxes = {'glc_e': 10, 'fru_e': 0, 'atpm': 8.39}, 
                         flux_bounds = (0, 100),
                         irr_reactions = None,
                         excluded_mb = None).solve()
    print(res.opt_fluxes)
    # res.opt_fluxes.save(OUTPUT_PATH+'/fba_fluxes.xlsx')

    # res = ecoli.optimize('fba', 
    #                      objective = {'biom': 1}, 
    #                      preset_fluxes = {'glc_e': 10, 'fru_e': 0, 'o2_e_f': 0, 'o2_e_b': 0, 'atpm': 8.39}, 
    #                      flux_bounds = (0, 100),
    #                      irr_reactions = None,
    #                      excluded_mb = None).solve()
    # print(res)
    # #print(res.opt_fluxes)
    # res.opt_fluxes.save(OUTPUT_PATH+'/fba_ana_fluxes.xlsx')
    

    # res = ecoli.optimize('tfba', 
    #                      objective = {'biom': 1}, 
    #                      preset_fluxes = {'glcpts': 10, 'frupts': 0, 'atpm': 8.39}, 
    #                      preset_concs = None,   # {'glc.e': 5.6} 1 g/L glucose
    #                      flux_bounds = (0, 100), conc_bounds = (0.001, 10),
    #                      irr_reactions = None,
    #                      excluded_mb = None).solve()
    # #print(res)
    # print(res.opt_fluxes)
    #print(res.opt_directions)
    #print(res.opt_gibbs_energy)

    # res = ecoli.optimize('tfba', 
    #                      objective = {'biom': 1}, 
    #                      preset_fluxes = {'glcpts': 10, 'frupts': 0, 'o2_e_f': 0, 'o2_e_b': 0, 'atpm': 8.39}, 
    #                      preset_concs = None,   # {'glc.e': 5.6} 1 g/L glucose
    #                      flux_bounds = (0, 100), conc_bounds = (0.001, 10),
    #                      irr_reactions = None,
    #                      excluded_mb = None).solve()
    # print(res)
    #print(res.opt_fluxes)
    #print(res.opt_directions)
    #print(res.opt_gibbs_energy)

    
    # res = ecoli.optimize('etfba', 
    #                      objective = {'biom': 1}, 
    #                      preset_fluxes = {'glcpts': 10, 'frupts': 0, 'atpm': 8.39},   #{'glcpts': 10, 'frupts': 0, 'o2_t': 0, 'atpm': 8.39} 
    #                      preset_concs = None, #{'glc.e': 5.6}
    #                      flux_bounds = (0, 100), conc_bounds = (0.001, 10),
    #                      excluded_mb = None,
    #                      included_epcs = ['pgi', 'pfk', 'fbp', 'fba', 'tpi', 'gap', 'pgk', 'gpm', 'eno', 'pyk', 'pps', 'pdh', 'zwf', 'pgl', 'gnd', 'rpi', 'rpe', 'tkt1', 'tal', 'tkt2', 'cs', 'acn1', 'acn2', 'icd', 'kgd', 'suc', 'sdh', 'fum', 'mdh', 'gs', 'gdh', 'gls', 'gogat', 'aldh', 'adh', 'pta', 'ak', 'ldh', 'pfl', 'icl', 'mals', 'me1', 'me2', 'ppc', 'ppck', 'nadh16', 'cytbd', 'atps4r', 'adk', 'atpm', 'nnt1', 'nnt2', 'glcpts', 'frupts', 'pi_t', 'co2_t', 'o2_t', 'nh4_t', 'acald_t', 'ac_t', 'akg_t', 'lac_t', 'etoh_t', 'for_t1', 'for_t2', 'fum_t', 'mal_t', 'pyr_t', 'suc_t1', 'suc_t2', 'gln_t', 'glu_t'],
    #                      enzyme_protein_ub = 0.1
    #                      ).solve()
    # print('optimization successful: ' + str(res.optimization_successful))
    # print(res.opt_fluxes)
    # print('\nopt enzyme cost')
    # print(res.opt_total_enzyme_cost)
    # print(res.opt_enzyme_costs)
    
    # res = ecoli.optimize('etfba', 
    #                      objective = {'biom': 1}, 
    #                      preset_fluxes = {'glcpts': 10, 'frupts': 0, 'o2_e_f': 0, 'o2_e_b': 0, 'atpm': 8.39},   #{'glcpts': 10, 'frupts': 0, 'o2_t': 0, 'atpm': 8.39} 
    #                      preset_concs = None, #{'glc.e': 5.6}
    #                      flux_bounds = (0, 100), conc_bounds = (0.001, 10),
    #                      excluded_mb = None,
    #                      included_epcs = ['pgi', 'pfk', 'fbp', 'fba', 'tpi', 'gap', 'pgk', 'gpm', 'eno', 'pyk', 'pps', 'pdh', 'zwf', 'pgl', 'gnd', 'rpi', 'rpe', 'tkt1', 'tal', 'tkt2', 'cs', 'acn1', 'acn2', 'icd', 'kgd', 'suc', 'sdh', 'fum', 'mdh', 'gs', 'gdh', 'gls', 'gogat', 'aldh', 'adh', 'pta', 'ak', 'ldh', 'pfl', 'icl', 'mals', 'me1', 'me2', 'ppc', 'ppck', 'nadh16', 'cytbd', 'atps4r', 'adk', 'atpm', 'nnt1', 'nnt2', 'glcpts', 'frupts', 'pi_t', 'co2_t', 'o2_t', 'nh4_t', 'acald_t', 'ac_t', 'akg_t', 'lac_t', 'etoh_t', 'for_t1', 'for_t2', 'fum_t', 'mal_t', 'pyr_t', 'suc_t1', 'suc_t2', 'gln_t', 'glu_t'],
    #                      enzyme_protein_ub = 0.04, 
    #                      ).solve()
    # print('optimization successful: ' + str(res.optimization_successful))
    # print(res.opt_fluxes)
    # res.opt_fluxes.save(OUTPUT_PATH+'/etfba_fluxes.xlsx')
    # #print(res.opt_directions)
    # # print('\nopt Gibbs energy')
    # # print(res.opt_gibbs_energy)
    # print('\nopt enzyme cost')
    # print(res.opt_total_enzyme_cost)
    # print(res.opt_enzyme_costs)


    # res = ecoli.optimize('epcm', 
    #                      preset_concs = None, #{'glc.e': 5.6}
    #                      conc_bounds = (0.001, 10),
    #                      excluded_thmd = ['nnt1', 'nnt2', 'pi_t', 'co2_t', 'o2_t', 'nh4_t', 'acald_t', 'ac_t', 'akg_t', 'lac_t', 'etoh_t', 'for_t1', 'for_t2', 'fum_t', 'mal_t', 'pyr_t', 'suc_t1', 'suc_t2', 'gln_t', 'glu_t'],   #! transport reactions and transhydrogenase excluded
    #                      included_epcs = ['pgi', 'pfk', 'fbp', 'fba', 'tpi', 'gap', 'pgk', 'gpm', 'eno', 'pyk', 'pps', 'pdh', 'zwf', 'pgl', 'gnd', 'rpi', 'rpe', 'tkt1', 'tal', 'tkt2', 'cs', 'acn1', 'acn2', 'icd', 'kgd', 'suc', 'sdh', 'fum', 'mdh', 'gs', 'gdh', 'gls', 'gogat', 'aldh', 'adh', 'pta', 'ak', 'ldh', 'pfl', 'icl', 'mals', 'me1', 'me2', 'ppc', 'ppck', 'nadh16', 'cytbd', 'atps4r', 'adk', 'atpm', 'glcpts', 'frupts'],   
    #                      use_fixed_fluxes = FIXED_FLUXES_FILE,
    #                      ).solve()
    # print('optimization successful: ' + str(res.optimization_successful))
    # # print('\nopt fluxes')
    # # print(res.opt_fluxes)
    # # print('\nopt concs')
    # # print(res.opt_concentrations)
    # # print('\nopt Gibbs energy')
    # # print(res.opt_gibbs_energy)
    # print('\nopt enzyme cost')
    # print(res.opt_total_enzyme_cost)
    # print(res.opt_enzyme_costs)
    # # print('\nopt enzyme base cost')
    # # print(res.opt_enzyme_base_costs)
    # # print('\nopt enzyme thermodynamic cost')
    # # print(res.opt_enzyme_thermodynamic_costs)
    # # print('\nopt enzyme kinetic cost')
    # # print(res.opt_enzyme_kinetic_costs)


def test_Ecoli_core_epcm():

    ECOLI_MODEL_FILE = r'C:\Users\cwu\Desktop\Software\ETFBA\Data\ecoli_core\E_coli_core.xlsx'
    #r'C:\Users\cwu\Desktop\Software\Papers\pH_effect\optimal_ph3\growth_simulation\E_coli_core_atpase.xlsx'
    FIXED_FLUXES_FILE = r'C:\Users\cwu\Desktop\Software\Papers\pH_effect\optimal_ph3\growth_simulation\etfba_kcat_adjusted\phout3-8\Q=0.04\0_fluxes.xlsx'
    #r'C:\Users\cwu\Desktop\Software\Papers\pH_effect\optimal_ph3\growth_simulation\etfba\phout3-8\Q=0.04\0_fluxes.xlsx'
    GIBBSENERGY_FILE = r'C:\Users\cwu\Desktop\Software\Papers\pH_effect\optimal_ph3\dgpm_estimation\dgpms_phout3-8_20_atpase.xlsx'

    import sys
    sys.path.append(ETFBA_PATH)
    from etfba import Model
    
    ecoli = Model('ecoli')
    ecoli.read_from_excel(ECOLI_MODEL_FILE)

    dgpmsInfo = pd.read_excel(GIBBSENERGY_FILE, header = 0, index_col = 0)
    for i, row in dgpmsInfo.iterrows():
        if i == 9:
            (phin, phout), dgpms = row[:2], row[2:]
            for rxnid, dgpm in dgpms.items():
                if rxnid in ecoli.reactions:
                    ecoli.reactions[rxnid].dgpm = dgpm

    res = ecoli.optimize('epcm', 
                         preset_concs = None, #{'glc.e': 5.6}
                         conc_bounds = (0.0005, 20), #(0.001, 10)
                         excluded_thmd = ['nnt1', 'nnt2', 'pi_t', 'co2_t', 'o2_t', 'nh4_t', 'acald_t', 'ac_t', 'akg_t', 'lac_t', 'etoh_t', 'for_t1', 'for_t2', 'fum_t', 'mal_t', 'pyr_t', 'suc_t1', 'suc_t2', 'gln_t', 'glu_t'],   #! transport reactions and transhydrogenase excluded
                         included_epcs = ['pgi', 'pfk', 'fbp', 'fba', 'tpi', 'gap', 'pgk', 'gpm', 'eno', 'pyk', 'pps', 'pdh', 'zwf', 'pgl', 'gnd', 'rpi', 'rpe', 'tkt1', 'tal', 'tkt2', 'cs', 'acn1', 'acn2', 'icd', 'kgd', 'suc', 'sdh', 'fum', 'mdh', 'gs', 'gdh', 'gls', 'gogat', 'aldh', 'adh', 'pta', 'ak', 'ldh', 'pfl', 'icl', 'mals', 'me1', 'me2', 'ppc', 'ppck', 'nadh16', 'cytbd', 'adk', 'atpm', 'glcpts', 'frupts'], # 'atps4r',   
                         use_fixed_fluxes = FIXED_FLUXES_FILE,
                         ).solve()
    print('optimization successful:', res.optimization_successful)
    print('\nopt enzyme cost')
    print(res.opt_total_enzyme_cost)
    print(res.opt_enzyme_costs)


def test_io():

    import sys
    sys.path.append(ETFBA_PATH)
    from etfba import Model
    
    ecoli = Model('ecoli')
    ecoli.read_from_excel(ECOLI_MODEL_FILE_ORI)

    ecoli.save(r'C:\Users\cwu\Desktop\ecoli_core.bin')

    ecoli = Model.load(r'C:\Users\cwu\Desktop\ecoli_core.bin')    
    print(ecoli.metabolites)


def test_build_model():

    import sys
    sys.path.append(ETFBA_PATH)
    from etfba import Model, Metabolite, Reaction
    
    ecoli = Model('ecoli')

    # altrn_c --> 2ddglcn_c + h2o_c
    altrn_c = Metabolite('altrn_c', 'altrn_c')
    # altrn_c2 = Metabolite('altrn_c', 'altrn_c')
    # altrn_c3 = Metabolite('altrn_c', 'altrn_c3')
    # altrn_c4 = Metabolite('altrn_c', 'altrn_c', is_biomass = True)
    # print(altrn_c is altrn_c2, altrn_c is altrn_c3, altrn_c is altrn_c4)
    
    ddglcn_c = Metabolite('2ddglcn_c', '2ddglcn_c')
    h2o_c = Metabolite('h2o_c', 'h2o_c', is_h2o = True)

    ALTRH = Reaction('ALTRH', 'Altronate hydrolase',
                     forward_kcat = 1000, backward_kcat = 1000, 
                     molecular_weight = 30, standard_gibbs_energy = -30, 
                     reversible = False)
    ALTRH.add_substrates(coes = {altrn_c: 1, ddglcn_c: 1},
                         kms = {altrn_c: 23})
    ALTRH.add_products(coes = {h2o_c: 1},
                       kms = {})

    # print(ALTRH.substrates['altrn_c'].km)
    # print(ALTRH.substrates['2ddglcn_c'].km)
    # print(ALTRH.products['h2o_c'].km)

    # ALTRH.remove_substrates([altrn_c])
    # print(ALTRH)
    # print(ALTRH.substrates['altrn_c'].km)
    
    ecoli.add_reactions([ALTRH])
    ecoli.remove_reactions([ALTRH])

    print(ecoli)
    



if __name__ == '__main__':
    
    #test_TCA()
    #test_TCA_etfba()
    #test_Ecoli_core_etfba()
    #test_Ecoli_core_etfba_ac()
    
    test_Ecoli_core_new()

    # test_Ecoli_core_epcm()

    #test_io()

    #test_build_model()