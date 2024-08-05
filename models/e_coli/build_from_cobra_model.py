'''Build an ETFBA model from a COBRA model (tested with COBRApy==0.21.0).
Default values are set for parameters such as kcat (catalytic rate constant), 
molecular weight (MW), and standard reaction Gibbs energy.
'''


import os
import re
from cobra.io import load_json_model
from etfba import Model, Metabolite, Reaction


COBRA_MODEL_FILE = './iML1515.json'

DEFAULT_MW = 40       # default enzyme molecular weight in kDa
DEFAULT_KCAT = 200    # default reaction catalytic rate constant in 1/s
DEFAULT_DGPM = 0      # default standard reaction Gibbs energy in KJ/mol


def main():
    cobra_model = load_json_model(COBRA_MODEL_FILE)

    filedir, filename = os.path.split(COBRA_MODEL_FILE)
    model_name = os.path.splitext(filename)[0]
    model = Model(model_name)
    for cobra_rxn in cobra_model.reactions:
        if re.search(r'biomass', cobra_rxn.name, flags = re.I):
            is_biomass_formation = True
        else:
            is_biomass_formation = False

        if re.search(r'exchange', cobra_rxn.subsystem, flags = re.I):
            is_exch_reaction = True
        else:
            is_exch_reaction = False

        if re.search(r'proton transport', cobra_rxn.name, flags = re.I):
            is_h_transport = True
        else:
            is_h_transport = False

        if re.search(r'H2O transport', cobra_rxn.name, flags = re.I):
            is_h2o_transport = True
        else:
            is_h2o_transport = False    

        if cobra_rxn.lower_bound >= 0 or cobra_rxn.upper_bound <= 0:
            if cobra_rxn.lower_bound == cobra_rxn.upper_bound == 0:
                print('\nzero flux reaction:', cobra_rxn.id)
            rev = False
        else:
            rev = True

        rxn = Reaction(
            cobra_rxn.id, 
            cobra_rxn.name, 
            cobra_rxn.subsystem, 
            forward_kcat = DEFAULT_KCAT, 
            backward_kcat = DEFAULT_KCAT if rev else None,
            molecular_weight = DEFAULT_MW, 
            standard_gibbs_energy = DEFAULT_DGPM, 
            reversible = rev,
            is_biomass_formation = is_biomass_formation, 
            is_exch_reaction = is_exch_reaction,
            is_h_transport = is_h_transport, 
            is_h2o_transport = is_h2o_transport
        )

        subs_cores = {}
        pros_coes = {}
        for cobra_metab, coe in cobra_rxn.metabolites.items():
            if cobra_metab.formula == 'H':
                is_h = True
            else:
                is_h = False

            if cobra_metab.formula == 'H2O':
                is_h2o = True
            else:
                is_h2o = False
            
            reac = Metabolite(
                cobra_metab.id, 
                cobra_metab.name, 
                cobra_metab.compartment, 
                is_h = is_h, 
                is_h2o = is_h2o
            )    
            if coe < 0:
                subs_cores[reac] = -coe
            elif coe > 0:
                pros_coes[reac] = coe

            rxn.add_substrates(coes = subs_cores)
            rxn.add_products(coes = pros_coes)

        model.add_reactions([rxn])

    print(model)
    model.save(f'{filedir}/etfba_{model_name}.bin')




if __name__ == '__main__':

    main()
