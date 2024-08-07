import numpy as np
from scipy.stats.mstats import gmean
from etfba import Reaction


def get_internal_ph(peri_ph):
    return 4.4 + 0.38*peri_ph
    

def H_leak_flux(cyto_ph, peri_ph):
    coes = [16.35021215, -7.10765306, 1.07854944, -0.05434789]
    g = coes@np.array([np.power(peri_ph, degree) for degree in range(4)])*3.6
    return g*(cyto_ph - peri_ph)


def pH2H(ph):
    return 10**(-ph)


def activity(ph, c, KaC, KdC):
    return c/(1 + pH2H(ph)/KaC + KdC/pH2H(ph))


def relative_activity(ph, KaC, KdC):
    return (1+2*(KdC/KaC)**0.5)/(1 + pH2H(ph)/KaC + KdC/pH2H(ph))


def add_H_leak_rxn(model):
    h_p = model.metabolites['h_p']
    h_c = model.metabolites['h_c']
    h_leak = Reaction(
        'H_leak', 
        'H_leak', 
        forward_kcat=None, 
        backward_kcat=None,
        molecular_weight=None, 
        standard_gibbs_energy=None, 
        reversible=True,
        is_biomass_formation=False, 
        is_exch_reaction=False,
        is_h_transport=True, 
        is_h2o_transport=False
    )
    h_leak.add_substrates(coes={h_p: 1})
    h_leak.add_products(coes={h_c: 1})
    model.add_reactions([h_leak])


def set_kcat_and_MW(model, kcat_data, mw_data):
    included_epcs = []
    for rxnid, rxn in model.reactions.items():
        if rxn.rev:
            if rxnid+'_f' in kcat_data.index and rxnid+'_b' in kcat_data.index:
                rxn.forward_kcat = kcat_data[rxnid+'_f']   
                rxn.backward_kcat = kcat_data[rxnid+'_b']   
                included_epcs.append(rxnid)
        else:
            if rxnid in kcat_data.index:
                rxn.forward_kcat = kcat_data[rxnid]   
                rxn.backward_kcat = None
                included_epcs.append(rxnid)
        if rxnid in mw_data.index:
            rxn.molecular_weight = mw_data[rxnid]

    return included_epcs


def get_default_KaC_and_KdC(kd_data):
    kd_data = kd_data[~(kd_data==0).any(axis=1)]   
    KaC_default = gmean(kd_data['KaC'])
    KdC_default = gmean(kd_data['KdC'])

    return KaC_default, KdC_default


def set_deltaGprimem(model, dgpm_data):
    for rxnid, dgpm in dgpm_data.items():
        if rxnid in model.reactions:
            model.reactions[rxnid].standard_gibbs_energy = dgpm


def get_H_leak_flux(cyto_ph, peri_ph):
    h_leak_flux = {}
    h_t = H_leak_flux(cyto_ph, peri_ph)
    h_leak_flux['H_leak_f'] = max(h_t, 0)
    h_leak_flux['H_leak_b'] = max(-h_t, 0)

    return h_leak_flux


def set_adjusted_kcat(model, kd_data, included_epcs, sep_rxns, cyto_ph, 
                      KaC_default, KdC_default):
    old_kcats = {}
    for rxnid in included_epcs:
        if rxnid in sep_rxns:
            KaC, KdC = kd_data.loc[rxnid+'_f', ['KaC', 'KdC']]
            corr_coe = relative_activity(cyto_ph, KaC, KdC)

            fkcat = model.reactions[rxnid].forward_kcat
            model.reactions[rxnid].forward_kcat = fkcat*corr_coe
            old_kcats[rxnid] = [fkcat]

            if model.reactions[rxnid].rev:
                KaC, KdC = kd_data.loc[rxnid+'_b', ['KaC', 'KdC']]   
                corr_coe = relative_activity(cyto_ph, KaC, KdC)
                bkcat = model.reactions[rxnid].backward_kcat
                model.reactions[rxnid].backward_kcat = bkcat*corr_coe
                old_kcats[rxnid].append(bkcat)
        else:
            if rxnid in kd_data.index:
                KaC, KdC = kd_data.loc[rxnid, ['KaC', 'KdC']]   
                corr_coe = relative_activity(cyto_ph, KaC, KdC)
            else:
                corr_coe = relative_activity(cyto_ph, KaC_default, KdC_default)

            fkcat = model.reactions[rxnid].forward_kcat
            model.reactions[rxnid].forward_kcat = fkcat*corr_coe
            old_kcats[rxnid] = [fkcat]
            
            if model.reactions[rxnid].rev:
                bkcat = model.reactions[rxnid].backward_kcat
                model.reactions[rxnid].backward_kcat = bkcat*corr_coe
                old_kcats[rxnid].append(bkcat)

    return old_kcats


def reset_kcat(model, old_kcats, included_epcs):
    for rxnid in included_epcs:
        model.reactions[rxnid].forward_kcat = old_kcats[rxnid][0]
        if model.reactions[rxnid].reversible:
            model.reactions[rxnid].backward_kcat = old_kcats[rxnid][1]


def set_cpu_affinity():
    import platform
    
    if platform.system() == 'Linux':
        import os
        os.sched_setaffinity(os.getpid(), range(os.cpu_count()))
