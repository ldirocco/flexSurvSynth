import numpy as np
import pandas as pd
import random
from collections import Counter

def syn_strata(data, strata=None, minstratumsize=None, tab_strataobs=True, tab_stratasyn=False, 
               method="cart", visit_sequence=None, predictor_matrix=None, m=1, k=None, proper=False, 
               minnumlevels=1, maxfaclevels=60, rules=None, rvalues=None, cont_na=None, semicont=None, 
               smoothing=None, event=None, denom=None, drop_not_used=False, drop_pred_only=False, 
               default_method=["normrank", "logreg", "polyreg", "polr"], numtocat=None, catgroups=None, 
               models=False, print_flag=True, seed="sample", **kwargs):

    # Set the random seed
    if seed != "sample":
        np.random.seed(seed)
    else:
        seed = random.randint(1, int(1e9))

    # Validation checks
    if strata is None:
        raise ValueError("Argument strata is missing.")

    # Convert strata to a factor if it's provided as variable names
    if isinstance(strata, list) and all(isinstance(i, str) for i in strata):
        varindex = [data.columns.get_loc(i) for i in strata]
        if any(i is None for i in varindex):
            raise ValueError(f"Unrecognized variable(s) in strata: {', '.join([str(i) for i in strata])}")
        
        strata_data = [data.iloc[:, i] for i in varindex]
        strata_lab = pd.factorize(pd.concat(strata_data, axis=1).astype(str).agg("_".join, axis=1))[0]

    else:
        if len(strata) != len(data):
            raise ValueError(f"The length of strata index ({len(strata)}) does not match the number of rows in the data ({len(data)}).")
        
        if any(pd.isna(strata)):
            raise ValueError("Strata indicator cannot have missing values.")
        
        strata_lab = pd.factorize(strata)[0]

    strata_levels = pd.unique(strata_lab)
    strata_n_obs = Counter(strata_lab)
    
    if tab_strataobs:
        print("Number of observations in strata (original data):")
        print(strata_n_obs)

    # Check min number of observations in strata
    if minstratumsize is None:
        minstratumsize = 10 + 10 * len(visit_sequence)
    
    smallstrata = sum(obs < minstratumsize for obs in strata_n_obs.values())
    if smallstrata > 5:
        raise ValueError(f"Multiple strata have fewer than the recommended number of observations. "
                         f"We advise that each should have at least {minstratumsize} observations.")
    
    if smallstrata > 0:
        print(f"CAUTION: In the original data some strata ({', '.join([str(level) for level in strata_levels if strata_n_obs[level] < minstratumsize])}) "
              f"have fewer than the recommended number of observations.")

    # Synthetic data creation
    syn_data = []
    for j in range(m):
        strata_syn = [random.choices(strata_levels, weights=[strata_n_obs[level] for level in strata_levels], k=k) for _ in range(m)]
        syn_data.append(strata_syn)

        if tab_stratasyn:
            print(f"\nNumber of observations in strata (synthetic data, m = {j}):")
            print(Counter(strata_syn))

    for j in range(m):
        strata_ind = []
        for i, level in enumerate(strata_levels):
            if print_flag:
                print(f"\nm = {j}, strata = {level}")
                print(f"-----------------------------------------------------")

            stratum_data = data[strata_lab == level]
            k_stratum = strata_n_obs[level]
            if k_stratum == 0 or m == 0:
                continue
            else:
                strata_ind.append({'data': stratum_data, 'k': k_stratum})

        syn_data[j].append(strata_ind)

    return syn_data
