import pandas as pd
import numpy as np

# Replicated Unique Units
def replicated_uniques(object, data, exclude=None):
    # Exclude specified columns
    if exclude:
        data = data.drop(columns=exclude)

    # Extract unique units from the real data
    uReal = data.drop_duplicates()

    no_uniques = len(uReal)

    if no_uniques == 0:
        no_duplicated = no_per_duplicated = np.zeros(object['m'])
        rm_syn = np.zeros(len(object['syn']), dtype=bool) if object['m'] == 1 else np.zeros((len(object['syn'][0]), object['m']), dtype=bool)
    else:
        if object['m'] == 1:
            Syn = object['syn'] if not exclude else object['syn'].drop(columns=exclude)
            rm_syn = np.zeros(len(Syn), dtype=bool)
            i_unique_Syn = Syn[~Syn.duplicated() & ~Syn.duplicated(keep='last')]
            uSyn = Syn.loc[i_unique_Syn]
            uAll = pd.concat([uReal, uSyn])
            dup_of_unique = uAll.duplicated().iloc[len(uReal):]
            rm_syn[i_unique_Syn.index] = dup_of_unique
            no_duplicated = rm_syn.sum()
        else:
            rm_syn = np.zeros((len(object['syn'][0]), object['m']), dtype=bool)
            for i in range(object['m']):
                Syn = object['syn'][i] if not exclude else object['syn'][i].drop(columns=exclude)
                i_unique_Syn = Syn[~Syn.duplicated() & ~Syn.duplicated(keep='last')]
                uSyn = Syn.loc[i_unique_Syn]
                uAll = pd.concat([uReal, uSyn])
                dup_of_unique = uAll.duplicated().iloc[len(uReal):]
                rm_syn[i_unique_Syn.index, i] = dup_of_unique
            no_duplicated = rm_syn.sum(axis=0)
        
        per_duplicated = no_duplicated / len(data) * 100

    return {
        'replications': rm_syn,
        'no_uniques': no_uniques,
        'no_replications': no_duplicated,
        'per_replications': per_duplicated
    }

# Statistical Disclosure Control (sdc)
def sdc(object, data, label=None, rm_replicated_uniques=False, uniques_exclude=None, 
        recode_vars=None, bottom_top_coding=None, recode_exclude=None, smooth_vars=None):
    if smooth_vars:
        if object['m'] == 1:
            if not all(var in object['syn'].columns for var in smooth_vars):
                raise ValueError("Some of smooth.vars not in the data")
            if not all(np.issubdtype(object['syn'][var], np.number) for var in smooth_vars):
                raise ValueError("Some of smooth.vars not numeric")
        else:
            if not all(var in object['syn'][0].columns for var in smooth_vars):
                raise ValueError("Some of smooth.vars not in the data")
            if not all(np.issubdtype(object['syn'][0][var], np.number) for var in smooth_vars):
                raise ValueError("Some of smooth.vars not numeric")

    if recode_vars:
        if bottom_top_coding and not isinstance(bottom_top_coding, list):
            bottom_top_coding = [bottom_top_c
