import pandas as pd
import numpy as np

def codebook_syn(data, maxlevs=10):
    if not isinstance(data, pd.DataFrame):
        raise ValueError("codebook.syn() requires a data frame as a parameter.")
    
    n, p = data.shape
    
    # Calculate number and % of missing and non-missing values
    nmiss = data.isna().sum()
    perctmiss = (nmiss / n) * 100
    nok = data.notna().sum()
    ndistinct = data.nunique()
    dfclass = data.dtypes
    
    details = [''] * len(nmiss)
    fortab2 = [''] * len(nmiss)
    
    for i in range(p):
        if data.iloc[:, i].dtype == 'object':  # Character type
            details[i] = f"Max length: {data.iloc[:, i].str.len().max()}"
        elif np.issubdtype(data.iloc[:, i].dtype, np.number):  # Numeric type
            min_val = data.iloc[:, i].min()
            max_val = data.iloc[:, i].max()
            details[i] = f"Range: {min_val} - {max_val}"
        elif pd.api.types.is_categorical_dtype(data.iloc[:, i]) and ndistinct[i] > maxlevs:  # Factor with too many levels
            details[i] = "See table in labs"
            fortab2[i] = " ".join([f"'{x}'" for x in data.iloc[:, i].cat.categories])
        elif pd.api.types.is_categorical_dtype(data.iloc[:, i]) and ndistinct[i] <= maxlevs:  # Factor with few levels
            details[i] = " ".join([f"'{x}'" for x in data.iloc[:, i].cat.categories])
    
    # Create the table of detailed factor levels if needed
    tabs2 = None
    if any(pd.api.types.is_categorical_dtype(data.iloc[:, i]) and ndistinct[i] > maxlevs for i in range(p)):
        vnum = [i for i in range(p) if pd.api.types.is_categorical_dtype(data.iloc[:, i]) and ndistinct[i] > maxlevs]
        tabs2 = {}
        for i in vnum:
            tabs2[data.columns[i]] = pd.DataFrame({'label': data.iloc[:, i].cat.categories})
    
    # Create the result dataframe
    result = pd.DataFrame({
        'variable': data.columns,
        'class': dfclass.astype(str),
        'nmiss': nmiss,
        'perctmiss': perctmiss,
        'ndistinct': ndistinct,
        'details': details
    })
    
    result.index = range(1, p + 1)
    
    return {'tab': result, 'labs': tabs2}

# Example usage (assuming you have a pandas DataFrame `data`)
# result = codebook_syn(data)
