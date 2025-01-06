import numpy as np
import pandas as pd
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
# Assuming that all methods (like 'cart', 'rf', etc.) will be defined in Python functions
from sklearn.linear_model import LogisticRegression
import random

def remove_lindep_syn(x, y, eps=0.00001, maxcor=0.99999, allow_na=False):
    """
    Remove linearly dependent variables from the data matrix
    """
    if x.shape[1] == 0:
        return None

    if eps <= 0:
        raise ValueError("Argument 'eps' must be positive.")

    x_obs = np.array(x).astype(float)
    y_obs = np.array(y).astype(float)
    keep = np.var(x_obs, axis=0) > eps
    keep[np.isnan(keep)] = False
    keep = keep & (np.abs(np.corrcoef(x_obs.T, y_obs)[-1, :-1]) < maxcor)

    if not np.any(keep):
        print("\nAll predictors are constant or have too high correlation.\n")
    
    ksum = np.sum(keep)
    cx = np.corrcoef(x_obs[:, keep].T)
    eig_values, eig_vectors = np.linalg.eig(cx)

    while eig_values[ksum - 1] / eig_values[0] < eps:
        j = np.argsort(np.abs(eig_vectors[:, ksum - 1]))[0]
        keep[keep][j] = False
        cx = cx[keep, :][:, keep]
        ksum -= 1
        eig_values, eig_vectors = np.linalg.eig(cx)
    
    return keep


def sampler_syn(p, data, m, syn, visit_sequence,
                rules, rvalues, event, proper,
                print_flag, k, pred_not_syn,
                models, numtocat, **kwargs):
    """
    This function controls the generation of conditional distributions
    """

    # --- Assign optional parameters (from kwargs) to appropriate synthesizing function
    meth_with_opt = ["cart", "cartbboot", "ctree", "survctree", "polyreg", 
                     "norm", "lognorm", "sqrtnorm", "cubertnorm", "normrank", "pmm",
                     "polr", "rf", "ranger", "bag", "ipf", "catall"]

    meth_check = [key for key in kwargs.keys() if any(m in key for m in meth_with_opt)]
    args_err = [key for key in kwargs.keys() if key not in meth_check]
    
    if args_err:
        raise ValueError(f"Unknown optional parameter(s): {', '.join(args_err)}")

    mth_args = None
    if kwargs:
        mth_args_dots = {key.split('.')[0]: value for key, value in kwargs.items()}
        mth_args = {key: {subkey: value for subkey, value in kwargs.items() if key in subkey} for key in mth_args_dots}

    fits = None
    if m > 0:
        if models:
            fits = [dict.fromkeys(p["method"].keys()) for _ in range(m)]
        
        for i in range(m):  # Synthesizing loop
            if print_flag and m > 1:
                print(f"\nSynthesis number {i+1}\n--------------------")
            if print_flag and m == 1:
                print("\nSynthesis\n-----------")

            rest_visit_sequence = p["visit_sequence"]
            
            if "catall" in p["method"] or "ipf" in p["method"]:
                ordmethod = [p["method"][index] for index in p["visit_sequence"]]
                grind = [index for index, method in enumerate(ordmethod) if method == ordmethod[0]]

                # Reorder dummies for grouped variables
                if any(name in p["visit_sequence"] for name in [f"{grind}1"]):
                    dumind = [index for index in p["visit_sequence"] if f"{grind}1" in name]
                    othind = [index for index in p["visit_sequence"] if index not in grind and index not in dumind]
                    p["visit_sequence"] = grind + dumind + othind
                    ordmethod = [p["method"][index] for index in p["visit_sequence"]]
                
                grouped = p["visit_sequence"][ordmethod[0]]
                if print_flag:
                    if len(rest_visit_sequence) > 0 and (len(data.columns) - len(numtocat)) > len(grouped):
                        print(f"First {len(grouped)} variables ({', '.join(grouped)}) synthesised together by method '{ordmethod[0]}'")

                # Call method
                x = p["data"][grouped]
                method_func = globals().get(f"syn_{ordmethod[0]}")
                if method_func:
                    synfun = method_func(x=x, k=k, proper=proper, **(mth_args.get(ordmethod[0], {})))
                    p["syn"][grouped] = synfun["res"]

                if models:
                    fits[i][grouped[0]] = synfun["fit"]
                    for j in range(1, len(grouped)):
                        fits[i][grouped[j]] = f"See first in group: {grouped[0]}"
                
                rest_visit_sequence = p["visit_sequence"][len(grouped):]
                
            # Process remaining variables
            if len(rest_visit_sequence) > 0:
                prcount = 0
                for j in rest_visit_sequence:
                    the_method = p["method"][j]
                    fun_args = mth_args.get(the_method) if the_method in mth_args else None
                    
                    vname = data.columns[j]
                    if print_flag and the_method != "dummy" and j <= len(data.columns) - len(numtocat):
                        print(vname, end=" ")
                        prcount += 1

                    if prcount % 10 == 0 and j <= len(data.columns) - len(numtocat):
                        print()

                    ya = np.arange(len(p["data"]))  # available rows
                    ypa = np.arange(k)  # predicted rows
                    
                    # Apply rules
                    if p["rules"][j] != "":
                        com_rules = " | ".join(p["rules"][j])
                        eval_rul_y = eval(com_rules)
                        ym = np.where(eval_rul_y == True)[0]
                        ya = np.setdiff1d(np.arange(len(p["data"])), ym)
                        
                        eval_rul_yp = eval(com_rules)
                        ypm = np.where(eval_rul_yp == True)[0]
                        ypa = np.setdiff1d(np.arange(k), ypm)
                    
                    if the_method != "" and not is_passive(the_method) and the_method != "dummy":
                        if the_method in ["sample", "sample.proper", "constant"]:
                            y = p["data"][ya, j]
                            if isinstance(y, pd.Series):
                                y = y.values
                            xp = len(ypa)
                            f = f"syn_{the_method}"
                            synfun = globals().get(f)(y=y, xp=xp, proper=proper, **fun_args)
                            p["syn"][ypa, j] = synfun["res"]
                            if models:
                                fits[i][j] = synfun["fit"]
                        else:
                            # Other synthesizing functions here...
                            pass

                    if p["rules"][j] != "":
                        for r in range(len(p["rules"][j])):
                            reval_rul_yp = eval(p["rules"][j][r])
                            rypm = np.where(reval_rul_yp == True)[0]
                            if len(rypm) > 0:
                                p["syn"][rypm, j] = p["rvalues"][j][r]
            
            syn[i] = p["syn"][:, :len(data.columns)]
            nms = data.columns.tolist()
            if np.sum(pred_not_syn) > 0:
                syn[i] = syn[i].drop(columns=[nms[idx] for idx in range(len(pred_not_syn)) if pred_not_syn[idx]])
                nms = [nms[idx] for idx in range(len(pred_not_syn)) if not pred_not_syn[idx]]

            # Prevent a single character column from being changed to a factor
            chgetochar = (np.sum(~pred_not_syn) == 1 and isinstance(syn[i].iloc[:, 0], str))
            syn[i] = pd.DataFrame(syn[i])
            if chgetochar:
                syn[i].iloc[:, 0] = syn[i].iloc[:, 0].astype(str)
                syn[i].columns = nms

            # Handling NA levels in factors or logicals
            for j in range(syn[i].shape[1]):
                if pd.api.types.is_categorical_dtype(syn[i].iloc[:, j]):
                    levels = syn[i].iloc[:, j].cat.categories
                    if "NAlogical" in levels:
                        syn[i].iloc[:, j].cat.categories = levels.str.replace("NAlogical", "NA")
                    else:
                        levels = levels.str.replace("NAtemp", "NA")
    return {"syn": syn, "fits": fits}

def is_passive(method):
    # Placeholder for passive check, based on the method name
    return method == "dummy"
