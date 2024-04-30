# the goal is to randomly pick 2 UT traces from each genotype 
#  plot them all in one graph

import csv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random
import os

COLOR_LIST = ['r', 'b', 'c', 'm', 'y', 'g']
def _random_select_traces(path: str):
    """Randomly pick 2 df/f from each batch and each genotype."""
    name_dict = {}
    all_fnames = {}
    path_dict = {}
    ut_fnames = []

    control = 0
    het = 0
    null = 0

    for folder, dirs, fnames in os.walk(path):
        for fname in fnames:
            if fname.endswith('.ome.tif') and "UT" in fname:
                ut_fnames.append(fname[:-8])
                path_dict[fname[:-8]] = os.path.join(folder, fname)

    for ut_fname in ut_fnames:
        fname_ele = ut_fname.split('_')
        geno = _genotype(fname_ele)
        all_fnames[ut_fname] = geno

    index = 0
    
    while index < 6:
        rand_fname = random.choice(list(all_fnames.keys()))
        geno = all_fnames[rand_fname]
        if geno == "control" and control < 2:
            path = path_dict[rand_fname]
            name_dict[index] = [rand_fname, path]
            index += 1
            control += 1
            all_fnames.pop(rand_fname)
        elif geno == "het" and het < 2:
            path = path_dict[rand_fname]
            name_dict[index] = [rand_fname, path]
            index += 1
            het += 1
            all_fnames.pop(rand_fname)
        elif geno == "null" and null < 2:
            path = path_dict[rand_fname]
            name_dict[index] = [rand_fname, path]
            index += 1
            null += 1
            all_fnames.pop(rand_fname)

    return name_dict

def _find_csv(ut_fnames: dict):
    """Find roi_data file for the selected recordings."""
    traces = {}
    dff_max = 0
    dff_min = 1

    for index, ut_fname in ut_fnames.items():
        dir = os.path.dirname(ut_fname[1])
        ut_file = os.path.basename(ut_fname[0])
        for folders, dirs, fnames in os.walk(dir):
            for fname in fnames:
                if ut_file in folders and fname == "dff.csv":
                    csv_file = os.path.join(folders, fname)
                    with open(csv_file, "r") as file:
                        dff_file = pd.read_csv(file)

                    dff = dff_file.sample(axis='columns')
                    dff_list, d_max, d_min = _turn_df_to_list(dff)
                    traces[index] = dff_list

                    if d_max > dff_max:
                        dff_max = d_max
                    if d_min < dff_min:
                        dff_min = d_min

    return traces, dff_max, dff_min

def _turn_df_to_list(df: pd.DataFrame):
    """Turn df to list"""
    data = []
    for col in df.columns:
        data = df[col].tolist()
        df_max = max(data)
        df_min = min(data)
    return data, df_max, df_min

def _normalize_dff(traces: dict, df_max, df_min):
    """Normalize dff."""
    dff_col = traces.keys()
    n_dff = {}
    for col in dff_col:
        dff = list(traces[col]) # roi_dff[col] is the dff of one ROI -> LIST
        df_min_array = np.ones(len(dff))*df_min
        nn_dff = (dff - df_min_array) / (df_max - df_min)
        n_dff[int(col)] = nn_dff
    
    return n_dff

def _plot(traces: dict, path: str):
    """Plot the traces."""
    fig, ax = plt.subplots()
    for index, dff in traces.items():
        ax.plot(dff, color=COLOR_LIST[index], label=index)
    ax.set_title("dff")
    ax.legend()
    save_path = os.path.join(path, "dff.png")
    plt.savefig(save_path)
    plt.close()
    
def _plot_ndff(traces: dict, path: str, df_max, df_min):
    fig, ax = plt.subplots()
    nn_dff = _normalize_dff(traces, df_max, df_min)
    for index, dff in nn_dff.items():
        ax.plot(dff, color=COLOR_LIST[index], label=index)
        ax.set_title("normalized dff")
    ax.legend()
    save_path = os.path.join(path, "nn_dff.png")
    plt.savefig(save_path)
    plt.close()

def _genotype(element_list: list[str]) -> str:
    """Define genotype."""
    neg = element_list.count('-')
    pos = element_list.count('+')

    if neg == 1 and pos == 1:
        return "het"
    elif neg == 2 and pos == 0:
        return "null"
    elif neg == 0 and pos == 2:
        return "control"

def _save_file(name_dict: dict, path: str):
    save_path = os.path.join(path, 'traces_name.txt')
    with open(save_path, 'w') as f:
        for index, fname in name_dict.items():
            fname = os.path.basename(fname[0])
            f.write(f"{index} is {fname}\n")


if __name__ == "__main__":
    folder = r"G:\Anna\TSC2\NC240123_Lam77_chronic\NC240123_240229_DIV37_Lam77_BatchA\spon"
    name_dict = _random_select_traces(folder)
    traces, df_max, df_min = _find_csv(name_dict)
    _plot(traces, folder)
    _plot_ndff(traces, folder, df_max, df_min)
    _save_file(name_dict, folder)