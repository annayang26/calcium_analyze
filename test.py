# the goal is to randomly pick 2 UT traces from each genotype 
#  plot them all in one graph

import csv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random
import os
import matplotlib.colors as mcolors

COLOR_LIST = ['xkcd:purple', 'xkcd:light purple', 'xkcd:blue', 'xkcd:light blue', 'xkcd:dark green', 'xkcd:soft green']
def _random_select_traces(path: str):
    """Randomly pick 2 df/f from each batch and each genotype."""
    name_dict = {}
    all_fnames = {}
    path_dict = {}
    ut_fnames = []
    geno_dict = {}

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
            geno_dict[index] = geno
            index += 1
            control += 1
            all_fnames.pop(rand_fname)
            
        elif geno == "het" and het < 2:
            path = path_dict[rand_fname]
            name_dict[index] = [rand_fname, path]
            geno_dict[index] = geno
            index += 1
            het += 1
            all_fnames.pop(rand_fname)
            
        elif geno == "null" and null < 2:
            path = path_dict[rand_fname]
            name_dict[index] = [rand_fname, path]
            geno_dict[index] = geno
            index += 1
            null += 1
            all_fnames.pop(rand_fname)

    return name_dict, geno_dict

def _find_csv(ut_fnames: dict):
    """Find roi_data file for the selected recordings."""
    dff_traces = {}
    raw_signal = {}
    dff_max = 0
    dff_min = 1

    for index, ut_fname in ut_fnames.items():
        dir = os.path.dirname(ut_fname[1])
        ut_file = os.path.basename(ut_fname[0])
        for folders, dirs, fnames in os.walk(dir):
            dff_file = pd.DataFrame([])
            raw_signal_file = pd.DataFrame([])            
            for fname in fnames:
                if ut_file in folders:
                    if fname == "dff.csv":
                        csv_file = os.path.join(folders, fname)
                        with open(csv_file, "r") as file:
                            dff_file = pd.read_csv(file)
                    if fname == "raw_signal.csv":
                        csv_file = os.path.join(folders, fname)
                        with open(csv_file, "r") as file:
                            raw_signal_file = pd.read_csv(file)

            if not dff_file.empty and not raw_signal_file.empty:
                dff = dff_file.sample(axis='columns')
                col = dff.columns.values.tolist()
                raw_signal_list = raw_signal_file.loc[:, col[0]].tolist()
                dff_list, _, _ = _turn_df_to_list(dff)

                dff_traces[index] = dff_list
                raw_signal[index] = raw_signal_list

    return dff_traces, raw_signal, dff_max, dff_min

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

def _plot(traces: dict, path: str, geno_dict: dict, csv_name: str):
    """Plot the traces."""
    fig, ax = plt.subplots()
    color_list = _select_color(geno_dict)
    for index, dff in traces.items():
        ax.plot(dff, color=color_list[index], label=geno_dict[index])
    ax.set_title(f"{csv_name} of untreated TSC2")
    ax.legend()
    save_path = os.path.join(path, f"{csv_name}.png")
    plt.savefig(save_path)
    plt.close()
    
def _select_color(geno_dict: dict):
    """Select color based on the genotype."""
    color_list = {}
    control_start = 0
    het_start = 2
    null_start = 4
    for index, geno in geno_dict.items():
        if geno == "control":
            color_list[index] = COLOR_LIST[control_start]
            control_start += 1
        elif geno == "het":
            color_list[index] = COLOR_LIST[het_start]
            het_start += 1
        elif geno == "null":
            color_list[index] = COLOR_LIST[null_start]
            null_start += 1
    return color_list
            

# def _plot_ndff(traces: dict, path: str, df_max, df_min):
#     fig, ax = plt.subplots()
#     nn_dff = _normalize_dff(traces, df_max, df_min)
#     for index, dff in nn_dff.items():
#         ax.plot(dff, color=COLOR_LIST[index], label=index)
#         ax.set_title("normalized dff")
#     ax.legend()
#     save_path = os.path.join(path, "nn_dff.png")
#     plt.savefig(save_path)
#     plt.close()

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

def _save_file(name_dict: dict, path: str, csv_name: str = "dff"):
    save_path = os.path.join(path, f'{csv_name}_traces_name.txt')
    with open(save_path, 'w') as f:
        for index, fname in name_dict.items():
            fname = os.path.basename(fname[0])
            f.write(f"{index} is {fname}\n")


if __name__ == "__main__":
    folder = r"G:\Anna\TSC2\NC240123_Lam77_chronic\spon\NC240123_240321_DIV58_Chronic_spon\BatchA"
    name_dict, geno_dict = _random_select_traces(folder)
    traces, raw_signal, df_max, df_min = _find_csv(name_dict)
    _plot(traces, folder, geno_dict, "dff")
    _plot(traces, folder, geno_dict, "raw_siganl")
    _save_file(name_dict, folder)