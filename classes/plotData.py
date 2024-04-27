import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os 
import csv

class PlotData():
    def __init__(self) -> None:
        pass

    def start_plotting(self, folder_path: str):
        """Start plotting."""
        # 1. find all compiled csv files for the same batch/condition/date

        # 1.a compile them into one giant csv file, with the original name saved in each csv file

        # 1.b save the giant csv file

        # 2. plot 
        
    def _find_csv(self, folder_path: str):
        # compiled_csv_list = []
        # find all the compiled.csv in the given folder
        for folders, dirs, fnames in os.walk(folder_path):
            for fname in fnames:
                if fname.endswith("_compiled.csv") and "._" not in fname:
                    file_path = os.path.join(folders, fname)
                    df = self._read_csv(file_path)
                    save_path = file_path[:-len("_compiled.csv")]
                    rec_name = df.loc[:, 'name']
                    genotypes, groups, _ = self._group_data(rec_name)

                    cols = df.columns
                    for col in cols:
                        if col == "name" or "Standard Deviation" in col:
                            continue
                        self._plot(genotypes, groups, df, col, save_path, None)

                elif fname.endswith("_compiled_st.csv") and "._" not in fname:
                    file_path = os.path.join(folders, fname)
                    df = self._read_csv(file_path)
                    save_path = file_path[:-len("_compiled_st.csv")] + "_ST"
                    rec_name = df.loc[:, 'name']
                    genotypes, groups, _ = self._group_data(rec_name)

                    cols = df.columns
                    for col in cols:
                        if col == "name" or "Standard Deviation" in col:
                            continue
                        self._plot(genotypes, groups, df, col, save_path, "ST")

                elif fname.endswith("_compiled_nst.csv") and "._" not in fname:
                    file_path = os.path.join(folders, fname)
                    df = self._read_csv(file_path)
                    save_path = file_path[:-len("_compiled_nst.csv")] + "_NST"
                    rec_name = df.loc[:, 'name']
                    genotypes, groups, _ = self._group_data(rec_name)

                    cols = df.columns
                    for col in cols:
                        if col == "name" or "Standard Deviation" in col:
                            continue
                        self._plot(genotypes, groups, df, col, save_path, "NST")

    def _read_csv(self, path: str) -> pd.DataFrame:
        """Read the csv file."""
        with open(path, "r") as file:
            dff_file = pd.read_csv(file)
        return dff_file

    def _group_data(self, fnames: list[str]):
        """Group the names."""
        genotypes={}
        groups = {}
        diff = []
        first_group = ""

        first_name = fnames[0].split('_')     
        for i, name in enumerate(fnames):
            elements = name.split('_')
            if i == 0:
                first_group = elements.index("MMStack") + 1

            genotype = self._genotype(elements)
            if not genotypes.get(genotype):
                genotypes[genotype] = []
            genotypes[genotype].append(genotype)

            diff_ele = [ele for ele in elements if ele not in first_name and len(ele)>1 and not ele.startswith("Pos")]
            if len(diff_ele) == 0:
                if not groups.get(elements[first_group]):
                    groups[elements[first_group]] = []
                groups[elements[first_group]].append(i)
            elif len(diff_ele) == 1:
                if not groups.get(diff_ele[0]):
                    groups[diff_ele[0]] = []
                groups[diff_ele[0]].append(i)
                diff.append(diff_ele)
        
        return genotypes, groups, diff

    def _genotype(self, element_list: list[str]) -> str:
        """Define genotype."""
        neg = element_list.count('-')
        pos = element_list.count('+')

        if neg == 1 and pos == 1:
            return "het"
        elif neg == 2 and pos == 0:
            return "null"
        elif neg == 0 and pos == 2:
            return "control"

    def _get_data(self, all_data: pd.DataFrame, metric: str, group_ind: list) -> list:
        """Get data for one group."""
        group_data = []
        for ind in group_ind:
            group_data.append(all_data.loc[ind, metric])
        
        return group_data

    def _plot(self, genotypes: dict, groups: dict, all_data: pd.DataFrame, 
              metric: str, path: str, evk: str | None):
        """Plot the metric """
        fig, ax = plt.subplots()
        start_x = 1

        for geno in genotypes:
            for group, index in groups.items():
                data = self._get_data(all_data, metric, index)

                x_range = np.ones(len(data)) * start_x
                ax.scatter(x_range, data, label=group)

                title = f"{geno}_{metric}"
                if evk:
                    title += f"_{evk}"
                ax.set_title(title)

                ax.set_xticks([])
                ax.legend()
                start_x += 1

        if "Average Frequency" in metric:
            metric = "Average Frequency"

        folder_path = os.path.join(path, f"{path}_{geno}_graphs")
        if not os.path.isdir(folder_path):
            os.mkdir(folder_path)

        save_path = os.path.join(folder_path, f"{metric}.png")
        plt.savefig(save_path)
        plt.close()
