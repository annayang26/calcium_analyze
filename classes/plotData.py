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
        compiled_csv_list = []
        # find all the compiled.csv in the given folder
        for folders, dirs, fnames in os.walk(folder_path):
            for fname in fnames:
                if fname.endswith(".csv") and "compiled" in fname:
                    df = self._read_csv(os.path.join(folders, fname))
                    rec_name = df.loc[:, 'name']
                    genotype_list, groups, diff = self._group_data(rec_name)
                    print(f"genotype list is: {genotype_list}")


    def _read_csv(self, path: str) -> pd.DataFrame:
        """Read the csv file."""
        with open(path, "r") as file:
            dff_file = pd.read_csv(file)
        return dff_file

    def _group_data(self, fnames: list[str]):
        """Group the names."""
        genotype_list=[]
        groups = {}
        diff = ['first_group']

        first_name = fnames[0].split('_')     
        for i, name in enumerate(fnames):
            elements = name.split('_')
            genotype_list.append(self._genotype(elements))
            diff_ele = [ele for ele in elements if ele not in first_name and len(ele)>1 and not ele.startswith("Pos")]
            if len(diff_ele) == 0:
                groups[elements[-3]].append(i)
            elif len(diff_ele) == 1:
                groups[diff_ele[0]].append(i)
                diff.append(diff_ele)
        
        return genotype_list, groups, diff

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