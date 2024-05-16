import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
import sys
from scipy import signal, stats
import numpy as np 
import random
from datetime import date
import pickle
import csv
from classes.cp_seg import SegmentNeurons
from classes.analyze import AnalyzeNeurons
from classes.plotData import PlotData

import tifffile as tff

from PyQt6.QtWidgets import (QApplication, QLabel, QDialog,
                             QGridLayout, QPushButton, QFileDialog,
                             QLineEdit, QCheckBox)

class Calcium(QDialog):
    """Test different ways to find spikes in Calcium Recordings."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.selection: str = None
        self.setWindowTitle("Calcium")

        # select the folder to analyze/reanalyze
        self._select_btn = QPushButton("Choose")
        self._select_btn.clicked.connect(self._select_folder)
        self.folder_c = QLabel("")

        self._check_seg = QCheckBox(text="Segment")
        self._check_ana = QCheckBox(text="Analyze")
        self._evk_seg = QCheckBox(text="Segment(evk)")
        self._evk_ana = QCheckBox(text="Analyze(evk)")

        # self.combobox = QComboBox()
        # self.combobox.addItem("mean")
        # self.combobox.addItem("median")
        # self.combobox.addItem("test both")
        # thres = QLabel("&select:")
        # thres.setBuddy(self.combobox)
        # self.combobox.activated.connect(self.mean_or_median)

        # frame = QLabel("Percentage of Frames to calculate mean/median: ")
        # self.ptg_line = QSpinBox()
        # self.ptg_line.setValue(90)

        # threshold = QHBoxLayout()
        # instruction = QLabel("Select the range of mean/median for the threshold (%): ")
        # from_lbl = QLabel("From ")
        # self.from_line = QLineEdit()
        # to_lbl = QLabel(" to ")
        # self.to_line = QLineEdit()

        # threshold.addWidget(instruction)
        # threshold.addWidget(from_lbl)
        # threshold.addWidget(self.from_line)
        # threshold.addWidget(to_lbl)
        # threshold.addWidget(self.to_line)

        filename = QLabel("Enter filename/date: ")
        self.fname = QLineEdit('TEST')

        group_name = QLabel("Enter group names for plot (separate by comma):")
        self.group_name = QLineEdit('optional')

        self._plot_btn = QPushButton("Plot")
        self._plot_btn.clicked.connect(self._plot_data)

        self.ok_btn = QPushButton("Run")
        self.ok_btn.clicked.connect(self._run)

        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.clicked.connect(self.close)

        self.layout = QGridLayout()
        self.layout.addWidget(QLabel("Select the folder: "), 0, 0)
        self.layout.addWidget(self._select_btn, 0, 1)
        self.layout.addWidget(self.folder_c, 1, 0)

        self.layout.addWidget(self._check_seg, 2, 0)
        self.layout.addWidget(self._check_ana, 2, 1)
        self.layout.addWidget(self._evk_seg, 3, 0)
        self.layout.addWidget(self._evk_ana, 3, 1)
        # self.layout.addWidget(thres, 2, 0)
        # self.layout.addWidget(self.combobox, 2, 1)
        # self.layout.addWidget(frame, 3, 0)
        # self.layout.addWidget(self.ptg_line, 3, 1)
        # self.layout.addLayout(threshold, 4, 0)
        self.layout.addWidget(filename, 5, 0)
        self.layout.addWidget(self.fname, 5, 1)
        self.layout.addWidget(self.ok_btn, 6, 0)
        self.layout.addWidget(self.cancel_btn, 6, 1)
        self.layout.addWidget(group_name, 7, 0)
        self.layout.addWidget(self.group_name, 7, 1)
        self.layout.addWidget(self._plot_btn, 8, 0)

        self.setLayout(self.layout)

        self.seg: SegmentNeurons = None
        self.analysis: AnalyzeNeurons = None
        self.plot_data: PlotData = PlotData()
        self.folder_list: list = []

    def _update_fname(self, file: str) -> None:
        '''Update the filename displayed on the selection window.'''
        _, fname = os.path.split(file)
        self.folder_c.setText(fname)
        self.layout.addWidget(self.folder_c, 1, 0)

    def _select_folder(self) -> None:
        """Select folder that contains dff.csv file."""
        dlg = QFileDialog(self)
        # dlg.setFileMode(ExistingFiles)
        self.folder_path = dlg.getExistingDirectory(None, "Select Folder")
        self._update_fname(self.folder_path)

    def _load_module(self):
        """Load Segmentation or Analyze."""
        if self._check_seg.isChecked() or self._evk_seg.isChecked():
            self.seg = SegmentNeurons()
            self.seg._load_model()

        if self._check_ana.isChecked() or self._evk_ana.isChecked():
            self.analysis = AnalyzeNeurons()

    def _run(self):
        """Run."""
        self._load_module()
        today = date.today().strftime("%y%m%d")
        additional_name = f"_{today}_{self.fname.text()}"

        # iterate through the folder, find ome.tif
        for (folder, _, filenames) in os.walk(self.folder_path):
            for file_name in filenames:
                if file_name.endswith(".ome.tif") and not file_name.startswith("._"):
                    self._record_folders(folder)
                    recording_name = file_name[:-8]
                    print(f"recording name: {recording_name}")
                    save_path = os.path.join(folder, recording_name)
                    save_path = save_path + additional_name
                    
                    # if segmentation, run segmentation
                    if self.seg and self.seg._check_model():
                        file_path = os.path.join(folder, file_name)             
                        img = tff.imread(file_path, is_mmstack=False)

                        self.roi_dict, self.raw_f = self.seg._run(img, save_path)

                        if self._evk_seg.isChecked():
                            self.seg._run_evk() ### TODO: to be implemented

                        # if also checked for analysis
                        if self.analysis:
                            if self._evk_ana.isChecked():
                                self.analysis.analyze_evk()

                            elif len(self.raw_f) > 0 and len(self.roi_dict) > 0:
                                mda_file = recording_name + "_metadata.txt"
                                mda_file = os.path.join(folder, mda_file)
                                self.analysis.analyze(self.roi_dict, None, self.raw_f, mda_file, save_path,
                                                      method="mean", frame_window_ptg=1, prom_pctg=0.25)
                            else:
                                print(f"No cells in {recording_name} to analyze. Check segmentation!")

                    # if only analysis, run analysis
                    if not self.seg and self.analysis:
                        if self._check_ana.isChecked():
                            self.analysis.reanalyze(folder, recording_name, save_path)
                        else:
                            self.analysis.analyze_evk(folder, recording_name, save_path)

            if len(self.folder_list) > 0:
                if self._evk_ana.isChecked():
                    groups = self.analysis.compile_files(self.folder_list[-1], "_compiled.csv", None, additional_name, "evk_summary.txt")
                    self._compile_plot(self.folder_list[-1], "_compiled.csv", None, groups)
                    groups = self.analysis.compile_files(self.folder_list[-1], "_compiled_st.csv", None, additional_name, "st_summary.txt")
                    self._compile_plot(self.folder_list[-1], "_compiled_st.csv", "_ST", groups)
                    groups = self.analysis.compile_files(self.folder_list[-1], "_compiled_nst.csv", None, additional_name, "nst_summary.txt")
                    self._compile_plot(self.folder_list[-1], "_compiled_nst.csv", "_NST", groups)
                else:
                    groups = self.analysis.compile_files(self.folder_list[-1], "_compiled.csv", None, additional_name)
                    self._compile_plot(self.folder_list[-1], "_compiled.csv", None, groups)

                del self.folder_list[-1]

        print("------------FINISHED-------------")
        self.analysis: AnalyzeNeurons = None
        self.seg = None
        self.folder_list: list = []

    def _record_folders(self, folder: str):
        """Record folder location for compilation."""
        if folder not in self.folder_list:
            self.folder_list.append(folder)

    def _plot_data(self):
        """Plot data."""
        # if self.group_name == "opional":
        #     self.group_name = ""
        
        # groups = self.group_name.text().split(",")
        print("Plotting")
        self.plot_data.just_plot(self.folder_path)

    def _compile_plot(self, base_folder: str, csv_name: str, evk: str | None, groups: list[str]):
        """To plot after compile."""
        for group in groups:
            compile_name = base_folder + group + csv_name
            csv_path = os.path.join(base_folder, compile_name)
            self.plot_data.ana_plot(csv_path, evk, group)

class InputGroup(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Input group name: ")

if __name__ == "__main__":
    sd_app = QApplication(sys.argv)
    sd = Calcium()
    sd.show()
    sys.exit(sd_app.exec())

