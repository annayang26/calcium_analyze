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
        self.fname = QLineEdit()

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
        # self.layout.addWidget(thres, 2, 0)
        # self.layout.addWidget(self.combobox, 2, 1)
        # self.layout.addWidget(frame, 3, 0)
        # self.layout.addWidget(self.ptg_line, 3, 1)
        # self.layout.addLayout(threshold, 4, 0)
        self.layout.addWidget(filename, 5, 0)
        self.layout.addWidget(self.fname, 5, 1)
        self.layout.addWidget(self.ok_btn, 6, 0)
        self.layout.addWidget(self.cancel_btn, 6, 1)

        self.setLayout(self.layout)

        self.seg: SegmentNeurons = None
        self.analysis: AnalyzeNeurons = None
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
        if self._check_seg.isChecked():
            self.seg = SegmentNeurons()
            self.seg._load_model()

        if self._check_ana.isChecked():
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

                        self.roi_dict, self.dff = self.seg._run(img, save_path)

                        # if also checked for analysis
                        if self.analysis:
                            if len(self.dff) > 0 and len(self.roi_dict) > 0:
                                mda_file = recording_name + "_metadata.txt"
                                mda_file = os.path.join(folder, mda_file)
                                self.analysis.analyze(self.roi_dict, None, self.dff, mda_file, save_path)
                            else:
                                print(f"No cells in {recording_name} to analyze. Check segmentation!")

                    # if only analysis, run analysis
                    if not self.seg and self.analysis:
                        # path = os.path.join(folder, recording_name)
                        self.analysis.reanalyze(folder, recording_name, save_path)

            if len(self.folder_list) > 0:

                self.analysis.compile_files(self.folder_list[-1], f"{additional_name}_compiled.csv", None)
                del self.folder_list[-1]

        print("------------FINISHED-------------")
        self.analysis: AnalyzeNeurons = None
        self.flder_list: list = []

    def _record_folders(self, folder: str):
        """Record folder location for compilation."""
        if not (folder in self.folder_list):
            self.folder_list.append(folder)

if __name__ == "__main__":
    sd_app = QApplication(sys.argv)
    sd = Calcium()
    sd.show()
    sys.exit(sd_app.exec())

