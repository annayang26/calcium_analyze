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
from .segment import SegmentNeurons
from .analyze import AnalyzeNeurons
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

    def _select_folder(self) -> None:
        """Select folder that contains dff.csv file."""
        dlg = QFileDialog(self)
        # dlg.setFileMode(QFileDialog.Directory)
        if dlg.exec():
            self.folder_path = dlg.selectedFiles()[0]
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
        # iterate through the folder, find ome.tif
        for (folder, _, filenames) in os.walk(self.folder_path):
            for file_name in filenames:
                if file_name.endswith(".ome.tif"):
                    recording_name = file_name[:-8]

                    # if segmentation, run segmentation
                    if self.seg and self.seg._check_model():
                        file_path = os.path.join(folder, file_name)
                        img = tff.imread(file_path, is_ome=False, is_mmstack=False)
                        print(f"the shape of img is {img.shape}")
                        self.roi_dict, self.dff = self.seg._run(img, file_path)

                        if self.analysis:
                            if len(self.dff) > 0 and len(self.roi_dict) > 0:
                                pass
                            else:
                                print(f"No cells in {recording_name} to analyze. Check segmentation!")
                            



                    # if only analysis, run analysis
                    if not self.seg and self.analysis:
                        path = os.path.join(folder, recording_name)
                        self.analysis._analyze(path)

if __name__ == "__main__":
    sd_app = QApplication(sys.argv)
    sd = Calcium()
    sd.show()
    sys.exit(sd_app.exec())

