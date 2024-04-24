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
from . import SegmentNeurons

from PyQt6.QtWidgets import (QApplication, QLabel, QDialog,
                             QGridLayout, QPushButton, QFileDialog,
                             QLineEdit)

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

        self._check_seg = QCheckBox()

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

        self.ok_btn = QPushButton("Analyze")
        self.ok_btn.clicked.connect(self._analyze)

        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.clicked.connect(self.close)

        self.layout = QGridLayout()
        self.layout.addWidget(QLabel("Select the folder: "), 0, 0)
        self.layout.addWidget(self._select_btn, 0, 1)
        self.layout.addWidget(self.folder_c, 1, 0)
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

        self._segment = SegmentNeurons()


if __name__ == "__main__":
    sd_app = QApplication(sys.argv)
    sd = Calcium()
    sd.show()
    sys.exit(sd_app.exec())

