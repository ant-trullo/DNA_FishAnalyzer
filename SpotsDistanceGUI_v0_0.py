"""This is the main window of the software to analyze spots distance in different smFiSH channels.
This is version 0.0, since January 2021.

"""


import sys
# import os.path
# import datetime
# from shutil import copyfile
from importlib import reload
import traceback
import numpy as np
import pyqtgraph as pg
from PyQt5.QtCore import Qt
from PyQt5 import QtGui, QtWidgets, QtCore
from scipy.ndimage import gaussian_laplace, gaussian_filter
from scipy.ndimage.morphology import binary_fill_holes
# from skimage.morphology import disk
from skimage.filters import gaussian

import RawDataLoader
import NucsSegmenter3
import NucsPileUp
import LabelsModify
import CrossInfo
import ModifCrossInfo
import RemoveBorderNuclei


class MainWindow(QtWidgets.QMainWindow):
    """Main windows: coordinates all the actions, algorithms, visualization tools and analysis tools"""
    def __init__(self, parent=None):

        QtWidgets.QMainWindow.__init__(self, parent)

        widget  =  QtWidgets.QWidget(self)
        self.setCentralWidget(widget)

        load_data_action  =  QtWidgets.QAction(QtGui.QIcon('Icons/load-hi.png'), "&Load data", self)
        load_data_action.setShortcut("Ctrl+L")
        load_data_action.setStatusTip("Load lsm files and the xls output of FQ")
        load_data_action.triggered.connect(self.load_data)

        save_analysis_action  =  QtWidgets.QAction(QtGui.QIcon('Icons/save-md.png'), "&Save Segmention", self)
        save_analysis_action.setShortcut("Ctrl+S")
        save_analysis_action.setStatusTip("Save segmentatyion analysis")
        save_analysis_action.triggered.connect(self.save_analysis)

        cross_info_action  =  QtWidgets.QAction(QtGui.QIcon('Icons/crossInfo.gif'), "&Cross Info", self)
        cross_info_action.setStatusTip("Cross info of nuclei and spots segmentation")
        cross_info_action.triggered.connect(self.cross_info)

        check_overlaps_action  =  QtWidgets.QAction(QtGui.QIcon('Icons/check_overlaps.jpg'), "&Check Overlaps", self)
        check_overlaps_action.setStatusTip("Check the randomness of the overlaps")
        check_overlaps_action.triggered.connect(self.check_overlaps)
        check_overlaps_action.setShortcut("Ctrl+P")

        change_sptsdet_action  =  QtWidgets.QAction(QtGui.QIcon('Icons/spots_detection.jpg'), "&Change Spots", self)
        # change_sptsdet_action_action.setShortcut("Ctrl+Q")
        change_sptsdet_action.setStatusTip("Change the spots filtering for the detection")
        change_sptsdet_action.triggered.connect(self.change_sptsdet)

        exit_action  =  QtWidgets.QAction(QtGui.QIcon('Icons/exit.png'), "&Exit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.setStatusTip("Exit application")
        exit_action.triggered.connect(self.close)

        menubar   =  self.menuBar()

        file_menu  =  menubar.addMenu("&File")
        file_menu.addAction(load_data_action)
        file_menu.addAction(save_analysis_action)
        file_menu.addAction(cross_info_action)
        file_menu.addAction(check_overlaps_action)
        file_menu.addAction(exit_action)

        modify_menu  =  menubar.addMenu("&Modify")
        modify_menu.addAction(change_sptsdet_action)

        fname_raw_lbl  =  QtWidgets.QLabel("File: ", self)
        fname_raw_lbl.setToolTip("Name of the file you are working on")

        pixsize_x_lbl  =  QtWidgets.QLabel("pix size XY =;")
        pixsize_z_lbl  =  QtWidgets.QLabel("Z step =")

        tabs_tot  =  QtWidgets.QTabWidget()
        tab_ch1   =  QtWidgets.QWidget()
        tab_ch2   =  QtWidgets.QWidget()
        tab_ch3   =  QtWidgets.QWidget()
        tab_dapi  =  QtWidgets.QWidget()

        tabs_tot.addTab(tab_ch1, "CH1")  # "Control")
        tabs_tot.addTab(tab_ch2, "CH2")  # "Escargo")
        tabs_tot.addTab(tab_ch3, "CH3")    # "Snail")
        tabs_tot.addTab(tab_dapi, "Dapi")

        # ###CONTROL START### #
        frame_ch1_raw  =  pg.ImageView(self, name="FrameControlRaw")
        frame_ch1_raw.ui.roiBtn.hide()
        frame_ch1_raw.ui.menuBtn.hide()
        frame_ch1_raw.timeLine.sigPositionChanged.connect(self.frame_ch1_raw_change)

        frame_ch1_flt  =  pg.ImageView()
        frame_ch1_flt.ui.roiBtn.hide()
        frame_ch1_flt.ui.menuBtn.hide()
        frame_ch1_flt.view.setXLink("FrameControlRaw")
        frame_ch1_flt.view.setYLink("FrameControlRaw")
        frame_ch1_flt.timeLine.sigPositionChanged.connect(self.frame_ch1_flt_change)

        frame_ch1_sgm  =  pg.ImageView()
        frame_ch1_sgm.ui.roiBtn.hide()
        frame_ch1_sgm.ui.menuBtn.hide()
        frame_ch1_sgm.view.setXLink("FrameControlRaw")
        frame_ch1_sgm.view.setYLink("FrameControlRaw")
        frame_ch1_sgm.timeLine.sigPositionChanged.connect(self.frame_ch1_sgm_change)

        tabs_ch1  =  QtWidgets.QTabWidget()
        tabs_ch1.addTab(frame_ch1_raw, "Raw Control (CH1)")
        tabs_ch1.addTab(frame_ch1_flt, "Filtered (CH1)")
        tabs_ch1.addTab(frame_ch1_sgm, "Segm (CH1)")

        sld_thr_ch1  =  QtWidgets.QScrollBar(QtCore.Qt.Vertical, self)
        sld_thr_ch1.valueChanged.connect(self.sld_thr_ch1_update)

        sld_ch1_lbl  =  QtWidgets.QLabel("0 %")
        sld_ch1_lbl.setFixedSize(30, 25)

        sld_ch1_box  =  QtWidgets.QVBoxLayout()
        sld_ch1_box.addWidget(sld_thr_ch1)
        sld_ch1_box.addWidget(sld_ch1_lbl)

        ch1_layout  =  QtWidgets.QHBoxLayout()
        ch1_layout.addWidget(tabs_ch1)
        ch1_layout.addLayout(sld_ch1_box)

        tab_ch1.setLayout(ch1_layout)
        # ###CONTROL END### #

        # ###ESCARGOL START### #
        frame_ch2_raw  =  pg.ImageView(self, name="FrameEscargoRaw")
        frame_ch2_raw.ui.roiBtn.hide()
        frame_ch2_raw.ui.menuBtn.hide()
        frame_ch2_raw.timeLine.sigPositionChanged.connect(self.frame_ch2_raw_change)

        frame_ch2_flt  =  pg.ImageView()
        frame_ch2_flt.ui.roiBtn.hide()
        frame_ch2_flt.ui.menuBtn.hide()
        frame_ch2_flt.view.setXLink("FrameEscargoRaw")
        frame_ch2_flt.view.setYLink("FrameEscargoRaw")
        frame_ch2_flt.timeLine.sigPositionChanged.connect(self.frame_ch2_flt_change)

        frame_ch2_sgm  =  pg.ImageView()
        frame_ch2_sgm.ui.roiBtn.hide()
        frame_ch2_sgm.ui.menuBtn.hide()
        frame_ch2_sgm.view.setXLink("FrameEscargoRaw")
        frame_ch2_sgm.view.setYLink("FrameEscargoRaw")
        frame_ch2_sgm.timeLine.sigPositionChanged.connect(self.frame_ch2_sgm_change)

        tabs_ch2  =  QtWidgets.QTabWidget()
        tabs_ch2.addTab(frame_ch2_raw, "Raw Control (CH2)")
        tabs_ch2.addTab(frame_ch2_flt, "Filtered (CH2)")
        tabs_ch2.addTab(frame_ch2_sgm, "Segm (CH2)")

        sld_thr_ch2  =  QtWidgets.QScrollBar(QtCore.Qt.Vertical, self)
        sld_thr_ch2.valueChanged.connect(self.sld_thr_ch2_update)

        sld_ch2_lbl  =  QtWidgets.QLabel("0 %")
        sld_ch2_lbl.setFixedSize(30, 25)

        sld_ch2_box  =  QtWidgets.QVBoxLayout()
        sld_ch2_box.addWidget(sld_thr_ch2)
        sld_ch2_box.addWidget(sld_ch2_lbl)

        ch2_layout  =  QtWidgets.QHBoxLayout()
        ch2_layout.addWidget(tabs_ch2)
        ch2_layout.addLayout(sld_ch2_box)

        tab_ch2.setLayout(ch2_layout)
        # ###ESCARGO END### #

        # ###SNAIL START### #
        frame_ch3_raw    =  pg.ImageView(self, name="FrameSnailRaw")
        frame_ch3_raw.ui.roiBtn.hide()
        frame_ch3_raw.ui.menuBtn.hide()
        frame_ch3_raw.timeLine.sigPositionChanged.connect(self.frame_ch3_raw_change)

        frame_ch3_flt  =  pg.ImageView()
        frame_ch3_flt.ui.roiBtn.hide()
        frame_ch3_flt.ui.menuBtn.hide()
        frame_ch3_flt.view.setXLink("FrameSnailRaw")
        frame_ch3_flt.view.setYLink("FrameSnailRaw")
        frame_ch3_flt.timeLine.sigPositionChanged.connect(self.frame_ch3_flt_change)

        frame_ch3_sgm  =  pg.ImageView()
        frame_ch3_sgm.ui.roiBtn.hide()
        frame_ch3_sgm.ui.menuBtn.hide()
        frame_ch3_sgm.view.setXLink("FrameSnailRaw")
        frame_ch3_sgm.view.setYLink("FrameSnailRaw")
        frame_ch3_sgm.timeLine.sigPositionChanged.connect(self.frame_ch3_sgm_change)

        tabs_ch3  =  QtWidgets.QTabWidget()
        tabs_ch3.addTab(frame_ch3_raw, "Raw Control (CH3)")
        tabs_ch3.addTab(frame_ch3_flt, "Filtered (CH3)")
        tabs_ch3.addTab(frame_ch3_sgm, "Segm (CH3)")

        sld_thr_ch3  =  QtWidgets.QScrollBar(QtCore.Qt.Vertical, self)
        sld_thr_ch3.valueChanged.connect(self.sld_thr_ch3_update)

        sld_ch3_lbl  =  QtWidgets.QLabel("0 %")
        sld_ch3_lbl.setFixedSize(30, 25)

        sld_ch3_box  =  QtWidgets.QVBoxLayout()
        sld_ch3_box.addWidget(sld_thr_ch3)
        sld_ch3_box.addWidget(sld_ch3_lbl)

        ch3_layout  =  QtWidgets.QHBoxLayout()
        ch3_layout.addWidget(tabs_ch3)
        ch3_layout.addLayout(sld_ch3_box)

        tab_ch3.setLayout(ch3_layout)
        # ###SNAIL END### #

        # ###DAPI START### #
        frame_dapi_raw    =  pg.ImageView(self, name="FrameDapiRaw")
        frame_dapi_raw.ui.roiBtn.hide()
        frame_dapi_raw.ui.menuBtn.hide()
        frame_dapi_raw.timeLine.sigPositionChanged.connect(self.frame_dapi_raw_change)

        frame_dapi_2D    =  pg.ImageView()
        frame_dapi_2D.ui.roiBtn.hide()
        frame_dapi_2D.ui.menuBtn.hide()
        frame_dapi_2D.view.setXLink("FrameDapiRaw")
        frame_dapi_2D.view.setYLink("FrameDapiRaw")
        frame_dapi_2D.timeLine.sigPositionChanged.connect(self.frame_dapi_2D_change)

        frame_dapi_2D_sgm    =  pg.ImageView()
        frame_dapi_2D_sgm.ui.roiBtn.hide()
        frame_dapi_2D_sgm.ui.menuBtn.hide()
        frame_dapi_2D_sgm.view.setXLink("FrameDapiRaw")
        frame_dapi_2D_sgm.view.setYLink("FrameDapiRaw")
        frame_dapi_2D_sgm.timeLine.sigPositionChanged.connect(self.frame_dapi_2D_sgm_change)
        frame_dapi_2D_sgm.getImageItem().mouseClickEvent  =  self.click

        frame_dapi_3D    =  pg.ImageView()
        frame_dapi_3D.ui.roiBtn.hide()
        frame_dapi_3D.ui.menuBtn.hide()
        frame_dapi_3D.view.setXLink("FrameDapiRaw")
        frame_dapi_3D.view.setYLink("FrameDapiRaw")
        frame_dapi_3D.timeLine.sigPositionChanged.connect(self.frame_dapi_3D_change)

        plus_dapi_thr_btn  =  QtWidgets.QPushButton("+", self)
        plus_dapi_thr_btn.setFixedSize(20, 25)
        plus_dapi_thr_btn.clicked.connect(self.plus_dapi_thr)

        minus_dapi_thr_btn  =  QtWidgets.QPushButton("-", self)
        minus_dapi_thr_btn.setFixedSize(20, 25)
        minus_dapi_thr_btn.clicked.connect(self.minus_dapi_thr)

        smooth_dapi_btn  =  QtWidgets.QPushButton("Smooth", self)
        smooth_dapi_btn.setFixedSize(110, 25)
        smooth_dapi_btn.clicked.connect(self.smooth_dapi)

        plusminus_thr_lbl  =  QtWidgets.QLabel("Thr")

        plusminus_box  =  QtWidgets.QHBoxLayout()
        plusminus_box.addWidget(plusminus_thr_lbl)
        plusminus_box.addWidget(plus_dapi_thr_btn)
        plusminus_box.addWidget(minus_dapi_thr_btn)

        tabs_dapi  =  QtWidgets.QTabWidget()
        tabs_dapi.addTab(frame_dapi_raw, "Raw Dapi")
        tabs_dapi.addTab(frame_dapi_2D, "Detected")
        tabs_dapi.addTab(frame_dapi_2D_sgm, "Segmented")
        tabs_dapi.addTab(frame_dapi_3D, "3D")

        nucs_2dsegm_first_lbl  =  QtWidgets.QLabel("first")
        nucs_2dsegm_first_lbl.setFixedSize(50, 25)

        nucs_2dsegm_first_btn  =  QtWidgets.QPushButton("First", self)
        nucs_2dsegm_first_btn.setFixedSize(50, 25)
        nucs_2dsegm_first_btn.clicked.connect(self.nucs_2dsegm_first)

        nucs_2dsegm_last_lbl  =  QtWidgets.QLabel("last")
        nucs_2dsegm_last_lbl.setFixedSize(50, 25)

        nucs_2dsegm_last_btn  =  QtWidgets.QPushButton("Last", self)
        nucs_2dsegm_last_btn.setFixedSize(50, 25)
        nucs_2dsegm_last_btn.clicked.connect(self.nucs_2dsegm_last)

        nucs_2dsegm_first_lbledt  =  QtWidgets.QHBoxLayout()
        nucs_2dsegm_first_lbledt.addWidget(nucs_2dsegm_first_lbl)
        nucs_2dsegm_first_lbledt.addWidget(nucs_2dsegm_first_btn)

        nucs_2dsegm_last_lbledt  =  QtWidgets.QHBoxLayout()
        nucs_2dsegm_last_lbledt.addWidget(nucs_2dsegm_last_lbl)
        nucs_2dsegm_last_lbledt.addWidget(nucs_2dsegm_last_btn)

        nucs_2dsegm_btn  =  QtWidgets.QPushButton("N Detect", self)
        nucs_2dsegm_btn.clicked.connect(self.nucs_2dsegm)
        nucs_2dsegm_btn.setToolTip("Segmentation of nuclei z by z")
        nucs_2dsegm_btn.setFixedSize(110, 25)

        nucs_2dsegm_ws_btn  =  QtWidgets.QPushButton("N Segment", self)
        nucs_2dsegm_ws_btn.clicked.connect(self.nucs_2dsegm_ws)
        nucs_2dsegm_ws_btn.setToolTip("Watershed of nuclei z by z")
        nucs_2dsegm_ws_btn.setFixedSize(75, 25)

        ws_mindist_edt  =  QtWidgets.QLineEdit(self)
        ws_mindist_edt.setFixedSize(25, 25)
        ws_mindist_edt.textChanged[str].connect(self.ws_mindist_var)
        ws_mindist_edt.setText("80")

        nucs_2dsegm_ws_mindist  =  QtWidgets.QHBoxLayout()
        nucs_2dsegm_ws_mindist.addWidget(nucs_2dsegm_ws_btn)
        nucs_2dsegm_ws_mindist.addWidget(ws_mindist_edt)

        nucs_2dpileup_btn  =  QtWidgets.QPushButton("Pile Up", self)
        nucs_2dpileup_btn.clicked.connect(self.nucs2D_pileup)
        nucs_2dpileup_btn.setToolTip("Pile up nuclei slices in all the z frames")
        nucs_2dpileup_btn.setFixedSize(110, 25)

        man_active_toggle  =  QtWidgets.QCheckBox("Hand Cut", self)
        man_active_toggle.setFixedSize(110, 25)
        man_active_toggle.stateChanged.connect(self.man_active)

        manual_cut_btn  =  QtWidgets.QPushButton("Modify", self)
        manual_cut_btn.clicked.connect(self.manual_cut)
        manual_cut_btn.setToolTip("Manual cutting of nuclei (Ctrl+Suppr)")
        manual_cut_btn.setFixedSize(110, 25)
        manual_cut_btn.setEnabled(False)

        dapi_keys  =  QtWidgets.QVBoxLayout()
        dapi_keys.addLayout(nucs_2dsegm_first_lbledt)
        dapi_keys.addLayout(nucs_2dsegm_last_lbledt)
        dapi_keys.addWidget(nucs_2dsegm_btn)
        dapi_keys.addLayout(plusminus_box)
        dapi_keys.addWidget(smooth_dapi_btn)
        dapi_keys.addLayout(nucs_2dsegm_ws_mindist)
        dapi_keys.addWidget(man_active_toggle)
        dapi_keys.addWidget(manual_cut_btn)
        dapi_keys.addWidget(nucs_2dpileup_btn)
        dapi_keys.addStretch()

        dapi_layout  =  QtWidgets.QHBoxLayout()
        dapi_layout.addWidget(tabs_dapi)
        dapi_layout.addLayout(dapi_keys)

        tab_dapi.setLayout(dapi_layout)

        # ###DAPI END### #

        busy_lbl  =  QtWidgets.QLabel("Ready")
        busy_lbl.setStyleSheet("color: green")

        bottom_labels_box  =  QtWidgets.QHBoxLayout()
        bottom_labels_box.addWidget(busy_lbl)
        bottom_labels_box.addStretch()
        bottom_labels_box.addWidget(pixsize_x_lbl)
        bottom_labels_box.addWidget(pixsize_z_lbl)

        layout  =  QtWidgets.QVBoxLayout(widget)
        layout.addWidget(fname_raw_lbl)
        layout.addWidget(tabs_tot)
        layout.addLayout(bottom_labels_box)

        mycmap           =  np.fromfile("mycmap.bin", "uint16").reshape((10000, 3))      # / 255.0
        self.colors4map  =  []
        for k in range(mycmap.shape[0]):
            self.colors4map.append(mycmap[k, :])
        self.colors4map     =  self.colors4map + self.colors4map + self.colors4map + self.colors4map + self.colors4map + self.colors4map
        self.colors4map[0]  =  np.array([0, 0, 0])

        self.soft_version             =  "DNA_FishAnalyzer_v0.0"
        self.fname_raw_lbl            =  fname_raw_lbl
        self.pixsize_x_lbl            =  pixsize_x_lbl
        self.pixsize_z_lbl            =  pixsize_z_lbl
        self.busy_lbl                 =  busy_lbl
        self.man_active_flag          =  0
        self.c_count                  =  0
        self.nucs_2dsegm_gkern_value  =  1.5
        self.nucs_2dsegm_thrf_value   =  1.1
        # self.spts_clusters_flag       =  "s"

        self.frame_ch1_raw          =  frame_ch1_raw
        self.frame_ch1_flt          =  frame_ch1_flt
        self.frame_ch1_sgm          =  frame_ch1_sgm
        self.frame_ch2_raw          =  frame_ch2_raw
        self.frame_ch2_flt          =  frame_ch2_flt
        self.frame_ch2_sgm          =  frame_ch2_sgm
        self.frame_ch3_raw          =  frame_ch3_raw
        self.frame_ch3_flt          =  frame_ch3_flt
        self.frame_ch3_sgm          =  frame_ch3_sgm
        self.frame_dapi_raw         =  frame_dapi_raw
        self.frame_dapi_2D          =  frame_dapi_2D
        self.frame_dapi_2D_sgm      =  frame_dapi_2D_sgm
        self.frame_dapi_3D          =  frame_dapi_3D
        self.sld_thr_ch1            =  sld_thr_ch1
        self.sld_thr_ch2            =  sld_thr_ch2
        self.sld_thr_ch3            =  sld_thr_ch3
        self.nucs_2dsegm_first_lbl  =  nucs_2dsegm_first_lbl
        self.nucs_2dsegm_last_lbl   =  nucs_2dsegm_last_lbl
        self.manual_cut_btn         =  manual_cut_btn
        self.sld_ch1_lbl            =  sld_ch1_lbl
        self.sld_ch2_lbl            =  sld_ch2_lbl
        self.sld_ch3_lbl            =  sld_ch3_lbl
        self.tabs_tot               =  tabs_tot

        self.ws_mindist_var("80")

        self.setGeometry(800, 100, 700, 500)
        self.setWindowTitle(self.soft_version)
        self.setWindowIcon(QtGui.QIcon('Icons/DrosophilaIcon.png'))
        self.show()

    def closeEvent(self, event):
        """Close the GUI, asking confirmation"""
        quit_msg  =  "Are you sure you want to exit the program?"
        reply     =  QtWidgets.QMessageBox.question(self, 'Message', quit_msg, QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)

        if reply == QtWidgets.QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

    def busy_indicator(self):
        """Write a red text (BUSY) as a label on the GUI (bottom left)"""
        self.busy_lbl.setText("Busy")
        self.busy_lbl.setStyleSheet('color: red')

    def ready_indicator(self):
        """Write a green text (READY) as a label on the GUI (bottom left)"""
        self.busy_lbl.setText("Ready")
        self.busy_lbl.setStyleSheet('color: green')

    def frame_dapi_raw_change(self):
        """Synchronize dapi frames when raw dapi is changed"""
        try:
            self.frame_dapi_2D.setCurrentIndex(self.frame_dapi_raw.currentIndex)
        except AttributeError:
            pass
        try:
            self.frame_dapi_2D_sgm.setCurrentIndex(self.frame_dapi_raw.currentIndex)
        except AttributeError:
            pass
        try:
            self.frame_dapi_3D.setCurrentIndex(self.frame_dapi_raw.currentIndex)
        except AttributeError:
            pass

    def frame_dapi_2D_change(self):
        """Synchronize dapi frames when 2D dapi is changed"""
        self.frame_dapi_raw.setCurrentIndex(self.frame_dapi_2D.currentIndex)
#         try:
#             self.frame_dapi_3D.setCurrentIndex(self.frame_dapi_2D.currentIndex)
#         except AttributeError:
#             pass

    def frame_dapi_2D_sgm_change(self):
        """Synchronize dapi frames when 3D dapi is changed"""
        self.frame_dapi_raw.setCurrentIndex(self.frame_dapi_2D_sgm.currentIndex)
#         self.frame_dapi_2D.setCurrentIndex(self.frame_dapi_3D.currentIndex)

    def frame_dapi_3D_change(self):
        """Synchronize dapi frames when 3D dapi is changed"""
        self.frame_dapi_raw.setCurrentIndex(self.frame_dapi_3D.currentIndex)
#         self.frame_dapi_2D.setCurrentIndex(self.frame_dapi_3D.currentIndex)

    def frame_ch1_raw_change(self):
        """Synchronize dapi frames when 3D dapi is changed"""
        self.frame_ch1_sgm.setCurrentIndex(self.frame_ch1_raw.currentIndex)
        self.frame_ch1_flt.setCurrentIndex(self.frame_ch1_raw.currentIndex)

    def frame_ch1_sgm_change(self):
        """Synchronize dapi frames when 3D dapi is changed"""
        self.frame_ch1_raw.setCurrentIndex(self.frame_ch1_sgm.currentIndex)
        self.frame_ch1_flt.setCurrentIndex(self.frame_ch1_sgm.currentIndex)

    def frame_ch1_flt_change(self):
        """Synchronize dapi frames when 3D dapi is changed"""
        self.frame_ch1_raw.setCurrentIndex(self.frame_ch1_flt.currentIndex)
        self.frame_ch1_sgm.setCurrentIndex(self.frame_ch1_flt.currentIndex)

    def frame_ch2_raw_change(self):
        """Synchronize dapi frames when 3D dapi is changed"""
        self.frame_ch2_sgm.setCurrentIndex(self.frame_ch2_raw.currentIndex)
        self.frame_ch2_flt.setCurrentIndex(self.frame_ch2_raw.currentIndex)

    def frame_ch2_sgm_change(self):
        """Synchronize dapi frames when 3D dapi is changed"""
        self.frame_ch2_raw.setCurrentIndex(self.frame_ch2_sgm.currentIndex)
        self.frame_ch2_flt.setCurrentIndex(self.frame_ch2_sgm.currentIndex)

    def frame_ch2_flt_change(self):
        """Synchronize dapi frames when 3D dapi is changed"""
        self.frame_ch2_raw.setCurrentIndex(self.frame_ch2_flt.currentIndex)
        self.frame_ch2_sgm.setCurrentIndex(self.frame_ch2_flt.currentIndex)

    def frame_ch3_raw_change(self):
        """Synchronize dapi frames when 3D dapi is changed"""
        self.frame_ch3_sgm.setCurrentIndex(self.frame_ch3_raw.currentIndex)
        self.frame_ch3_flt.setCurrentIndex(self.frame_ch3_raw.currentIndex)

    def frame_ch3_sgm_change(self):
        """Synchronize dapi frames when 3D dapi is changed"""
        self.frame_ch3_raw.setCurrentIndex(self.frame_ch3_sgm.currentIndex)
        self.frame_ch3_flt.setCurrentIndex(self.frame_ch3_sgm.currentIndex)

    def frame_ch3_flt_change(self):
        """Synchronize dapi frames when 3D dapi is changed"""
        self.frame_ch3_raw.setCurrentIndex(self.frame_ch3_flt.currentIndex)
        self.frame_ch3_sgm.setCurrentIndex(self.frame_ch3_flt.currentIndex)

    def nucs_2dsegm_first(self):
        """Set the first of the frames on which calculate the threshold"""
        self.nucs_2dsegm_first_value  =  self.frame_dapi_raw.currentIndex
        self.nucs_2dsegm_first_lbl.setText(str(self.nucs_2dsegm_first_value))

    def nucs_2dsegm_last(self):
        """Set the last of the frames on which calculate the threshold"""
        self.nucs_2dsegm_last_value  =  self.frame_dapi_raw.currentIndex
        self.nucs_2dsegm_last_lbl.setText(str(self.nucs_2dsegm_last_value))

    def sld_thr_ch1_update(self):
        """Set manual threshold for ch1 spots segmentation"""
        cif                =  self.frame_ch1_raw.currentIndex
        hh                 =  self.frame_ch1_sgm.view.viewRange()
        self.spts_ch1_sgm  =  (self.spts_ch1_flt > (self.sld_thr_ch1.value() * (self.spts_ch1_flt.max() - self.spts_ch1_flt.min()) / 100)).astype(np.uint8)
        self.frame_ch1_sgm.setImage(self.spts_ch1_sgm)
        self.frame_ch1_sgm.setCurrentIndex(cif)
        self.frame_ch1_sgm.view.setXRange(hh[0][0], hh[0][1], padding=.0002)
        self.frame_ch1_sgm.view.setYRange(hh[1][0], hh[1][1], padding=.0002)

        self.sld_ch1_lbl.setText(str(self.sld_thr_ch1.value()) + "%")

    def sld_thr_ch2_update(self):
        """Set manual threshold for ch2 spots segmentation"""
        cif                =  self.frame_ch2_raw.currentIndex
        hh                 =  self.frame_ch2_raw.view.viewRange()
        self.spts_ch2_sgm  =  (self.spts_ch2_flt > (self.sld_thr_ch2.value() * (self.spts_ch2_flt.max() - self.spts_ch2_flt.min()) / 100)).astype(np.uint8)
        self.frame_ch2_sgm.setImage(self.spts_ch2_sgm)
        self.frame_ch2_sgm.setCurrentIndex(cif)
        self.frame_ch2_sgm.view.setXRange(hh[0][0], hh[0][1], padding=.0002)
        self.frame_ch2_sgm.view.setYRange(hh[1][0], hh[1][1], padding=.0002)
        self.sld_ch2_lbl.setText(str(self.sld_thr_ch2.value()) + "%")

    def sld_thr_ch3_update(self):
        """Set manual threshold for ch3 spots segmentation"""
        cif                =  self.frame_ch3_raw.currentIndex
        hh                 =  self.frame_ch2_raw.view.viewRange()
        self.spts_ch3_sgm  =  (self.spts_ch3_flt > (self.sld_thr_ch3.value() * (self.spts_ch3_flt.max() - self.spts_ch3_flt.min()) / 100)).astype(np.uint8)
        self.frame_ch3_sgm.setImage(self.spts_ch3_sgm)
        self.frame_ch3_sgm.setCurrentIndex(cif)
        self.frame_ch3_sgm.view.setXRange(hh[0][0], hh[0][1], padding=.0002)
        self.frame_ch3_sgm.view.setYRange(hh[1][0], hh[1][1], padding=.0002)
        self.sld_ch3_lbl.setText(str(self.sld_thr_ch3.value()) + "%")

    def plus_dapi_thr(self):
        """Manual threshold increase for a single dapi z-frame"""
        cif  =  self.frame_dapi_2D.currentIndex
        hh   =  self.frame_dapi_2D.view.viewRange()

        self.nucs_2d_det.vals[self.frame_dapi_2D.currentIndex]     +=  self.nucs_2d_det.val_step
        self.nucs_2d_det.nucs_sgm[self.frame_dapi_2D.currentIndex]  =  (self.nucs_m[self.frame_dapi_2D.currentIndex] > self.nucs_2d_det.vals[self.frame_dapi_2D.currentIndex]).astype(np.uint8)
        self.nucs_2d_det.nucs_sgm[self.frame_dapi_2D.currentIndex]  =  binary_fill_holes(self.nucs_2d_det.nucs_sgm[self.frame_dapi_2D.currentIndex])
        # self.frame_dapi_2D.updateImage()
        self.frame_dapi_2D.setImage(self.nucs_2d_det.nucs_sgm)
        self.frame_dapi_2D.setCurrentIndex(cif)
        self.frame_dapi_2D.view.setXRange(hh[0][0], hh[0][1], padding=.0002)
        self.frame_dapi_2D.view.setYRange(hh[1][0], hh[1][1], padding=.0002)

    def minus_dapi_thr(self):
        """Manual threshold decrease for a single dapi z-frame"""
        cif  =  self.frame_dapi_2D.currentIndex
        hh   =  self.frame_dapi_2D.view.viewRange()
        self.nucs_2d_det.vals[self.frame_dapi_2D.currentIndex]     -=  self.nucs_2d_det.val_step
        self.nucs_2d_det.nucs_sgm[self.frame_dapi_2D.currentIndex]  =  (self.nucs_m[self.frame_dapi_2D.currentIndex] > self.nucs_2d_det.vals[self.frame_dapi_2D.currentIndex]).astype(np.uint8)
        self.nucs_2d_det.nucs_sgm[self.frame_dapi_2D.currentIndex]  =  binary_fill_holes(self.nucs_2d_det.nucs_sgm[self.frame_dapi_2D.currentIndex])
        # self.frame_dapi_2D.updateImage()
        self.frame_dapi_2D.setImage(self.nucs_2d_det.nucs_sgm)
        self.frame_dapi_2D.setCurrentIndex(cif)
        self.frame_dapi_2D.view.setXRange(hh[0][0], hh[0][1], padding=.0002)
        self.frame_dapi_2D.view.setYRange(hh[1][0], hh[1][1], padding=.0002)

    def ws_mindist_var(self, text):
        """Set the watershed parameter"""
        self.ws_mindist_value  =  np.int64(text)

    def load_data(self):
        """Loads raw data into the GUI, with pixel size info"""
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()
        try:
            self.spts_clusters_flag  =  SpotsOrClusters.getFlag()
            self.raw_data_fname      =  str(QtWidgets.QFileDialog.getOpenFileName(None, "Select raw data file to analyze", filter='*.czi')[0])
            self.raw_data            =  RawDataLoader.RawDataLoader(self.raw_data_fname)
            self.fname_raw_lbl.setText(self.raw_data_fname)

            self.frame_ch1_raw.setImage(self.raw_data.ch1)
            self.frame_ch2_raw.setImage(self.raw_data.ch2)
            self.frame_ch3_raw.setImage(self.raw_data.ch3)
            self.frame_dapi_raw.setImage(self.raw_data.dapi)

            self.pixsize_x_lbl.setText("pix size XY = " + str(np.round(self.raw_data.pix_sizeX, decimals=4)) + "µm;")
            self.pixsize_z_lbl.setText("Z step = " + str(np.round(self.raw_data.pix_sizeZ, decimals=4)) + "µm")

            if self.spts_clusters_flag  == "s":
                self.spts_ch1_flt  =  np.abs(gaussian_laplace(self.raw_data.ch1.astype(float), 1))
            elif self.spts_clusters_flag == "c":
                self.spts_ch1_flt  =  gaussian(self.raw_data.ch1.astype(float), 1)
            self.spts_ch2_flt  =  np.abs(gaussian_laplace(self.raw_data.ch2.astype(float), 1))
            self.spts_ch3_flt  =  np.abs(gaussian_laplace(self.raw_data.ch3.astype(float), 1))

            self.frame_ch1_flt.setImage(self.spts_ch1_flt)
            self.sld_thr_ch1_update()

            self.frame_ch2_flt.setImage(self.spts_ch2_flt)
            self.sld_thr_ch2_update()

            self.frame_ch3_flt.setImage(self.spts_ch3_flt)
            self.sld_thr_ch3_update()

            self.nucs_m  =  NucsSegmenter3.NucsSegmenter3Pre(self.raw_data.dapi).nucs_m

        except Exception:
            traceback.print_exc()

        self.ready_indicator()

    def change_sptsdet(self):
        if self.tabs_tot.currentIndex() == 0:
            self.spts_ch1_flt  =  np.abs(gaussian_filter(self.raw_data.ch1.astype(float), 1.2) - gaussian_filter(self.raw_data.ch1.astype(float), 2.4))
            self.sld_thr_ch1_update()
        elif self.tabs_tot.currentIndex() == 1:
            self.spts_ch2_flt  =  np.abs(gaussian_filter(self.raw_data.ch2.astype(float), 1.2) - gaussian_filter(self.raw_data.ch2.astype(float), 2.4))
            self.sld_thr_ch2_update()
        elif self.tabs_tot.currentIndex() == 2:
            self.spts_ch3_flt  =  np.abs(gaussian_filter(self.raw_data.ch3.astype(float), 1.2) - gaussian_filter(self.raw_data.ch3.astype(float), 2.4))
            self.sld_thr_ch3_update()

    def nucs_2dsegm(self):
        """Nuclei detection frame by frame"""
        reload(NucsSegmenter3)
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()
        try:
            self.nucs_2d_det  =  NucsSegmenter3.NucsSegmenter3(self.nucs_m, np.int64(self.nucs_2dsegm_first_value), np.int64(self.nucs_2dsegm_last_value))
            self.frame_dapi_2D.setImage(self.nucs_2d_det.nucs_sgm)
            self.nucs_ws_flag  =  0

        except Exception:
            traceback.print_exc()

        self.ready_indicator()

    def smooth_dapi(self):
        """Smooth the nuclei detection"""
        reload(NucsSegmenter3)
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()
        try:
            self.nucs_2d_smooth  =  NucsSegmenter3.NucsSegmenter3Smooth(self.nucs_2d_det.nucs_sgm, self.nucs_2dsegm_first_value, self.nucs_2dsegm_last_value).nucs_smooth
            cif                  =  self.frame_dapi_2D.currentIndex
            hh                   =  self.frame_dapi_2D.view.viewRange()
            self.frame_dapi_2D.setImage(self.nucs_2d_smooth)
            self.frame_dapi_2D.setCurrentIndex(cif)
            self.frame_dapi_2D.view.setXRange(hh[0][0], hh[0][1], padding=.0002)
            self.frame_dapi_2D.view.setYRange(hh[1][0], hh[1][1], padding=.0002)
        except Exception:
            traceback.print_exc()

        self.ready_indicator()

    def nucs_2dsegm_ws(self):
        """Nuclei detection frame by frame"""
        reload(NucsSegmenter3)
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()
        try:
            self.nucs_2d_sgm  =  NucsSegmenter3.NucsSegmenter3Split(self.nucs_2d_smooth, self.ws_mindist_value).nucs
            self.rnd_cmap     =  pg.ColorMap(np.linspace(0, 1, self.nucs_2d_sgm.max()), color=self.colors4map)
            self.frame_dapi_2D_sgm.setImage(self.nucs_2d_sgm)
            self.frame_dapi_2D_sgm.setColorMap(self.rnd_cmap)

        except Exception:
            traceback.print_exc()

        self.ready_indicator()

    def man_active(self, state):
        """Activate manual corrections for segmented nuclei"""
        if state == QtCore.Qt.Checked:
            self.man_active_flag  =  1
            self.manual_cut_btn.setEnabled(True)
        else:
            self.man_active_flag  =  0
            self.manual_cut_btn.setEnabled(False)

    def keyPressEvent(self, event):
        """Undo command for nuclei manual modification"""
        if event.key() == (QtCore.Qt.ControlModifier and Qt.Key_Z):
            cif                          =  self.frame_dapi_2D_sgm.currentIndex
            hh                           =  self.frame_dapi_2D_sgm.view.viewRange()
            self.nucs_2d_sgm[cif, :, :]  =  self.bufframe
            self.frame_dapi_2D_sgm.setImage(self.nucs_2d_sgm)
            self.frame_dapi_2D_sgm.setCurrentIndex(cif)
            self.frame_dapi_2D_sgm.view.setXRange(hh[0][0], hh[0][1], padding=.0002)
            self.frame_dapi_2D_sgm.view.setYRange(hh[1][0], hh[1][1], padding=.0002)

        if event.key() == (QtCore.Qt.ShiftModifier and Qt.Key_Delete):
            self.manual_cut()

    def click(self, event):
        """Put a segment on the nuclei image to cut or merge"""
        event.accept()
        pos        =  event.pos()
        modifiers  =  QtWidgets.QApplication.keyboardModifiers()

        if self.man_active_flag == 1:
            if modifiers  ==  QtCore.Qt.ShiftModifier:
                if self.c_count - 2 * np.int64(self.c_count / 2) == 0:
                    self.pos1  =  pos
                else:
                    try:
                        self.frame_dapi_2D_sgm.removeItem(self.roi_m)
                    except AttributeError:
                        pass

                    self.roi  =  pg.LineSegmentROI([self.pos1, pos], pen='r')
                    self.frame_dapi_2D_sgm.addItem(self.roi)

                self.c_count  +=  1

    def manual_cut(self):
        """Manual nuclei cutting"""
        cif      =  self.frame_dapi_2D_sgm.currentIndex
        hh       =  self.frame_dapi_2D_sgm.view.viewRange()
        pp       =  self.roi.getHandles()
        pp       =  [self.roi.mapToItem(self.frame_dapi_2D_sgm.imageItem, p.pos()) for p in pp]
        end_pts  =  np.array([[np.int64(pp[0].x()), np.int64(pp[0].y())], [np.int64(pp[1].x()), np.int64(pp[1].y())]])

        self.bufframe                =  np.copy(self.nucs_2d_sgm[cif, :, :])
        self.nucs_2d_sgm[cif, :, :]  =  LabelsModify.LabelsModify(self.nucs_2d_sgm[cif, :, :], end_pts).labels_fin
        self.frame_dapi_2D_sgm.setImage(self.nucs_2d_sgm)
        self.frame_dapi_2D_sgm.setCurrentIndex(cif)
        self.frame_dapi_2D_sgm.view.setXRange(hh[0][0], hh[0][1], padding=.0002)
        self.frame_dapi_2D_sgm.view.setYRange(hh[1][0], hh[1][1], padding=.0002)
        # self.frame_dapi_3D.removeItem(self.roi)

    def nucs2D_pileup(self):
        """Pile up segmented 2D nuclei"""
        reload(NucsPileUp)
        self.busy_indicator()
        QtWidgets.QApplication.processEvents()
        QtWidgets.QApplication.processEvents()
        try:

            self.nucs_3d_det  =  NucsPileUp.NucsPileUp(self.nucs_2d_sgm).nucs_lbls_piled
            self.nucs_3d_det  =  RemoveBorderNuclei.RemoveBorderNuclei(self.nucs_3d_det).nucs_dil
            self.frame_dapi_3D.setImage(self.nucs_3d_det)
            self.rnd_cmap     =  pg.ColorMap(np.linspace(0, 1, self.nucs_3d_det.max()), color=self.colors4map)
            self.frame_dapi_3D.setColorMap(self.rnd_cmap)

        except Exception:
            traceback.print_exc()

        self.ready_indicator()

    def save_analysis(self):
        """Save segmentation"""
        analysis_folder   =  QtWidgets.QFileDialog.getExistingDirectory(None, "Define Directory to save analysis")
        np.save(analysis_folder + '/nucs_dapi.npy', self.nucs_3d_det)
        np.save(analysis_folder + '/spts_ch1.npy', self.spts_ch1_sgm)
        np.save(analysis_folder + '/spts_ch2.npy', self.spts_ch2_sgm)
        np.save(analysis_folder + '/spts_ch3.npy', self.spts_ch3_sgm)
        np.save(analysis_folder + '/pix_sizes.npy', np.array([self.raw_data.pix_sizeZ, self.raw_data.pix_sizeX]))
        file  =  open(analysis_folder + '/spts_params_text.txt', "w")
        file.write("Control = " + str(self.sld_thr_ch1.value()))
        file.write("\n" + "Escargo = " + str(self.sld_thr_ch2.value()))
        file.write("\n" + "Snail = " + str(self.sld_thr_ch3.value()))
        file.close()

        if self.spts_clusters_flag == "s":
            np.save(analysis_folder + '/spts_clstrs_flag.npy', np.array([1]))
        elif self.spts_clusters_flag == "c":
            np.save(analysis_folder + '/spts_clstrs_flag.npy', np.array([2]))

    def cross_info(self):
        """Cross information of spots from all the channels and nuclei"""
        analysis_folder  =  QtWidgets.QFileDialog.getExistingDirectory(None, "Define Directory to save analysis")
        crsinfo          =  CrossInfo.CrossInfo(analysis_folder)
        self.pp          =  SingleNuc(crsinfo.mtx2show, crsinfo.nucs_dapi)
        self.pp.show()

    def check_overlaps(self):
        """Check overlapping by rotating and flipping channels one respect to the others"""
        analysis_folder  =  QtWidgets.QFileDialog.getExistingDirectory(None, "Define Directory to save analysis")
        # analysis_folder  =  '/home/atrullo/Desktop/Helene_Update'
        self.mpp1        =  CheckOverlapps(analysis_folder)
        self.mpp1.show()


class SingleNuc(QtWidgets.QWidget):
    """Pop up tool to show just one nucleus"""
    def __init__(self, mtx2show, nucs_dapi):
        """Widget for results visualization"""
        QtWidgets.QWidget.__init__(self)

        self.mtx2show      =  mtx2show
        self.nucs_dapi     =  nucs_dapi
        self.mtx_mdf2show  =  np.copy(self.mtx2show)

        frame  =  pg.ImageView()
        frame.setImage(self.mtx_mdf2show)
        frame.ui.roiBtn.hide()
        frame.ui.menuBtn.hide()

        nucs_lbl  =  QtWidgets.QLabel("Nuc Tag")
        nucs_lbl.setFixedSize(55, 25)

        nucs_edt  =  QtWidgets.QLineEdit(self)
        nucs_edt.setFixedSize(35, 25)
        nucs_edt.returnPressed.connect(self.nucs_show)
        nucs_edt.textChanged.connect(self.nucs_var)

        nucs_box  =  QtWidgets.QHBoxLayout()
        nucs_box.addStretch()
        nucs_box.addWidget(nucs_lbl)
        nucs_box.addWidget(nucs_edt)

        layout  =  QtWidgets.QVBoxLayout()
        layout.addWidget(frame)
        layout.addLayout(nucs_box)

        self.frame  =  frame

        self.setLayout(layout)
        self.setGeometry(300, 300, 800, 600)
        self.setWindowTitle("Visualize Spots and Nuclei")

    def nucs_show(self):
        """Show a single nucleus"""
        cif  =  self.frame.currentIndex
        hh   =  self.frame.view.viewRange()
        self.mtx_mdf2show  = np.copy(self.mtx2show)
        self.mtx_mdf2show[:, :, :, 0]  +=  (self.nucs_dapi == self.nucs_value) * np.uint8(50) * np.sign(self.mtx2show[:, :, :, 0]) * np.sign(self.mtx2show[:, :, :, 1]) * np.sign(self.mtx2show[:, :, :, 2])
        self.mtx_mdf2show[:, :, :, 1]  +=  (self.nucs_dapi == self.nucs_value) * np.uint8(50) * np.sign(self.mtx2show[:, :, :, 0]) * np.sign(self.mtx2show[:, :, :, 1]) * np.sign(self.mtx2show[:, :, :, 2])
        self.frame.setImage(self.mtx_mdf2show)
        self.frame.setCurrentIndex(cif)
        self.frame.view.setXRange(hh[0][0], hh[0][1], padding=.0002)
        self.frame.view.setYRange(hh[1][0], hh[1][1], padding=.0002)

    def nucs_var(self, text):
        """Define the tag of the nucleus to show"""
        self.nucs_value  =  np.int64(text)


class SpotsOrClusters(QtWidgets.QDialog):
    """Choose the spots channel to analyse"""
    def __init__(self, parent=None):
        super(SpotsOrClusters, self).__init__(parent)

        choose_lbl  =  QtWidgets.QLabel("Spots or Clusters?", self)
        choose_lbl.setFixedSize(120, 22)

        spots_btn  =  QtWidgets.QPushButton("Spots", self)
        spots_btn.setToolTip("Work with spots detection")
        spots_btn.setFixedSize(60, 25)
        spots_btn.clicked.connect(self.spots)

        clusters_btn  =  QtWidgets.QPushButton("Clusters", self)
        clusters_btn.setToolTip("Work with clusters detection")
        clusters_btn.setFixedSize(60, 25)
        clusters_btn.clicked.connect(self.clusters)

        choose_box  =  QtWidgets.QHBoxLayout()
        choose_box.addWidget(spots_btn)
        choose_box.addStretch()
        choose_box.addWidget(clusters_btn)

        layout  =  QtWidgets.QVBoxLayout()
        layout.addWidget(choose_lbl)
        layout.addLayout(choose_box)

        self.setWindowModality(Qt.ApplicationModal)
        self.setLayout(layout)
        self.setGeometry(300, 300, 200, 25)
        self.setWindowTitle("Spots or Clusters")

    def spots(self):
        """Choose flag a"""
        self.flag_s_c  =  "s"
        self.close()

    def clusters(self):
        """Choose flag b"""
        self.flag_s_c  =  "c"
        self.close()

    def params(self):
        """Function to send choice"""
        return self.flag_s_c

    @staticmethod
    def getFlag(parent=None):
        """Send choice"""
        dialog  =  SpotsOrClusters(parent)
        result  =  dialog.exec_()
        flag    =  dialog.params()
        return flag


class CheckOverlapps(QtWidgets.QWidget):
    """Pop up tool to show just one nucleus"""
    def __init__(self, analysis_folder):
        """Widget for results visualization"""
        QtWidgets.QWidget.__init__(self)

        spts_ch1              =  np.load(analysis_folder + '/spts_ch1.npy')                                     # load analysis results
        spts_ch2              =  np.load(analysis_folder + '/spts_ch2.npy')
        spts_ch3              =  np.load(analysis_folder + '/spts_ch3.npy')

        spts_chs_orig              =  np.zeros(spts_ch1.shape + (3,), dtype=np.uint16)
        spts_chs_orig[:, :, :, 0]  =  spts_ch1
        spts_chs_orig[:, :, :, 1]  =  spts_ch2
        spts_chs_orig[:, :, :, 2]  =  spts_ch3

        spts_chs_modif  =  np.copy(spts_chs_orig)

        fldr_name  =  QtWidgets.QLabel("Folder: " + analysis_folder, self)

        frame_numb  =  QtWidgets.QLabel("Frame: " + str(0), self)
        frame_numb.setFixedSize(90, 25)

        frame_orig  =  pg.ImageView(self, name='Original')
        frame_orig.ui.roiBtn.hide()
        frame_orig.ui.menuBtn.hide()
        frame_orig.setImage(spts_chs_orig)
        frame_orig.timeLine.sigPositionChanged.connect(self.frame_update_from_orig)

        frame_modif  =  pg.ImageView()
        frame_modif.ui.roiBtn.hide()
        frame_modif.ui.menuBtn.hide()
        frame_modif.setImage(spts_chs_modif)
        frame_modif.timeLine.sigPositionChanged.connect(self.frame_update_from_modif)
        frame_modif.view.setXLink("Original")
        frame_modif.view.setYLink("Original")

        rot_ch1_btn  =  QtWidgets.QPushButton("Rotate 90°", self)
        rot_ch1_btn.clicked.connect(self.rot_ch1)
        rot_ch1_btn.setToolTip("Rotate channel 1 90 degrees clockwise")
        rot_ch1_btn.setFixedSize(110, 25)

        flip_vert_ch1_btn  =  QtWidgets.QPushButton("Flip Vert", self)
        flip_vert_ch1_btn.clicked.connect(self.flip_vert_ch1)
        flip_vert_ch1_btn.setToolTip("Pile up nuclei slices in all the z frames")
        flip_vert_ch1_btn.setFixedSize(110, 25)

        flip_hor_ch1_btn  =  QtWidgets.QPushButton("Flip Hor", self)
        flip_hor_ch1_btn.clicked.connect(self.flip_hor_ch1)
        flip_hor_ch1_btn.setToolTip("Pile up nuclei slices in all the z frames")
        flip_hor_ch1_btn.setFixedSize(110, 25)

        spts_ch1_transf  =  QtWidgets.QVBoxLayout()
        spts_ch1_transf.addWidget(rot_ch1_btn)
        spts_ch1_transf.addWidget(flip_hor_ch1_btn)
        spts_ch1_transf.addWidget(flip_vert_ch1_btn)

        spts_ch1_trabsf_group  =  QtWidgets.QGroupBox("CH1")
        spts_ch1_trabsf_group.setLayout(spts_ch1_transf)
        spts_ch1_trabsf_group.setFixedSize(130, 120)

        rot_ch2_btn  =  QtWidgets.QPushButton("Rotate 90°", self)
        rot_ch2_btn.clicked.connect(self.rot_ch2)
        rot_ch2_btn.setToolTip("Rotate channel 2 90 degrees clockwise")
        rot_ch2_btn.setFixedSize(110, 25)

        flip_vert_ch2_btn  =  QtWidgets.QPushButton("Flip Vert", self)
        flip_vert_ch2_btn.clicked.connect(self.flip_vert_ch2)
        flip_vert_ch2_btn.setToolTip("Pile up nuclei slices in all the z frames")
        flip_vert_ch2_btn.setFixedSize(110, 25)

        flip_hor_ch2_btn  =  QtWidgets.QPushButton("Flip Hor", self)
        flip_hor_ch2_btn.clicked.connect(self.flip_hor_ch2)
        flip_hor_ch2_btn.setToolTip("Pile up nuclei slices in all the z frames")
        flip_hor_ch2_btn.setFixedSize(110, 25)

        spts_ch2_transf  =  QtWidgets.QVBoxLayout()
        spts_ch2_transf.addWidget(rot_ch2_btn)
        spts_ch2_transf.addWidget(flip_hor_ch2_btn)
        spts_ch2_transf.addWidget(flip_vert_ch2_btn)

        spts_ch2_transf_group  =  QtWidgets.QGroupBox("CH2")
        spts_ch2_transf_group.setLayout(spts_ch2_transf)
        spts_ch2_transf_group.setFixedSize(130, 120)

        rot_ch3_btn  =  QtWidgets.QPushButton("Rotate 90°", self)
        rot_ch3_btn.clicked.connect(self.rot_ch3)
        rot_ch3_btn.setToolTip("Rotate channel 3 90 degrees clockwise")
        rot_ch3_btn.setFixedSize(110, 25)

        flip_vert_ch3_btn  =  QtWidgets.QPushButton("Flip Vert", self)
        flip_vert_ch3_btn.clicked.connect(self.flip_vert_ch3)
        flip_vert_ch3_btn.setToolTip("Pile up nuclei slices in all the z frames")
        flip_vert_ch3_btn.setFixedSize(110, 25)

        flip_hor_ch3_btn  =  QtWidgets.QPushButton("Flip Hor", self)
        flip_hor_ch3_btn.clicked.connect(self.flip_hor_ch3)
        flip_hor_ch3_btn.setToolTip("Pile up nuclei slices in all the z frames")
        flip_hor_ch3_btn.setFixedSize(110, 25)

        spts_ch3_transf  =  QtWidgets.QVBoxLayout()
        spts_ch3_transf.addWidget(rot_ch3_btn)
        spts_ch3_transf.addWidget(flip_hor_ch3_btn)
        spts_ch3_transf.addWidget(flip_vert_ch3_btn)

        spts_ch3_transf_group  =  QtWidgets.QGroupBox("CH3")
        spts_ch3_transf_group.setLayout(spts_ch3_transf)
        spts_ch3_transf_group.setFixedSize(130, 120)

        reset_btn  =  QtWidgets.QPushButton("Reset", self)
        reset_btn.clicked.connect(self.reset)
        reset_btn.setToolTip("Reset all the changes")
        reset_btn.setFixedSize(130, 25)

        run_btn  =  QtWidgets.QPushButton("RUN", self)
        run_btn.clicked.connect(self.run)
        run_btn.setToolTip("Run the analysis on modified spots matrix")
        run_btn.setFixedSize(130, 25)

        commands_box  =  QtWidgets.QVBoxLayout()
        commands_box.addWidget(spts_ch1_trabsf_group)
        commands_box.addWidget(spts_ch2_transf_group)
        commands_box.addWidget(spts_ch3_transf_group)
        commands_box.addWidget(reset_btn)
        commands_box.addStretch()
        commands_box.addWidget(run_btn)
        commands_box.addWidget(frame_numb)

        layout_frames_cmmds  =  QtWidgets.QHBoxLayout()
        layout_frames_cmmds.addWidget(frame_orig)
        layout_frames_cmmds.addWidget(frame_modif)
        layout_frames_cmmds.addLayout(commands_box)

        layout  =  QtWidgets.QVBoxLayout()
        layout.addWidget(fldr_name)
        layout.addLayout(layout_frames_cmmds)

        self.frame_orig   =  frame_orig
        self.frame_modif  =  frame_modif
        self.frame_numb      =  frame_numb

        self.spts_chs_orig    =  spts_chs_orig
        self.spts_chs_modif   =  spts_chs_modif
        self.analysis_folder  =  analysis_folder

        self.setLayout(layout)
        self.setGeometry(300, 300, 800, 600)
        self.setWindowTitle("Visualize Spots and Nuclei")

    def frame_update_from_orig(self):
        """Update frame modif index when frame orig changes index"""
        self.frame_modif.setCurrentIndex(self.frame_orig.currentIndex)
        self.frame_numb.setText("Frame: " + str(self.frame_orig.currentIndex))

    def frame_update_from_modif(self):
        """Update frame modif index when frame orig changes index"""
        self.frame_orig.setCurrentIndex(self.frame_modif.currentIndex)

    def rot_ch1(self):
        """"""
        self.spts_chs_modif[:, :, :, 0]  =  np.rot90(self.spts_chs_modif[:, :, :, 0], axes=(1, 2))
        cif  =  self.frame_modif.currentIndex
        hh   =  self.frame_modif.view.viewRange()
        self.frame_modif.setImage(self.spts_chs_modif)
        self.frame_updater(cif, hh)

    def flip_vert_ch1(self):
        self.spts_chs_modif[:, :, :, 0]  =  self.spts_chs_modif[:, :, ::-1, 0]
        cif  =  self.frame_modif.currentIndex
        hh   =  self.frame_modif.view.viewRange()
        self.frame_modif.setImage(self.spts_chs_modif)
        self.frame_updater(cif, hh)

    def flip_hor_ch1(self):
        self.spts_chs_modif[:, :, :, 0]  =  self.spts_chs_modif[:, ::-1, :, 0]
        cif  =  self.frame_modif.currentIndex
        hh   =  self.frame_modif.view.viewRange()
        self.frame_modif.setImage(self.spts_chs_modif)
        self.frame_updater(cif, hh)

    def rot_ch2(self):
        self.spts_chs_modif[:, :, :, 1]  =  np.rot90(self.spts_chs_modif[:, :, :, 1], axes=(1, 2))
        cif  =  self.frame_modif.currentIndex
        hh   =  self.frame_modif.view.viewRange()
        self.frame_modif.setImage(self.spts_chs_modif)
        self.frame_updater(cif, hh)

    def flip_vert_ch2(self):
        self.spts_chs_modif[:, :, :, 1]  =  self.spts_chs_modif[:, :, ::-1, 1]
        cif  =  self.frame_modif.currentIndex
        hh   =  self.frame_modif.view.viewRange()
        self.frame_modif.setImage(self.spts_chs_modif)
        self.frame_updater(cif, hh)

    def flip_hor_ch2(self):
        self.spts_chs_modif[:, :, :, 1]  =  self.spts_chs_modif[:, ::-1, :, 1]
        cif  =  self.frame_modif.currentIndex
        hh   =  self.frame_modif.view.viewRange()
        self.frame_modif.setImage(self.spts_chs_modif)
        self.frame_updater(cif, hh)

    def rot_ch3(self):
        self.spts_chs_modif[:, :, :, 2]  =  np.rot90(self.spts_chs_modif[:, :, :, 2], axes=(1, 2))
        cif  =  self.frame_modif.currentIndex
        hh   =  self.frame_modif.view.viewRange()
        self.frame_modif.setImage(self.spts_chs_modif)
        self.frame_updater(cif, hh)

    def flip_vert_ch3(self):
        self.spts_chs_modif[:, :, :, 2]  =  self.spts_chs_modif[:, :, ::-1, 2]
        cif  =  self.frame_modif.currentIndex
        hh   =  self.frame_modif.view.viewRange()
        self.frame_modif.setImage(self.spts_chs_modif)
        self.frame_updater(cif, hh)

    def flip_hor_ch3(self):
        self.spts_chs_modif[:, :, :, 2]  =  self.spts_chs_modif[:, ::-1, :, 2]
        cif  =  self.frame_modif.currentIndex
        hh   =  self.frame_modif.view.viewRange()
        self.frame_modif.setImage(self.spts_chs_modif)
        self.frame_updater(cif, hh)

    def reset(self):
        self.spts_chs_modif  =  np.copy(self.spts_chs_orig)
        cif  =  self.frame_modif.currentIndex
        hh   =  self.frame_modif.view.viewRange()
        self.frame_modif.setImage(self.spts_chs_modif)
        self.frame_updater(cif, hh)

    def frame_updater(self, cif, hh):
        self.frame_modif.setCurrentIndex(cif)
        self.frame_modif.view.setXRange(hh[0][0], hh[0][1], padding=.0002)
        self.frame_modif.view.setYRange(hh[1][0], hh[1][1], padding=.0002)

    def run(self):
        crsinfo  =  ModifCrossInfo.CrossInfo(self.analysis_folder, self.spts_chs_modif)
        self.pp  =  SingleNuc(crsinfo.mtx2show, crsinfo.nucs_dapi)
        self.pp.show()


def main():
    app         =  QtWidgets.QApplication(sys.argv)
    splash_pix  =  QtGui.QPixmap('Icons/van.png')
    splash      =  QtWidgets.QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
    splash.setMask(splash_pix.mask())
    splash.show()
    app.processEvents()
    ex   =  MainWindow()
    splash.finish(ex)
    sys.exit(app.exec_())


def except_hook(cls, exception, traceback):
    sys.__excepthook__(cls, exception, traceback)


if __name__ == '__main__':

    main()
    sys.excepthook = except_hook
