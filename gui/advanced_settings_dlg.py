#!/usr/bin/python
# -*- coding: utf-8 -*-
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *

import os
import sys
from functools import partial

import xml.etree.ElementTree as ET
import lxml.etree as lxmlET
import yaml 

from .advanced_settings_dlg_ui import * 

from PyQt5 import QtCore, QtGui, QtWidgets

if hasattr(QtCore.Qt, 'AA_EnableHighDpiScaling'):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)

if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)

class AdvancedSettingsDialog(QDialog):

    def __init__(self, *args, **kwargs):
        super(AdvancedSettingsDialog, self).__init__(*args, **kwargs)

        self.ui = Ui_Dialog()
        self.ui.setupUi(self)

        # Handle findXML buttons
        self.ui.find488.clicked.connect(partial(self.findXML,
                self.ui.pathLbl488, self.ui.thresh488, 0, None))
        self.ui.find561.clicked.connect(partial(self.findXML,
                self.ui.pathLbl561, self.ui.thresh561, 1, None))
        self.ui.find647.clicked.connect(partial(self.findXML,
                self.ui.pathLbl647, self.ui.thresh647, 2, None))
        self.ui.find750.clicked.connect(partial(self.findXML,
                self.ui.pathLbl750, self.ui.thresh750, 3, None))

        self.ui.updateBtn.clicked.connect(self.updateXML)

        self.labels = [self.ui.pathLbl488, self.ui.pathLbl561,
                       self.ui.pathLbl647, self.ui.pathLbl750]
        self.threshs = [self.ui.thresh488, self.ui.thresh561,
                        self.ui.thresh647, self.ui.thresh750]
        self.sfs = [self.ui.sf488, self.ui.sf561, self.ui.sf647, self.ui.sf750]
        self.efs = [self.ui.ef488, self.ui.ef561, self.ui.ef647, self.ui.ef750]
        self.bss = [self.ui.bs488, self.ui.bs561, self.ui.bs647, self.ui.bs750]
        self.fss = [self.ui.fs488, self.ui.fs561, self.ui.fs647, self.ui.fs750]
        self.psfs = [self.ui.psf488, self.ui.psf561, self.ui.psf647, self.ui.psf750]
        self.checkBoxes = [self.ui.checkBox488, self.ui.checkBox561, self.ui.checkBox647, self.ui.checkBox750]
        
        # Auto-Load 
        with open(r'C:\Users\Vatsal\QT_Projects\STORM_GUI\configs\XMLs.yml') as f:
            self.config = yaml.load(f)
        self.findXML(self.ui.pathLbl488, self.ui.thresh488, 0, self.config[488])
        self.findXML(self.ui.pathLbl561, self.ui.thresh561, 1, self.config[561])
        self.findXML(self.ui.pathLbl647, self.ui.thresh647, 2, self.config[647])
        self.findXML(self.ui.pathLbl750, self.ui.thresh750, 3, self.config[750])
          

    def findXML(self, label, thresh, channel, fileName=None):
        
        print('pressed') 
        if fileName is None:
            print('NONE') 
            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            (fileName, _) = QFileDialog.getOpenFileName(self,
                    'QFileDialog.getOpenFileName()', '', 'XML Files (*.xml)'
                    , options=options)

        if fileName:
            print(fileName)
            label.setText(fileName)
            
            # Load in the XML tree 
            if os.path.exists(fileName):
                tree = lxmlET.parse(fileName)
                root = tree.getroot()
            else:
                self.error_dialog = QErrorMessage()
                self.error_dialog.setWindowTitle('Error')
                self.error_dialog.showMessage('Please enter valid file path!')
                return

            sf = self.sfs[channel]
            ef = self.efs[channel]
            bs = self.bss[channel]
            fs = self.fss[channel]
            psf = self.psfs[channel] 
            cb = self.checkBoxes[channel]
            
            for child in root:
                if child.tag == "start_frame":
                    sf.setValue(float(child.text))
                if child.tag == "max_frame":
                    ef.setValue(float(child.text))
                if child.tag == "background_sigma":
                    bs.setValue(float(child.text))
                if child.tag == "foreground_sigma":
                    fs.setValue(float(child.text))
                if child.tag == "threshold":
                    thresh.setValue(float(child.text))
                if child.tag == "sigma":
                    psf.setValue(float(child.text))
                if child.tag == "drift_correction":
                    if child.text == str(1):
                        cb.setChecked(True)

    def updateXML(self):

        tree = None
        gui_objs = zip(self.labels, self.threshs, self.sfs, self.efs, self.bss, self.fss, self.psfs, self.checkBoxes)

        counter = 0

        for (label, thresh, sf, ef, bs, fs, psf, cb) in gui_objs:
            file_path = label.text()
            print(file_path)
            if not file_path.endswith('xml'):
                pass
            else:
                # Load in the XML file
                if os.path.exists(file_path):
                    counter += 1
                    tree = lxmlET.parse(file_path)
                    root = tree.getroot()
                else:
                    self.error_dialog = QErrorMessage()
                    self.error_dialog.setWindowTitle('Error')
                    self.error_dialog.showMessage('Please enter valid file path!')
                    return

                # Second need to find necessary variables of interest that will be edited for each XML file
                for child in root:
                    if child.tag == "start_frame":
                        child.text = str(sf.value())
                    if child.tag == "max_frame":
                        child.text = str(ef.value())
                    if child.tag == "background_sigma":
                        child.text = str(bs.value())
                    if child.tag == "foreground_sigma":
                        child.text = str(fs.value())
                    if child.tag == "threshold":
                        child.text = str(thresh.value())
                    if child.tag == "sigma":
                        child.text = str(psf.value())
                    if child.tag == "drift_correction":
                        if cb.isChecked():
                            child.text = str(1)
                        else:
                            child.text = str(0)
                if tree:
                    tree.write(file_path)

        msg = QMessageBox()
        msg.setWindowTitle("Success")
        msg.setText('{} XMLs were updated!'.format(counter))
        msg.exec_()
        
        self.close()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = AdvancedSettingsDialog()
    ex.show()
    app.exec_()
    sys.exit()
