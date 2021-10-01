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

from gui import input_parameters_dlg_ui
from gui.advanced_settings_dlg import AdvancedSettingsDialog

from PyQt5 import QtCore, QtGui, QtWidgets

if hasattr(QtCore.Qt, 'AA_EnableHighDpiScaling'):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)

if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)

class InputParametersDialog(QDialog):

    def __init__(self, *args, **kwargs):
        super(InputParametersDialog, self).__init__(*args, **kwargs)

        self.ui = input_parameters_dlg_ui.Ui_Dialog()
        self.ui.setupUi(self)

        # Advanced Dialog
        advanced_settings_dlg = AdvancedSettingsDialog()

        # Handle findXML buttons
        self.ui.find488.clicked.connect(partial(self.findXML,
                self.ui.pathLbl488))
        self.ui.find647.clicked.connect(partial(self.findXML,
                self.ui.pathLbl647))
        self.ui.find750.clicked.connect(partial(self.findXML,
                self.ui.pathLbl750))
        self.ui.find561.clicked.connect(partial(self.findXML,
                self.ui.pathLbl561))

        # Update button
        self.ui.updateBtn.clicked.connect(self.updateXML)

        # Creating dialog
        self.ui.advancedBtn.clicked.connect(self.open_advanced_settings_dlg)

        # Lists for updating all XMLs
        self.labels = [self.ui.pathLbl488, self.ui.pathLbl561,
                       self.ui.pathLbl647, self.ui.pathLbl750]
        self.threshs = [self.ui.thresh488, self.ui.thresh561,
                        self.ui.thresh647, self.ui.thresh750]

    def open_advanced_settings_dlg(self):
        advanced_settings_dlg = AdvancedSettingsDialog()
        advanced_settings_dlg.show()

    def findXML(self, label):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        (fileName, _) = QFileDialog.getOpenFileName(self,
                'QFileDialog.getOpenFileName()', '', 'XML Files (*.xml)'
                , options=options)
        if fileName:
            label.setText(fileName)

    def updateXML(self):

        tree = None
        counter = 0

        for (label, thresh) in zip(self.labels, self.threshs):
            file_path = label.text()
            print(file_path)
            if not file_path.endswith('xml'):
                pass
            else:
                # Load in the XML file
                if os.path.exists(file_path):
                    tree = lxmlET.parse(file_path)
                    root = tree.getroot()
                    counter += 1
                else:
                    self.error_dialog = QErrorMessage()
                    self.error_dialog.setWindowTitle('Error')
                    self.error_dialog.showMessage('Please enter valid file path!')
                    return

                # Second need to find necessary variables of interest that will be edited for each XML file
                for child in root:
                    if child.tag == 'start_frame':
                        child.text = str(self.ui.sf.value())

                    if child.tag == 'max_frame':
                        child.text = str(self.ui.ef.value())

                    if child.tag == 'threshold':
                        child.text = str(thresh.value())
                if tree:
                    tree.write(file_path)

        msg = QMessageBox()
        msg.setWindowTitle("Success")
        msg.setText('{} XMLs were updated!'.format(counter))
        msg.exec_()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = InputParametersDialog()
    ex.show()
    app.exec_()
    sys.exit()
