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

from .frame_range_dlg_ui import * 
# import frame_range_dlg_ui 

class FrameRangeDialog(QDialog):

    def __init__(self, *args, **kwargs):
        super(FrameRangeDialog, self).__init__(*args, **kwargs)

        self.ui = Ui_Dialog()
        self.ui.setupUi(self)

        # Handle findXML buttons
        self.sfs = [self.ui.sf_spinBox_488, self.ui.sf_spinBox_561,
                       self.ui.sf_spinBox_647, self.ui.sf_spinBox_750]
        self.efs = [self.ui.ef_spinBox_488, self.ui.ef_spinBox_561,
                        self.ui.ef_spinBox_647, self.ui.ef_spinBox_750]

        self.ui.okBtn.clicked.connect(self.send_params)
        self.ui.cancelBtn.clicked.connect(self.cancel)
        
        with open(r'C:\Users\Vatsal\QT_Projects\STORM_GUI\configs\XMLs.yml') as f:
            self.config = yaml.load(f)

  
    def send_params(self):

        sfs = []
        efs = []

        for sf, ef, channel in zip(self.sfs, self.efs, [488, 561, 647, 750]):

            file_path = self.config[channel]
            if os.path.exists(file_path):
                tree = lxmlET.parse(file_path)
                root = tree.getroot()
            else:
                self.error_dialog = QErrorMessage()
                self.error_dialog.setWindowTitle('Error')
                self.error_dialog.showMessage('Please enter valid XML path for channel {}!'.format(channel))
                return

            # Second need to find necessary variables of interest that will be edited for each XML file
            for child in root:
                if child.tag == "start_frame":
                    child.text = str(sf.value())
                if child.tag == "max_frame":
                    child.text = str(ef.value())
            if tree:
                tree.write(file_path)

        msg = QMessageBox()
        msg.setWindowTitle("Success")
        msg.setText('Frame ranges are set!')
        msg.exec_()
        
        self.close() 

        return (sfs, efs)

    def cancel(self):
        self.close()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = FrameRangeDialog()
    ex.show()
    app.exec_()
    sys.exit()
