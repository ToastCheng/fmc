#!/usr/bin/python3
# ~*~ coding: utf-8 ~*~

import sys
import os
import re
from pathlib import Path
import webbrowser
#webbrowser.open('http://example.com')

from PyQt5.QtWidgets import QMainWindow, QApplication, QLineEdit,\
    QPushButton,QLabel,QSlider, QFileDialog, QMessageBox
from PyQt5.QtCore import Qt

from parameter import calculate_mua_mus

PARA=['upA','upk','downA','downk','upg','downg',
      'StO2','CHb','Ccol','thick']
#format : [min,max,scale,precision]
LIMIT = {'upA':[0,100000,1,'%.f'],
         'upk':[0,1000,0.01,'%.2f'],
         'downA':[0,100000,1,'%.f'],
         'downk':[0,1000,0.01,'%.2f'],
         'upg':[0,100,0.01,'%.2f'],
         'downg':[0,100,0.01,'%.2f'],
         'StO2':[0,100,0.01,'%.2f'],
         'CHb':[0,10000,0.01,'%.2f'],
         'Ccol':[0,10000,0.01,'%.2f'],
         'thick':[0,1000,0.0001,'%.4f']
         }
#e.g.\n\
# \tsample_name upA upk downA downk upg downg StO2 CHb Ccol thick'

ERROR_TEXT = {
    'fileNameError':'Please enter file name',
    'autoFileNotImportError':'Auto Run needs to import data file',
    'zeroGError':'anisotropy factor g should greater than zero',



}



class Window(QMainWindow):

    def __init__(self,default=None):
        super().__init__()
        if default != None and len(default)==10:
            self.default = default
        else:
            self.default = ['0']*10
        self.para=[]
        self.paraSlide={}
        self.paraValue={}
        self.fname = None
        self.initUI()

    def initUI(self):

        statusbar = self.statusBar()
        statusbar.showMessage('created by TC,  BIOMEDICAL OPTICAL SPECTROSCOPY AND IMAGING LAB')


        #MENUBAR

        menubar = self.menuBar()
        menubar.setNativeMenuBar(False)
        menuFile = menubar.addMenu('&File')
        menuEdit = menubar.addMenu('&Edit')
        menuSimulate = menubar.addMenu('&Simulate')
        menuAbout = menubar.addMenu('&About')

        file_new = menuFile.addAction('New')
        file_new.triggered.connect(self.new)

        file_open = menuFile.addAction('Auto')
        file_open.triggered.connect(self.auto)

        file_close = menuFile.addAction('Close')
        file_close.triggered.connect(QApplication.instance().quit)

        menuEdit.addAction('^3^')

        simulate_run = menuSimulate.addAction('Run')
        simulate_run.triggered.connect(self.run)

        about_website = menuAbout.addAction('Website')
        about_website.triggered.connect(self.website)

        self.fileNameLabel = QLabel('output file name',self)
        self.fileNameLabel.setGeometry(250,33,100,30)

        self.fileName = QLineEdit('example',self)
        self.fileName.setGeometry(370,40,100,20)

        for i,j in enumerate(PARA):
            self.para.append(QLabel(j,self))
            self.para[i].setGeometry(10,33+25*i,100,30)

            self.paraSlide[j] = QSlider(Qt.Horizontal,self)
            self.paraSlide[j].setGeometry(60,40+25*i,100,30)
            self.paraSlide[j].valueChanged.connect(self.slideChange)
            self.paraSlide[j].setObjectName(j)
            self.paraSlide[j].setMinimum(LIMIT[j][0])
            self.paraSlide[j].setMaximum(LIMIT[j][1])

            self.paraValue[j] = QLineEdit(self.default[i],self)
            self.paraValue[j].setGeometry(170, 40+25*i, 50, 20)
            self.paraValue[j].textChanged.connect(self.textChange)
            self.paraValue[j].setObjectName(j)

        btn = QPushButton('Run',self)
        btn.move(390,290)
        btn.setStatusTip('simulate')
        btn.clicked.connect(self.run)

        auto_btn = QPushButton('Auto Run',self)
        auto_btn.move(390,265)
        auto_btn.setStatusTip('simulate using data in file')
        auto_btn.clicked.connect(self.autoRun)

        help_btn = QPushButton('Help',self)
        help_btn.move(390,240)
        help_btn.setStatusTip('help message')
        help_btn.clicked.connect(self.help)



        self.setFixedSize(500,350)
        self.setWindowTitle('Fluorescence Monte Carlo')
        self.show()

    def slideChange(self,value):

        key = self.sender().objectName()
        value = value*LIMIT[key][2]
        self.paraValue[key].setText(LIMIT[key][3] % value)

    def textChange(self,value):

        try:
            key = self.sender().objectName()
            value = float(value)/LIMIT[key][2]
            if value>=LIMIT[key][0] and value<=LIMIT[key][1]:
                self.paraSlide[key].setValue(value)
        except:
            pass

    def new(self):
        self.__init__()

    def help(self):
        msg = QMessageBox()
        msg.setFixedSize(800,400)
        msg.setInformativeText(
"""Fluorescence Monte Carlo(FMC) is an GPU based simulation program,
Nvidia graphic card, and CUDA environment is needed for simulation.

To run simulation, you can:

1. set the parameters by GUI, and press the 'Run' button.

2. store the parameters in the form given below, and open with 
File->Auto, then press the 'Auto Run' button. It will iterate 
all the case in the input file.

[format]
sample_name1, upA, upk, downA, downk, upg, downg, StO2, CHb, Ccol, thick
sample_name2, upA, upk, downA, downk, upg, downg, StO2, CHb, Ccol, thick
"""
        )
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()

    def auto(self):
        home = str(Path.home())
        self.fname = QFileDialog.getOpenFileName(self, 'Open file',
                                            home)
    def autoRun(self):
        if self.fname == None:
            self.warning('autoFileNotImportError')
            return
        with open(self.fname[0],'r') as f:
            for line in f.readlines():
                line = re.split('\t|\n|,',line)[:-1]
                self.fileName.setText(line.pop(0))

                for key,value in zip(PARA,line):
                    self.paraValue[key].setText(value)

                self.run()

    def website(self):
        webbrowser.open('http://bosi.ee.ntu.edu.tw/')


    def warning(self,errorValue):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Warning)
        msg.setText(ERROR_TEXT[errorValue])
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()


    def run(self):
        name = self.fileName.text()
        # print name
        if name == '':
            self.warning('fileNameError')
            return

        para = [float(j.text()) for i,j in self.paraValue.items()]
        if para[4]==0 or para[5]==0:
            self.warning('zeroGError')
            return

        wavelength,parameters = calculate_mua_mus(para[0],para[1],para[2],
                                      para[3],para[4],para[5],
                                      para[6],para[7],para[8],
                                      para[9])

        with open('DRSresult/'+name + '.txt', 'w+') as output:
            output.write('%f\n' % parameters[4])
            for i in range(len(wavelength)):
                output.write('%d\t%f\t%f\t%f\t%f\n' % (wavelength[i],
                parameters[0][i],parameters[1][i],parameters[2][i],
                parameters[3][i]))

        os.system('mkdir ' + name)
        os.system('./mcf ' + name)





if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = Window()
    sys.exit(app.exec_())
