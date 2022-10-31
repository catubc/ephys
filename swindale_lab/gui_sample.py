import sys
from PyQt4 import QtGui, QtCore

class Window(QtGui.QMainWindow):

    def __init__(self):
        super(Window, self).__init__()
        self.setGeometry(50, 50, 500, 300)
        self.setWindowTitle("OpenNeuron")
        self.setWindowIcon(QtGui.QIcon('pythonlogo.png'))

        toolMenu = QtGui.QMenuBar()
        toolMenu.setNativeMenuBar(False) # <--Sets the menu with the widget; otherwise shows up as global (i.e. at top desktop screen)
        self.setMenuBar(toolMenu)

        loadExperiment = QtGui.QAction("&Load Experiment", self)
        loadExperiment.setStatusTip('Load Experiment')
        loadExperiment.triggered.connect(self.load_experiment)
        
        processExperiment = QtGui.QAction("&Process Experiment", self)
        processExperiment.setStatusTip('Process Experiment')
        processExperiment.triggered.connect(self.process_experiment)

        exitApplication = QtGui.QAction("&Exit Application", self)
        exitApplication.setStatusTip('Exit')
        exitApplication.triggered.connect(self.close_application)

        self.statusBar()

        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu('&Main')
        fileMenu.addAction(exitApplication)

        fileMenu = mainMenu.addMenu('&Experiment')
        fileMenu.addAction(loadExperiment)

        fileMenu = mainMenu.addMenu('&Analyze')
        #fileMenu.addAction(loadFile)
        
        self.home()

    def home(self):
        btn = QtGui.QPushButton("Quit", self)
        btn.clicked.connect(self.close_application)
        btn.resize(btn.minimumSizeHint())
        btn.move(0,100)
        self.show()


    def load_experiment(self):
        print("....loading experiment ...")
        

    def process_experiment(self):
        print("....processing experiment ...")
        
        
    def close_application(self):
        print("whooaaaa so custom!!!")
        sys.exit()


def run():
    app = QtGui.QApplication(sys.argv)
    GUI = Window()
    sys.exit(app.exec_())


run()
