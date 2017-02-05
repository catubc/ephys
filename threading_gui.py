import sys
from PyQt4 import QtGui, QtCore
import time
import random


class MyThread(QtCore.QThread):
    trigger = QtCore.pyqtSignal(int)

    def __init__(self, parent=None):
        super(MyThread, self).__init__(parent)

    def setup(self, thread_no):
        self.thread_no = thread_no

    def run(self):
        time.sleep(random.random()*5)  # random sleep to imitate working
        self.trigger.emit(self.thread_no)


class Main(QtGui.QMainWindow):
    def __init__(self, parent=None):
        super(Main, self).__init__(parent)
        self.text_area = QtGui.QTextBrowser()
        self.thread_button = QtGui.QPushButton('Start threads')
        self.thread_button.clicked.connect(self.start_threads)

        central_widget = QtGui.QWidget()
        central_layout = QtGui.QHBoxLayout()
        central_layout.addWidget(self.text_area)
        central_layout.addWidget(self.thread_button)
        central_widget.setLayout(central_layout)
        self.setCentralWidget(central_widget)

    def start_threads(self):
        self.threads = []              # this will keep a reference to threads
        for i in range(10):
            thread = MyThread(self)    # create a thread
            thread.trigger.connect(self.update_text)  # connect to it's signal
            thread.setup(i)            # just setting up a parameter
            thread.start()             # start the thread
            self.threads.append(thread) # keep a reference

    def update_text(self, thread_no):
        self.text_area.append('thread # %d finished' % thread_no)

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)

    mainwindow = Main()
    mainwindow.show()

    sys.exit(app.exec_())
