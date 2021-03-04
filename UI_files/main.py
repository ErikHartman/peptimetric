import sys
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5 import QtGui, QtCore
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton

from ui_main import Ui_MainWindow


class MainWindow(QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        self.ui.button_menu_1.clicked.connect(lambda: self.ui.stackedWidget.setCurrentWidget(self.ui.page_1))
        self.ui.button_menu_2.clicked.connect(lambda: self.ui.stackedWidget.setCurrentWidget(self.ui.page_2))
        self.ui.button_menu_3.clicked.connect(lambda: self.ui.stackedWidget.setCurrentWidget(self.ui.page_3))
        self.ui.button_menu_4.clicked.connect(lambda: self.ui.stackedWidget.setCurrentWidget(self.ui.page_4))

        self.show()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    sys.exit(app.exec_())
