# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'UI_file.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1125, 644)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.left_menu_frame = QtWidgets.QFrame(self.centralwidget)
        self.left_menu_frame.setGeometry(QtCore.QRect(20, 0, 170, 610))
        self.left_menu_frame.setAutoFillBackground(False)
        self.left_menu_frame.setStyleSheet("background-color: rgb(170, 170, 127)")
        self.left_menu_frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.left_menu_frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.left_menu_frame.setObjectName("left_menu_frame")
        self.widget = QtWidgets.QWidget(self.left_menu_frame)
        self.widget.setGeometry(QtCore.QRect(10, 40, 141, 135))
        self.widget.setObjectName("widget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.widget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setSpacing(0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.button_menu1 = QtWidgets.QPushButton(self.widget)
        self.button_menu1.setStyleSheet("QPushButton:hover{\n"
"    background-color:rgb(253, 255, 247)\n"
"}")
        self.button_menu1.setObjectName("button_menu1")
        self.verticalLayout.addWidget(self.button_menu1)
        self.button_menu_2 = QtWidgets.QPushButton(self.widget)
        self.button_menu_2.setStyleSheet("QPushButton:hover{\n"
"    background-color:rgb(253, 255, 247)\n"
"}")
        self.button_menu_2.setObjectName("button_menu_2")
        self.verticalLayout.addWidget(self.button_menu_2)
        self.button_menu_3 = QtWidgets.QPushButton(self.widget)
        self.button_menu_3.setStyleSheet("QPushButton:hover{\n"
"    background-color:rgb(253, 255, 247)\n"
"}")
        self.button_menu_3.setObjectName("button_menu_3")
        self.verticalLayout.addWidget(self.button_menu_3)
        self.button_menu_4 = QtWidgets.QPushButton(self.widget)
        self.button_menu_4.setStyleSheet("QPushButton:hover{\n"
"    background-color:rgb(253, 255, 247)\n"
"}")
        self.button_menu_4.setObjectName("button_menu_4")
        self.verticalLayout.addWidget(self.button_menu_4)
        self.stackedWidget = QtWidgets.QStackedWidget(self.left_menu_frame)
        self.stackedWidget.setGeometry(QtCore.QRect(10, 190, 141, 381))
        self.stackedWidget.setObjectName("stackedWidget")
        self.page_1 = QtWidgets.QWidget()
        self.page_1.setObjectName("page_1")
        self.label_1 = QtWidgets.QLabel(self.page_1)
        self.label_1.setGeometry(QtCore.QRect(60, 150, 55, 16))
        self.label_1.setObjectName("label_1")
        self.stackedWidget.addWidget(self.page_1)
        self.page_2 = QtWidgets.QWidget()
        self.page_2.setObjectName("page_2")
        self.label_2 = QtWidgets.QLabel(self.page_2)
        self.label_2.setGeometry(QtCore.QRect(40, 100, 61, 71))
        self.label_2.setObjectName("label_2")
        self.stackedWidget.addWidget(self.page_2)
        self.page_3 = QtWidgets.QWidget()
        self.page_3.setObjectName("page_3")
        self.label_3 = QtWidgets.QLabel(self.page_3)
        self.label_3.setGeometry(QtCore.QRect(40, 160, 55, 16))
        self.label_3.setObjectName("label_3")
        self.stackedWidget.addWidget(self.page_3)
        self.page_4 = QtWidgets.QWidget()
        self.page_4.setObjectName("page_4")
        self.label_4 = QtWidgets.QLabel(self.page_4)
        self.label_4.setGeometry(QtCore.QRect(40, 190, 55, 16))
        self.label_4.setObjectName("label_4")
        self.stackedWidget.addWidget(self.page_4)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1125, 26))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.button_menu1.setText(_translate("MainWindow", "Menu 1"))
        self.button_menu_2.setText(_translate("MainWindow", "Menu 2"))
        self.button_menu_3.setText(_translate("MainWindow", "Menu 3"))
        self.button_menu_4.setText(_translate("MainWindow", "Menu 4"))
        self.label_1.setText(_translate("MainWindow", "Page 1"))
        self.label_2.setText(_translate("MainWindow", "Page 2"))
        self.label_3.setText(_translate("MainWindow", "Page 3"))
        self.label_4.setText(_translate("MainWindow", "Page 4"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

