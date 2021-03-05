# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui_file.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1607, 898)
        MainWindow.setStyleSheet("QPushButton {\n"
"    color:rgb(0, 0, 0);\n"
"    background-color:rgb(190, 190, 190);\n"
"    border: 0px solid red\n"
"}\n"
"QPushButton:hover{\n"
"    background-color:rgb(212, 232, 236)\n"
"}\n"
"")
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.left_menu_frame = QtWidgets.QFrame(self.centralwidget)
        self.left_menu_frame.setGeometry(QtCore.QRect(0, 0, 301, 871))
        self.left_menu_frame.setAutoFillBackground(False)
        self.left_menu_frame.setStyleSheet("background-color: rgb(226, 226, 223);")
        self.left_menu_frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.left_menu_frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.left_menu_frame.setObjectName("left_menu_frame")
        self.widget = QtWidgets.QWidget(self.left_menu_frame)
        self.widget.setGeometry(QtCore.QRect(0, 60, 301, 151))
        self.widget.setObjectName("widget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.widget)
        self.verticalLayout.setContentsMargins(10, 0, 10, 0)
        self.verticalLayout.setSpacing(0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.button_menu_1 = QtWidgets.QPushButton(self.widget)
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(12)
        self.button_menu_1.setFont(font)
        self.button_menu_1.setStyleSheet("QPushButton {\n"
"    color:rgb(0, 0, 0);\n"
"    background-color: rgb(210, 240, 209);;\n"
"    border: 0px solid red\n"
"}\n"
"QPushButton:hover{\n"
"    color:rbg(264,264,264);\n"
"    background-color:rgb(212, 232, 236)\n"
"}\n"
"")
        self.button_menu_1.setObjectName("button_menu_1")
        self.verticalLayout.addWidget(self.button_menu_1)
        self.button_menu_2 = QtWidgets.QPushButton(self.widget)
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(12)
        self.button_menu_2.setFont(font)
        self.button_menu_2.setStyleSheet("QPushButton {\n"
"    color:rgb(0, 0, 0);\n"
"    background-color: rgb(210, 240, 209);;\n"
"    border: 0px solid red\n"
"}\n"
"QPushButton:hover{\n"
"    color:rbg(264,264,264);\n"
"    background-color:rgb(212, 232, 236)\n"
"}\n"
"")
        self.button_menu_2.setObjectName("button_menu_2")
        self.verticalLayout.addWidget(self.button_menu_2)
        self.button_menu_3 = QtWidgets.QPushButton(self.widget)
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(12)
        self.button_menu_3.setFont(font)
        self.button_menu_3.setStyleSheet("QPushButton {\n"
"    color:rgb(0, 0, 0);\n"
"    background-color: rgb(210, 240, 209);;\n"
"    border: 0px solid red\n"
"}\n"
"QPushButton:hover{\n"
"    color:rbg(264,264,264);\n"
"    background-color:rgb(212, 232, 236)\n"
"}\n"
"")
        self.button_menu_3.setObjectName("button_menu_3")
        self.verticalLayout.addWidget(self.button_menu_3)
        self.button_menu_4 = QtWidgets.QPushButton(self.widget)
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(12)
        self.button_menu_4.setFont(font)
        self.button_menu_4.setStyleSheet("QPushButton {\n"
"    color:rgb(0, 0, 0);\n"
"    background-color: rgb(210, 240, 209);;\n"
"    border: 0px solid red\n"
"}\n"
"QPushButton:hover{\n"
"    color:rbg(264,264,264);\n"
"    background-color:rgb(212, 232, 236)\n"
"}\n"
"")
        self.button_menu_4.setObjectName("button_menu_4")
        self.verticalLayout.addWidget(self.button_menu_4)
        self.stackedWidget = QtWidgets.QStackedWidget(self.left_menu_frame)
        self.stackedWidget.setGeometry(QtCore.QRect(10, 220, 281, 621))
        self.stackedWidget.setStyleSheet("background-color: rgb(220, 220, 220);")
        self.stackedWidget.setObjectName("stackedWidget")
        self.page_1 = QtWidgets.QWidget()
        self.page_1.setObjectName("page_1")
        self.label_1 = QtWidgets.QLabel(self.page_1)
        self.label_1.setGeometry(QtCore.QRect(90, 300, 55, 16))
        self.label_1.setObjectName("label_1")
        self.stackedWidget.addWidget(self.page_1)
        self.page_2 = QtWidgets.QWidget()
        self.page_2.setObjectName("page_2")
        self.label_2 = QtWidgets.QLabel(self.page_2)
        self.label_2.setGeometry(QtCore.QRect(80, 270, 131, 91))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_2.sizePolicy().hasHeightForWidth())
        self.label_2.setSizePolicy(sizePolicy)
        self.label_2.setObjectName("label_2")
        self.stackedWidget.addWidget(self.page_2)
        self.page_3 = QtWidgets.QWidget()
        self.page_3.setObjectName("page_3")
        self.label_3 = QtWidgets.QLabel(self.page_3)
        self.label_3.setGeometry(QtCore.QRect(80, 310, 55, 16))
        self.label_3.setObjectName("label_3")
        self.stackedWidget.addWidget(self.page_3)
        self.page_4 = QtWidgets.QWidget()
        self.page_4.setObjectName("page_4")
        self.label_4 = QtWidgets.QLabel(self.page_4)
        self.label_4.setGeometry(QtCore.QRect(80, 350, 55, 16))
        self.label_4.setObjectName("label_4")
        self.stackedWidget.addWidget(self.page_4)
        self.top_bar = QtWidgets.QFrame(self.centralwidget)
        self.top_bar.setGeometry(QtCore.QRect(-1, -1, 1611, 45))
        self.top_bar.setStyleSheet("background-color: rgb(226, 226, 223);")
        self.top_bar.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.top_bar.setFrameShadow(QtWidgets.QFrame.Raised)
        self.top_bar.setObjectName("top_bar")
        self.button_FAQ = QtWidgets.QPushButton(self.top_bar)
        self.button_FAQ.setEnabled(True)
        self.button_FAQ.setGeometry(QtCore.QRect(1230, 10, 70, 30))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(12)
        self.button_FAQ.setFont(font)
        self.button_FAQ.setStyleSheet("QPushButton {\n"
"    color:rgb(0, 0, 0);\n"
"    background-color:rgb(200, 200, 200);\n"
"    border: 0px solid red\n"
"}\n"
"QPushButton:hover{\n"
"    background-color:rgb(212, 232, 236)\n"
"}\n"
"")
        self.button_FAQ.setObjectName("button_FAQ")
        self.button_documentation = QtWidgets.QPushButton(self.top_bar)
        self.button_documentation.setGeometry(QtCore.QRect(1310, 10, 160, 30))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(12)
        self.button_documentation.setFont(font)
        self.button_documentation.setStyleSheet("QPushButton {\n"
"    color:rgb(0, 0, 0);\n"
"    background-color:rgb(200, 200, 200);\n"
"    border: 0px solid red\n"
"}\n"
"QPushButton:hover{\n"
"    background-color:rgb(212, 232, 236)\n"
"}\n"
"")
        self.button_documentation.setObjectName("button_documentation")
        self.button_feedback = QtWidgets.QPushButton(self.top_bar)
        self.button_feedback.setGeometry(QtCore.QRect(1480, 10, 111, 30))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(12)
        self.button_feedback.setFont(font)
        self.button_feedback.setStyleSheet("QPushButton {\n"
"    color:rgb(0, 0, 0);\n"
"    background-color:rgb(200, 200, 200);\n"
"    border: 0px solid red\n"
"}\n"
"QPushButton:hover{\n"
"    background-color:rgb(212, 232, 236)\n"
"}\n"
"")
        self.button_feedback.setObjectName("button_feedback")
        self.peptide_table = QtWidgets.QTableWidget(self.centralwidget)
        self.peptide_table.setGeometry(QtCore.QRect(1250, 530, 320, 300))
        self.peptide_table.setObjectName("peptide_table")
        self.peptide_table.setColumnCount(0)
        self.peptide_table.setRowCount(0)
        self.protein_widget = QtWidgets.QWidget(self.centralwidget)
        self.protein_widget.setGeometry(QtCore.QRect(380, 120, 850, 350))
        self.protein_widget.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.protein_widget.setObjectName("protein_widget")
        self.peptide_widget = QtWidgets.QWidget(self.centralwidget)
        self.peptide_widget.setGeometry(QtCore.QRect(380, 530, 850, 300))
        self.peptide_widget.setStyleSheet("background-color:rgb(255, 255, 255)")
        self.peptide_widget.setObjectName("peptide_widget")
        self.protein_table = QtWidgets.QTableWidget(self.centralwidget)
        self.protein_table.setGeometry(QtCore.QRect(1250, 120, 320, 300))
        self.protein_table.setObjectName("protein_table")
        self.protein_table.setColumnCount(0)
        self.protein_table.setRowCount(0)
        self.search_bar = QtWidgets.QLineEdit(self.centralwidget)
        self.search_bar.setGeometry(QtCore.QRect(1030, 80, 200, 30))
        self.search_bar.setObjectName("search_bar")
        self.label_protein_info = QtWidgets.QLabel(self.centralwidget)
        self.label_protein_info.setGeometry(QtCore.QRect(1340, 90, 130, 30))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(12)
        self.label_protein_info.setFont(font)
        self.label_protein_info.setStyleSheet("color:rgb(0, 0, 0);\n"
"background-color: rgb(210, 240, 209);;\n"
"border: 0px solid red\n"
"\n"
"\n"
"")
        self.label_protein_info.setAlignment(QtCore.Qt.AlignCenter)
        self.label_protein_info.setObjectName("label_protein_info")
        self.label_peptide_table = QtWidgets.QLabel(self.centralwidget)
        self.label_peptide_table.setGeometry(QtCore.QRect(1330, 500, 130, 30))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(12)
        self.label_peptide_table.setFont(font)
        self.label_peptide_table.setStyleSheet("color:rgb(0, 0, 0);\n"
"background-color: rgb(210, 240, 209);;\n"
"border: 0px solid red\n"
"")
        self.label_peptide_table.setAlignment(QtCore.Qt.AlignCenter)
        self.label_peptide_table.setObjectName("label_peptide_table")
        self.label_protein_view = QtWidgets.QLabel(self.centralwidget)
        self.label_protein_view.setGeometry(QtCore.QRect(380, 90, 150, 30))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(14)
        font.setBold(True)
        font.setWeight(75)
        self.label_protein_view.setFont(font)
        self.label_protein_view.setStyleSheet("background-color: rgb(210, 240, 209);\n"
"border: 0px solid red;")
        self.label_protein_view.setAlignment(QtCore.Qt.AlignCenter)
        self.label_protein_view.setObjectName("label_protein_view")
        self.label_peptide_view = QtWidgets.QLabel(self.centralwidget)
        self.label_peptide_view.setGeometry(QtCore.QRect(380, 500, 150, 30))
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(14)
        font.setBold(True)
        font.setWeight(75)
        self.label_peptide_view.setFont(font)
        self.label_peptide_view.setStyleSheet("background-color: rgb(210, 240, 209);\n"
"border: 0px solid red;")
        self.label_peptide_view.setAlignment(QtCore.Qt.AlignCenter)
        self.label_peptide_view.setObjectName("label_peptide_view")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1607, 26))
        self.menubar.setObjectName("menubar")
        self.menuFiles = QtWidgets.QMenu(self.menubar)
        self.menuFiles.setObjectName("menuFiles")
        self.menuExport_data = QtWidgets.QMenu(self.menubar)
        self.menuExport_data.setObjectName("menuExport_data")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionAdd_files = QtWidgets.QAction(MainWindow)
        self.actionAdd_files.setObjectName("actionAdd_files")
        self.actionExport_as = QtWidgets.QAction(MainWindow)
        self.actionExport_as.setObjectName("actionExport_as")
        self.actionView_files = QtWidgets.QAction(MainWindow)
        self.actionView_files.setObjectName("actionView_files")
        self.menuFiles.addAction(self.actionAdd_files)
        self.menuFiles.addAction(self.actionView_files)
        self.menuExport_data.addAction(self.actionExport_as)
        self.menubar.addAction(self.menuFiles.menuAction())
        self.menubar.addAction(self.menuExport_data.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.button_menu_1.setText(_translate("MainWindow", "Menu 1"))
        self.button_menu_2.setText(_translate("MainWindow", "Menu 2"))
        self.button_menu_3.setText(_translate("MainWindow", "Menu 3"))
        self.button_menu_4.setText(_translate("MainWindow", "Menu 4"))
        self.label_1.setText(_translate("MainWindow", "Page 1"))
        self.label_2.setText(_translate("MainWindow", "Page 2"))
        self.label_3.setText(_translate("MainWindow", "Page 3"))
        self.label_4.setText(_translate("MainWindow", "Page 4"))
        self.button_FAQ.setText(_translate("MainWindow", "FAQ"))
        self.button_documentation.setText(_translate("MainWindow", "Documentation"))
        self.button_feedback.setText(_translate("MainWindow", "Feedback"))
        self.search_bar.setText(_translate("MainWindow", "Search for proteins..."))
        self.label_protein_info.setText(_translate("MainWindow", "Protein info"))
        self.label_peptide_table.setText(_translate("MainWindow", "Peptide table"))
        self.label_protein_view.setText(_translate("MainWindow", "Protein View"))
        self.label_peptide_view.setText(_translate("MainWindow", "Peptide View"))
        self.menuFiles.setTitle(_translate("MainWindow", "Files"))
        self.menuExport_data.setTitle(_translate("MainWindow", "Export data"))
        self.actionAdd_files.setText(_translate("MainWindow", "Add files"))
        self.actionExport_as.setText(_translate("MainWindow", "Export as..."))
        self.actionView_files.setText(_translate("MainWindow", "View files"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

