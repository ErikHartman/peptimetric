from PyQt5.QtWidgets import QWidget, QApplication, QVBoxLayout
import sys
import graphistest
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from methods import*
from lists import *


class MatplotlibWindow(QWidget):
    def __init__(self, parent=None):
        super(MatplotlibWindow, self).__init__(parent)
        self.figure = Figure()
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.axis = self.figure.add_subplot(111)

        self.layoutvertical = QVBoxLayout(self)
        self.layoutvertical.addWidget(self.canvas)


class MainWidget(QWidget, graphistest.Ui_Form):
    def __init__(self):
        super(MainWidget, self).__init__()
        self.setupUi(self)
        self.init_widget()
        self.pushButton.clicked.connect(self.plot_widget())

    def init_widget(self):
        self.matplotlibwidget = MatplotlibWindow()
        self.layoutvertical = QVBoxLayout(self.widget)
        self.layoutvertical.addWidget(self.matplotlibwidget)

    def plot_widget(self):
        g1 = concatenate_dataframes(read_files_gui())
        g2 = concatenate_dataframes(read_files_gui())
        master = g1.merge(g2, on=['Peptide', 'Accession'], how='outer', suffixes=['_g1', '_g2'])
        protein_list = create_protein_list(master)
        self.matplotlibwidget.canvas.draw(create_graphic(protein_list))



if __name__ == '__main__':
    app = QApplication(sys.argv)
    w = MainWidget()
    w.show()
    sys.exit(app.exec_())
