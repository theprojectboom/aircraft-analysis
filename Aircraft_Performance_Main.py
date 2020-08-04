# Project Boom:  Supersonic UAV
# Turbine-Powered UAV Performance Assessment
# Created & Maintained by Johnathan M. Burgess
# Summer 2020


import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import QDialog, QApplication, QFileDialog, QMessageBox
from PyQt5.QtGui import QCursor
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from Aircraft_Performance_GUI import Ui_Dialog
from Boom_Library import BoomLibrary


class PlotCanvas(FigureCanvas):
    def __init__(self, parent, width=None, height=None, dpi=100):
        if width is None:
            width = parent.width()/100
        if height is None:
            height = parent.height()/100
        fig = Figure(figsize=(width, height), dpi=dpi)
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)
        self.ax = None

    def plot_thrust(self, x, x_sub, x_ts, x_ss, y_sub, y_ts, y_ss, z, a, b):
        x_min, x_max, x_step = 0, 1312*0.3048, 25
        y_min, y_max, y_step = 0, 600, 100
        self.ax = self.figure.add_subplot(111)
        self.ax.plot(x_sub*0.3048, y_sub*4.448222, 'blue', label='Thrust Required')
        self.ax.plot(x_ts*0.3048, y_ts*4.448222, 'blue')
        self.ax.plot(x_ss*0.3048, y_ss*4.448222, 'blue')
        self.ax.plot(a*0.3048, b*4.448222, 'black', marker='o', label='Mach 1.0')
        self.ax.plot(x*0.3048, z*4.448222, 'red', linestyle='--', label='Target Thrust')
        self.ax.set_title('Thrust Required vs. Thrust Available\n')
        self.ax.set_xlabel('Velocity - [ m/s ]')
        self.ax.set_ylabel('Thrust - [ N ]')
        self.ax.set_xlim(x_min, x_max)
        self.ax.xaxis.set_ticks(np.arange(x_min, x_max, x_step))
        self.ax.set_ylim(y_min, y_max)
        self.ax.yaxis.set_ticks(np.arange(y_min, y_max, y_step))
        self.ax.grid(True)
        self.ax.legend(shadow=True, ncol=1)
        self.draw()

    def plot_polar(self, x, y):
        x_min, x_max, x_step = 0, 0.28, 0.02
        y_min, y_max, y_step = 0, 2.2, 0.2
        self.ax = self.figure.add_subplot(111)
        self.ax.plot(x, y, 'limegreen', label='CD = CD_0 + K*CL^2')
        self.ax.set_title('Drag Polar < Mach 0.9:  CD vs. CL\n')
        self.ax.set_xlabel('CD - [ - ]')
        self.ax.set_ylabel('CL - [ - ]')
        self.ax.set_xlim(x_min, x_max)
        self.ax.xaxis.set_ticks(np.arange(x_min, x_max, x_step))
        self.ax.set_ylim(y_min, y_max)
        self.ax.yaxis.set_ticks(np.arange(y_min, y_max, y_step))
        self.ax.grid(True)
        self.ax.legend(bbox_to_anchor=(0.6, 0.75), shadow=True, ncol=1)
        self.draw()


class AircraftPerformance(QDialog):
    def __init__(self):
        super(AircraftPerformance, self).__init__()
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)
        self.assign_widgets()
        plot_window_1 = self.ui.graphicsView_1
        plot_window_2 = self.ui.graphicsView_2
        self.m_1 = PlotCanvas(plot_window_1)
        self.m_2 = PlotCanvas(plot_window_2)
        self.show()
        self.i = BoomLibrary()

        self.rho = None
        self.Temp = None
        self.Pressure = None
        self.a = None
        self.a_18 = None
        self.length = None
        self.b_wing = None
        self.S_wing = None
        self.AR_wing = None
        self.W_max = None
        self.V_stream = np.linspace(0.1, 1312, 10**4)
        self.V_sub = np.linspace(0.1, 850, 10**4)
        self.V_ts = np.linspace(850, 1175, 10**4)
        self.V_ss = np.linspace(1175, 1312, 10**4)
        self.V_stall = None
        self.V_max_ss = None
        self.V_max_ts = None
        self.V_TO = None
        self.Mach_ts = None
        self.Mach_ss = None
        self.q = None
        self.q_sub = None
        self.q_ts = None
        self.q_ss = None
        self.q_a = None
        self.CL_max = None
        self.CL = None
        self.CL_sub = None
        self.CL_ts = None
        self.CL_ss = None
        self.CD = None
        self.CD_sub = None
        self.CD_ts = None
        self.CD_ss = None
        self.CD_0 = None
        self.CD_0_sub = None
        self.CD_0_ts = None
        self.CD_0_ss = None
        self.CL_a = None
        self.CD_a = None
        self.T_req_a = None
        self.e = None
        self.e_1 = None
        self.e_2 = None
        self.K = None
        self.K_r = 2.0
        self.K_s = 0.087
        self.K_p = 1.035
        self.K_t = 0.975
        self.R_air = 1716                                   # Gas Constant R for Ambient Air - [ (ft*lbm)/(lbm*R) ]
        self.g_c = 32.174                                   # Gravitational Prop. Constant - [ (ft*lbm)/(lbf*s^2) ]
        self.eta_mechanical = None
        self.dist_TO = None
        self.T_req = None
        self.T_req_sub = None
        self.T_req_ts = None
        self.T_req_ss = None
        self.T_req_a = None
        self.T_avl = None
        self.T_avl_input = None
        self.T_target = None
        self.T_target_SL = None
        self.T_target_18 = None
        self.T_TO_min = None
        self.F_TO_min = None
        self.T_ratio = None
        self.LE_sweep = None
        self.altitude = None
        self.mod_altitude = None                           # Altitude for Initial Start-Up Calculations - [ ft ]
        self.delta_P = None
        self.delta_t = None
        self.T_2_W = None
        self.W_S = None
        self.L_2_D_sub = None
        self.L_2_D_ts = None
        self.L_2_D_ss = None
        self.L_2_D_sub_max = None
        self.L_2_D_ts_max = None
        self.L_2_D_ss_max = None
        self.V_T_req_min = None
        self.T_req_2_W_min = None
        self.T_req_min = None
        self.V_range = None
        self.theta_ROC = None
        self.V_ROC_max = None
        self.ROC_max = None
        self.Z = None
        self.P_SL = 101324.204
        self.P_18k = 50633.374
        self.h_18k = 5486.4

    def assign_widgets(self):
        self.ui.pushButton_2.clicked.connect(exit_app)
        self.ui.pushButton.clicked.connect(self.data_retrieval)
        self.ui.radioButton.clicked.connect(self.Standard_Day)
        self.ui.radioButton_2.clicked.connect(self.Hot_Day)
        self.ui.radioButton_3.clicked.connect(self.Cold_Day)
        self.ui.radioButton_4.clicked.connect(self.Tropic_Day)
        self.ui.horizontalSlider.valueChanged.connect(self.gl_zoom_slider)

    def gl_zoom_slider(self):                           # Slider to Control GL Zooming
        self.mod_altitude = float(self.ui.horizontalSlider.value())*1000
        self.data_retrieval()

    def Standard_Day(self):
        if self.altitude is not None:
            self.Temp = self.i.T_SL-0.00356*self.altitude
            self.Pressure = 2216*(self.Temp/self.i.T_SL)**5.256
            self.rho = self.Pressure/(self.i.R_air*self.Temp)
        self.data_retrieval()

    def Hot_Day(self):
        if self.altitude is not None:
            Temp = self.i.T_SL-0.00356*self.altitude
            self.Temp = Temp*self.i.theta_break
            self.Pressure = 2216*(Temp/self.i.T_SL)**5.256
            self.rho = self.Pressure/(self.i.R_air*self.Temp)
        self.data_retrieval()

    def Cold_Day(self):
        if self.altitude is not None:
            Temp = self.i.T_SL-0.00356*self.altitude
            self.Temp = Temp/self.i.theta_break
            self.Pressure = 2216*(Temp/self.i.T_SL)**5.256
            self.rho = self.Pressure/(self.i.R_air*self.Temp)
        self.data_retrieval()

    def Tropic_Day(self):
        self.data_retrieval()

    def data_retrieval(self):
        # Retrieve Input Values from GUI
        self.b_wing = (float(self.ui.lineEdit_1.text())*3.28084)
        self.S_wing = (float(self.ui.lineEdit_2.text())*10.76391)
        self.W_max = (float(self.ui.lineEdit_3.text())*0.224809)
        self.dist_TO = (float(self.ui.lineEdit_26.text())*3.28084)
        self.length = (float(self.ui.lineEdit_5.text())*3.28084)
        self.T_avl_input = (float(self.ui.lineEdit_21.text())*0.224809)
        self.eta_mechanical = (float(self.ui.lineEdit_22.text())/100)
        self.CD_0 = float(self.ui.lineEdit_23.text())
        self.CL_max = float(self.ui.lineEdit_32.text())
        self.T_ratio = (1+float(self.ui.lineEdit_33.text())/100)
        self.LE_sweep = (float(self.ui.lineEdit_31.text())*(np.pi / 180))
        self.altitude = float(self.ui.lineEdit.text())*1000
        if self.mod_altitude is not None:
            self.altitude = self.mod_altitude
        if self.Temp is None:
            self.Cold_Day()
        self.library()

    def library(self):
        self.q = (1/2)*self.V_stream**2*self.rho
        self.q_sub = (1/2)*self.V_sub**2*self.rho
        self.q_ts = (1/2)*self.V_ts**2*self.rho
        self.q_ss = (1/2)*self.V_ss**2*self.rho
        self.T_target_SL = (self.T_avl_input*self.T_ratio)*self.eta_mechanical
        self.T_target = self.T_target_SL*(self.rho/self.i.rho_SL)
        self.T_target_18 = self.T_target_SL*(self.i.rho_18_HD/self.i.rho_SL)
        self.T_2_W = self.T_target/self.W_max
        self.W_S = self.W_max/self.S_wing
        self.a = np.sqrt(self.i.cp*self.i.R_air*self.Temp)
        self.a_18 = np.sqrt(self.i.cp*self.i.R_air*self.i.T_18_CD)
        self.q_a = (1/2)*self.a**2*self.rho
        self.aircraft_parameters()

    def aircraft_parameters(self):
        self.CD_0_sub = self.CD_0
        self.CD_0_ts = self.CD_0*2.2
        self.CD_0_ss = self.CD_0*1.2
        self.AR_wing = self.b_wing**2/self.S_wing
        self.e_1 = 4.61*(1-0.045*self.AR_wing**0.68)*(np.cos(self.LE_sweep)**0.15)-3.1
        self.e_2 = (2/(2-self.AR_wing+np.sqrt(4+(self.AR_wing**2)*(1+np.tan(self.LE_sweep)**2))))
        self.e = (self.e_1+self.e_2)/2
        self.K = 1/(np.pi*self.AR_wing*self.e)
        self.V_stall = np.sqrt((2/self.i.rho_SL_HD)*(self.W_max/self.S_wing)*(1/self.CL_max))
        self.V_TO = 1.25*self.V_stall
        self.V_max_ss = np.sqrt(((self.T_target / self.W_max) * (self.W_max / self.S_wing) + (self.W_max / self.S_wing) * np.sqrt((
                self.T_target/self.W_max) ** 2 - 4 * self.CD_0_ss * self.K)) / (self.rho*self.CD_0_ss))
        self.V_max_ts = np.sqrt(((self.T_target / self.W_max) * (self.W_max / self.S_wing) + (self.W_max / self.S_wing) * np.sqrt((
                self.T_target/self.W_max) ** 2 - 4 * self.CD_0_ts * self.K)) / (self.rho*self.CD_0_ts))
        self.Mach_ts = self.V_max_ts / self.a_18
        self.Mach_ss = self.V_max_ss / self.a_18
        self.thrust_required()

    def thrust_required(self):
        self.CL = self.W_max / (self.q*self.S_wing)
        self.CL_sub = self.W_max / (self.q_sub*self.S_wing)
        self.CL_ts = self.W_max / (self.q_ts*self.S_wing)
        self.CL_ss = self.W_max / (self.q_ss*self.S_wing)
        self.CD = self.CD_0+self.K*self.CL**2
        self.CD_sub = self.CD_0_sub+self.K*self.CL_sub**2
        self.CD_ts = self.CD_0_ts+self.K*self.CL_ts**2
        self.CD_ss = self.CD_0_ss+self.K*self.CL_ss**2
        self.L_2_D_sub = self.CL_sub/self.CD_sub
        self.L_2_D_ts = self.CL_ts/self.CD_ts
        self.L_2_D_ss = self.CL_ss/self.CD_ss
        self.L_2_D_sub_max = 1/(np.sqrt(4*self.CD_0_sub*self.K))
        self.L_2_D_ts_max = 1/(np.sqrt(4*self.CD_0_ts*self.K))
        self.L_2_D_ss_max = 1/(np.sqrt(4*self.CD_0_ss*self.K))
        self.T_req_2_W_min = np.sqrt(4*self.CD_0_sub*self.K)
        self.V_T_req_min = np.sqrt((2/self.rho)*np.sqrt(self.K/self.CD_0_sub)*self.W_S)*0.3048
        self.V_range = 1.32*self.V_T_req_min
        self.Z = (1+(np.sqrt(1+(3/((self.L_2_D_sub_max**2)*(self.T_2_W**2))))))
        self.V_ROC_max = np.sqrt(((self.T_2_W*self.W_S)/(3*self.rho*self.CD_0_sub))*self.Z)*0.3048
        # self.ROC_max = np.sqrt((self.W_S*self.Z)/(3*self.rho*self.CD_0_sub))*(self.T_2_W**(3/2))*(1-(self.Z/6)-(3/(2*(self.T_2_W**2)*(
        #     self.L_2_D_sub_max**2)*self.Z)))
        # self.theta_ROC = np.arcsin((self.T_2_W-np.sqrt(4*self.CD_0_sub*self.K))*(np.pi / 180))
        # print(self.V_ROC_max)
        self.T_req_min = (self.T_req_2_W_min*self.W_max)/0.225
        self.T_req = self.q*self.S_wing*self.CD
        self.T_req_sub = self.q_sub*self.S_wing*self.CD_sub
        self.T_req_ts = self.q_ts*self.S_wing*self.CD_ts
        self.T_req_ss = self.q_ss*self.S_wing*self.CD_ss
        self.T_TO_min = ((1.69*self.W_max**2)/(self.i.g_c*self.i.rho_SL_HD*self.dist_TO*self.S_wing*self.CL_max))/0.225
        self.F_TO_min = self.T_TO_min/self.eta_mechanical
        self.CL_a = self.W_max / (self.q_a*self.S_wing)
        self.CD_a = self.CD_0_ts + self.K*self.CL_a**2
        self.T_req_a = self.q_a * self.S_wing * self.CD_a

        if self.Mach_ss < 1.0:
            self.V_max_ss = 0
            self.Mach_ss = 0
            self.delta_P = 0
            self.delta_t = 0
            self.upload_values()
        else:
            self.sonic_boom()

    def sonic_boom(self):
        self.delta_P = self.K_p * self.K_r * np.sqrt(self.P_18k*self.P_SL) * (self.Mach_ss ** 2 - 1) ** (1 / 8) * self.h_18k ** (-3 / 4) * self.length ** (
                3/4) * self.K_s
        self.delta_t = self.K_t * (3.42 / (self.V_max_ss * 0.3048)) * (self.Mach_ss / ((self.Mach_ss ** 2 - 1) ** (3 / 8))) * self.h_18k ** (
                1/4) * self.length ** (3/4) * self.K_s

        self.upload_values()

    def upload_values(self):
        try:
            self.ui.lineEdit_6.setText('{:.1f}'.format(self.AR_wing))
            self.ui.lineEdit_7.setText('{:.2f}'.format(self.e))
            self.ui.lineEdit_8.setText('{:.3f}'.format(self.K))
            self.ui.lineEdit_9.setText('{:.2f}'.format(self.V_stall*0.3048))
            self.ui.lineEdit_10.setText('{:.2f}'.format(self.V_TO*0.3048))
            self.ui.lineEdit_24.setText('{:.2f}'.format(self.T_TO_min))
            self.ui.lineEdit_25.setText('{:.2f}'.format(self.F_TO_min))
            self.ui.lineEdit_27.setText('{:.2f}'.format(self.delta_P))
            self.ui.lineEdit_28.setText('{:.2f}'.format(self.Mach_ts))
            self.ui.lineEdit_30.setText('{:.4f}'.format(self.delta_t))
            self.ui.lineEdit_29.setText('{:.2f}'.format(self.V_max_ts * 0.3048))
            self.ui.lineEdit_35.setText('{:.2f}'.format(self.V_max_ss * 0.3048))
            self.ui.lineEdit_34.setText('{:.2f}'.format(self.Mach_ss))
            self.ui.lineEdit.setText('{:.1f}'.format(self.altitude/1000))
        except Exception as err:
            print(err)
            QApplication.restoreOverrideCursor()
            bad_file()
        self.thrust_plot_initiate()

    def thrust_plot_initiate(self):
        self.m_1.figure.clf()
        self.T_avl = self.T_avl_input
        self.m_1.plot_thrust(self.V_stream, self.V_sub, self.V_ts, self.V_ss, self.T_req_sub, self.T_req_ts, self.T_req_ss, self.T_target*self.V_stream**0, self.a, self.T_req_a)
        self.polar_plot_initiate()

    def polar_plot_initiate(self):
        self.m_2.figure.clf()
        CL = np.linspace(0.001, 2.5, 10**4)
        CD_plot = self.CD_0+self.K*CL**2

        self.m_2.plot_polar(CD_plot, CL)


def exit_app():
    app.exit()


def no_file():
    msg = QMessageBox()
    msg.setText('There was No File Selected.')
    msg.setWindowTitle("No File")
    msg.exec_()
    return None


def bad_file():
    msg = QMessageBox()
    msg.setText('Unable to Process the Selected File.')
    msg.setWindowTitle("Bad File")
    msg.exec_()
    return None


if __name__ == "__main__":
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)
    app.aboutToQuit.connect(app.deleteLater)
    main_win = AircraftPerformance()
    sys.exit(app.exec_())
