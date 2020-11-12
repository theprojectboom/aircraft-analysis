# Mann Mansy, Jason Conway, and Malorie Travis
# Future work: email mtravis024@gmail.com or abdalrahm1998@gmail.com with questions about current program
# 4-2-2020
# MAE 3403 and MAE 4374

# This GUI is designed to take user inputs about an aircraft and its mission and output useful information.
# I can write more about it, but probably the most helpful way to catch up is by running the gui!

# Inputs: vehicle, engine, battery details, weights, and leg details
# Leg details are output to the text box under "Current Mission Profile"
    # Presently this cannot be removed! Will need to be added in the future. Will need to restart the program for a new mission
# Outputs: after clicking run, program will output a profile graphic and a report featuring (for each leg):
    # - Range
    # - TSFC
    # - Best Velocity

# Most of these equations come from the Anderson textbook, Aircraft Performance and Design (applied aero textbook)

# This file is the main program which runs the gui and calls necessary programs: OpenGL_2D_class and missionProfileRevC
# also sorry it's very long! I (malorie) just learned python this semester.

# Import general capabilities
import numpy as np
from scipy.interpolate import griddata
from matplotlib.figure import Figure

import sys
from PyQt5.QtWidgets import QDialog, QApplication
from PyQt5.QtWidgets import QFileDialog, QMessageBox
from PyQt5.QtGui import QCursor
from PyQt5.QtCore import Qt
from copy import deepcopy

# Figure in GUI
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# Import the UI
from RevE import Ui_Dialog
# To create the error messages. Might not work?

# Import drawing capabilities for the Mission Profile graphic from Delahoussaye file OpenGL_2D_class
# from OpenGL_2D_class import gl2D, gl2DText, gl2DCircle, gl2DArc, gl2DArrow


class PlotCanvas(FigureCanvas):
    def __init__(self, parent, width=None, height=None, dpi=100):
        if width == None: width = parent.width()/100
        if height == None: height = parent.height()/100
        fig = Figure(figsize=(width, height), dpi=dpi)
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

    def BarGraphs(self,NumberOfLegs, x, y, titles):

        self.figure.clf()

        title = titles[0]
        xlabel = titles[1]
        ylabel = titles[2]

        Legs = [None]*(NumberOfLegs)
        ax = self.figure.add_subplot(111)
        ax.set_title(title)

        for i in range(NumberOfLegs):
            Legs[i] = x[i]
        pos = np.arange(len(Legs))
        ax.bar(pos, y)
        # ax.xticks(pos, Legs)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        self.show()

class main_window(QDialog):
    def __init__(self):
        # Boilerplate code
        super(main_window, self).__init__()
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)
        self.checkClicked()

        plotwin1 = self.ui.graphicsViewFuel     # Initialize Fuel plot
        self.Fplot = PlotCanvas(plotwin1)

        plotwin2 = self.ui.graphicsViewRange    # Initialize Range plot
        self.Rplot = PlotCanvas(plotwin2)

        plotwin3 = self.ui.graphicsViewTime     # Initialize Time plot
        self.Tplot = PlotCanvas(plotwin3)
        # self.setupGLWindows()

        # Plotting initialize
        graphWindow = self.ui.openGLWidget
        self.show()

    # initialize any necessary variables
        # User Mission Profile input values
        self.MissionLegs = ['','','','','','','','']    # Mission legs types: string array
        self.begalt = np.zeros(8)        # Beginning altitude: float array
        self.endalt = np.zeros_like(self.begalt)        # Ending altitude: float array
        self.MissionReport = ' '  # The stuff that goes into the "current mission profile box"
        self.dist = np.zeros_like(self.begalt)

        # Cruise eg functionality
        self.cruiseLegs = 0
        
        # Mugin aircraft specifications
        self.b_wing = 14.76     # ft
        self.S_wing = 19.39     # ft^2
        self.V_stall = 41       # ft/s
        self.dist_TO = 50       # ft
        self.CD_0 = 0.045
        self.CL_max = 1.7
        self.TSFCMax = 0
        
        # System specs
        self.eta_elec = 0.70
        self.eta_prop = 0.80
        self.f = 0.025
        self.loss_install = 0.25

        # Derived parameters
        self.AR_wing = (self.b_wing ** 2) / self.S_wing   # Aspect Ratio

        self.sweep = np.pi / 180
        self.e = 0.9
        self.K = 1 / (10.08 * np.pi)
        self.CL_calc_HD = 0      # Hot day CL

        self.CL_CDmax = np.sqrt(self.CD_0 / self.K)      # CL for max lift to drag ratio
        self.CD_CLmax = 2 * self.CD_0                   # CD for max lift to drag ratio
        self.LoD_max = self.CL_CDmax  / self.CD_CLmax   # Max lift to drag ratio

        # Power
        self.P_a = 0     # Actual Power
        self.P_req = 0   # Required Power
        
        #Thrust
        self.T_a = 0     # Actual thrust
        self.T_req = 0   # Required thrust
        
        # weights
        self.EW = 0  # Empty weight
        self.PW = 0  # Payload weight
        self.FW = 0  # Fuel weight
        self.TW = 0  # Total weight

        # configuration
        self.Config = ''
        self.Engine = ''
        self.NEngine = ''
        self.NumberOfLegs = 0

        # Atmospheric values
        self.altitude = 0
        self.temperature = 0
        # self.pressure = np.zeros([6])
        self.density = 0
        self.density_SL = 2.1546e-3 #slug/ft^3
        self.Rconst = 1716  # (ft lbf) / (slug R)
        self.g_c = 32.2     # (lbm ft)/ (lbf s^2)

        # Outputs

        self.RangeCruise = 0
        self.EnduranceTotal = 0
        self.RangeTotal = 0


    # Check if all the buttons are clicked
    def checkClicked(self):
        # check syntax
        self.ui.pushButton_Run.clicked.connect(self.main)               # Run button
        self.ui.pushButton_Exit.clicked.connect(app.exit)               # Exit button
        self.ui.pushButton_AddLeg.clicked.connect(self.AddLeg)          # Add leg button

    # def setupGLWindows(self):
    #     # Make drawing widget
    #     self.glwindow1 = gl2D(self.ui.openGLWidget, self.DrawProfile)
    #     self.glwindow1.setViewSize(-10, 500, -10, 200, allowDistortion=False)
        

    # Add a leg to the mission
    def AddLeg(self):
        # Everytime this runs, keep track of it by NumberOfLegs variable
        self.NumberOfLegs += 1
        # Pull the altitude values
        self.begalt[self.NumberOfLegs - 1] = float(self.ui.lineEdit_AS.text())
        self.endalt[self.NumberOfLegs - 1] = float(self.ui.lineEdit_AF.text())

        # Pull the leg type: takeoff, cruise, etc
        self.MissionLegs[self.NumberOfLegs - 1] = (self.ui.comboBox_Leg.currentText())
        
        # Check to see if the first leg is takeoff
        if self.MissionLegs[0] != 'Takeoff roll':
            self.ui.textEdit_Report.setText("Dude! First leg must be takeoff roll.")

        if self.MissionLegs[self.NumberOfLegs - 1] == 'Loiter':
            self.loitertime = float(self.ui.lineEdit_LT.text())

        # Add this leg to the info box
        self.GenerateMissionProfile()

        # Clear boxes so they're ready for a new input
        self.ui.lineEdit_AS.setText('0')
        self.ui.lineEdit_AF.setText('0')
        self.ui.lineEdit_LT.setText('0')

    # As the user inputs a leg and adds it, add the info to the box next door
    def GenerateMissionProfile(self):
        # Add the leg info to the GUI
        self.MissionReport += '\n Leg ' + str(self.NumberOfLegs) + ':'
        self.MissionReport += '\t Leg type: ' + str(self.MissionLegs[self.NumberOfLegs - 1])
        self.MissionReport += '\t Beginning Altitude: ' + str(self.begalt[self.NumberOfLegs - 1]) + 'ft'
        self.MissionReport += '\t End Altitude: ' + str(self.endalt[self.NumberOfLegs - 1]) + 'ft' + '\n'

        self.ui.textEdit_Current.setText(self.MissionReport)

    # Draws the mission profile - requires OpenGL
    # def DrawProfile(self):
    #
    #     NumberOfLegs = self.NumberOfLegs
    #     altitudes = self.altitude
    #
    #     # glColor3f(0, 0, 0)
    #     # glLineWidth(1.5)
    #
    #     # Cycle through the altitudes and draw accordingly
    #     for i in range(NumberOfLegs):
    #         glBegin(GL_LINE_STRIP)  # begin drawing connected lines
    #         glVertex2f(100 * i, self.begalt[i])
    #         glVertex2f(100 * (i + 1), self.endalt[i])
    #         glEnd()
    #
    #     # Find the highest altitude and set it as the y maximum
    #     for i in range(1, len(altitudes)):
    #
    #         if altitudes[i - 1] < altitudes[i]:
    #             maxaltitude = altitudes[i]
    #
    #         else:
    #             maxaltitude = altitudes[i]
    #
    #     # set bounds
    #     xmax = NumberOfLegs * 100 + 20
    #     xmin = -20
    #     ymin = -10
    #     ymax = maxaltitude + 10
    #
    #     self.glwindow1.setViewSize(xmin, xmax, ymin, ymax, allowDistortion=False)
    #
    #     self.glwindow1.glUpdate()
    
    # Pull atmospheric values from table
    def getAtmVals(self):
    
        altitudecol, tempcol, presscol = np.loadtxt("atmvalues.txt", skiprows=1, unpack=True)
    
        self.temperature = np.zeros_like(self.altitude)
        self.pressure = np.zeros_like(self.altitude)
        self.density = np.zeros_like(self.altitude)
            
        for i in range(len(self.altitude)):
            self.temperature[i] = float(griddata(altitudecol, tempcol, self.altitude[i])) + 460 # Convert to Rankine from Fahrenheit
            self.pressure[i] = float(griddata(altitudecol, presscol, self.altitude[i]))
            self.density[i] = self.pressure[i] / (self.temperature[i] * self.Rconst) * 144 #Slug/ft^3
    
    
    # Step through each of the user-defined legs: direct the program to the appropriate function
    def stepThruLegs(self):
            
        NumberOfLegs = self.NumberOfLegs

        # Initialize outputs
        self.V = np.zeros(NumberOfLegs)
        self.Range = np.zeros(NumberOfLegs)
        self.Endurance = np.zeros(NumberOfLegs)
        self.RunningTotalTime = np.zeros(NumberOfLegs)
        self.distS = np.zeros(NumberOfLegs)
        self.distF = np.zeros(NumberOfLegs)
        self.FB = np.zeros(NumberOfLegs)
        self.TotW = np.zeros(NumberOfLegs)

        # Error Messages
        if self.MissionLegs[0] != 'Takeoff roll':
            self.ui.textEdit_Report.setText("Error: First leg must be takeoff roll")

        if self.MissionLegs[NumberOfLegs-1] != 'Landing':
            self.ui.textEdit_Report.setText("Error: Last leg must be landing")
        
        # Need this for report generating
        self.allalt = np.zeros([NumberOfLegs, 100])


        
        # loop through each leg
        for i in range(NumberOfLegs):
            # Compute total weight
            self.TW = self.PW + self.EW + self.FW
        
            # call getAtmVals and create 100 pts of altitude IF the leg is changing altitude
            # if self.begalt[i] != self.endalt[i]:

            self.altitude = np.linspace(self.begalt[i], self.endalt[i], 100)
            self.getAtmVals()
                        
            for j in range(0, 99):
                self.allalt[i][j] = self.altitude[j]
        
            # Takeoff leg
            if self.MissionLegs[i] == 'Takeoff roll':
                Leg = i
                self.TakeoffCalc(Leg)
        
            # Climb Leg
            if self.MissionLegs[i] == 'Constant speed climb':
                Leg = i
                self.ClimbCalc(Leg)
        
            # Cruise
            if self.MissionLegs[i] == 'Subsonic cruise':
                Leg = i
                self.cruiseLegs += 1
        
            # Payload drop
            if self.MissionLegs[i] == 'Payload delivery':
                Leg = i
                self.PayloadCalc(Leg)

            # Loiter Leg
            if self.MissionLegs[i] == 'Loiter':
                Leg = i
                self.LoiterCalc(Leg)

            # Descend
            if self.MissionLegs[i] == 'Constant speed descent':
                Leg = i
                self.DescendCalc(Leg)
        
            # Land
            if self.MissionLegs[i] == 'Landing':
                Leg = i
                self.LandingCalc(Leg)
        
            # update fuel weight
            self.FW -= self.FB[i]
            self.TW = self.PW + self.EW + self.FW
            self.TotW[i] = self.TW
            # add to the totals for range and endurance
            self.RangeTotal += self.Range[i]
            self.EnduranceTotal += self.Endurance[i]

        
        # self.cruiseLeg = np.zeros([1, self.cruiseLegs])
             
        for i in range(NumberOfLegs):
            # Compute total weight
            self.TW = self.PW + self.EW + self.FW
        
            self.altitude = np.linspace(self.begalt[i], self.endalt[i], 100)
            self.getAtmVals()

            
            if self.MissionLegs[i] == 'Subsonic cruise':
                Leg = i
                self.CruiseCalc(Leg)

        # Draw Mission Profile
        np.reshape(self.allalt, [1, NumberOfLegs * 100])
    
          
    def TakeoffCalc(self, Leg):

        # Takeoff thrust
        self.T = (1.44 * self.TW**2) / (32.2 * self.density[1] * self.dist_TO * self.S_wing * self.CL_max)

        # Load factor
        n = .5 * self.density_SL * (1.15 * self.V_stall) ** 2 * self.S_wing * 0.9 * self.CL_max / self.TW

        # Radius of takeoff
        Rad = (1.15 * self.V_stall) ** 2 / (32.2 * (n - 1))

        # Assume height of obstacle is 30 ft
        hOB = 30    # ft
        
        # Calculate angle made by the: point of takeoff (a straight, vertical line) and 
        thet = np.arccos(1 - hOB / Rad)

        self.sA = Rad * np.sin(thet)

        self.distF[Leg] = self.dist_TO + self.sA

        self.distS[Leg] = 0

        # Range and endurance start at the beginning of the next leg
        self.Endurance[Leg] = 0
        self.Range[Leg] = self.distF[Leg]

        # Find fuel burnt. Max throttle at Takeoff velocity
        V_TO = 1.2 * self.V_stall
        self.V[Leg] = V_TO
        t = self.dist_TO / V_TO         # seconds
        self.FB[Leg] = self.mdotfuelmax / self.g_c * t
        self.Endurance[Leg] = self.distF[Leg] / V_TO
        self.RunningTotalTime[Leg] += self.Endurance[Leg] + self.RunningTotalTime[Leg - 1]

    def ClimbCalc(self, Leg):
        # Power calculation
        self.P_a = self.P_max

        P = self.P_a * 0.00134102 * 550 #Watts to W to HP to lb-ft / s ....

        t1 = (P / self.TW) #Ft/s

        RoC_max_SL = t1  - np.sqrt((2 / self.density_SL)) * np.sqrt(self.K / (3 * self.CD_0)) * np.sqrt(self.TW / self.S_wing) * (1.155 / self.LoD_max)
        
        # Rate of climb slope calculation (arbitrary points for slope)
        RoC_max1 = t1 - np.sqrt((2 / self.density[33])) * np.sqrt(self.K / (3 * self.CD_0)) * np.sqrt(self.TW/self.S_wing) * (1.155 / self.LoD_max)

        RoC_max2 = t1 - np.sqrt((2 / self.density[66])) * np.sqrt(self.K / (3 * self.CD_0)) * np.sqrt(self.TW/self.S_wing) * (1.155 / self.LoD_max)

        slope_RoC = (RoC_max2 - RoC_max1) / (self.altitude[33] - self.altitude[66])
        
        # Time to climb equation. Integral of Rate of Climb equation between the start and end altitude
        tclimb_min = 1 / slope_RoC * (np.log(RoC_max_SL + slope_RoC * self.endalt[Leg]) - np.log(RoC_max_SL))
        
        # # Average max Climb Velocity and Thrust calculation
        V_RoC_max1 = np.sqrt(2 / self.density_SL * np.sqrt(self.K / (3 * self.CD_0))*self.TW/self.S_wing)
        V_RoC_max2 = np.sqrt(2 / self.density_SL * np.sqrt(self.K / (3 * self.CD_0))*self.TW/self.S_wing)

        # V_RoC_max_avg = (4 * self.TW / self.S_wing * self.K) / (self.density_SL * P / self.TW)

        V_RoC_max_avg = (V_RoC_max1 + V_RoC_max2) / 2
        self.V[Leg] = V_RoC_max_avg
        T_RoC_max_avg = P / V_RoC_max_avg
        F_RoC_max_avg = T_RoC_max_avg / (1 - self.loss_install)
        
        # Average TSFC calculation
        mdot_0 = self.mdotfuelmax / self.f
        self.SFC = self.f / (T_RoC_max_avg / mdot_0)
        self.TSFC = self.SFC / (1 - self.loss_install)

        # Check for max TSFC
        if self.TSFC > self.TSFCMax:
            self.TSFCMax = deepcopy(self.TSFC)
        
        # Fuel Burn
        self.FB[Leg] = self.mdotfuelmax / self.g_c * tclimb_min
        
        # V_RoC_max1 = V_RoC_max_avg

        # distance traveled
        A = P / (V_RoC_max1 * self.TW)
        B = 0.5 * self.density_SL * V_RoC_max1**2 *(self.TW / self.S_wing)**-1 * self.CD_0
        C = (self.TW / self.S_wing) * (2 * self.K) / (self.density_SL * V_RoC_max1**2)
        
        theta = np.arcsin(A - B - C)
        
        alt_change = self.endalt[Leg] - self.begalt[Leg]
        dist_travel = alt_change / np.arctan(theta)
        
        self.distS[Leg] = self.distF[Leg-1]
        self.Range[Leg] = self.distS[Leg] + dist_travel

        self.Endurance[Leg] = tclimb_min
        self.RunningTotalTime[Leg] += self.Endurance[Leg] + self.RunningTotalTime[Leg - 1]

    def CruiseCalc(self, leg):
        # Power calculation
        self.P_a = self.eta_prop * self.eta_elec * self.P_max

        # set L = w, solve 
        self.CL_HD = np.sqrt(self.CD_0 / self.K)

        self.V[leg] = np.sqrt(2 * self.TW / (self.density[99] * self.S_wing * self.CL_HD))
        
        # dynamic pressure
        self.q_SL_HD = (1/2) * self.V[leg]**2 * self.density[99] /32.2   # dynamic pressure, units of psf

        # Calculate aero coefficients, required thrust, and required power
            # Assume hot day, standard day is commented out
        
        self.CD_HD = self.CD_0 + self.K*self.CL_HD**2
        self.T_req_HD = self.q_SL_HD * self.S_wing * self.CD_HD
        self.P_req_HD = self.T_req_HD * self.V[leg] * 1000
        
        # Average TSFC calculation
        mdotfuel = self.mdotfuelmax * 2/3        # Assume 2/3 throttle
        mdot_0 = mdotfuel / self.f 
        self.SFC = self.f / (self.T_req_HD / mdot_0)
        self.TSFC = self.SFC / (1 + self.loss_install)

        # Check for max TSFC
        if self.TSFC > self.TSFCMax:
            self.TSFCMax = deepcopy(self.TSFC)

        # Use battery- this function only calculates the MAX power the battery can put out (if demanded)
        # self.UseBattery()
        # if self.P_req_HD > self.BattPower:
        #     self.ui.textEdit_Report.setText("Error: Battery power not sufficient")

        # Endurance of fuel divided equally among the cruise legs

        # print(self.FW / (mdotfuel / self.g_c))
        # units: minutes
        # burn remaining fuel (split between cruise legs) & find time that it can run for
        self.Endurance[leg] = (self.FW*0.8 / (mdotfuel / self.g_c)) / (self.cruiseLegs) / 60
        self.RunningTotalTime[leg] += self.Endurance[leg] + self.RunningTotalTime[leg - 1]
        self.Range[leg] = self.Endurance[leg] * self.V[leg]
        self.FB[leg] = self.FW/self.cruiseLegs

        self.Vcruise = deepcopy(self.V[leg])
    
    def PayloadCalc(self, leg):
        self.TW = self.TW - self.PW
        self.FB[leg] = 0

        self.Range[leg] = 0
        self.Endurance[leg] = 0
        self.RunningTotalTime[leg] += self.Endurance[leg] + self.RunningTotalTime[leg - 1]

    def LoiterCalc(self, leg):
        mdotfuel = self.mdotfuelmax * 2/3
        self.FB[leg] = mdotfuel/self.g_c * self.loitertime
        self.Range[leg] = 0
        self.Endurance[leg] = self.loitertime
        self.RunningTotalTime[leg] += self.Endurance[leg] + self.RunningTotalTime[leg - 1]

    def DescendCalc(self, leg):
        # Assume glide!
        
        # solve for glide angle
        self.GlideAngle = np.arctan(1 / (self.LoD_max)) # units of radians
        
        # find the glide range **pass altitude of previous leg**
        self.Range[leg] = self.begalt[leg] / np.tan(self.GlideAngle)

        # equilibrium glide velocity
        self.CL_HD = np.sqrt(self.CD_0 / self.K)
        self.V[leg] = np.sqrt((2*np.cos(self.GlideAngle)*self.TW)/(self.density[99]*self.CL_HD*self.S_wing))

        # horizontal velocity
        self.V_horiz = self.V[leg] * np.cos(self.GlideAngle)

        # no fuel burnt: assumed glide
        self.FB[leg] = 0

        # Time to descent
        self.Endurance[leg] = self.Range[leg] / self.V_horiz

        self.RunningTotalTime[leg] += self.Endurance[leg] + self.RunningTotalTime[leg-1]

    def LandingCalc(self, leg):
        
        # Total landing distance is approach distance, flare distance, free roll, and ground roll.
        # Using simple brake distance equation
        V_land = 1.3 * self.V_stall
        self.V[leg] = V_land
        self.dist_Land = V_land**2 / (32.2)
        self.fuelremaining = self.FW
        self.FB[leg] = 0
        self.RunningTotalTime[leg] += self.Endurance[leg] + self.RunningTotalTime[leg - 1]
    
    def main(self):
    
        # Get weights
        self.EW = float(self.ui.lineEdit_EW.text())
        self.PW = float(self.ui.lineEdit_PW.text())
        self.FW = float(self.ui.lineEdit_FW.text())

    
        # Get configuration
        self.Config = (self.ui.comboBox_Config.currentText())
        self.Engine = (self.ui.comboBox_Engine.currentText())
    
    # Set settings according to configuration
        if self.Config == 'Subsonic':

            self.ui.textEdit_Report.setText("Good")

        elif self.Config == 'Supersonic':
            # output error message
            self.ui.textEdit_Report.setText("Error: Configuration not established. Please choose a different setting.")
    
    # Set settings according to engine
        f = 0.025                       # Fuel-air ratio assumed to be a constant 2.5%
        if self.Engine == 'FT - 250':
            # Need to change these values
            self.mdotfuelmax = 0.265    # lbm/s
            self.P_max = 5000                # W
            # output error message
            self.ui.textEdit_Report.setText("Error: Configuration not established. Please choose a different setting.")
    
        elif self.Engine == 'P - 400':
            # Need to change these values

            self.mdotfuelmax = 0.353   # lbm/s
            self.P_max = 7000                # W
    

    
        # Power available equation
        self.P_a =  self.P_max
        
        # Step through legs
        self.stepThruLegs()

        self.Fplot.BarGraphs(self.NumberOfLegs, self.MissionLegs, self.FB,titles = ['Fuel Burned in each Leg', 'Legs', 'Fuel Burned (lbf)'])
        self.Rplot.BarGraphs(self.NumberOfLegs, self.MissionLegs, self.Range/5280,titles = ['Distance for each Leg', 'Legs', 'miles'])
        self.Tplot.BarGraphs(self.NumberOfLegs, self.MissionLegs, self.Endurance,titles = ['Time for each Leg', 'Legs', 'minutes'])

        # Generate report and put into proper place in the GUI
        # print(self.V, '\n',self.Range,'\n', self.Endurance,'\n',self.RunningTotalTime,'\n', self.cruiseLegs)
        self.report = self.GenerateReport()
        self.ui.textEdit_Report.setText(self.report)

        # self.setupGLWindows()
        # Python is being dumb and not drawing for some reason
        # self.DrawProfile()

    # Create report
    def GenerateReport(self):
        report = 'Output Report'
        report += '\n\nNumber of legs:\t'+ str(self.NumberOfLegs)
        report += '\nOverall range:\t'+ '{:.2f}'.format(self.RangeTotal/6076) + ' nm'+ ' ({:.2f}'.format(self.RangeTotal/5280) + ' miles)'
        report += '\nOverall Mission Time:\t'+ '{:.2f}'.format(self.EnduranceTotal/60) + ' minutes'
        # report += '\nMaximum TSFC:\t'+ '{:.2f}'.format(self.TSFCMax)
    
        # Loop through each leg and tell us about it
        for i in range(self.NumberOfLegs):     # check bounds
            report += '\n\nLeg '+ str(i+1) + ':\t' +str(self.MissionLegs[i])
            report += '\n\tAltitude from '+ str(self.begalt[i]) + ' ft to ' + str(self.endalt[i]) + ' ft'
            if str(self.MissionLegs[i]) == 'Subsonic cruise':
                report += '\n\tElapsed time: ' + '{:.2f}'.format(self.RunningTotalTime[i]/60) + ' minutes'
            report += '\n\tLeg velocity: '+ '{:.2f}'.format(self.V[i]/6076 *3600) + ' knots' + ' ({:.2f}'.format(self.V[i]) + ' ft/s)'    # should convert to mph and knots
            report += '\n\tRange: '+ '{:.2f}'.format(self.Range[i]) + ' ft' + ' ({:.2f}'.format(self.Range[i]/5280) + ' miles)'
            report += '\n\tLeg Time: ' + '{:.2f}'.format(self.Endurance[i]) + ' seconds' + ' ({:.2f}'.format(self.Endurance[i]/60) + ' minutes)'
            report += '\n\tElapsed time: ' + '{:.2f}'.format(self.RunningTotalTime[i]/60) + ' minutes'
            report += '\n\tFuel Burned: ' + '{:.2f}'.format(self.FB[i]) + ' lbf'
            report += '\n\tUpdated Aircraft Weight: ' + '{:.2f}'.format(self.TotW[i]) + ' lbf'
            if str(self.MissionLegs[i]) == 'Landing':
                report += '\n\t Remaining Fuel: {:.2f} lbf'.format(self.fuelremaining)
        report += '\n\nSeveral key assumptions were made during calculations:'
        report += '\nDuring cruise, the battery is fully discharged at its specified capability, split equally along all cruise legs'
        report += '\nTakeoff roll is 50 ft, cleared obstacle is 30 ft. These are variables, but not user inputs at this time.'
        report += '\n\nPlease note: the battery recharge calculations are not included.'
        report += '\nThis is left as an exercise for the user.'
    
        report += '\n\n'
        return report

# Boilerplate code: main statement for GUIs
if __name__ == "__main__":                                          
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)
    app.aboutToQuit.connect(app.deleteLater)
    main_window = main_window()
    sys.exit(app.exec_())
