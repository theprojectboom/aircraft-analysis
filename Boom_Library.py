
class BoomLibrary:
    def __init__(self):

        self.theta_break = 1.086

        # Atmospheric Conditions - Density - [ slugs/ft^3 ]
        self.rho_SL = 2.3769e-3
        self.rho_SL_HD = 2.19086e-3
        # self.rho_SL_HD = 3.08346e-3
        self.rho_SL_CD = 2.56705e-3

        self.rho_1 = 2.3081e-3
        self.rho_1_HD = 2.12728e-3
        self.rho_1_CD = 2.87502e-3

        self.rho_5 = 2.0482e-3
        self.rho_5_HD = 1.8873e-3
        self.rho_5_CD = 2.30612e-3

        self.rho_10 = 1.7556e-3
        self.rho_10_HD = 1.6173e-3
        self.rho_10_CD = 1.90855e-3

        self.rho_12 = 1.6480e-3
        self.rho_12_HD = 1.6173e-3
        self.rho_12_CD = 1.8873e-3

        self.rho_14 = 1.5455e-3
        self.rho_14_HD = 1.6173e-3
        self.rho_14_CD = 1.8873e-3

        self.rho_16 = 1.4480e-3
        self.rho_16_HD = 1.6173e-3
        self.rho_16_CD = 1.8873e-3

        self.rho_18 = 1.3553e-3
        self.rho_18_HD = 1.2481e-3
        self.rho_18_CD = 1.4739e-3

        # Atmospheric Conditions - Temperature - [ Rankine ]
        self.T_SL = 518.67
        self.T_SL_HD = 562.68
        # self.T_SL_HD = 399.78
        self.T_SL_CD = 450.25

        self.T_5 = 500.84
        self.T_5_HD = 543.48
        self.T_5_CD = 444.78

        self.T_10 = 483.04
        self.T_10_HD = 524.28
        self.T_10_CD = 444.26

        self.T_18 = 454.55
        self.T_18_HD = 493.55
        self.T_18_CD = 417.93

        self.R_air = 1716                                   # Gas Constant R for Ambient Air - [ (ft*lbm)/(lbm*R) ]
        # self.a_0 = None                                     # Ambient Speed of Sound - [ ft/s ]
        # self.V_0 = None                                     # Free-Stream Velocity at Engine Inlet - [ ft/s ]
        self.g_c = 32.174

        # Atmospheric Conditions - Pressure - [ lbf/ft^2 ]
        self.P_SL = 2116.2
        self.P_SL_HD = None
        self.P_SL_CD = None

        self.P_18k = 1057.5
        self.P_18k_HD = None
        self.P_18k_CD = None
