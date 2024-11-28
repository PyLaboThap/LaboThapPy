# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 14:09:18 2024

@author: Elise Neven
@email: elise.neven@uliege.be

"""


class WorkConnector:
    """
    A class to handle transfer of work power.

    **Attributes**:

        W_dot : float, optional
            Work power in W.
        N : float, optional
            Speed in rpm.
        C : float, optional
            Torque in Nm.
        variables_input : list of lists
            A list of the variables used to define the work connector. Each entry is a list of [variable_name, value].

    **Methods**:

        __init__(self):
            Initializes the WorkConnector object with.

        set_W_dot(self, value):
            Sets the work power and updates the list of known variables.

        set_N(self, value):
            Sets the speed and updates the list of known variables.
        
        set_C(self, value):
            Sets the torque and updates the list of known variables.

        print_resume(self):
            Print a summary of the work connector properties.
    """

    def __init__(self):

        self.variables_input = []
        
        self.W_dot = None          # Work power [W]
        self.N = None              # Speed [rpm]
        self.C = None            # Torque [Nm]


    def calculate_properties(self):
        pass

    def set_W_dot(self, value):
        self.W_dot = value
        self.variables_input = self.variables_input + [['W_dot', value]]
        self.calculate_properties()

    def set_N(self, value):
        self.N = value
        self.variables_input = self.variables_input + [['N', value]]
        self.calculate_properties()

    def set_C(self, value):
        self.C = value
        self.variables_input = self.variables_input + [['C', value]]
        self.calculate_properties()


    def print_resume(self):
        """
        Print a summary of the work connector properties
        """
        print("Work power: " + str(self.W_dot) + " [W]")
        print("Speed: " + str(self.N) + " [rpm]")
        print("Torque: " + str(self.C) + " [Nm]")
