# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 14:09:18 2024

@author: Elise Neven
@email: elise.neven@uliege.be

"""

from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

class BaseComponent:
    """

    **Attributes**:

        calculable : bool
            Indicates whether the component has enough inputs to perform calculations.
        parametrized : bool
            Indicates whether the component has all required parameters set.
        solved : bool
            Indicates whether the component has been successfully solved (i.e., its state has been computed).
        inputs : dict
            A dictionary holding the input variables for the component.
        params : dict
            A dictionary holding the parameters required for the component.
        guesses : dict
            A dictionary holding initial guesses for solving the component.

    **Methods**:

        set_inputs(inputs):
            Sets the input values for the component.

        sync_inputs():
            Synchronizes the inputs dictionary with the current state of the component's connectors.
            
        set_parameters(parameters):
            Sets the parameter values for the component and checks if it is fully parametrized.
            
        set_guesses(guesses):
            Sets initial guesses for variables to be solved.
            
        check_calculable():
            Checks if the component has all the required inputs to perform calculations.
            
        check_parametrized():
            Checks if the component has all the required parameters set.
            
        get_required_inputs():
            Returns a list of required input variables for the component. Meant to be overridden in derived classes.
            
        get_required_parameters():
            Returns a list of required parameters for the component. Meant to be overridden in derived classes.
        
        get_required_guesses():
            Returns a list of required guesses for the component.
            
        solve():
            Solves the component's state, to be implemented in derived classes.
            
    **Notes**:

    - This is a base class and should be extended for specific types of components (e.g., heat exchangers, pumps, turbines).
    - The `solve` method is not implemented here and must be defined in derived classes for actual computation.
    """

    def __init__(self):
        self.calculable = False
        self.parametrized = False
        self.solved = False
        self.inputs = {}
        self.params = {}
        self.guesses = {}
        self.print_flag = 1

    def set_inputs(self, **kwargs):
        """Set inputs directly through a dictionary and update connector properties."""
        self.inputs.update(kwargs)

        # Define mappings from input keys to methods
        input_methods = {
            # su connector inputs
            'fluid':     lambda val: self.su.set_fluid(val),
            'x_su':      lambda val: self.su.set_x(val),
            'T_su':      lambda val: self.su.set_T(val),
            'h_su':      lambda val: self.su.set_h(val),
            'P_su':      lambda val: self.su.set_p(val),
            'm_dot':     lambda val: self.su.set_m_dot(val),

            # su_1 connector inputs
            'fluid_su_1': lambda val: self.su_1.set_fluid(val),
            'T_su_1':    lambda val: self.su_1.set_T(val),
            'h_su_1':    lambda val: self.su_1.set_h(val),
            'P_su_1':    lambda val: self.su_1.set_p(val),
            'm_dot_su_1': lambda val: self.su_1.set_m_dot(val),

            # su_2 connector inputs
            'fluid_su_2': lambda val: self.su_2.set_fluid(val),
            'T_su_2':    lambda val: self.su_2.set_T(val),
            'h_su_2':    lambda val: self.su_2.set_h(val),
            'P_su_2':    lambda val: self.su_2.set_p(val),
            'm_dot_su_2': lambda val: self.su_2.set_m_dot(val),

            # su_H connector inputs
            'fluid_H': lambda val: self.su_H.set_fluid(val),
            'T_su_H':    lambda val: self.su_H.set_T(val),
            'h_su_H':    lambda val: self.su_H.set_h(val),
            'P_su_H':    lambda val: self.su_H.set_p(val),
            'm_dot_H':   lambda val: self.su_H.set_m_dot(val),

            # su_C connector inputs
            'fluid_C': lambda val: self.su_C.set_fluid(val),
            'T_su_C':    lambda val: self.su_C.set_T(val),
            'h_su_C':    lambda val: self.su_C.set_h(val),
            'P_su_C':    lambda val: self.su_C.set_p(val),
            'm_dot_C':   lambda val: self.su_C.set_m_dot(val),

            # ex connector inputs
            'P_ex':      lambda val: self.ex.set_p(val),
            'T_ex':      lambda val: self.ex.set_T(val),
            'h_ex':      lambda val: self.ex.set_h(val),
            'x_ex':      lambda val: self.ex.set_x(val),

            # ex_1 connector inputs
            'P_ex_1':    lambda val: self.ex_1.set_p(val),
            'T_ex_1':    lambda val: self.ex_1.set_T(val),
            'h_ex_1':    lambda val: self.ex_1.set_h(val),

            # ex_2 connector inputs
            'P_ex_2':    lambda val: self.ex_2.set_p(val),
            'T_ex_2':    lambda val: self.ex_2.set_T(val),
            'h_ex_2':    lambda val: self.ex_2.set_h(val),         

            # ex_H connector inputs
            'P_ex_H':    lambda val: self.ex_H.set_p(val),
            'T_ex_H':    lambda val: self.ex_H.set_T(val),
            'h_ex_H':    lambda val: self.ex_H.set_h(val),

            # ex_C connector inputs
            'P_ex_C':    lambda val: self.ex_C.set_p(val),
            'T_ex_C':    lambda val: self.ex_C.set_T(val),
            'h_ex_C':    lambda val: self.ex_C.set_h(val),

            # Storage connector related inputs
            'sto_fluid': lambda val: self.sto_fluid.set_fluid(val),

            # W connector inputs
            'N_rot':     lambda val: self.W.set_N(val),

            # Q_amb connector inputs
            'T_amb':     lambda val: self.Q_amb.set_T_cold(val),
        }

        unknown_keys = []  # To collect any keys that do not match the input methods

        for key, value in self.inputs.items():
            method = input_methods.get(key)
            if method:
                try:
                    method(value)
                except Exception as e:
                    # Optionally log the exception or raise with more context
                    pass  # Replace with logging if desired
            else:
                unknown_keys.append(key)

        if unknown_keys:
            raise ValueError(f"Unrecognized input keys: {', '.join(unknown_keys)}")
        return

    def sync_inputs(self):
        """Synchronize the inputs dictionary with the connector states."""

        # Lazy getters: only access if the connector exists
        attribute_map = {
            # su connectors
            'fluid':     lambda: self.su.fluid if hasattr(self,'su') else None,
            'T_su':      lambda: self.su.T if hasattr(self,'su') else None,
            'h_su':      lambda: self.su.h if hasattr(self,'su') else None,
            'P_su':      lambda: self.su.p if hasattr(self,'su') else None,
            'm_dot':     lambda: self.su.m_dot if hasattr(self, 'su') else None,

            # su_H connector
            'fluid_H':   lambda: self.su_H.fluid if hasattr(self,'su_H') else None,
            'T_su_H':    lambda: self.su_H.T if hasattr(self,'su_H') else None,
            'h_su_H':    lambda: self.su_H.h if hasattr(self,'su_H') else None,
            'P_su_H':    lambda: self.su_H.p if hasattr(self,'su_H') else None,
            'm_dot_H':   lambda: self.su_H.m_dot if hasattr(self,'su_H') else None,

            # su_C connector
            'fluid_C':   lambda: self.su_C.fluid if hasattr(self,'su_C') else None,
            'T_su_C':    lambda: self.su_C.T if hasattr(self,'su_C') else None,
            'h_su_C':    lambda: self.su_C.h if hasattr(self,'su_C') else None,
            'P_su_C':    lambda: self.su_C.p if hasattr(self,'su_C') else None,
            'm_dot_C':   lambda: self.su_C.m_dot if hasattr(self,'su_C') else None,

            # ex connector
            'P_ex':      lambda: self.ex.p if hasattr(self,'ex') else None,
            'T_ex':      lambda: self.ex.T if hasattr(self,'ex') else None,
            'h_ex':      lambda: self.ex.h if hasattr(self,'ex') else None,

            # ex_C connector
            'P_ex_C':    lambda: self.ex_C.p if hasattr(self,'ex_C') else None,
            'T_ex_C':    lambda: self.ex_C.T if hasattr(self,'ex_C') else None,
            'h_ex_C':    lambda: self.ex_C.h if hasattr(self,'ex_C') else None,

            # ex_H connector
            'P_ex_H':    lambda: self.ex_H.p if hasattr(self,'ex_H') else None,
            'T_ex_H':    lambda: self.ex_H.T if hasattr(self,'ex_H') else None,
            'h_ex_H':    lambda: self.ex_H.h if hasattr(self,'ex_H') else None,

            # sto_fluid connector
            'sto_fluid': lambda: self.sto_fluid.fluid if hasattr(self, 'sto_fluid') else None,

            # W connector
            'N_rot':     lambda: self.W.N if hasattr(self,'W') else None,

            # Q_amb connector
            'T_amb':     lambda: self.Q_amb.T_cold if hasattr(self,'Q_amb') else None,
        }

        self.inputs = getattr(self,'inputs',{})

        for key, getter in attribute_map.items():
            try:
                value = getter()
                if value is not None:
                    self.inputs[key] = value
            except Exception:
                pass  # Optional: add logging for debugging

        return

    def mute_print(self):
        self.print_flag = 0
        return

    def print_setup(self):
        self.sync_inputs()
        print("=== Component Setup ===")
        print("\nInputs:")
        for input in self.get_required_inputs():
            if input in self.inputs:
                print(f"  - {input}: {self.inputs[input]}")
            else:
                print(f"  - {input}: Not set")


        print("\nParameters:")
        for param in self.get_required_parameters():
            if param in self.params:
                print(f"  - {param}: {self.params[param]}")
            else:
                print(f"  - {param}: Not set")

        print("======================")

    def set_parameters(self, **parameters):
        for key, value in parameters.items():
            self.params[key] = value

    def set_guesses(self, **guesses):
        for key, value in guesses.items():
            self.guesses[key] = value

    def check_calculable(self):
        self.sync_inputs()
        required_inputs = self.get_required_inputs() 
        
        self.calculable = all(self.inputs.get(inp) is not None for inp in required_inputs) # check if all required inputs are set
        if not self.calculable:
            if self.print_flag:
                print(f"Component {self.__class__.__name__} is not calculable. Missing inputs: {', '.join([inp for inp in required_inputs if self.inputs.get(inp) is None])}")
        return self.calculable

    def check_parametrized(self):
        required_params = self.get_required_parameters()
        self.parametrized = all(self.params.get(param) is not None for param in required_params) # check if all required parameters are set
        if not self.parametrized:
            if self.print_flag:
                print(f"Component {self.__class__.__name__} is not parametrized. Missing parameters: {', '.join([param for param in required_params if self.params.get(param) is None])}")
        return self.parametrized

    def get_required_inputs(self):
        # This method should be overridden in derived classes
        return []

    def get_required_parameters(self):
        # This method should be overridden in derived classes
        return []
    
    def get_required_guesses(self):

        return []

    def solve(self):
        # This method should be overridden in derived classes
        raise NotImplementedError("The 'solve' method should be implemented in derived classes.")
    

    # def plot_component(self, inputs_names=None, outputs_names=None, parameters_names=None):
    #     """
    #     Plot a visual representation of the component with inputs/outputs and parameters.
    #     """

    #     # Enable LaTeX rendering in Matplotlib
    #     plt.rcParams['text.usetex'] = True

    #     # Create figure and axis
    #     fig, ax = plt.subplots(figsize=(8, 8))

    #     # Draw the block (component)
    #     block = patches.FancyBboxPatch((0.4, 0.6), 0.3, 0.1, boxstyle="round,pad=0.1", 
    #                                     edgecolor="black", facecolor="#D4E9C7", zorder=2)
    #     ax.add_patch(block)
    #     ax.text(0.55, 0.65, self.__class__.__name__, horizontalalignment='center', 
    #             verticalalignment='center', fontsize=20, fontweight='bold')

    #     # Define positions for inputs, outputs, and parameters
    #     input_pos = [0.2, 0.65]  # Position for input
    #     output_pos = [0.8, 0.65]  # Position for output
    #     param_pos = [0.55, 0.4]  # Position for parameters

    #     # Inputs
    #     if inputs_names:
    #         ax.text(input_pos[0] - 0.15, input_pos[1] + 0.2, r'\textbf{Inputs}', fontsize=17, fontweight='bold', color='black')
    #         ax.annotate('', xy=(input_pos[0] + 0.1, input_pos[1]), xytext=(input_pos[0] - 0.05, input_pos[1]),
    #                     fontsize=12, color='black',
    #                     arrowprops=dict(facecolor='black', edgecolor='black', width=0.5, headwidth=8, shrink=0.05))

    #          # Calculate spacing for inputs
    #         num_inputs = len(inputs_names)
    #         vertical_spacing = 0.1
    #         start_y = input_pos[1] - (num_inputs - 1) * vertical_spacing / 2

    #         # List all inputs next to the arrow
    #         for i, input_name in enumerate(inputs_names):
    #             y = start_y + i * vertical_spacing
    #             ax.text(input_pos[0] - 0.15, y, f'$\\mathbf{{{input_name}}}$', fontsize=17, color='black', verticalalignment='center')


    #     # Outputs
    #     if outputs_names:
    #         ax.text(output_pos[0] + 0.15, output_pos[1] + 0.2, r'\textbf{Outputs}', fontsize=17, fontweight='bold', color='black')
    #         ax.annotate('', xy=(output_pos[0] + 0.15, output_pos[1]), xytext=(output_pos[0], output_pos[1]),
    #                     fontsize=12, color='black',
    #                     arrowprops=dict(facecolor='black', edgecolor='black', width=0.5, headwidth=8, shrink=0.05))

    #         # Calculate spacing for outputs
    #         num_outputs = len(outputs_names)
    #         vertical_spacing = 0.1
    #         start_y = output_pos[1] - (num_outputs - 1) * vertical_spacing / 2

    #         # List all outputs next to the arrow
    #         for i, output_name in enumerate(outputs_names):
    #             y = start_y + i * vertical_spacing
    #             ax.text(output_pos[0] + 0.2, y, f'$\\mathbf{{{output_name}}}$', fontsize=17, color='black', verticalalignment='center')


    #     # Parameters
    #     if parameters_names:
    #         ax.text(param_pos[0], param_pos[1], r'\textbf{Parameters}', fontsize=17, fontweight='bold', color='black', horizontalalignment='center')
    #         ax.annotate('', xy=(param_pos[0], param_pos[1] + 0.1), xytext=(param_pos[0], param_pos[1] + 0.05),
    #                     fontsize=12, color='black',
    #                     arrowprops=dict(facecolor='black', edgecolor='black', width=0.5, headwidth=8, shrink=0.05))

    #         # Calculate spacing for parameters
    #         num_params = len(parameters_names)
    #         vertical_spacing = 0.1
    #         start_y = param_pos[1] - (num_params - 1) * vertical_spacing / 2

    #         # List all parameters below "Parameters"
    #         for i, param_name in enumerate(parameters_names):
    #             y = start_y - i * vertical_spacing
    #             ax.text(param_pos[0], y-0.05, f'$\\mathbf{{{param_name}}}$', fontsize=17, color='black', verticalalignment='center', horizontalalignment='center')

    #     # Plot formatting
    #     ax.set_xlim(0, 1)
    #     ax.set_ylim(0, 1)
    #     ax.axis('off')  # Hide axis


    # def plot_connectors(self, supply_connectors_names=None, exhaust_connectors_names=None):

    #     """
    #     Plot a visual representation of the component with inputs/outputs and parameters.
    #     """

    #     # Enable LaTeX rendering in Matplotlib
    #     plt.rcParams['text.usetex'] = True

    #     # Create figure and axis
    #     fig, ax = plt.subplots(figsize=(8, 8))

    #     # Draw the block (component)
    #     block = patches.FancyBboxPatch((0.4, 0.6), 0.3, 0.1, boxstyle="round,pad=0.1", 
    #                                     edgecolor="black", facecolor="skyblue", zorder=2)
    #     ax.add_patch(block)
    #     ax.text(0.55, 0.65, self.__class__.__name__, horizontalalignment='center', 
    #             verticalalignment='center', fontsize=20, fontweight='bold')

    #     # Define positions for inputs, outputs, and parameters
    #     input_pos = [0.2, 0.65]  # Position for input
    #     output_pos = [0.8, 0.65]  # Position for output

    #     # Inputs
    #     if supply_connectors_names:
    #         ax.annotate('', xy=(input_pos[0] + 0.1, input_pos[1]), xytext=(input_pos[0] - 0.05, input_pos[1]),
    #                     fontsize=12, color='black',
    #                     arrowprops=dict(facecolor='black', edgecolor='black', width=0.5, headwidth=8, shrink=0.05))

    #          # Calculate spacing for inputs
    #         num_inputs = len(supply_connectors_names)
    #         vertical_spacing = 0.1
    #         start_y = input_pos[1] - (num_inputs - 1) * vertical_spacing / 2

    #         # List all inputs next to the arrow
    #         for i, input_name in enumerate(supply_connectors_names):
    #             y = start_y + i * vertical_spacing
    #             ax.text(input_pos[0] - 0.15, y, f'$\\mathbf{{{input_name}}}$', fontsize=17, color='black', verticalalignment='center')


    #     # Outputs
    #     if exhaust_connectors_names:
    #         ax.annotate('', xy=(output_pos[0] + 0.15, output_pos[1]), xytext=(output_pos[0], output_pos[1]),
    #                     fontsize=12, color='black',
    #                     arrowprops=dict(facecolor='black', edgecolor='black', width=0.5, headwidth=8, shrink=0.05))

    #         # Calculate spacing for outputs
    #         num_outputs = len(exhaust_connectors_names)
    #         vertical_spacing = 0.1
    #         start_y = output_pos[1] - (num_outputs - 1) * vertical_spacing / 2

    #         # List all outputs next to the arrow
    #         for i, output_name in enumerate(exhaust_connectors_names):
    #             y = start_y + i * vertical_spacing
    #             ax.text(output_pos[0] + 0.2, y, f'$\\mathbf{{{output_name}}}$', fontsize=17, color='black', verticalalignment='center')


    #     # Plot formatting
    #     ax.set_xlim(0, 1)
    #     ax.set_ylim(0, 1)
    #     ax.axis('off')  # Hide axis




