from component.base_component import BaseComponent
from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector


class Template(BaseComponent):
    def __init__(self):
        super().__init__() # Call the constructor of the parent class BaseComponent
        # Define connectors here
        self.su = MassConnector()
        self.ex = MassConnector() # Mass_connector
        self.work_su = WorkConnector()
        self.heat_ex = HeatConnector()

    def get_required_inputs(self): # Used in check_calculablle to see if all of the required inputs are set
        self.sync_inputs()
        # Return a list of required inputs
        return ['su_p', 'su_T', 'ex_p', 'su_fluid'] 
    
    def sync_inputs(self):
        """Synchronize the inputs dictionary with the connector states."""
        "This function is used to update the inputs dictionary with the current state of the connectors when the inputs are set through the connectors"
        if self.su.fluid is not None:
            self.inputs['su_fluid'] = self.su.fluid
        if self.su.T is not None:
            self.inputs['su_T'] = self.su.T
        elif self.su.h is not None:
            self.inputs['su_h'] = self.su.h
        if self.su.p is not None:
            self.inputs['su_p'] = self.su.p
        if self.ex.p is not None:
            self.inputs['ex_p'] = self.ex.p

    def set_inputs(self, **kwargs):
        """Set inputs directly through a dictionary and update connector properties."""
        self.inputs.update(kwargs) # This line merges the keyword arguments ('kwargs') passed to the 'set_inputs()' method into the eisting 'self.inputs' dictionary.

        # Update the connectors based on the new inputs
        if 'su_fluid' in self.inputs:
            self.su.set_fluid(self.inputs['su_fluid'])
        if 'su_T' in self.inputs:
            self.su.set_T(self.inputs['su_T'])
        elif 'su_h' in self.inputs:
            self.su.set_h(self.inputs['su_h'])
        if 'su_p' in self.inputs:
            self.su.set_p(self.inputs['su_p'])
        if 'ex_p' in self.inputs:
            self.ex.set_p(self.inputs['ex_p'])

    def get_required_parameters(self): # Used in check_parametrized to see if all of the required parameters are set
        return ['param1', 'param2', 'param3']
    
    def print_setup(self):
        # Print the setup of the component (the connectors and the inputs)

        print("=== Component Setup ===")
        print("Connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T}")
        print(f"  - ex: ...")
        print(f"  - W_dot: ...")
        print(f"  - Q_dot_amb: ...")

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

    def System(self):
        # Define the system of equations to solve
        # This function should return the residuals of the equations
        # The residuals should be stored in self.res
        pass

    def solve(self):
        # Perform checks to ensure the model can be calculated and has parameters
        self.check_calculable()
        self.check_parametrized()

        if not (self.calculable and self.parametrized):
            self.solved = False # If the model is not calculable or parametrized, set the solved flag to False and return
            return

        try:
            """COMPONENT MODEL"""
            # Here is the model in itself (to put in another function)

            """Update connectors after the calculations"""
            self.update_connectors()

            # Mark the model as solved if successful
            self.solved = True
        except Exception as e:
            # Handle any errors that occur during solving
            self.solved = False
            print(f"Convergence problem in pump model: {e}")

        
    def update_connectors(self, h_ex, w_pp):
        """Update the connectors with the calculated values."""
        self.ex.set_h(h_ex)
        self.ex.set_fluid(self.su.fluid)
        self.ex.set_m_dot(self.su.m_dot)
        self.W_pp.set_w(w_pp)
    
    def print_results(self):
        print("=== Pump Results ===")
        print(f"  - h_ex: {self.ex.h} [J/kg]")
        print(f"  - T_ex: {self.ex.T} [K]")
        print(f"  - w_pp: {self.W_pp.w} [J/kg]")
        print("=========================")

    def print_states_connectors(self):
        print("=== Pump States ===")
        print("Mass connectors:")
        print(f"  - su: fluid={self.su.fluid}, T={self.su.T} [K], p={self.su.p} [Pa], h={self.su.h} [J/kg], s={self.su.s} [J/K.kg], m_dot={self.su.m_dot} [kg/s]")
        print(f"  - ex: fluid={self.ex.fluid}, T={self.ex.T} [K], p={self.ex.p} [Pa], h={self.ex.h} [J/kg], s={self.ex.s} [J/K.kg], m_dot={self.ex.m_dot} [kg/s]")
        print("=========================")
        print("Work connector:")
        print(f"  - W_pp: w={self.W_pp.w} [J/kg]")
        print("=========================")


#-------------------------------------------------------------
if __name__ == '__main__':
    # Example usage
    template = Template()

    #Print the set up
    template.print_setup()

    "If the inputs are not set directly BUT throught the connectors"
    # Set fluid states for connectors
    template.su.set_fluid('R134a')
    template.ex.set_fluid('R134a')

    # Set properties for su connector
    template.su.set_T(300)

    # Set properties for ex connector
    
    # Set properties for work connector

    # Set properties for heat connector

    "If the inputs are set directlys"
    # Setting inputs
    template.set_inputs(
        su_T=300,
        su_fluid='R134a'  # Make sure to include fluid information
    )
    
    # Setting parameters
    template.set_parameters(
        param1=1, param2=2, param3=3
    )

    # Solve the component
    template.solve()
    template.print_results()
    template.print_states_connectors()

