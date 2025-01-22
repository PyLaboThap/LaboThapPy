# -*- coding: utf-8 -*-
"""
Created on Wed Jul 07 11:47:52 2024
    
@author: elise neven
"""

from connector.mass_connector import MassConnector
from connector.work_connector import WorkConnector
from connector.heat_connector import HeatConnector

from CoolProp.CoolProp import PropsSI

class Circuit:
    class Component:
        def __init__(self, name, model, fluid=None):
            self.name = name
            self.model = model
            self.previous = {}  # Dictionary to store preceding components by connector
            self.next = {}  # Dictionary to store succeeding components by connector
            self.fluid = fluid

        def add_previous(self, connector, component):
            self.previous[connector] = component 

        def add_next(self, connector, component):
            self.next[connector] = component

        def link(self, output_connector, target_component, input_connector):
            # Determine the type of connector based on the output connector
            connector_type = output_connector.split('-')[0]

            if connector_type == "m":  # Mass connector
                connector = MassConnector(fluid=self.fluid)
            elif connector_type == "q":  # Heat connector
                connector = HeatConnector()
            elif connector_type == "w":  # Work connector
                connector = WorkConnector()
            else:
                raise ValueError(f"Unknown connector type: {connector_type}")

            # Set the connector in both the source and target components
            setattr(self.model, output_connector.split('-')[1], connector)
            setattr(target_component.model, input_connector.split('-')[1], connector)

            # Add the connection to the tracking dictionaries
            self.add_next(output_connector, target_component)
            target_component.add_previous(input_connector, self)

            # print(f"Linked {self.name}.{output_connector} to {target_component.name}.{input_connector}")

        def set_properties(self, connector_name, **kwargs):
            # Set properties for a specific connector
            connector = getattr(self.model, connector_name)
            connector.set_properties(**kwargs)

        def solve(self):
            # Solve the component if it is calculable
            self.model.check_calculable()
            if self.model.calculable:
                self.model.solve()
            else:
                print(f"{self.name} is not calculable")
                self.model.print_states_connectors()

    def __init__(self, fluid=None):
        self.components = {}  # Store components using a dictionary for easy access
        self.fixed_properties = {}
        self.guesses = {}
        self.parameters = {}
        self.fluid = fluid
        self.residuals_var = []  # Store the residuals to check
        self.n_it = 0
        self.tolerance = 1e-6  # Convergence tolerance for residuals

    def add_component(self, model, name):
        # Add a component to the cycle
        component = Circuit.Component(name, model, self.fluid)
        self.components[name] = component

    def get_component(self, name):
        # Retrieve a component by name
        # print(f"Name: {name}")
        if name in self.components:
            return self.components[name]
        raise ValueError(f"Component '{name}' not found")

    def link_components(self, component1_name, output_connector, component2_name, input_connector):
        # Link two components through specified connectors
        component1 = self.get_component(component1_name)
        component2 = self.get_component(component2_name)
        component1.link(output_connector, component2, input_connector)

    def set_cycle_properties(self, **kwargs):
        # Set properties for a specific component and connector
        target = kwargs.pop('target')
        component_name, connector_name = target.split(':')
        component = self.get_component(component_name)
        component.set_properties(connector_name, **kwargs)

        # Update the fixed properties for the cycle
        if target not in self.fixed_properties:
            self.fixed_properties[target] = {}
        self.fixed_properties[target].update(kwargs)

    def set_cycle_guess(self, **kwargs):
        # Set initial guesses for a specific component and connector
        target = kwargs.pop('target')
        component_name, connector_name = target.split(':')
        component = self.get_component(component_name)
        component.set_properties(connector_name, **kwargs)
        
        # Update the guesses for the cycle
        if target not in self.guesses:
            self.guesses[target] = {}
            # print("problem here?")
        self.guesses[target].update(kwargs)
        # print(f"Guesses: {self.guesses}")

    def set_cycle_parameters(self, **kwargs):
        # Set parameters for the cycle
        self.parameters.update(kwargs)

    def set_start_key(self, start_key):
        # Set the starting component for the cycle
        self.start_key = start_key

    def solve(self, start_key=None):
        # self.set_cycle_guesses()

        max_iterations = 100
        converged = False
        
        while not converged and self.n_it < max_iterations:
            print(f"Iteration {self.n_it}: Solving the cycle.")

            # Calculate residuals before solving
            self.prev_residuals = self.get_residuals()

            # Recursively solve the cycle starting from the pump
            visited = set()  # Reset visited components for each iteration
            self.recursive_solve(self.get_component(self.start_key), visited)

            # Calculate residuals after solving
            final_residuals = self.get_residuals()

            # Update guesses for the next iteration
            self.update_guesses()

            # Check if the residuals are within the tolerance
            if self.n_it > 0:
                converged = self.check_residuals(final_residuals)

            self.n_it += 1

        if converged:
            print("Cycle solved successfully.")
        else:
            print("Cycle did not converge within the maximum number of iterations.")

    def update_guesses(self):
        """
        Update the guesses dynamically based on the current state of components,
        considering the nested structure of the 'self.guesses' dictionary.
        """
        for target, properties in self.guesses.items():
            # Parse the target to identify component and connector
            component_name, connector_name = target.split(':')
            # Retrieve the component and connector
            component = self.get_component(component_name)
            connector = getattr(component.model, connector_name)

            # Iterate over each property in the nested dictionary
            for property_name in properties.keys():
                # Retrieve the current value of the property from the connector
                new_value = getattr(connector, property_name)
                
                # Update the property in self.guesses
                self.set_cycle_guess(target=target, **{property_name: new_value})

    def recursive_solve(self, component, visited):
        # Prevent infinite loops by skipping already visited components
        if component in visited:
            return

        # Mark the component as visited and solve it
        visited.add(component)
        print('COUOU')
        print(f"Solving {component.name}")

        if isinstance(component, Circuit.Component):
            component.solve()
            # component.model.print_results()

        # Recursively solve connected components
        for next_component in component.next.values():
            self.recursive_solve(next_component, visited)


    def check_residuals(self, final_residuals):
        """
        Check if the difference between the previous and current residuals is within the tolerance.
        """
        if not self.prev_residuals:
            return False  # No previous residuals to compare to
    
        # Output residuals for debugging
        
        # for i, (f, p) in enumerate(zip(final_residuals, self.prev_residuals)):
        #     key = self.residuals_var[i]
        #     print(f"Key: {key}, Final Residual (f): {f}, Previous Residual (p): {p}")

        residual_diff = [abs(f - p) for f, p in zip(final_residuals, self.prev_residuals)]
        
        # Output residuals for debugging
        for i, diff in enumerate(residual_diff):
            print(f"Residual {self.residuals_var[i]}: {diff}")

        return all(diff < self.tolerance for diff in residual_diff)

    def get_residuals(self):
        """
        Calculate the residuals based on the specified variables.
        """
        
        residuals = [] # Store the calculated residuals
        for residual_target in self.residuals_var: # Iterate over the residuals to calculate
            component_name, connector_prop = residual_target.split(':') # Split the target into component and connector
            connector_name, prop_name = connector_prop.split('-') # Split the connector into name and property
            component = self.get_component(component_name) # Get the component object
            residual_value = getattr(getattr(component.model, connector_name), prop_name) # Get the value of the property from the connector
            residuals.append(residual_value) # Append the calculated residual to the list
           

        return residuals
    
    def set_cycle_guesses_residuals(self, guesses, residuals_var):
        """
        Set the initial guesses and define the variables used to compute the residuals.
        """
        self.residuals_var = residuals_var
        # Set initial guesses for the cycle based on provided values
        for guess, value in guesses.items():
            component_name, prop = guess.split(':')
            connector_name, prop_name = prop.split('-')
            self.set_cycle_guess(target=f"{component_name}:{connector_name}", **{prop_name: value})

