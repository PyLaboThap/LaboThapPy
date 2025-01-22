from CoolProp.CoolProp import PropsSI

from labothappy.machine.circuit import Circuit
from labothappy.machine.boundary_conditions.mass_source import MassSource
from labothappy.machine.boundary_conditions.mass_sink import MassSink

from labothappy.component.heat_exchanger.steady_state.pinch_cst.simulation_model import HXPinchCst
from labothappy.component.volumetric_machine.expander.steady_state.constant_isentropic_efficiency.simulation_model import ExpanderCstEff
from labothappy.component.pump.steady_state.constant_efficiency.simulation_model import PumpCstEff

class System:
    def __init__(self, tolerance=1e-6, max_iterations=100):
        self.circuits = {}  # Store all circuits in the system
        self.main_circuit = None  # Main circuit to start solving
        self.guesses = {}  # Initial guesses for the entire system
        self.tolerance = tolerance  # Convergence tolerance
        self.n_it = 0  # Iteration counter
        self.max_iterations = max_iterations  # Maximum number of iterations
    
    def add_circuit(self, name, circuit):
        """Add a circuit to the system with a unique name."""
        if name in self.circuits:
            raise ValueError(f"Circuit name '{name}' is already used.")
        self.circuits[name] = circuit
    
    def set_main_circuit(self, circuit_name):
        """Set the main circuit to start solving."""
        print(f"Setting main circuit to '{circuit_name}'.")
        print(f"Available circuits: {list(self.circuits)}")
        if circuit_name not in self.circuits:
            raise ValueError(f"Circuit '{circuit_name}' not found.")
        self.main_circuit = self.circuits[circuit_name]
    
    def set_system_guesses(self, guesses):
        """Set initial guesses for the entire system."""
        self.guesses.update(guesses)
    
    def solve(self):
        """Solve the entire system of interconnected circuits."""
        if not self.main_circuit:
            raise ValueError("Main circuit is not set.")
        
        converged = False
        
        while not converged and self.n_it < self.max_iterations:
            print(f"Iteration {self.n_it}: Solving the system.")
            
            # Step 1: Set initial guesses
            self.apply_guesses()
            
            # Step 2: Start solving from the main circuit
            visited_circuits = set()  # Track solved circuits
            self.recursive_solve(self.main_circuit, visited_circuits)
            
            # Step 5: Check system residuals
            final_residuals = self.get_residuals()
            if self.n_it > 0:
                converged = self.check_residuals(final_residuals)
            
            # Step 6: Update guesses
            self.update_guesses()
            
            self.n_it += 1
        
        if converged:
            print("System solved successfully.")
        else:
            print("System did not converge within the maximum number of iterations.")
    
    def apply_guesses(self):
        """Apply initial guesses to all circuits and components."""
        for target, properties in self.guesses.items():
            circuit_name, component_prop = target.split(':')
            component_name, prop_name = component_prop.split('-')
            circuit = self.circuits[circuit_name]
            circuit.set_cycle_guess(target=f"{component_name}:{prop_name}", **properties)
    
    def recursive_solve(self, circuit, visited_circuits):
        """Recursively solve interconnected circuits."""
        if circuit in visited_circuits:
            return  # Avoid redundant solving
        
        visited_circuits.add(circuit)
        print(f"Solving circuit {list(self.circuits.keys())[list(self.circuits.values()).index(circuit)]}")
        
        for component in circuit.components.values():
            component.solve()
            
            # Check for shared circuits and recursively solve them
            for next_component in component.next.values():
                for shared_circuit_name, shared_circuit in self.circuits.items():
                    if shared_circuit is not circuit and next_component in shared_circuit.components.values():
                        self.recursive_solve(shared_circuit, visited_circuits)
    
    def get_residuals(self):
        """Calculate residuals for the entire system."""
        residuals = []
        for circuit in self.circuits.values():
            residuals.extend(circuit.get_residuals())
        return residuals
    
    def check_residuals(self, final_residuals):
        """Check if the system's residuals are within tolerance."""
        if not hasattr(self, "prev_residuals"):
            self.prev_residuals = [0] * len(final_residuals)
        
        residual_diff = [abs(f - p) for f, p in zip(final_residuals, self.prev_residuals)]
        self.prev_residuals = final_residuals  # Update for the next iteration
        
        return all(diff < self.tolerance for diff in residual_diff)
    
    def update_guesses(self):
        """Update guesses dynamically for the system."""
        for circuit in self.circuits.values():
            circuit.update_guesses()

                

if __name__ == "__main__":

    #-------------------------------------------------
    # Define the ORC cycle
    orc_cycle = Circuit(fluid='R1233zd(E)')

    # Add components
    PUMP = PumpCstEff()
    EVAP = HXPinchCst()
    EXP = ExpanderCstEff()
    COND = HXPinchCst()

    # Set component parameters
    PUMP.set_parameters(eta_is=0.6)
    EVAP.set_parameters(Pinch=5, Delta_T_sh_sc=5, type_HX='evaporator')
    COND.set_parameters(Pinch=5, Delta_T_sh_sc=5, type_HX='condenser')
    EXP.set_parameters(eta_is=0.8)

    # Add components to the cycles
    orc_cycle.add_component(PUMP, "Pump")
    orc_cycle.add_component(EVAP, "Evaporator")
    orc_cycle.add_component(EXP, "Expander")
    orc_cycle.add_component(COND, "Condenser")

    # Link components
    orc_cycle.link_components("Pump", "m-ex", "Evaporator", "m-su_C")
    orc_cycle.link_components("Evaporator", "m-ex_C", "Expander", "m-su")
    orc_cycle.link_components("Expander", "m-ex", "Condenser", "m-su_H")
    orc_cycle.link_components("Condenser", "m-ex_H", "Pump", "m-su")

    # Set the cycle properties
    orc_cycle.set_cycle_properties(m_dot=0.4, target='Pump:su')

    # orc_cycle.set_cycle_properties(T=15 + 273.15, fluid='Water', m_dot=2, target='Condenser:su_C', p = 4e5)
    # orc_cycle.set_cycle_properties(cp=4186, target='Condenser:su_C')
    # orc_cycle.set_cycle_properties(fluid='Water', target='Condenser:ex_C')

    orc_cycle.set_cycle_properties(T=90 + 273.15, fluid='Water', m_dot=1, target='Evaporator:su_H', p = 4e5)
    orc_cycle.set_cycle_properties(cp=4186, target='Evaporator:su_H')
    orc_cycle.set_cycle_properties(fluid='Water', target='Evaporator:ex_H')

    # Set parameters for the cycle
    SC_cd = 5
    SH_ev = 5
    orc_cycle.set_cycle_parameters(SC_cd=SC_cd, SH_ev=SH_ev)

    # Initial guesses for pressures
    T_ev_guess = 90 + 273.15
    P_ev_guess = PropsSI('P', 'T', T_ev_guess, 'Q', 0.5, 'R1233zd(E)')

    T_cd_guess = 15 + 273.15
    P_cd_guess = PropsSI('P', 'T', T_cd_guess, 'Q', 0.5, 'R1233zd(E)')

    # Define guesses and residuals
    guesses = {
        "Evaporator:ex_C-p": P_ev_guess,
        "Condenser:ex_H-p": P_cd_guess,
        "Pump:su-T": T_cd_guess - SC_cd,
        "Pump:ex-p": P_ev_guess,
        "Expander:ex-p": P_cd_guess, 
        "Evaporator:ex_C-T": T_ev_guess + SH_ev
    }

    residuals_var = [
        "Evaporator:ex_C-h",
        "Condenser:ex_H-h"
    ]

    # orc_cycle.set_cycle_guesses_residuals(guesses, residuals_var)

    # Solve the cycle
    orc_cycle.set_start_key("Pump")
    # orc_cycle.solve(start_key="Pump")

    #-------------------------------------------------
    # Define the cold water loop

    cold_water_loop = Circuit(fluid='Water')

    # Add components to the cycle
    cold_water_loop.add_component(COND, "Condenser")

    # Add sink and source components
    mass_sink = MassSink() 
    mass_source = MassSource()
    # I couldn't set the properties direclty in the mass source because then all informations are erased when we link the components together.

    # Add components to the cycle
    cold_water_loop.add_component(mass_source, "ColdSource")
    cold_water_loop.add_component(mass_sink, "ColdSink")

    # Link components
    cold_water_loop.link_components("ColdSource", "m-ex", "Condenser", "m-su_C")
    cold_water_loop.link_components("Condenser", "m-ex_C", "ColdSink", "m-su")

    # Set the properties of the mass source
    mass_source.set_properties(T=15 + 273.15, fluid='Water', m_dot=2, p=1e5)
    print("Condenser supply temp: ", COND.su_C.T)
    print("Mass source temp: ", mass_source.ex.T)
    # Set the cycle properties
    # cold_water_loop.set_cycle_properties(m_dot=2, target='Pump:su')
    # cold_water_loop.set_cycle_properties(T=20 + 273.15, fluid='Water', target='Dry Cooler:su_C', P = 1e5)

    # FAIRE AVEC MASS SINK ET MASS SOURCE
    cold_water_loop.set_start_key("Condenser")
    # When we set the properties of the mass source, the properties don't pass to the inlet of the condenser
    # SYSTEM DEFINITION
    ORC_system = System()
    ORC_system.add_circuit(orc_cycle, "ORC")
    ORC_system.add_circuit(cold_water_loop, "ColdWater")

    ORC_system.set_main_circuit("ORC")
    ORC_system.set_system_guesses(guesses)

    # Solve the system
    ORC_system.solve()



# Right now the residuals are on the circuit ORC, they should be on the system
# SHould they be several solve for each circuits or only one solve for the entire system?
# Indeed, the residuals are on the system not on the circuit level, as well as the guesses.
# Should the circuit just be a way to connect components while the real recursivity and iteration be on 
# the system level?
# Also, it's once we solve the entire system that we can update the guesses,
# not after each circuit.
# System should also be able to work with only one circuit.
# EN PARLER AVEC BASILE DEMAIN!!!!