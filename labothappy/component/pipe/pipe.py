
from labothappy.component.base_component import BaseComponent
from labothappy.connector.mass_connector import MassConnector
# from labothappy.connector.heat_connector import HeatConnector

from labothappy.component.pipe.straight_pipe import StraightPipe
from labothappy.component.pipe.curved_elbow import CurvedElbow
import CoolProp.CoolProp as CP
import numpy as np

class Pipe(BaseComponent):
    """
    Composite pipe: chain of straight sections and bends.

    Internally uses StraightPipe and CurvedElbow components,
    automatically connects them and solves sequentially.
    
    Computes total refrigerant charge inventory by summing
    charges from all segments.
    """

    def __init__(self):
        super().__init__()
        self.su = MassConnector()  # Mass connector for the suction side
        self.ex = MassConnector()  # Mass connector for the exhaust side
        self.segments = []
        self.components = []  # Store component instances for diagnostics

    def add_straight(self, D, L, K=0.0, theta=0.0):
        """Add a straight pipe section."""
        self.segments.append({
            'type': 'straight',
            'D': D, 'L': L, 'K': K, 'theta': theta
        })
        return self  # ← Returns the Pipe object itself for chaining

    def add_curved_elbow(self, D, delta, R0_D):
        """Add a curved elbow."""
        self.segments.append({
            'type': 'curved_elbow',
            'D': D, 'delta': delta, 'R0_D': R0_D
        })
        return self
    
    def solve(self):
        """Solve all segments sequentially and accumulate charge and pressure drop."""
        current_state = self.su
        self.components = []  # Reset component list
        self.m_charge = 0.0  # Initialize total charge
        self.deltaP_total = 0.0  # Initialize total pressure drop
        
        for i, seg in enumerate(self.segments):
            if seg['type'] == 'straight':
                component = StraightPipe()
                component.set_parameters(D=seg['D'], L=seg['L'], K=seg['K'], theta=seg['theta'])
            elif seg['type'] == 'curved_elbow':
                component = CurvedElbow()
                component.set_parameters(D=seg['D'], delta=seg['delta'], R0_D=seg['R0_D'])

            # Connect and solve
            component.su = current_state
            component.solve()
            current_state = component.ex
            
            # Accumulate charge and pressure drop
            self.m_charge += component.m_charge
            self.deltaP_total += component.deltaP
            
            # Store component for diagnostics
            self.components.append(component)

        # Set outlet
        self.ex = current_state
        self.solved = True


    def print_results(self):
        """Print summary of all segments and total charge."""
        print("=" * 70)
        print("COMPOSITE PIPE RESULTS")
        print("=" * 70)
        print(f"Number of segments: {len(self.components)}")
        print(f"\n{'Segment':<10} {'Type':<12} {'ΔP (Pa)':<15} {'m_charge (kg)':<15}")
        print("-" * 70)
        
        for i, comp in enumerate(self.components):
            seg_type = "Straight" if isinstance(comp, StraightPipe) else "Curved Elbow"
            dp = comp.deltaP
            charge = comp.m_charge
            
            print(f"{i+1:<10} {seg_type:<12} {dp:<15.2f} {charge:<15.6f}")
        
        print("-" * 70)
        print(f"{'TOTAL':<10} {'':<12} {self.deltaP_total:<15.2f} {self.m_charge:<15.6f}")
        print("=" * 70)



