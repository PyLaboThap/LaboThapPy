
from labothappy.component.base_component import BaseComponent
from labothappy.connector.mass_connector import MassConnector
from labothappy.component.pipe.straight_pipe import StraightPipe
from labothappy.component.pipe.curved_elbow import CurvedElbow

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

    def get_required_inputs(self):
        return ['P_su', 'h_su', 'm_dot', 'fluid']

    def add_straight(self, D, L, K=0.0, theta=0.0):
        """Add a straight pipe section."""
        self.segments.append({
            'type': 'straight',
            'D': D, 'L': L, 'K': K, 'theta': theta
        })
        return self  # ← Returns the Pipe object itself for chaining

    def add_curved_elbow(self, D, delta, R0):
        """Add a curved elbow."""
        self.segments.append({
            'type': 'curved_elbow',
            'D': D, 'delta': delta, 'R0': R0
        })
        return self


    def solve(self):
        """Solve all segments sequentially and accumulate charge and pressure drop."""

        self.check_calculable()
        
        current_state = self.su
        # print(f"Starting Pipe solve with inlet state: p={current_state.p:.2f} Pa, h={current_state.h:.2f} J/kg, m_dot={current_state.m_dot:.4f} kg/s")
        self.components = []
        self.m_charge = 0.0
        self.dP_total = 0.0
        
        for i, seg in enumerate(self.segments):
            if seg['type'] == 'straight':
                component = StraightPipe()
                component.set_parameters(D=seg['D'], L=seg['L'], K=seg['K'], theta=seg['theta'])
            elif seg['type'] == 'curved_elbow':
                component = CurvedElbow()
                component.set_parameters(D=seg['D'], delta=seg['delta'], R0=seg['R0'])

            # Connect inlet
            component.su = current_state
            
            # Solve
            component.solve()
            current_state = component.ex
            
            # Accumulate
            self.m_charge += component.m_charge
            self.dP_total += component.dP
            self.components.append(component)

        # Set outlet
        self.ex.set_fluid(self.su.fluid)
        self.ex.set_m_dot(self.su.m_dot)
        self.ex.set_h(self.su.h)
        self.ex.set_p(self.su.p - self.dP_total)
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
            dp = comp.dP
            charge = comp.m_charge
            
            print(f"{i+1:<10} {seg_type:<12} {dp:<15.2f} {charge:<15.6f}")
        
        print("-" * 70)
        print(f"{'TOTAL':<10} {'':<12} {self.deltaP_total:<15.2f} {self.m_charge:<15.6f}")
        print("=" * 70)



