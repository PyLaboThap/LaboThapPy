from labothappy.connector.mass_connector import MassConnector

class MassSink:
    def __init__(self):
        self.su = MassConnector()

    def check_calculable(self):
        self.calculable = True # Needs a calculable to be treated as a component even though it is not really one.

    
    def solve(self):
        # Process or "accept" incoming mass flow; optionally check fluid properties
        pass