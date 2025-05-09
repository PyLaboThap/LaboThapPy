from labothappy.connector.mass_connector import MassConnector

class MassSource:
    def __init__(self):
        self.ex = MassConnector()

    def set_properties(self, **properties):
        # Set the properties of the source
        self.ex.set_properties(**properties)

    def check_calculable(self):
        self.calculable = True # We need to set it here as it is not really a component but still needs a 'calculable' to be treated as one

    def solve(self):
        pass
  
    # def get_properties(self):
    #     # Return the properties of the sink
    #     return self.ex.get_properties()
    