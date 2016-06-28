import numpy as np
import convenience_functions as cf
from functools import wraps


class Population:
    def __init__(self, *args):
        self.populations = []
        self.members = []
    
        for arg in args:
            pass
        
    def __iter__(self):
        for pop in self.populations:
            yield pop
        
        for member in self.members:
            yield member
            
a 
    
    