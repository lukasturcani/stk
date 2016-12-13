from functools import partial

class Normalization:
    def __init__(self, func_data):
        self.scaling_func = partial(getattr(self, func_data.name),
                                    **func_data.params)
                                    
    def __call__(self, population):
        self.scaling_func(population)        
        
    @staticmethod
    def cage(population):
        pass