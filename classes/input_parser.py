from .ga import FunctionData

class GAInput:
    def __init__(self, input_file):
        self.input_file = input_file
        self.extract_data()
        
    def extract_data(self):
        with open(self.input_file, 'r') as input_file:
            for line in input_file:
                line = line.split()
                if not line:
                    continue
                if 'generational_select_func' in line:
                    params = {key : value for key, value in line[:]}
                    self.gen_select_func = FunctionData(line[1], )