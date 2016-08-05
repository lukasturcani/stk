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
                
                kw = line[0]
                
                if '_select_func' in kw:
                                     
                    name = line[1]                    
                    params = {}
                    for i, word in enumerate(line):
                        if i >= 2:
                            word = word.split("=")
                            params[word[0]] = word[1]
                    if 'gen' in kw:
                        self.gen_select_func = FunctionData(name, **params)
                    if 'parent' in kw:
                        self.mating_select_func = FunctionData(name, **params)
                    if 'mutant' in kw:
                        self.mutant_select_func = FunctionData(name, **params)
                        
                    