import matplotlib.pyplot as plt

class GAProgress:
    def __init__(self):
        self.gens = []
        self.means = []
        self.mins = []
        self.maxs = []

    def update(self, population):
        self.gens.append(len(self.gens))
        self.means.append(population.mean('fitness'))
        self.maxs.append(max(x.fitness for x in population))
        self.mins.append(min(x.fitness for x in population))
        
    def epp(self, plot_name):
        fig = plt.figure()
        plt.plot(self.gens, self.means, color='green')
        plt.plot(self.gens, self.mins, color='blue')
        plt.plot(self.gens, self.maxs, color='red')
        fig.savefig(plot_name, dpi=fig.dpi)
        plt.close('all')