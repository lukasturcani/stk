import os
import shutil
import sys

from .classes import (Population, GATools, Selection, Mutation, Mating, 
                      FunctionData, FourPlusSix, EightPlusTwelve, 
                      GAInput)
                      
a = GAInput(sys.argv[1])
print(a.gen_select_func.name, a.gen_select_func.params)
print(a.mating_select_func.name, a.mating_select_func.params)
print(a.mutant_select_func.name, a.mutant_select_func.params) 
#for name in os.listdir():
#    if name == 'output':
#        shutil.rmtree('output')
#
#print(sys.argv)
#
#os.mkdir('output')
#os.chdir('output')
#root_dir = os.getcwd()
#os.mkdir('initial')
#os.chdir('initial')
#gen_select_func = FunctionData('fittest', size=10)
#mat_select_func = FunctionData('all_combinations')
#mut_select_func = FunctionData('fittest', size=3)
#selector = Selection(gen_select_func, mat_select_func, mut_select_func)
#
#mat_func = FunctionData('bb_lk_exchange')
#mator = Mating(mat_func)
#
#bb_db = r'C:\Users\lukas\Projects\MMEA\Database_prec\amines3f'
#lk_db = r'C:\Users\lukas\Projects\MMEA\Database_prec\aldehydes2f'
#
#mut_func = FunctionData('random_bb', database=bb_db)
#mutator = Mutation(mut_func)
#ga_tools = GATools(selector, mator, mutator)
#topologies = [EightPlusTwelve, FourPlusSix]
#size = 20
#pop = Population.init_random_cages(bb_db, lk_db, topologies, size, ga_tools)
#
#num_gens = 10
#num_muts = 2
#
#
#for x in range(0, num_gens):
#    os.chdir(root_dir)
#    os.mkdir(str(x))
#    os.chdir(str(x))
#    print(0)
#    offspring = pop.gen_offspring()
#    print(1, len(pop))
#    mutants = pop.gen_mutants()
#    print(2)
#    pop += offspring + mutants
#    print(3)
#    pop = pop.select('generational')