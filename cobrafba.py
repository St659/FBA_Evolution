import unittest
import array
import numpy as np
import cobra
import random
import os
from deap import base
from deap import creator
from deap import tools
from deap import algorithms
from deap.benchmarks.tools import diversity, convergence, hypervolume
from jamieexcel_to_json import cobrajsonfromexcel
from fba_gene_min import evaluate_essential_reactions
import matplotlib.pyplot as plt
import csv

class TestFBAEvolver(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._evolver = CobraFBA('E:\\Chrome Download\\FBA\\FBA\\iJO1366.json')

    def test_biomass_amount(self):
        solution = self._evolver.solution
        self.assertAlmostEqual(solution, 0.9865144469529762)

    def test_reaction_values(self):
        model = self._evolver.model
        self.assertEqual(len(model.reactions), 2583)
        self.assertEqual(model.reactions[10].lower_bound, -1000)
        self.assertEqual(model.reactions[10].upper_bound, 1000)
        self.assertEqual(model.reactions[10].id, '12PPDStpp')

    def test_reaction_values(self):
        init_values = self._evolver.initial_reaction_bounds
        self.assertEqual(len(init_values), 2255)

    def test_Ex_removed_from_reactions(self):
        init_values = self._evolver.initial_reaction_bounds
        for reaction_id, reaction_bounds in init_values.items():
            self.assertFalse('EX' in reaction_id)

    def test_set_reactions(self):
        zero_individual = np.zeros(len(self._evolver.initial_reaction_bounds))
        self._evolver.set_reaction_bounds(zero_individual)

        for reaction in self._evolver.model.reactions:
            if 'EX' in reaction.id:
                pass
            else:
                self.assertEqual(reaction.lower_bound, 0 )
                self.assertEqual(reaction.upper_bound, 0)

        ones_individual = np.ones(len(self._evolver.initial_reaction_bounds))
        self._evolver.set_reaction_bounds(ones_individual)
        for reaction in self._evolver.model.reactions:
            if 'EX' in reaction.id:
                pass
            else:
                self.assertEqual(reaction.lower_bound, self._evolver.initial_reaction_bounds[reaction.id][0])
                self.assertEqual(reaction.upper_bound, self._evolver.initial_reaction_bounds[reaction.id][1])



class CobraFBA():
    def __init__(self, input):

        # import iaf1260
        self.model = cobra.io.load_json_model(input)
        self.model.objective = "Biomass"
        self.non_essential_reactions, self.essential_reactions = evaluate_essential_reactions(self.model)
        self.reaction_names, self.initial_reaction_bounds = self.find_reaction_intial_bounds(self.non_essential_reactions)
        print(len(self.non_essential_reactions))
        # aerobic growth under  18.5 glucose
        #self.model.reactions.get_by_id("EX_glc__D_e").lower_bound = -18.5
        #self.model.reactions.get_by_id("EX_o2_e").lower_bound = -20


    def find_reaction_intial_bounds(self, reactions):
        reaction_names = list()
        reaction_bounds = list()
        for reaction in reactions:
            if 'EX' not in reaction.id[:2]:
                reaction_names.append(reaction.id)
                reaction_bounds.append([reaction.lower_bound, reaction.upper_bound])
            else:
                reaction.lower_bound = -1000
                print(reaction.id + ' Lower Bound: ' + str(reaction.lower_bound))


        return reaction_names, dict(zip(reaction_names, reaction_bounds))
        #with open("solution.json", 'w') as f:
            #.write(solution.fluxes.to_json())
    def set_reaction_bounds(self, individual):
        for reaction_active, reaction in zip(individual, self.reaction_names):
            if reaction_active:
                bounds = self.initial_reaction_bounds[reaction]
                self.model.reactions.get_by_id(reaction).lower_bound = bounds[0]
                self.model.reactions.get_by_id(reaction).upper_bound = bounds[1]
            else:
                self.model.reactions.get_by_id(reaction).lower_bound = 0
                self.model.reactions.get_by_id(reaction).upper_bound = 0
    def run_fba(self):
        return self.model.slim_optimize()

class CobraFBAEvolver():
    def __init__(self, input):
        self.fba = CobraFBA(input)
        self.create_evo(len(self.fba.non_essential_reactions))
        ones_individual = np.ones(len(self.fba.initial_reaction_bounds))



    def fitness_function(self,individual):
        self.fba.set_reaction_bounds(individual)
        growth = self.fba.run_fba()
        if growth > 0.001:
            reactions = sum(individual)
        else:
            reactions = len(self.fba.non_essential_reactions)
        return growth,reactions

    def run_nsga2evo(self, seed =None):
        random.seed(seed)

        NGEN = 500
        MU = 100


        stats = tools.Statistics(lambda ind: ind.fitness.values)
        # stats.register("avg", numpy.mean, axis=0)
        # stats.register("std", numpy.std, axis=0)
        stats.register("min", np.min, axis=0)
        stats.register("max", np.max, axis=0)

        logbook = tools.Logbook()
        logbook.header = "gen", "evals", "std", "min", "avg", "max"

        try:
            pop = self.toolbox.population_restart()
        except FileNotFoundError:
            pop = self.toolbox.population(n=MU)


        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in pop if not ind.fitness.valid]
        fitnesses = self.toolbox.map(self.toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        # This is just to assign the crowding distance to the individuals
        # no actual selection is done
        pop = self.toolbox.select(pop, len(pop))

        record = stats.compile(pop)
        logbook.record(gen=0, evals=len(invalid_ind), **record)
        print(logbook.stream)

        # Begin the generational process
        for gen in range(1, NGEN):
            # Vary the population
            offspring = tools.selTournamentDCD(pop, len(pop))
            offspring = [self.toolbox.clone(ind) for ind in offspring]

            for ind1 in offspring:
                self.toolbox.mutate(ind1)
                del ind1.fitness.values

            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            fitnesses = self.toolbox.map(self.toolbox.evaluate, invalid_ind)
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit

            # Select the next generation population
            pop = self.toolbox.select(pop + offspring, MU)
            record = stats.compile(pop)
            logbook.record(gen=gen, evals=len(invalid_ind), **record)
            print(logbook.stream)
            with open('currentpop.csv', 'w') as csvfile:
                writer = csv.writer(csvfile, delimiter=',')
                for p in pop:
                    writer.writerow(p)


        print("Final population hypervolume is %f" % hypervolume(pop, [11.0, 11.0]))

        return pop, logbook

    def run_evo(self, iteration):
        random.seed()
        min_fits = list()


        pop = self.toolbox.population(n=50)

        #pop = [np.ones(len(self.fba.reaction_names)) for p in pop]

        CXPB, MUTPB, NGEN = 0, 1, 2000
        MU = 50
        LAMBDA = 100

        print("Start of evolution")

        # Evaluate the entire population
        fitnesses = list(map(self.toolbox.evaluate, pop))

        for ind, fit in zip(pop, fitnesses):
            ind.fitness.values = fit

        print("  Evaluated %i individuals" % len(pop))
        filename = str(iteration) + '.csv'
        final_filename = str(iteration) + '_final.csv'
        # Begin the evolution

        hof = tools.ParetoFront()
        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("avg", np.mean, axis=0)
        stats.register("std", np.std, axis=0)
        stats.register("min", np.min, axis=0)
        stats.register("max", np.max, axis=0)

        algorithms.eaMuPlusLambda(pop, self.toolbox, MU, LAMBDA, CXPB, MUTPB, NGEN, stats,
                                  halloffame=hof)


    def initInd(self, icls,content):
        return icls(content)
    def initPop(self,pcls,ind_init, filename):
        with open(filename,'r') as csv_file:
            initial_list = list()
            reader = csv.reader(csv_file, delimiter=',')
            for ind in reader:
                if ind:
                    initial_list.append([int(i) for i in ind])
            return pcls(ind_init(ind) for ind in initial_list)

    def create_evo(self, num_reactions):
        creator.create("FitnessMin", base.Fitness, weights=(1,-1))
        creator.create("Individual", list, fitness=creator.FitnessMin)
        self.toolbox = base.Toolbox()
        # Attribute generator
        self.toolbox.register("attr_bool", random.randint, 1, 1)
        # Structure initializers
        self.toolbox.register("individual", tools.initRepeat, creator.Individual,
                         self.toolbox.attr_bool, num_reactions)

        self.toolbox.register("population", tools.initRepeat, list, self.toolbox.individual)
        self.toolbox.register("individual_restart", self.initInd, creator.Individual)
        self.toolbox.register("population_restart", self.initPop, list, self.toolbox.individual_restart, 'currentpop.csv')
        self.toolbox.register("evaluate", self.fitness_function)
        self.toolbox.register("mate", tools.cxTwoPoint)
        self.toolbox.register("mutate", tools.mutFlipBit, indpb=0.0005)
        self.toolbox.register("select", tools.selNSGA2)
        self.toolbox.register("HallOfFame", tools.HallOfFame)

if __name__ == "__main__":
    directory = 'E:\\Git\\FBA_Evolution'
    #json_file = os.path.join(directory, 'iPRAE34_edit.json')
    #if not os.path.isfile(json_file):
    #    json_file = cobrajsonfromexcel(os.path.join(directory,'iPRAE34_edit.xlsx'),directory)
    evolver = CobraFBAEvolver(os.path.join(directory, 'iRH826.json'))
    pop, log = evolver.run_nsga2evo()

    fig, ax = plt.subplots()
    print(len(pop))
    for p in pop:
        growth,reactions = evolver.fitness_function(p)
        ax.plot(reactions,growth,'o')

    ax.set_xlabel('Reactions')
    ax.set_ylabel('Growth')
    plt.show()
