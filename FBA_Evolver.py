
from py4j.java_gateway import JavaGateway
import xlrd
from xlutils.copy import copy
import random
import subprocess
import numpy

from deap import base
from deap import creator
from deap import tools
from deap import algorithms
import csv
from orderedset import OrderedSet


class FBA_Base():

    def __init__(self, input):
        self.model= input
        print(self.model)
        self.start = 0
        self.o2 = 0
        self.num_exchange_reactions = 0
        self.total_reactions = 0
        self.reaction_flux = list()
        workbook = xlrd.open_workbook(self.model)
        self.num_exchange_reactions, self.total_reactions, self.o2, self.reaction_names, self.reactions = self.find_number_of_exchanges(workbook.sheet_by_name('Reactions'))
        print(self.reaction_names)
        self.total_compounds, self.compound_names, self.compounds = self.find_number_of_compounds(workbook.sheet_by_name('Compounds'))

    def find_number_of_compounds(self, ws):
        nrows = ws.nrows
        name_rows = ws.col_slice(0, end_rowx=nrows)
        compound_rows = ws.col_slice(1, end_rowx=nrows)
        compound_names = list()
        compounds = list()
        for name_row, compound_row in zip(name_rows, compound_rows):
            compound_names.append(name_row.value)
            compounds.append(str(int(compound_row.value)))
        return nrows, compound_names, compounds


    def find_number_of_exchanges(self, ws):
        nrows = ws.nrows
        reactions_name_rows = ws.col_slice(0,end_rowx=nrows)
        reactions_rows = ws.col_slice(1,end_rowx=nrows)
        num_exchange_reactions = 0
        reaction_names = list()
        reactions = list()

        exchange_start_num = False
        i = 0
        for  reaction_name_row, reaction_row in zip(reactions_name_rows, reactions_rows):
            reaction_names.append(reaction_name_row.value)
            reactions.append(reaction_row.value)
            self.reaction_flux.append([ws.cell(i,2).value,ws.cell(i,3).value])
            if 'EX_' in reaction_name_row.value:

                if not exchange_start_num:
                    self.start = i + 1

                    exchange_start_num = True
                if 'EX_o2(e)' in reaction_name_row.value:

                    o2_pos = i + 1
                num_exchange_reactions += 1
            i+=1
        return num_exchange_reactions, nrows, o2_pos, reaction_names, reactions

    def set_exchanges(self, values):

        rb = xlrd.open_workbook(self.model, formatting_info=True)
        r_sheet = rb.sheet_by_name('Reactions')  # read only copy to introspect the file
        wb = copy(rb)  # a writable copy (I can't read values out of this, only write to it)
        w_sheet = wb.get_sheet(1)  # the sheet to write to within the writable copy
        names = r_sheet.col_slice(0)
        for i, set_exchange in enumerate(values):
            if int(set_exchange):

                w_sheet.write(i + self.start-1  , 2, -1000)
            else:
                w_sheet.write(i + self.start-1, 2, 0)

        if values[self.o2  - self.start]:
            w_sheet.write(self.o2 -1, 2, -20)

        wb.save(self.model)

class FBAEvolver(FBA_Base):

    def __init__(self, input):
        FBA_Base.__init__(self, input)
        self.toolbox = 0


    def create_evo(self, num_reactions):
        creator.create("FitnessMin", base.Fitness, weights=(-1,))
        creator.create("Individual", list, fitness=creator.FitnessMin)
        self.toolbox = base.Toolbox()
        # Attribute generator
        self.toolbox.register("attr_bool", random.randint, 0, 1)
        # Structure initializers
        self.toolbox.register("individual", tools.initRepeat, creator.Individual,
                         self.toolbox.attr_bool, num_reactions)
        self.toolbox.register("population", tools.initRepeat, list, self.toolbox.individual)
        self.toolbox.register("evaluate", self.fitness_function)
        self.toolbox.register("mate", tools.cxTwoPoint)
        self.toolbox.register("mutate", tools.mutFlipBit, indpb=0.1)
        self.toolbox.register("select", tools.selBest)
        self.toolbox.register("HallOfFame", tools.HallOfFame)
        #self.toolbox.register("clone", self.clone_mu)

    def set_prob(self, prob):
        self.toolbox.register("mutate", tools.mutFlipBit, indpb=prob)

    def write_constant_data(self):
        with open('data.csv', 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter = ',')
            writer.writerow(self.compound_names)
            writer.writerow(self.compounds)
            writer.writerow(self.reaction_names)
            writer.writerow(self.reactions)





    def clone_mu(self, mu, pop_len):
        offspring = list()
        for i in range(0, pop_len -1):
            offspring.append(copy.deep_copy(mu))
        return offspring

    def fitness_function(self,individual):
        self.set_exchanges(individual)
        growth = self.run_fba(self.model, './out.xls')

        if growth > 10:
            fitness = sum(individual)
        else:
            for i, val in enumerate(individual):
                individual[i] = 1
            fitness = sum(individual)

        return fitness,
    def run_fba(self, input, output):
        print(len(self.reaction_names))
        print(self.total_reactions)
        p1 = subprocess.Popen(['java', '-cp',
                               "E:\\Chrome Download\\FBA\\FBA\\bin;E:\\Chrome Download\\jexcelapi_2_6_12\\jexcelapi\\jxl.jar;E:\\Chrome Download\\winglpk-4.61\\glpk-4.61\\w64\\glpk-java.jar;E:\\Chrome Download\\winglpk-4.61\\glpk-4.61\\w64\\glpk-java-sources.jar;E:\\Chrome Download\\winglpk-4.61\\glpk-4.61\\w64",
                               "fba.FBAEvolution", "E:\\Chrome Download\\FBA\\FBA\\iPRAE34.xls", "output.xls", str(self.total_compounds), str(self.total_reactions),
                               "./data.csv"],
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        growth_raw = str(p1.communicate()[0])


        growth_first= growth_raw.split('\'')[1]
        growth = growth_first.split('\\')[0]
        print(growth)
        return float(growth)

    def run_evo(self, iteration):
        random.seed()
        min_fits = list()

        pop = self.toolbox.population(n=10)
        CXPB, MUTPB, NGEN = 0.3, 0.8, 2000

        print("Start of evolution")

        # Evaluate the entire population
        fitnesses = list(map(self.toolbox.evaluate, pop))
        print(len(pop))
        for ind, fit in zip(pop, fitnesses):
            ind.fitness.values = fit

        print("  Evaluated %i individuals" % len(pop))
        filename = str(iteration) + '.csv'
        final_filename = str(iteration) + '_final.csv'
        # Begin the evolution
        for g in range(NGEN):
            print("-- Generation %i --" % g)

            # Select the next generation individuals
            mu = self.toolbox.select(pop, 1)
            print(len(mu))
            # Clone the selected individuals
            lamda = [self.toolbox.clone(mu[0]) for x in range(0, len(pop)-1)]
            print(lamda[0].fitness.values[0])
            #Apply crossover and mutation on the offspring
            # for child1, child2 in zip(offspring[::2], offspring[1::2]):
            #     if random.random() < CXPB:
            #         self.toolbox.mate(child1, child2)
            #         del child1.fitness.values
            #         del child2.fitness.values

            for mutant in lamda:
                if random.random() < MUTPB:
                    self.toolbox.mutate(mutant)
                    del mutant.fitness.values

            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in lamda if not ind.fitness.valid]
            print(len(invalid_ind))
            fitnesses = map(self.toolbox.evaluate, invalid_ind)
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit

            print("  Evaluated %i individuals" % len(invalid_ind))

            # The population is entirely replaced by the offspring
            pop[:] = [item for sublist in [mu,lamda] for item in sublist]
            

            # Gather all the fitnesses in one list and print the stats
            fits = [ind.fitness.values[0] for ind in pop]

            length = len(pop)
            mean = sum(fits) / length
            sum2 = sum(x * x for x in fits)
            std = abs(sum2 / length - mean ** 2) ** 0.5
            min_fits.append(min(fits))
            print("  Min %s" % min(fits))
            print("  Max %s" % max(fits))
            print("  Avg %s" % mean)
            print("  Std %s" % std)


        print("-- End of (successful) evolution --")


        with open(final_filename, 'w') as csvfile:
            final_writer = csv.writer(csvfile, delimiter=',')
            best_ind = tools.selBest(pop, 1)[0]
            print("Best individual is %s, %s" % (best_ind, best_ind.fitness.values))
            final_writer.writerow(min_fits)
            final_writer.writerow(best_ind)

    def run_mu_lambda(self):

        MU, LAMBDA = 2, 10
        pop = self.toolbox.population(n=MU)
        hof = tools.ParetoFront()
        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("avg", numpy.mean, axis=0)
        stats.register("std", numpy.std, axis=0)
        stats.register("min", numpy.min, axis=0)
        stats.register("max", numpy.max, axis=0)
        algorithms.eaMuPlusLambda(pop, self.toolbox, mu=MU, lambda_=LAMBDA,
                                  cxpb=0.2, mutpb=0.6, ngen=1000,
                                  stats=stats, halloffame=hof)
        print("Min:" + str(stats.min))

class FBAPlotter(FBA_Base):
    def __init__(self, input, model):
        FBA_Base.__init__(self,model)
        self.active_exchange = list()
        self.active_exchange_count = list()
        self.fitnesses = list()
        self.phenotypes = list()

        if isinstance(input, str):
            self.num_files = 1
        else:
            self.read_files(input)
            self.num_files = len(input)
            self.unique = self.get_unique_reactions(self.phenotypes)
            self.get_num_of_used_reactions()

    def read_files(self, files):
        for file in files:
            with open(file, 'r') as csvfile:

                reader = csv.reader(csvfile, delimiter=',')
                row_count = sum(1 for row in reader)
                csvfile.seek(0)

                if row_count > 2:
                   self.fitnesses.append(next(reader))
                   #This is to skip the blank line between the rows
                   next(reader)

                   self.phenotypes.append(next(reader))
                else:
                    self.phenotypes.append(next(reader))


    def filter_reactions(self, num):
        filtered_name = list()
        filtered_count = list()
        for name, count in zip(OrderedSet(self.active_exchange), self.active_exchange_count):
            if int(count) > num:

                filtered_count.append(count)
                filtered_name.append(name)
        return filtered_name, filtered_count


    def get_unique_reactions(self, input):

        for reactions in input:

            reaction_ints = [int(react) for react in reactions]
            self.set_exchanges(reaction_ints)
            workbook = xlrd.open_workbook(self.model)
            ws = workbook.sheet_by_name('Reactions')

            reaction_values = ws.col_slice(2, start_rowx=self.start -1, end_rowx=self.start + self.num_exchange_reactions -1)
            reaction_names = ws.col_slice(0, start_rowx=self.start -1, end_rowx=self.start + self.num_exchange_reactions -1)

            for reaction, name in zip(reaction_values, reaction_names):
                if int(reaction.value):
                    self.active_exchange.append(name.value)

        return OrderedSet(self.active_exchange)

    def get_num_of_used_reactions(self):
        for name in OrderedSet(self.active_exchange):
            self.active_exchange_count.append(self.active_exchange.count(name))

if __name__ == "__main__":
    input = './iPRAE34_ppc_knockout.xls'
    evolver = FBAEvolver(input)
    num_reactions, start, o2 = evolver.find_number_of_exchanges()
    #best = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0]
    #evolver.set_exchanges(best)
    evolver.create_evo(num_reactions)
    bit_fit_prob = [0.01, 0.01,0.01]
    for i, prob in enumerate(bit_fit_prob):

        evolver.set_prob(prob)
        evolver.run_evo(i)
