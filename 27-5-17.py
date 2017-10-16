import os
import matplotlib.pyplot as plt
import numpy as np
import csv
from FBA_Evolver import FBAPlotter, FBAEvolver

filenames = list()
for f_file in os.listdir('./27-5-17'):
    filenames.append('./27-5-17/' + f_file)
plotter= FBAPlotter(filenames, 'E:\\Chrome Download\\FBA\\FBA\\iPRAE34_ppc_knockout.xls')

evolver = FBAEvolver('E:\\Chrome Download\\FBA\\FBA\\iPRAE34_ppc_knockout.xls')

filtered_name, filtered_count = plotter.filter_reactions(0)


for phenotype in plotter.phenotypes:
    print(phenotype)
    evolver.set_exchanges(phenotype)
    growth = evolver.run_fba('E:\\Chrome Download\\FBA\\FBA\\iPRAE34_ppc_knockout.xls', 'E:\\Chrome Download\\FBA\\FBA\\iPRAE34_ppc_knockout.xls')
    print(growth)

fig = plt.figure()
ax = fig.add_subplot(111)

ind = np.arange(len(filtered_name))
width = 0.5
ax.bar(ind, filtered_count, width)
ax.set_xticks(ind + width / 2)
ax.set_xticklabels(filtered_name, rotation='vertical')
plt.gcf().subplots_adjust(bottom=0.3)
plt.savefig(os.path.join('./27-5-17/','media metabolites 26-5-17.png'), dpi=300)
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
for fitness in plotter.fitnesses:
    ax2.plot(fitness)

ax2.set_xlabel('Generation')
ax2.set_ylabel('Fitness')
plt.tick_params(axis='both', which='major', labelsize=8)
plt.show()