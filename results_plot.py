import os
import matplotlib.pyplot as plt
import numpy as np
import csv
from FBA_Evolver import FBAPlotter

#igure = plt.figure()
#ax = figure.add_subplot(111)
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
phenotype_filenames = list()
fitness_filenames = list()
for p_file, f_file in zip( os.listdir('./Results_21_5_17/Phenotype'), os.listdir('./Results_21_5_17/Fitness')):

    phenotype_filenames.append('./Results_21_5_17/Phenotype/' + p_file)
    fitness_filenames.append('./Results_21_5_17/Fitness/' + f_file)

plotter = FBAPlotter(phenotype_filenames, 'E:\\Chrome Download\\FBA\\FBA\\iPRAE34_plot_test.xls')
filtered_name, filtered_count = plotter.filter_reactions(1)

for file in fitness_filenames:
    with open(file, 'r') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            ax2.plot(np.asarray(next(reader)))

fig = plt.figure()
ax = fig.add_subplot(111)
ind = np.arange(len(filtered_name))
width = 0.5
ax.bar(ind, filtered_count, width)
ax.set_xticks(ind + width / 2)
ax.set_xticklabels(filtered_name, rotation='vertical')
ax2.set_xlabel('Generation')
ax2.set_ylabel('Fitness')
plt.tick_params(axis='both', which='major', labelsize=8)
plt.show()