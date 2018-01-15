import csv
import matplotlib.pyplot as plt
from cobrafba import CobraFBAEvolver
from collections import Counter
import numpy as np
import os
import collections


def pareto_pop_plot(filename, ax,color,label):
    with open(filename, 'r') as csv_file:
        reader = csv.reader(csv_file, delimiter=',')
        removed_reactions = list()
        for row in reader:
            int_list = [int(num) for num in row]
            growth, reactions = evolver.fitness_function(int_list)
            for reaction, ind in zip(evolver.fba.non_essential_reactions, int_list):
                if ind:
                    pass
                else:
                    removed_reactions.append(reaction.id)

            handle = ax.plot(reactions, growth, 'o',c=color,label=label)
        return removed_reactions

filename = '1_450.csv'

directory = os.getcwd()
evolver = CobraFBAEvolver(os.path.join(directory, 'iRH826.json'))
for reaction in evolver.fba.non_essential_reactions:
    if 'PPC' in reaction.id:
        print(reaction.id)
fig, ax = plt.subplots()
fig2, ax2 = plt.subplots()

gens = ['100','200','300','400','500']
colors = {'100':'xkcd:blue','200':'xkcd:red','300':'xkcd:green','400':'xkcd:black','500':'xkcd:violet'}
removed_reactions = list()
for directory in ['blood_1','blood_2','blood_3','blood_4','blood_5']:

    print(directory)
    files = os.listdir(directory)
    for file in files:
        gen_csv = file.split('_')
        gen_file = gen_csv[-1].split('.')
        if gen_file[0] in gens:
            print(gen_file)
            color = colors[gen_file[0]]
            if directory == 'blood_5':
                label = gen_file[0]
            else:
                label = ""

            removed_reactions.append(pareto_pop_plot(os.path.join(directory, file), ax, color, label))

flat_removed_reactions= [item for sublist in removed_reactions for item in sublist]

counted_reactions = Counter(flat_removed_reactions)

for k,v in counted_reactions.items():
    print(str(k) + " " + str(v))

common_reactions = dict((k,v) for k,v in counted_reactions.items() if v>=50)
print("Total Reactions Removed: " + str(len(counted_reactions)))

print("PPC Removed: " + str(counted_reactions["PPC"]))
ind = np.arange(len(common_reactions))
width = 1
ax2.bar(ind, common_reactions.values(), width)
ax2.set_xticks(ind + width / 2)
ax2.set_xticklabels(common_reactions.keys(), rotation='vertical')


#ax.set_ylim([0.7,1])
#ax.set_xlim([500,850])

handles, labels =ax.get_legend_handles_labels()
by_label = collections.OrderedDict(zip(labels, handles))
ax.legend(by_label.values(), by_label.keys())
ax.set_xlabel('Reactions')
ax.set_ylabel('Growth')
# print(by_label.values())

plt.show()