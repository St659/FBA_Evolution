import csv
import matplotlib.pyplot as plt
from cobrafba import CobraFBAEvolver
from collections import Counter
import numpy as np
import os



filename = '1_450.csv'

directory = os.getcwd()
evolver = CobraFBAEvolver(os.path.join(directory, 'iRH826.json'))
fig, ax = plt.subplots()
fig2, ax2 = plt.subplots()

with open(filename, 'r') as csv_file:
    reader = csv.reader(csv_file, delimiter=',')
    removed_reactions= list()
    for row in reader:
        int_list = [int(num) for num in row]
        growth, reactions = evolver.fitness_function(int_list)
        for reaction, ind in zip(evolver.fba.non_essential_reactions, int_list):
            if ind:
                pass
            else:
                removed_reactions.append(reaction.id)

        ax.plot(reactions, growth, 'o')

counted_reactions = Counter(removed_reactions)

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
ax.set_xlabel('Reactions')
ax.set_ylabel('Growth')
plt.show()