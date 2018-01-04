import cobra
import os
from jamieexcel_to_json import cobrajsonfromexcel

# file = 'iGLOSS22_1_edit.xlsx'
# directory = 'E:\Chrome Download'
# cobrajsonfromexcel(file)
# directory = 'E:\\Git\\FBA_Evolution\\'
# json_file = 'iGLOSS22_1_edit.json'
# json_file_2 = 'iJO1366.json'

def evaluate_essential_reactions(model):
    #model = cobra.io.load_json_model(os.path.join(directory,json_file))
    #model.objective = objective

    print(len(model.metabolites))
    print(len(model.reactions))
    reactions_list = list()
    reaction_bounds = list()

    for reaction in model.reactions:
        if 'EX' in reaction.id[:2]:
            print('Passed')
            pass
        else:
            print('Added')
            reactions_list.append(reaction)
            reaction_bounds.append([reaction.lower_bound,reaction.upper_bound])


    non_essential_reactions = list()
    essential_reactions = list()

    for reaction, reaction_bound in zip(reactions_list, reaction_bounds):

        growth = model.slim_optimize()
        print('Growth: ' + str(growth))
        model.reactions.get_by_id(reaction.id).lower_bound = 0
        model.reactions.get_by_id(reaction.id).upper_bound = 0
        growth = model.slim_optimize()
        print('Knockout Growth: ' + str(growth))
        print(reaction.id)
        if growth > 0.001:
            non_essential_reactions.append(reaction)
        else:
            essential_reactions.append(reaction)

        model.reactions.get_by_id(reaction.id).lower_bound = reaction_bound[0]
        model.reactions.get_by_id(reaction.id).upper_bound = reaction_bound[1]

    print('Number of essential genes: ' + str(len(essential_reactions)))
    print('Numer of non essential genes' + str(len(non_essential_reactions)))

    return non_essential_reactions, essential_reactions