import numpy as np
import eco_module as eco
from random import shuffle
import matplotlib.pyplot as plt

#specifies the directory and a small plant pollinator network
#for the network used in generating most figures in the paper, use filename = 'M_PL044.csv'
directory = 'large_poll_nets'
filename = 'M_PL_009.csv'

#constructs the network and calculates predicted survival probabilities
simpleGraph,primarySpecies,secondarySpecies = eco.graph_maker(directory, filename)
randomSurvival = eco.random_extinctions(simpleGraph,primarySpecies,secondarySpecies,sensitivityRatio = 1)
survivalXSpecies = [x/len(primarySpecies) for x in range(len(primarySpecies)+1)]

#simulates extinctions on the network for 1000 iterations
simulationArray = []
for iteration in range(1000):
    removalGraph = simpleGraph.copy()
    shuffle(primarySpecies)
    localSurvival = [1-[d for v,d in removalGraph.degree(secondarySpecies)].count(0)/len(secondarySpecies)]
    for primary in primarySpecies:
        removalGraph.remove_node(primary)
        localSurvival.append(1-[d for v,d in removalGraph.degree(secondarySpecies)].count(0)/len(secondarySpecies))
    simulationArray.append(localSurvival)
    
simulatedSurvival = np.average(simulationArray,axis = 0)

#plots the predictions and simulated results for random extinctions, similar to Figure 2(a)
plt.figure()
plt.plot(survivalXSpecies,simulatedSurvival,markersize=0,linewidth=2,label = 'Simulated Extinctions, 1000 iterations')
plt.plot(survivalXSpecies,randomSurvival,markersize=0,linewidth=2,color = 'cyan', linestyle = ':', label = 'Predicted Extinctions')
plt.xlim(0,1)
plt.ylim(0,1)
plt.legend(loc='lower left')
plt.xlabel(r'$\frac{\varphi}{N_p}$',fontsize=16)
plt.ylabel('P(survive)',fontsize=12)

#predicts extinctions targeted according to degree value
ascendTargSurvival = eco.targeted_extinctions(simpleGraph,primarySpecies,secondarySpecies,sensitivityRatio = 1,descending = False)
descendTargSurvival = eco.targeted_extinctions(simpleGraph,primarySpecies,secondarySpecies,sensitivityRatio = 1,descending = True)

#plots targeted extinctions, similar to Figure 4
plt.figure()
plt.plot(survivalXSpecies,ascendTargSurvival,markersize=0,linewidth=2,label = 'Low to high removal')
plt.plot(survivalXSpecies,descendTargSurvival,markersize=0,linewidth=2,color = 'orange',label = 'High to low removal')
plt.xlim(0,1)
plt.ylim(0,1.01)
plt.legend(loc='lower left')
plt.xlabel(r'$\frac{\varphi}{N_p}$',fontsize=16)
plt.ylabel('P(survive)',fontsize=12)

#predicts random extinctions on the network with 100 false positive and 100 false negative edges
falsePositiveSurvival = eco.random_extinctions(simpleGraph,primarySpecies,secondarySpecies,sensitivityRatio = 1,falsePositives = True, falseEdges = 100)
falseNegativeSurvival = eco.random_extinctions(simpleGraph,primarySpecies,secondarySpecies,sensitivityRatio = 1,falseNegatives = True, falseEdges = 100)

#plots predictions for the original network and the network with false positives and negatives, similar to Figure 8(a)
plt.figure()
plt.plot(survivalXSpecies,randomSurvival,markersize=0,linewidth=2, label = 'Original Network')
plt.plot(survivalXSpecies,falsePositiveSurvival,markersize=0,linewidth=2,color = 'orange',label = '100 False Positive Edges')
plt.plot(survivalXSpecies,falseNegativeSurvival,markersize=0,linewidth=2,color = 'green',label = '100 False Negative Edges')
plt.xlim(0,1)
plt.ylim(0,1.01)
plt.legend(loc='lower left')
plt.xlabel(r'$\frac{\varphi}{N_p}$',fontsize=16)
plt.ylabel('P(survive)',fontsize=12)

#constructs the network with multi-edges to reflect variable interaction strength and predicts interaction strength sensitive extinctions at depth = 2
multiGraph,primarySpecies,secondarySpecies = eco.graph_maker(directory, filename,multiEdges = True)
interactionSpeciesSurvival = eco.interaction_strength_species_extinctions(multiGraph,primarySpecies,secondarySpecies,sensitivityRatio = 0.5,depth = 2)
survivalXInteractions = [x/multiGraph.number_of_edges() for x in range(multiGraph.number_of_edges()+1)]

#plots predictions for interaction strength sensitive extinctions, similar to Figure 9(a)
plt.figure()
plt.plot(survivalXSpecies,interactionSpeciesSurvival,markersize=0,linewidth=2,color = 'orange',label = 'T = 0.5 Predicted, Depth = 2')
plt.xlim(0,1)
plt.ylim(0,1)
plt.legend(loc='lower left')
plt.xlabel(r'$\frac{\varphi}{N_p}$',fontsize=16)
plt.ylabel('P(survive)',fontsize=12)

#predicts survival probabilities for the scenario in which interaction strength is lost, similar to Figure (11)
interactionLossSurvival = eco.interaction_loss_extinctions(multiGraph,primarySpecies,secondarySpecies,sensitivityRatio = 0.5)

#plots predictions for interaction loss extinctions
plt.figure()
plt.plot(survivalXInteractions,interactionLossSurvival,markersize=0,linewidth=2,color = 'orange',label = 'T = 0.5')
plt.xlim(0,1)
plt.ylim(0,1)
plt.legend(loc='lower left')
plt.xlabel(r'$\frac{\varphi}{E}$',fontsize=16)
plt.ylabel('P(survive)',fontsize=12)