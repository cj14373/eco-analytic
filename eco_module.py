import numpy as np
import networkx as nx
import os
from scipy.stats import hypergeom
import math
from collections import Counter

#constructs the ecological network based on an interaction matrix
def graph_maker(directory,filename,multiEdges = False):
    
    #the data file is opened and the matrix is prepared
    with open(os.path.join(directory, filename)) as inputfile:
        lines = inputfile.readlines()
    matrix = [[int(n) for n in list(l.split(","))] for l in lines]
    matrix = np.array(matrix)
    
    #secondary and primary species are labelled numerically
    secondarySpecies = [k for k in range(len(matrix))]
    primarySpecies = [n+len(secondarySpecies) for n in range(len(matrix[0]))]
    
    #an edge list is constructed from the matrix. If multiEdges is True, multi-edges are allowed
    edgeList = []
    for secondary in secondarySpecies:
        for primary in primarySpecies:
            if multiEdges == True:
                for k in range(matrix[secondary][primary-len(matrix)]):
                    edgeList.append((secondary,primary))
            else:
                if matrix[secondary][primary-len(matrix)] != 0:
                    edgeList.append((secondary,primary))
    
    #the network is constructed, with species and interactions being added
    if multiEdges == True:
        graph = nx.MultiGraph()
    else:
        graph = nx.Graph()
        
    graph.add_nodes_from(secondarySpecies,bipartite=0)
    graph.add_nodes_from(primarySpecies,bipartite=1)
    graph.add_edges_from(edgeList)

    return graph,primarySpecies,secondarySpecies

#predicts random extinctions, and can include false positives or negatives, and can also calculate extinctions via interaction loss
def random_extinctions(graph,primarySpecies,secondarySpecies,sensitivityRatio,falsePositives = False,falseNegatives = False, falseEdges = 0):
        
    #measures the degree sequence and degree distribution of the secondary species
    secondaryDegreeSeq = [graph.degree(secondary) for secondary in secondarySpecies]
    degreeDist = [secondaryDegreeSeq.count(k)/len(secondarySpecies) for k in range(max(secondaryDegreeSeq)+1)]
    
    #includes false positives or negatives if specified
    if falsePositives == True and falseNegatives == True:
        exit('Error: Inclusion of both false positives and false negatives is not supported')
    elif falsePositives == True:
        degreeDist = false_positives(graph,primarySpecies,secondarySpecies,falseEdges)
    elif falseNegatives == True:
        degreeDist = false_negatives(graph,primarySpecies,secondarySpecies,falseEdges)
            
    predictSurvival = []
    
    #calculates survival probability as primary species are removed
    for removedSpecies in range(len(primarySpecies)+1):
        localSurviveProb = 0
        
        #iterates over the possible degree values for secondary species
        for degree in range(len(degreeDist)):
            
            #calculates the number of neighbours that must be lost for extinction to occur
            thresholdValue = math.ceil(sensitivityRatio*degree)
            
            if degreeDist[degree] > 0:
                
                #extinction probability is calculated using the hypergeometric distirbution and weighted by the degree distribution
                localSurviveProb += degreeDist[degree]*(hypergeom.cdf(thresholdValue-1,len(primarySpecies),degree,removedSpecies))
                    
        #survival probability is recorded
        predictSurvival.append(localSurviveProb)
    
    return predictSurvival
    
#calculates survival probabilities for targeted extinctions, here they are targeted according to degree
def targeted_extinctions(graph,primarySpecies,secondarySpecies,sensitivityRatio,descending = True):
    
    #measures the degree of primary species, counting the degree of each secondary species' neighbours and calculates extinction thresholds
    primaryDegreeSeq = [graph.degree(p) for p in primarySpecies]
    secondaryNeighbours = {s:[graph.degree(n) for n in graph.neighbors(s)] for s in secondarySpecies}
    secondaryNeighDegrees = {s:[secondaryNeighbours[s].count(d) for d in range(max(primaryDegreeSeq)+1)] for s in secondarySpecies}
    thresholdValues = {s:math.ceil(sensitivityRatio*graph.degree(s)) for s in secondarySpecies}
    
    #sorts removal sequence and counts species per group
    if descending == True:
        sortedPrimarySeq = sorted(primaryDegreeSeq,reverse=True)
    elif descending == False:
        sortedPrimarySeq = sorted(primaryDegreeSeq)
    primaryDegreeCount = [sortedPrimarySeq.count(degree) for degree in range(max(primaryDegreeSeq)+1)]
    
    predictSurvival = [1-[d for v,d in graph.degree(secondarySpecies)].count(0)/len(secondarySpecies)]
    removed = 0       
    
    #iterates over each primary degree value
    for degree in sortedPrimarySeq:
        
        #specifies the current degree value being removed
        removedDegree = sortedPrimarySeq[:removed+1].count(degree)
        localSurviveProb = 0
        
        #iterates over secondary species
        for secondary in secondarySpecies:
            threshold = thresholdValues[secondary]
            
            #calculates survival probabilities
            if descending == True:
                upperRemoved = sum(secondaryNeighDegrees[secondary][degree:])
                lowerRemoved = sum(secondaryNeighDegrees[secondary][degree+1:])
            elif descending == False:
                upperRemoved = sum(secondaryNeighDegrees[secondary][:degree+1])
                lowerRemoved = sum(secondaryNeighDegrees[secondary][:degree])
            if upperRemoved < threshold or removedDegree < threshold - lowerRemoved:
                localSurviveProb += float(1)/len(secondarySpecies)
            elif upperRemoved >= threshold and lowerRemoved < threshold and removedDegree >= threshold - lowerRemoved:
                localSurviveProb += hypergeom.cdf(threshold - lowerRemoved - 1,primaryDegreeCount[degree],secondaryNeighDegrees[secondary][degree],removedDegree)/len(secondarySpecies)
                
        #records survival probabilities
        predictSurvival.append(localSurviveProb)
        removed += 1
    
    return predictSurvival

#adjusts the degree distribution to include false positives
def false_positives(graph,primarySpecies,secondarySpecies,falseEdges):
    
    secondaryDegreeSeq = [graph.degree(secondary) for secondary in secondarySpecies]
    degreeDist = [secondaryDegreeSeq.count(k)/len(secondarySpecies) for k in range(max(secondaryDegreeSeq)+1)]

    #recurs a specified number of times, depending on the number of false edges
    for e in range(falseEdges):
        edgeAddProb = [0]
        degree = 1
        while degree < len(degreeDist):
            
            #calculates the probability of an edge being added to a secondary species of a certain degree
            edgeAddProb.append(degreeDist[degree]*(len(primarySpecies)-degree)/(len(primarySpecies)*len(secondarySpecies)-sum(secondaryDegreeSeq)))
            degree += 1
            
        #adjusts the degree distribution
        degreeDist = [degreeDist[0]] + [degreeDist[k]-edgeAddProb[k]+edgeAddProb[k-1] for k in range(1,len(degreeDist))] + [edgeAddProb[-1]]

    return degreeDist

def false_negatives(graph,primarySpecies,secondarySpecies,falseEdges):
    
    secondaryDegreeSeq = [graph.degree(secondary) for secondary in secondarySpecies]
    degreeDist = [secondaryDegreeSeq.count(k)/len(secondarySpecies) for k in range(max(secondaryDegreeSeq)+1)]

    #recurs a specified number of times, depending on the number of false edges
    for e in range(falseEdges):
        edgeRemoveProb = [0]
        degree = 1
        while degree < len(degreeDist):
            
            #calculates the probability of an edge being removed from a secondary species of a certain degree
            edgeRemoveProb.append(degreeDist[degree]*(degree/sum(secondaryDegreeSeq)))
            degree += 1
            
        #adjusts the degree distribution
        degreeDist = [degreeDist[k]-edgeRemoveProb[k]+edgeRemoveProb[k+1] for k in range(0,len(degreeDist)-1)] + [degreeDist[-1]-edgeRemoveProb[-1]]

    return degreeDist

#finds the ``positition'' of an ordering of indices
def position_finder(degrees,index):
    count = 1
    previousIndex = 0
    position = 0
    
    #operates recursively until the position of the specified combination is found
    for i in index[:-1]:
        for x in range(previousIndex,i):
            position += math.comb(len(degrees)-1-x,len(index)-count)
        previousIndex = i + 1
        count += 1
    
    position += index[-1] - previousIndex
    
    return position

#finds the upper and lower weight sums of a bracket. This function is specified in the Supplementary Materials as Algorithm 1
def upper_lower_combos(samples,prefix,indices,choices,threshold,counter,depth):
    
    #checks that there are still available choices
    if len(samples) >= prefix[-1]+1+choices-len(prefix):
        
        #calculates the upper and lower values of the bracket
        upperIndex = prefix + [r for r in range(prefix[-1]+1,prefix[-1]+1+choices-len(prefix))]
        lowerIndex = prefix + [r for r in range(len(indices)-choices+len(prefix),len(indices))]
        upper = [samples[p] for p in prefix] + [samples[r] for r in range(prefix[-1]+1,prefix[-1]+1+choices-len(prefix))]
        lower = [samples[p] for p in prefix] + [samples[r] for r in range(len(indices)-choices+len(prefix),len(indices))]
        
        #if the threshold is within the bracket, the counter is updated if maximum depth is reached, otherwise the prefix is extended
        if sum(upper) >= threshold and sum(lower) < threshold:
            if len(prefix) > (depth-1):
                counter += (position_finder(indices,lowerIndex) - position_finder(indices,upperIndex) + 1)*(sum(upper)-threshold)/(sum(upper)-sum(lower))
                return (prefix[:-1] + [prefix[-1]+1],counter)
            else:
                for d in indices[prefix[-1]+1:]:
                    return upper_lower_combos(samples,prefix + [d],indices,choices,threshold,counter,depth)
        
        #if the bracket is above the threshold, the counter is updated with the size of the bracket and the prefix is extended
        elif sum(upper) >= threshold and sum(lower) >= threshold:
            counter += position_finder(indices,lowerIndex) - position_finder(indices,upperIndex) + 1
            return (prefix[:-1] + [prefix[-1]+1],counter)
        
        #if the bracket is below the threshold, the prefix is shortened and updated
        else:
            if len(prefix) > 1:
                return (prefix[:-2] + [prefix[-2] + 1],counter)
            else:
                return ([],counter)
    
    #if there are no remaining free choices, the empty prefix and final count are returned
    else:
        return ([],counter)

#predicts survival probabilities when extinctions occur due to loss of interaction strength
def interaction_strength_species_extinctions(graph,primarySpecies,secondarySpecies,sensitivityRatio,depth):
    
    #determines the distribution of unique neighbours for secondary species
    degreeSequence = [len([n for n in graph.neighbors(secondary)]) for secondary in secondarySpecies]
    neighbourFrequency = Counter(degreeSequence)
    neighbourFrequency = {f:neighbourFrequency[f]/len(degreeSequence) for f in neighbourFrequency}
    degreeValues = list(set(degreeSequence))
    extinctByDegree = {d:[0 for x in range(d)] for d in degreeValues}

    #iterates over secondary species
    for secondary in secondarySpecies:

        #initialises various parameters, including the neighbours of the secondary species and the species' extinction sensitivity
        neighs = [n for n in graph.neighbors(secondary)]
        degrees = [graph.number_of_edges(secondary,n) for n in neighs]
        degrees = sorted(degrees,reverse=True)
        samples = [d for d in degrees]
        thresholdValue = math.ceil(sum(samples)*sensitivityRatio)
        indices = [r for r in range(len(samples))]
        localExtinct = []
        counter = 0
        prefix = [0]
        
        #iterates over number of neighbours chosen
        for choices in range(1,len(samples)+1):
            
            #runs the bracket checking algorithm recursively down to the specified depth
            while prefix != []:
                prefix,counter = upper_lower_combos(samples,prefix,indices,choices,thresholdValue,counter,depth)
                
            #records the proportion of combinations which result in extinction
            localExtinct.append(counter/math.comb(len(samples),choices))
            
            #resets parameters
            prefix = [0]
            counter = 0
        
        #records the probability of a secondary species of a given degree going extinct
        extinctByDegree[len(localExtinct)] = [d+l for d,l in zip(extinctByDegree[len(localExtinct)],localExtinct)]
        
    #normalises probability distributions
    for degree in degreeValues:
        norm = extinctByDegree[degree][-1]
        extinctByDegree[degree] = [n/norm for n in extinctByDegree[degree]]
            
    numberPrimary = len(primarySpecies)
    predictSurvival = [1-[d for v,d in graph.degree(secondarySpecies)].count(0)/len(secondarySpecies)]
    
    #iterates over removed primary species
    for removed in range(numberPrimary):
        extinctProb = 0
        
        #iterates over degree values
        for degree in degreeValues:
            for k in range(1,degree+1):
                if removed >= k:
                    
                    #updates the extinction probability, using the hypergeometric distribution weighted by 
                    #the probability of exceeding the threshold and the distribution of secondary species
                    extinctProb += extinctByDegree[degree][k-1]*neighbourFrequency[degree]*hypergeom.pmf(k,numberPrimary,degree,removed)
        
        #records the survival probability
        predictSurvival.append(1-extinctProb)

    return predictSurvival

#predicts extinctions when interactions are lost
def interaction_loss_extinctions(graph,primarySpecies,secondarySpecies,sensitivityRatio):
        
    #measures the degree sequence and degree distribution of the secondary species
    secondaryDegreeSeq = [graph.degree(secondary) for secondary in secondarySpecies]
    degreeDist = [secondaryDegreeSeq.count(k)/len(secondarySpecies) for k in range(max(secondaryDegreeSeq)+1)]
    interactions = [e for e in graph.edges()]
            
    predictSurvival = []
    
    #calculates survival probability as primary species are removed
    for removedInteract in range(len(interactions)+1):
        localSurviveProb = 0
        
        #iterates over the possible degree values for secondary species
        for degree in range(len(degreeDist)):
            
            #calculates the number of interactions that must be lost for extinction to occur
            thresholdValue = math.ceil(sensitivityRatio*degree)
            
            if degreeDist[degree] > 0:
                
                #extinction probability is calculated using the hypergeometric distirbution and weighted by the degree distribution
                localSurviveProb += degreeDist[degree]*(hypergeom.cdf(thresholdValue-1,len(interactions),degree,removedInteract))
                    
        #survival probability is recorded
        predictSurvival.append(localSurviveProb)
    
    return predictSurvival