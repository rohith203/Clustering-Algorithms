import numpy as np
import random
from collections import Counter
import matplotlib.pyplot as plt
from dissimilarity import *
import os
# print(os.listdir("./"))

#Reading the Amino acid data as a text file
amino = open('./Data/aminoacids.txt')

#Data- preprocessing:-

"""There is a new line character at the end of each line of the text file
    and -1 removes it"""
data = []
for line in amino:
    data.append(line[:-1])
#data[:10]

"""labels(dict) stores the name of the amino acid as key
    and the sequence as a value both as strings"""
labels = {}
for i in range(len(data)):
    if(data[i][0]=='>'):
        key = data[i][1:]
        val = ""
        j = i+1
        while(j<len(data) and data[j][0]!='>'):
            val = val + data[j]
            if j<len(data)-1: val = val[:-1]
            j+=1
        labels[key] = val

N = len(labels)
n = len(labels)
sim_mat = [[0 for i in range(n)] for j in range(n)]
keys = list(labels.keys())
key_dict = {}
for i in range(len(keys)):
    key_dict[keys[i]] = i

"""The following commented out code has been used to find the dissimlarity value between every sequence 
    with every other sequence and the values have been stored in a text file(similarity_mat)"""
# for i in range(n):
#     for j in range(i,n):
#         sim_mat[i][j] = get_dissimilarity(labels[keys[i]], labels[keys[j]])
#         sim_mat[j][i] = sim_mat[i][j]


'''
    Taking the dissimilarity matrix from a file which contains the values that have been calculated before.
'''
simfile = open("./Data/similarity_mat.txt")
simstring = simfile.readlines()
sim_mat = eval(simstring[0])


def assign_clusters(n, k, clusters, mediods):
    """Assigns cluster-centers to all the clusters"""
    for i in range(n):
        mincenter = 0
        min_dist = 99999
        for j in range(k):
            if min_dist>sim_mat[key_dict[mediods[j]]][i]:
                min_dist = sim_mat[key_dict[mediods[j]]][i]
                mincenter = j
        clusters[i] = mincenter
        
def get_total_cost(n, k, clusters, mediods):
    """Function that calculates the total cost of assigning a point to a cluster"""
    cost = [0]*k
    for i in range(k):
        for j in range(n):
            if clusters[j]==i:
                cost[i] += sim_mat[key_dict[mediods[i]]][j]
    return sum(cost)
    
def kmeans(clusters, k = 3):
    """Input:k the number of clusters to be formed and the clusters list.
        This function makes use of all the above functions and 
        forms the required no of clusters"""
    # initialization
    n = len(keys)
    random.seed(10)
    mediods = random.sample(keys, k)
    mediods_copy = mediods.copy()
#    clusters = [0 for i in range(n)]
    
    # assigning clusters to all data points
    assign_clusters(n, k, clusters, mediods)
    prev_cost = get_total_cost(n, k, clusters, mediods)

    # print(prev_cost)
    
    new_clusters = [0]*n
    for _ in range(100):
        for i in range(k):
            min_cost = prev_cost
            old_med = key_dict[mediods[i]]
            best_mediods = mediods.copy()
            for j in range(n):
                if clusters[j]==i and j!=old_med:
                    new_med = j
                    mediods[i] = keys[new_med]
                    assign_clusters(n, k, new_clusters, mediods)
                    cur_cost = get_total_cost(n, k, new_clusters, mediods)
                    if cur_cost < min_cost:
                        min_cost = cur_cost
                        best_mediods = mediods.copy()
                    mediods[i] = keys[old_med]

            mediods = best_mediods.copy()
            assign_clusters(n, k, clusters, mediods)
            prev_cost = get_total_cost(n, k, clusters, mediods)
        if(mediods_copy==mediods): break
        mediods_copy = mediods.copy()
    # print(clusters)
    
    
k = 1
clusters = [0 for i in range(N)]  # Each element is assigned a cluster number in this clusters.

'''
    Call to the kmeans function.
'''
kmeans(clusters, k)


cnt = Counter(clusters)
# print(len(cnt))


# print(clusters)
cluster_dict = {}
for i in range(k): cluster_dict[i] = []
for i in range(N):
    cluster_dict[clusters[i]].append(keys[i])


for ks,vs in cluster_dict.items():
    print(ks,":",vs)
    print()