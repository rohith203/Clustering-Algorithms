import numpy as np
from scipy.cluster.hierarchy import dendrogram,fcluster
import random
import matplotlib.pyplot as plt
from collections import Counter
from dissimilarity import *

import os

# Reading the Amino acid data as a text file 
amino = open('./Data/aminoacids.txt')

# Data- preprocessing:-

"""There is a new line character at the end of each line of the text file
    and -1 removes it"""
data = []
for line in amino:
    data.append(line[:-1])

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
# # calculate similarity
# for i in range(n):
#     for j in range(i,n):
#         sim_mat[i][j] = get_similarity(labels[keys[i]], labels[keys[j]])
#         sim_mat[j][i] = sim_mat[i][j]


simfile = open("./Data/similarity_mat.txt")
simstring = simfile.readlines()
sim_mat = eval(simstring[0])
# print(len(sim_mat))


sim_mat = np.array(sim_mat)
N = len(keys)
linkage = 'single'


def find_max_cluster(cluster_dict):
    '''
        This function finds the largest cluster and returns it.
    '''
    maxkey = frozenset([])
    maxd = 0
    for x in cluster_dict:
        if len(x) > maxd:
            maxkey = x
            maxd = len(x)
    return maxkey


def get_min_dist(a, b):
    '''
        for Minimum/Single Linkage.
        returns the minimum distance among all the pairs
        in the clusters a and b.
    '''
    mind = 999999
    for x in list(a):
        for y in list(b):
            dist = sim_mat[key_dict[x]][key_dict[y]]
            if mind>dist:
                mind = dist
    return mind

def get_max_dist(a, b):
    '''
        for Maximum/Complete Linkage.
        returns the maximum distance among all the pairs
        of the clusters a and b.
    '''
    maxd = 0
    for x in a:
        for y in b:
            dist = sim_mat[key_dict[x]][key_dict[y]]
            if maxd<dist:
                maxd = dist
    return maxd

def get_ave_dist(a,b):
    '''
        for Average Linkage.
        returns the average distance between all the pairs
        of data points of the clusters a and b.
    '''
    sumd = 0
    for x in a:
        for y in b:
            dist = sim_mat[key_dict[x]][key_dict[y]]
            sumd += dist
    return sumd/(len(a)*len(b))


def divide_cluster(cluster):
    '''
        This function takes a cluster as input,
        divides it into two and returns them.
    '''

    # if cluster has only two elements return them.
    if len(cluster)==2:
        # print(cluster)
        return [cluster[0]],[cluster[1]]
    

    splinter_list = []  #  a list that holds the elements of new cluster
    clust_list = list(cluster)
    
    # initializing the splinter cluster with the max distant object in the first cluster
    max_aved = 0  # maximum average distance
    jump_obj = -1  # a variable that holds the object to be transferred to splinter_cluster 
    for i in range(len(clust_list)):
        aved = 0
        for j in range(len(clust_list)):
            if i!=j:
                aved += sim_mat[key_dict[clust_list[i]]][key_dict[clust_list[j]]]
        if max_aved<aved:
            max_aved = aved
            jump_obj = i
    if jump_obj != -1:
        splinter_list.append(clust_list[jump_obj])
        clust_list.remove(clust_list[jump_obj])
    
    # transferring objects of clust_list with more dissimilarity to splinter
    larger_clust = clust_list
    smaller_clust = splinter_list
    if len(splinter_list)>len(clust_list): 
        larger_clust = splinter_list
        smaller_clust = clust_list
    
    maxdiff = 0
    
    while jump_obj!=-1 and len(larger_clust)>1 and len(smaller_clust)>0:
        jump_obj = -1
        for i in range(len(larger_clust)):
            aved1 = 0
            aved2 = 0
            for j in range(len(larger_clust)):
                if i!=j: aved1 += sim_mat[key_dict[larger_clust[i]]][key_dict[larger_clust[j]]]
            aved1 = aved1/(len(larger_clust)-1)
            for j in range(len(smaller_clust)):
                aved2 += sim_mat[key_dict[larger_clust[i]]][key_dict[smaller_clust[j]]]
            aved2 = aved2/(len(smaller_clust))
            diff = aved1-aved2
            if diff>0 and diff>maxdiff:
                maxdiff = diff
                jump_obj = i

        #  the object with maximum positive difference is moved from the larger cluster to the smaller cluster
        if jump_obj!=-1:
            smaller_clust.append(larger_clust[jump_obj])
            larger_clust.remove(larger_clust[jump_obj])

    return larger_clust,smaller_clust



iteration = 1
def Divisive(clust):
    '''
        a function that initially takes the cluster with
        all objects and recursively divides it into clusters
        until each object is in its own cluster.
    '''

    if len(clust)<2: return
    global iteration
    x,y = divide_cluster(clust)
    
    if len(x)>0: 
        cluster_dict[frozenset(x)] = iteration
    if len(y)>0: 
        cluster_dict[frozenset(y)] = iteration
    iteration+=1

    if len(y)>=1: Divisive(x)
    if len(x)>=1: Divisive(y)


cluster_dict = {}  # A dictionary the stores all clusters as keys. Clusters to be merged have same values.
cluster_dict[frozenset(keys)] = 0
clust = find_max_cluster(cluster_dict) # the cluster with all the data objects.

'''
    call to the divisive function.
'''
Divisive(clust) 

link_mat = np.zeros((N-1,4))   # matrix that is required by dendrogram function as input to plot dendrogram.


cluster_indices = {}   # A dictionary the stores all clusters as keys and node indices as values.
for ks in keys:
    cluster_indices[frozenset([ks])] = key_dict[ks]


'''
    The code below calculates the indices for clusters other than the singleton clusters.
'''
ite = 0
a = -1
b = -1
items = list(cluster_dict.items())
for j in range(len(items)-1,1,-2):
    a = items[j][0]
    b = items[j-1][0]
    a = frozenset(a)
    b = frozenset(b)
    cluster_indices[a.union(b)] = N+ite
    ite+=1

# print(cluster_dict)


'''
    This code creates the linkage matrix which is required to plot dendrogram.
'''
ite = 0
a = -1
b = -1
items = list(cluster_dict.items())
for j in range(len(items)-1,1,-2):
    a = items[j][0]
    b = items[j-1][0]
    # print(a,b)
    a = frozenset(a)
    b = frozenset(b)        
    link_mat[ite][0] = cluster_indices[a]
    link_mat[ite][1] = cluster_indices[b]
    if linkage=='single': link_mat[ite][2] = get_min_dist(a,b)
    elif linkage=='complete': link_mat[ite][2] = get_max_dist(a,b)
    elif linkage=='average': link_mat[ite][2] = get_ave_dist(a,b)
    link_mat[ite][3] = len(a.union(b))
    ite+=1


'''
    Plotting dendrogram using scipy.cluster.hierarchy.dendrogram which
    linkage matrix constructed earlier as input.
'''
fig = plt.figure(figsize=(20,10))
plt.ylabel('Dissimilarity', fontsize=15)
plt.xlabel('Amino Acids',fontsize=15)
plt.title("Dendrogram of Divisive Clustering", fontsize=20)

dendro = dendrogram(link_mat, labels=keys)

plt.show()


''' This code prints the clusters obtained at a particular value of k.'''
k = 150
assignments = fcluster(link_mat, k,'maxclust')

cnt = Counter(assignments)

clusters_final = {}
for i in range(len(cnt)): clusters_final[i+1] = []
for i in range(N):
    clusters_final[assignments[i]].append(keys[i])

for ks,vs in clusters_final.items():
    print(ks,":",vs)
    print()