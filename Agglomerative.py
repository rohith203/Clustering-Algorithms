import numpy as np
import random
from scipy.cluster.hierarchy import dendrogram,fcluster
import matplotlib.pyplot as plt
from collections import Counter
from dissimilarity import *
import matplotlib as mpl

import os
# print(os.listdir("./"))

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

sim_mat = [[0 for i in range(N)] for j in range(N)]
keys = list(labels.keys())
key_dict = {}
for i in range(len(keys)):
    key_dict[keys[i]] = i
    
"""The following commented out code has been used to find the dissimlarity value between every sequence 
    with every other sequence and the values have been stored in a text file(similarity_mat)"""

# # calculate dissimilarity
# for i in range(N):
#     for j in range(i,N):
#         sim_mat[i][j] = get_dissimilarity(labels[keys[i]], labels[keys[j]])
#         sim_mat[j][i] = sim_mat[i][j]


simfile = open("./Data/similarity_mat.txt")
simstring = simfile.readlines()
sim_mat = eval(simstring[0])



def get_min(s_mat):
    """Input: A dissimilarity matrix
        Returns a tuple containing indices of the minimum dissimilar object and its dissimilarity"""
    m = 999999
    m_i = 0
    m_j = 0
    n = len(s_mat)
    for i in range(0,n):
        for j in range(i,n):
            if i!=j and m>s_mat[i][j]:
                m = s_mat[i][j]
                m_i = i
                m_j = j
    return (m_i, m_j, m)

def get_min_dist(a, i, centroids):
    """Input: two clusters
        Returns minimum distance between two clusters"""
    mind = 999999
    for x in centroids[a]:
        for y in centroids[i]:
            dist = sim_mat[key_dict[x]][key_dict[y]]
            if mind>dist:
                mind = dist
    return mind

def get_max_dist(a, i, centroids):
    """Input: two clusters
        Returns maximum distance between two clusters"""
    maxd = 0
    for x in centroids[a]:
        for y in centroids[i]:
            dist = sim_mat[key_dict[x]][key_dict[y]]
            if maxd<dist:
                maxd = dist
    return maxd

def get_ave_dist(a, i, centroids):
    """Input:two clusters
        Returns the average distance between two clusters"""
    sumd = 0
    for x in centroids[a]:
        for y in centroids[i]:
            dist = sim_mat[key_dict[x]][key_dict[y]]
            sumd += dist
    return sumd/(len(centroids[a])*len(centroids[i]))


def merge(a,b,mat,centroids,i,linkage):
    '''
        This function takes two clusters as input,
         combines them and calculates dissimilarity with other
         objects according to the linkage type.
    '''
    for x in centroids[b]:
        centroids[a].append(x)
    del centroids[b]
    mat = np.delete(mat, b, 1)
    mat = np.delete(mat, b, 0)
    for i in range(len(centroids)):
        if a!=i:
            if linkage=='single': mat[a][i] = get_min_dist(a,i,centroids)
            elif linkage=='complete': mat[a][i] = get_max_dist(a,i,centroids)
            elif linkage=='average': mat[a][i] = get_ave_dist(a,i,centroids)
            mat[i][a] = mat[a][i]
    return mat


mat = np.array(sim_mat.copy())
centroids = keys.copy()
cluster_dict = {}
for k,v in key_dict.items():
    cluster_dict[frozenset([k])] = v

for i in range(len(keys)): centroids[i] = [centroids[i]]
    


N = len(keys)
link_mat =  np.zeros((N-1,4))



def Agglomerate(mat, centroids, link_mat, cluster_dict,linkage='single'):
    '''
        This function uses the merge function to combine clusters
        using bottom-up strategy until a single cluster with all data
        is formed.
        This function also constructs the linkage matrix used as input
        to plot the dendrogram.
    '''
    for i in range(N-1):
        if mat.shape == (1,1):break
        a,b,dist = get_min(mat)
        link_mat[i][0] = cluster_dict[frozenset(centroids[a])]
        link_mat[i][1] = cluster_dict[frozenset(centroids[b])]
        link_mat[i][2] = dist
        mat = merge(a,b,mat,centroids,i,linkage)
        cluster_dict[frozenset(centroids[a])] = N + i
        link_mat[i][3] = len(centroids[a])


linkage_type = 'complete'
Agglomerate(mat, centroids, link_mat, cluster_dict,linkage_type)


"""The following code prints the dendrogram of the linkage matrix"""

fig = plt.figure(figsize=(20,10))
plt.ylabel('Dissimilarity', fontsize=15)
plt.xlabel('Amino Acids',fontsize=15)
plt.title("Dendrogram of Agglomerative Clustering - "+linkage_type, fontsize=20)
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