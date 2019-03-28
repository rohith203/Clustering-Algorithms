
gap = 2
substitution = 1
match = 0

def get_dissimilarity(sequenceA, sequenceB):
    opt = [[0 for i in range(len(sequenceB) + 1)] for j in range(len(sequenceA)+1)]
    for i in range(1,len(sequenceA) + 1):
        opt[i][0] = opt[i - 1][0] + gap
    for j in range(1,len(sequenceB) + 1):
        opt[0][j] = opt[0][j - 1] + gap

    for i in range(1,len(sequenceA) + 1):
        for j in range(1,len(sequenceB) + 1):
            scoreDiag = opt[i - 1][j - 1]
            if sequenceA[i-1] == sequenceB[j-1]: scoreDiag+=match 
            else: scoreDiag+=substitution
            scoreLeft = opt[i][j - 1] + gap
            scoreUp = opt[i - 1][j] + gap

            opt[i][j] = min(min(scoreDiag, scoreLeft), scoreUp)

    return opt[len(sequenceA)][len(sequenceB)]