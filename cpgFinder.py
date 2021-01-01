import sys
import math

states = ['A','C','G','T','a','c','g','t']
emits = {'A':0,'C':1,'G':2,'T':3}

def viterbi(nucs,tMat, eMat, initProbs):
    nucCnt = len(nucs)
    stateCnt = len(states)
    trellis = [[0.0 for i in range(stateCnt)] for j in range(nucCnt)]
    backtrack = [[0.0 for i in range(stateCnt)] for j in range(nucCnt)]
    tMat = [[math.log(val) for val in row] for row in tMat]
    eMat = [[math.log(val) for val in row] for row in eMat]
    for i in range(stateCnt):
        trellis[0][i] = eMat[i][emits[nucs[0]]] + math.log(initProbs[i]) #weights for initial state

    #fill trellis
    for i in range(1,nucCnt):
        for j in range(stateCnt):
            heavy = -math.inf
            heavy_idx = 0
            for k in range(stateCnt):
                temp = trellis[i-1][k] + tMat[k][j] #accessing
                if temp > heavy:
                    heavy = temp
                    heavy_idx = k
            trellis[i][j] = heavy + eMat[j][emits[nucs[i]]] #setting
            backtrack[i][j] = heavy_idx

    #walk back through heaviest trellis
    idxPrediction = [0]*nucCnt
    heaviest = -math.inf
    heaviest_idx = 0
    for i in range(stateCnt):
        if trellis[nucCnt-1][i] > heaviest:
            heaviest = trellis[nucCnt - 1][i]
            heaviest_idx = i
    idxPrediction[nucCnt-1] = heaviest_idx
    for i in range(nucCnt-1,0,-1):
        idxPrediction[i-1] = backtrack[i][idxPrediction[i]]
    statePrediction = ''.join(states[i] for i in idxPrediction)
    return statePrediction


def main(argv):
    #not particularly important for percentages to sum correctly. Relative to weighting
    #only allowed to go +->- in G state (end of dinucleotide) and -->+ in c state (start of dinucleotide)
            #A         C       G        T        a        c         g        t
    tMat_default =[
         [0.18000, 0.27400, 0.42600, 0.12000, 1.0e-37, 1.0e-37, 1.0e-37, 1.0e-37], #A
         [0.17100, 0.36800, 0.27400, 0.18800, 1.0e-37, 1.0e-37, 1.0e-37, 1.0e-37], #C
         [0.16100, 0.33900, 0.37500, 0.12500, 0.00100, 0.00100, 0.00100, 0.00100], #G
         [0.07900, 0.35500, 0.38400, 0.18200, 1.0e-37, 1.0e-37, 1.0e-37, 1.0e-37], #T
         [1.0e-37, 0.00025, 1.0e-37, 1.0e-37, 0.30000, 0.20500, 0.28500, 0.21000], #a
         [1.0e-37, 0.00025, 1.0e-37, 1.0e-37, 0.32200, 0.29800, 0.07800, 0.30200], #c
         [1.0e-37, 0.00025, 1.0e-37, 1.0e-37, 0.24800, 0.24600, 0.29800, 0.20800], #g
         [1.0e-37, 0.00025, 1.0e-37, 1.0e-37, 0.17700, 0.23900, 0.29200, 0.29200]  #t
        ]
    #transition matrix values for 8-state look-back is very hard to find information on. Above are placeholders that perform relatively well but need adjusting.
    #difficult to perform traditional EM to generate as next emission probabilities are abstracted into transition matrix.
    #intra-state transition weights taken
    #from https://ocw.mit.edu/courses/electrical-engineering-and-computer-science/6-047-computational-biology-fall-2015/lectures_slides/MIT6_047F15_Lecture05.pdf
    #cross-state weights estimated from textbook information and could use some tuning


    #emission matrix is sort of silly for 8 state look-back, but keeping it for clarity
    #0s->lower end of 4 bytes (1e-37-ish) for log-space
    eMat_default = [
          [1.00000, 1.0e-37, 1.0e-37, 1.0e-37],
          [1.0e-37, 1.00000, 1.0e-37, 1.0e-37],
          [1.0e-37, 1.0e-37, 1.00000, 1.0e-37],
          [1.0e-37, 1.0e-37, 1.0e-37, 1.00000],
          [1.00000, 1.0e-37, 1.0e-37, 1.0e-37],
          [1.0e-37, 1.00000, 1.0e-37, 1.0e-37],
          [1.0e-37, 1.0e-37, 1.00000, 1.0e-37],
          [1.0e-37, 1.0e-37, 1.0e-37, 1.00000]
        ]

    #placeholders again. Perform well enough. Estimated from book.
    #                       A        C        G        T        a        c        g        t
    initProbs_default = [0.01500, 0.03000, 0.03000, 0.01500, 0.25000, 0.21000, 0.21000, 0.25000]

    try:
        with open(argv[0],'r') as f:
            gName = f.readline().strip()
            nucs = "".join(f.read().splitlines())
    except:
        print("usage: bad file or filepath")
        return
    statePrediction = viterbi(nucs,tMat_default,eMat_default,initProbs_default)

    #finds beginning and end of islands from state list
    print(gName)
    if statePrediction[0].isupper():
        island_start = 0
    else:
        island_start = None
    for i in range(1, len(statePrediction)):
        if statePrediction[i].isupper():
            if island_start == None:
                island_start = i
        elif island_start != None:
            print(island_start,"to",i)
            island_start = None

if __name__ == '__main__':
    main(sys.argv[1:])