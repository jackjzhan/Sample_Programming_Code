#!/usr/bin/env python
#Computational Genomics
#Solution to Global Alignment with Scoring Matrix and Affine Gap
#http://rosalind.info/problems/gaff/
#Written by Jack Zhan

from numpy import zeros

def BLOSUM62(xc,yc):
	map = {
			"A": {"A":-4, "C":0, "D":2, "E":1, "F":2, "G":0, "H":2, "I":1, "K":1, "L":1, "M":1, "N":2, "P":1, "Q":1, "R":1, "S":-1, "T":0, "V":0, "W":3, "Y":2}, 
			"C": {"A":0, "C":-9, "D":3, "E":4, "F":2, "G":3, "H":3, "I":1, "K":3, "L":1, "M":1, "N":3, "P":3, "Q":3, "R":3, "S":1, "T":1, "V":1, "W":2, "Y":2},
			"D": {"A":2, "C":3, "D":-6, "E":-2, "F":3, "G":1, "H":1, "I":3, "K":1, "L":4, "M":3, "N":-1, "P":1, "Q":0, "R":2, "S":0, "T":1, "V":3, "W":4, "Y":3},
			"E": {"A":1, "C":4, "D":-2, "E":-5, "F":3, "G":2, "H":0, "I":3, "K":-1, "L":3, "M":2, "N":0, "P":1, "Q":-2, "R":0, "S":0, "T":1, "V":2, "W":3, "Y":2},
			"F": {"A":2, "C":2, "D":3, "E":3, "F":-6, "G":3, "H":1, "I":0, "K":3, "L":0, "M":0, "N":3, "P":4, "Q":3, "R":3, "S":2, "T":2, "V":1, "W":-1, "Y":-3},
			"G": {"A":0, "C":3, "D":1, "E":2, "F":3, "G":-6, "H":2, "I":4, "K":2, "L":4, "M":3, "N":0, "P":2, "Q":2, "R":2, "S":0, "T":2, "V":3, "W":2, "Y":3},
			"H": {"A":2, "C":3, "D":1, "E":0, "F":1, "G":2, "H":-8, "I":3, "K":1, "L":3, "M":2, "N":-1, "P":2, "Q":0, "R":0, "S":1, "T":2, "V":3, "W":2, "Y":-2},
			"I": {"A":1, "C":1, "D":3, "E":3, "F":0, "G":4, "H":3, "I":-4, "K":3, "L":-2, "M":-1, "N":3, "P":3, "Q":3, "R":3, "S":2, "T":1, "V":-3, "W":3, "Y":1},
			"K": {"A":1, "C":3, "D":1, "E":-1, "F":3, "G":2, "H":1, "I":3, "K":-5, "L":2, "M":1, "N":0, "P":1, "Q":-1, "R":-2, "S":0, "T":1, "V":2, "W":3, "Y":2},
			"L": {"A":1, "C":1, "D":4, "E":3, "F":0, "G":4, "H":3, "I":-2, "K":2, "L":-4, "M":-2, "N":3, "P":3, "Q":2, "R":2, "S":2, "T":1, "V":-1, "W":2, "Y":1},
			"M": {"A":1, "C":1, "D":3, "E":2, "F":0, "G":3, "H":2, "I":-1, "K":1, "L":-2, "M":-5, "N":2, "P":2, "Q":0, "R":1, "S":1, "T":1, "V":-1, "W":1, "Y":1},
			"N": {"A":2, "C":3, "D":-1, "E":0, "F":3, "G":0, "H":-1, "I":3, "K":0, "L":3, "M":2, "N":-6, "P":2, "Q":0, "R":0, "S":-1, "T":0, "V":3, "W":4, "Y":2},
			"P": {"A":1, "C":3, "D":1, "E":1, "F":4, "G":2, "H":2, "I":3, "K":1, "L":3, "M":2, "N":2, "P":-7, "Q":1, "R":2, "S":1, "T":1, "V":2, "W":4, "Y":3},
			"Q": {"A":1, "C":3, "D":0, "E":-2, "F":3, "G":2, "H":0, "I":3, "K":-1, "L":2, "M":0, "N":0, "P":1, "Q":-5, "R":-1, "S":0, "T":1, "V":2, "W":2, "Y":1},
			"R": {"A":1, "C":3, "D":2, "E":0, "F":3, "G":2, "H":0, "I":3, "K":-2, "L":2, "M":1, "N":0, "P":2, "Q":-1, "R":-5, "S":1, "T":1, "V":3, "W":3, "Y":2},
			"S": {"A":-1, "C":1, "D":0, "E":0, "F":2, "G":0, "H":1, "I":2, "K":0, "L":2, "M":1, "N":-1, "P":1, "Q":0, "R":1, "S":-4, "T":-1, "V":2, "W":3, "Y":2},
			"T": {"A":0, "C":1, "D":1, "E":1, "F":2, "G":2, "H":2, "I":1, "K":1, "L":1, "M":1, "N":0, "P":1, "Q":1, "R":1, "S":-1, "T":-5, "V":0, "W":2, "Y":2},
			"V": {"A":0, "C":1, "D":3, "E":2, "F":1, "G":3, "H":3, "I":-3, "K":2, "L":-1, "M":-1, "N":3, "P":2, "Q":2, "R":3, "S":2, "T":0, "V":-4, "W":3, "Y":1},
			"W": {"A":3, "C":2, "D":4, "E":3, "F":-1, "G":2, "H":2, "I":3, "K":3, "L":2, "M":1, "N":4, "P":4, "Q":2, "R":3, "S":3, "T":2, "V":3, "W":-11, "Y":-2},
			"Y": {"A":2, "C":2, "D":3, "E":2, "F":-3, "G":3, "H":-2, "I":1, "K":2, "L":1, "M":1, "N":2, "P":3, "Q":1, "R":2, "S":2, "T":2, "V":1, "W":-2, "Y":-7},
			}
	return(map[xc][yc])

def bBLOSUM62(xc,yc):
    map = {
            "A": {"A":4, "C":0, "D":-2, "E":-1, "F":-2, "G":0, "H":-2, "I":-1, "K":-1, "L":-1, "M":-1, "N":-2, "P":-1, "Q":-1, "R":-1, "S":1, "T":0, "V":0, "W":-3, "Y":-2}, 
            "C": {"A":0, "C":9, "D":-3, "E":-4, "F":-2, "G":-3, "H":-3, "I":-1, "K":-3, "L":-1, "M":-1, "N":-3, "P":-3, "Q":-3, "R":-3, "S":-1, "T":-1, "V":-1, "W":-2, "Y":-2},
            "D": {"A":-2, "C":-3, "D":6, "E":2, "F":-3, "G":-1, "H":-1, "I":-3, "K":-1, "L":-4, "M":-3, "N":1, "P":-1, "Q":0, "R":-2, "S":0, "T":-1, "V":-3, "W":-4, "Y":-3},
            "E": {"A":-1, "C":-4, "D":2, "E":5, "F":-3, "G":-2, "H":0, "I":-3, "K":1, "L":-3, "M":-2, "N":0, "P":-1, "Q":2, "R":0, "S":0, "T":-1, "V":-2, "W":-3, "Y":-2},
            "F": {"A":-2, "C":-2, "D":-3, "E":-3, "F":6, "G":-3, "H":-1, "I":0, "K":-3, "L":0, "M":0, "N":-3, "P":-4, "Q":-3, "R":-3, "S":-2, "T":-2, "V":-1, "W":1, "Y":3},
            "G": {"A":0, "C":-3, "D":-1, "E":-2, "F":-3, "G":6, "H":-2, "I":-4, "K":-2, "L":-4, "M":-3, "N":0, "P":-2, "Q":-2, "R":-2, "S":0, "T":-2, "V":-3, "W":-2, "Y":-3},
            "H": {"A":-2, "C":-3, "D":-1, "E":0, "F":-1, "G":-2, "H":8, "I":-3, "K":-1, "L":-3, "M":-2, "N":1, "P":-2, "Q":0, "R":0, "S":-1, "T":-2, "V":-3, "W":-2, "Y":2},
            "I": {"A":-1, "C":-1, "D":-3, "E":-3, "F":0, "G":-4, "H":-3, "I":4, "K":-3, "L":2, "M":1, "N":-3, "P":-3, "Q":-3, "R":-3, "S":-2, "T":-1, "V":3, "W":-3, "Y":-1},
            "K": {"A":-1, "C":-3, "D":-1, "E":1, "F":-3, "G":-2, "H":-1, "I":-3, "K":5, "L":-2, "M":-1, "N":0, "P":-1, "Q":1, "R":2, "S":0, "T":-1, "V":-2, "W":-3, "Y":-2},
            "L": {"A":-1, "C":-1, "D":-4, "E":-3, "F":0, "G":-4, "H":-3, "I":2, "K":-2, "L":4, "M":2, "N":-3, "P":-3, "Q":-2, "R":-2, "S":-2, "T":-1, "V":1, "W":-2, "Y":-1},
            "M": {"A":-1, "C":-1, "D":-3, "E":-2, "F":0, "G":-3, "H":-2, "I":1, "K":-1, "L":2, "M":5, "N":-2, "P":-2, "Q":0, "R":-1, "S":-1, "T":-1, "V":1, "W":-1, "Y":-1},
            "N": {"A":-2, "C":-3, "D":1, "E":0, "F":-3, "G":0, "H":1, "I":-3, "K":0, "L":-3, "M":-2, "N":6, "P":-2, "Q":0, "R":0, "S":1, "T":0, "V":-3, "W":-4, "Y":-2},
            "P": {"A":-1, "C":-3, "D":-1, "E":-1, "F":-4, "G":-2, "H":-2, "I":-3, "K":-1, "L":-3, "M":-2, "N":-2, "P":7, "Q":-1, "R":-2, "S":-1, "T":-1, "V":-2, "W":-4, "Y":-3},
            "Q": {"A":-1, "C":-3, "D":0, "E":2, "F":-3, "G":-2, "H":0, "I":-3, "K":1, "L":-2, "M":0, "N":0, "P":-1, "Q":5, "R":1, "S":0, "T":-1, "V":-2, "W":-2, "Y":-1},
            "R": {"A":-1, "C":-3, "D":-2, "E":0, "F":-3, "G":-2, "H":0, "I":-3, "K":2, "L":-2, "M":-1, "N":0, "P":-2, "Q":1, "R":5, "S":-1, "T":-1, "V":-3, "W":-3, "Y":-2},
            "S": {"A":1, "C":-1, "D":0, "E":0, "F":-2, "G":0, "H":-1, "I":-2, "K":0, "L":-2, "M":-1, "N":1, "P":-1, "Q":0, "R":-1, "S":4, "T":1, "V":-2, "W":-3, "Y":-2},
            "T": {"A":0, "C":-1, "D":-1, "E":-1, "F":-2, "G":-2, "H":-2, "I":-1, "K":-1, "L":-1, "M":-1, "N":0, "P":-1, "Q":-1, "R":-1, "S":1, "T":5, "V":0, "W":-2, "Y":-2},
            "V": {"A":0, "C":-1, "D":-3, "E":-2, "F":-1, "G":-3, "H":-3, "I":3, "K":-2, "L":1, "M":1, "N":-3, "P":-2, "Q":-2, "R":-3, "S":-2, "T":0, "V":4, "W":-3, "Y":-1},
            "W": {"A":-3, "C":-2, "D":-4, "E":-3, "F":1, "G":-2, "H":-2, "I":-3, "K":-3, "L":-2, "M":-1, "N":-4, "P":-4, "Q":-2, "R":-3, "S":-3, "T":-2, "V":-3, "W":11, "Y":2},
            "Y": {"A":-2, "C":-2, "D":-3, "E":-2, "F":3, "G":-3, "H":2, "I":-1, "K":-2, "L":-1, "M":-1, "N":-2, "P":-3, "Q":-1, "R":-2, "S":-2, "T":-2, "V":-1, "W":2, "Y":7},
            }
    return(map[xc][yc])

def align(x, y):

# Initialize the matrices.
	S_lower = zeros((len(x)+1, len(y)+1), dtype=int)
	S_middle = zeros((len(x)+1, len(y)+1), dtype=int)
	S_upper = zeros((len(x)+1, len(y)+1), dtype=int)
	backtrack = zeros((len(x)+1, len(y)+1), dtype=int)

	for i in range(1,len(x)+1):
		S_lower[i][0] = 11 +(i-1)
		S_upper[i][0] = 100
		S_middle[i][0] = 100
	for j in range(1,len(y)+1):
		S_lower[0][j] = 100
		S_upper[0][j] = 11 +(j-1)
		S_middle[0][j] = 100
	print(S_lower)
    
	# Fill in the Score and Backtrack matrices.
	for i in range(1, len(x)+1):
		for j in range(1, len(y)+1):
			S_lower[i][j] = min([S_lower[i-1][j] + 1, S_middle[i-1][j] + 11, S_upper[i-1][j] + 11])
			S_upper[i][j] = min([S_upper[i][j-1] + 1, S_middle[i][j-1] + 11, S_lower[i][j-1] + 11])           
			middle_scores = [S_lower[i][j], S_middle[i-1][j-1] + BLOSUM62(x[i-1], y[j-1]), S_upper[i][j]]
			S_middle[i][j] = min(middle_scores)
			backtrack[i][j] = middle_scores.index(S_middle[i][j])
	print(S_lower)
# Initialize the aligned strings as the input strings up to the position of the high score.
	x_aligned= ""
	y_aligned = ""

	# Backtrack to start of the local alignment starting at the highest scoring cell.
	# Note: the solution format specifically asks for substrings, so no indel insertion necessary.
	score = 0
	while i*j != 0:
		if backtrack[i][j] == 0:
			i -= 1
			if y_aligned != "" and y_aligned[0] == "-":
				score += 1
			else:
				score += 11
			x_aligned = x[i] + x_aligned
			y_aligned = "-" + y_aligned
		elif backtrack[i][j] == 1:
			i -= 1
			j -= 1
			score += BLOSUM62(x[i],y[j])
			x_aligned = x[i] + x_aligned
			y_aligned = y[j] + y_aligned
		elif backtrack[i][j] == 2:
			j -= 1
			if x_aligned != "" and x_aligned[0] == "-":
				score += 1
			else:
				score += 11
			y_aligned = y[j] + y_aligned
			x_aligned = "-" + x_aligned   
	score = -score
	print(score)
	print(x_aligned)
	print(y_aligned)
	return str(score), x_aligned, y_aligned

alignment = align("IIAIAQRNDPAKVQIEWRTSSHERLPFNYKDNWKKRGCCPQMKDAHQLKVETPPQINHSVSQYQTKGCMAQEFWNV","IIARAQFNDDLQSDHFYHAKMAQQIEWRTSSHERLPFHMPSFYYKYWNMEELPHQLDCPQINHDISPSGYQTKGCMAKGYWNV")


def global_alignment_affine_gap_penalty(v, w):
    '''Returns the global alignment score of v and w with constant gap peantaly 11 subject to the scoring_matrix.'''
    # Initialize the matrices.
    S = [[[0 for j in range(len(w)+1)] for i in range(len(v)+1)] for k in range(3)]
    backtrack = [[[0 for j in range(len(w)+1)] for i in range(len(v)+1)] for k in range(3)]

    # Initialize the edges with the given penalties.
    for i in range(1, len(v)+1):
        S[0][i][0] = -11 - (i-1)*1
        S[1][i][0] = -11 - (i-1)*1
        S[2][i][0] = -10*11
    for j in range(1, len(w)+1):
        S[2][0][j] = -11 - (j-1)*1
        S[1][0][j] = -11 - (j-1)*1
        S[0][0][j] = -10*11

    # Fill in the scores for the lower, middle, upper, and backtrack matrices.
    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            lower_scores = [S[0][i-1][j] - 1, S[1][i-1][j] - 11]
            S[0][i][j] = max(lower_scores)
            backtrack[0][i][j] = lower_scores.index(S[0][i][j])

            upper_scores = [S[2][i][j-1] - 1, S[1][i][j-1] - 11]
            S[2][i][j] = max(upper_scores)
            backtrack[2][i][j] = upper_scores.index(S[2][i][j])

            middle_scores = [S[0][i][j], S[1][i-1][j-1] + bBLOSUM62(v[i-1], w[j-1]), S[2][i][j]]
            S[1][i][j] = max(middle_scores)
            backtrack[1][i][j] = middle_scores.index(S[1][i][j])

   # Initialize the values of i, j and the aligned sequences.
    i,j = len(v), len(w)
    v_aligned, w_aligned = v, w

    # Get the maximum score, and the corresponding backtrack starting position.
    matrix_scores = [S[0][i][j], S[1][i][j], S[2][i][j]]
    max_score = max(matrix_scores)
    backtrack_matrix = matrix_scores.index(max_score)

    # Quick lambda function to insert indels.
    insert_indel = lambda word, i: word[:i] + '-' + word[i:]

    # Backtrack to the edge of the matrix starting bottom right.
    while i*j != 0:
        if backtrack_matrix == 0:  # Lower backtrack matrix conditions.
            if backtrack[0][i][j] == 1:
                backtrack_matrix = 1
            i -= 1
            w_aligned = insert_indel(w_aligned, j)

        elif backtrack_matrix == 1:  # Middle backtrack matrix conditions.
            if backtrack[1][i][j] == 0:
                backtrack_matrix = 0
            elif backtrack[1][i][j] == 2:
                backtrack_matrix = 2
            else:
                i -= 1
                j -= 1

        else:  # Upper backtrack matrix conditions.
            if backtrack[2][i][j] == 1:
                backtrack_matrix = 1
            j -= 1
            v_aligned = insert_indel(v_aligned, i)

    # Prepend the necessary preceeding indels to get to (0,0).
    for _ in range(i):
        w_aligned = insert_indel(w_aligned, 0)
    for _ in range(j):
        v_aligned = insert_indel(v_aligned, 0)

    return str(max_score), v_aligned, w_aligned


global_alignment_affine_gap_penalty("QSSNHMAIRSMSAEFIWSLSTFMEHGQSHHDGQLSNVAQWHMSQKKLAAILSAGDLAWVAIFNMAYAQPRQKSVARE","QSSIHMPIRSMSAEAPQKYNIWSLSTFMEHAWMYMQIYHCGDHHDGQLSNVAMSQKGLAILSAGDLVLVAYAQPRQKSVARE")

from numpy import zeros

def Cost(xc, yc):
    #Cost function assigning 0 to match, 2 to transition, 4 to transversion, and 8 to a gap
    if xc == yc: return -5 # match
    if xc == '-' or yc == '-': return 4 # gap
    minc, maxc = min(xc, yc), max(xc, yc)
    if minc == 'A' and maxc == 'G': return 2 # transition
    elif minc == 'C' and maxc == 'T': return 1 # transition
    return 3 # transversion

def globalAlignment(x, y, s):
    #Calculate global alignment value of sequences x and y using dynamic programming.  Return global alignment value.
    D = zeros((len(x)+1, len(y)+1), dtype=int)
    for j in range(1, len(y)+1):
        D[0, j] = D[0, j-1] + s('-', y[j-1])
    for i in range(1, len(x)+1):
        D[i, 0] = D[i-1, 0] + s(x[i-1], '-')
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            D[i, j] = min(D[i-1, j-1] + s(x[i-1], y[j-1]), # diagonal
                          D[i-1, j  ] + s(x[i-1], '-'),    # vertical
                          D[i  , j-1] + s('-',    y[j-1])) # horizontal
    return D, D[len(x), len(y)]


D, globalAlignmentValue = globalAlignment('TACAG', 'TATAG', Cost)
globalAlignmentValue