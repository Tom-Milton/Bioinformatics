import time
import sys

# cd \Users\Tom\Documents\GitHub\Bioinformatics-Coursework
# py Question2.py length500_A.txt length500_B.txt

def matrices(seq1, seq2):
    l1 = len(seq1) + 1
    l2 = len(seq2) + 1
    matrix_1 = [[0 for i in range(l2)] for i in range(l1)] #initialises score matrix
    matrix_2 = [[0 for i in range(l2)] for i in range(l1)] #initialises direction matrix
    maximum = 0
    maximum_index = ""
    for i in range(l1): #iterates through matrices and calculates cell data
        for j in range(l2):
            if i == 0 or j == 0: #if on first row or column, can only be 0 and end
                matrix_2[i][j] = "E"
            else:
                option1 = score(seq1[i-1], seq2[j-1]) + matrix_1[i-1][j-1] #diagonal option
                option2 = matrix_1[i-1][j] - 4 #up option
                option3 = matrix_1[i][j-1] - 4 #left option
                matrix_1[i][j] = max([option1, option2, option3, 0]) #chooses biggest such that it is non-negative
                m1 = matrix_1[i][j]
                if m1 > maximum: #checks against current highest score
                    maximum = m1 #stores score
                    maximum_index = [i, j] #stores index of score
                if m1 == 0:
                    matrix_2[i][j] = "E"
                elif m1 == option1:
                    matrix_2[i][j] = "D"
                elif m1 == option2:
                    matrix_2[i][j] = "U"
                elif m1 == option3:
                    matrix_2[i][j] = "L"
    best_score.append(maximum)
    i, j = maximum_index[0], maximum_index[1] #index of highest score
    seq1 = seq1[:i] #reduces sequences so they end at the index of the highest score
    seq2 = seq2[:j]
    return alignment(seq1, seq2, matrix_2)


def alignment(seq1, seq2, matrix_2):
    element = matrix_2[len(seq1)][len(seq2)]
    while element != "E": #backtracks until end is reached
        element = matrix_2[len(seq1)][len(seq2)] #current element letter
        if element == "D": #if diagonal adds last letter from both sequences
            best_alignment1.append(seq1[-1])
            best_alignment2.append(seq2[-1])
            seq1 = seq1[:-1]
            seq2 = seq2[:-1]
        elif element == "L": #if left adds last element from seq2 and "-" from seq1
            best_alignment1.append("-")
            best_alignment2.append(seq2[-1])
            seq2 = seq2[:-1]
        elif element == "U": #if left adds last element from seq1 and "-" from seq2
            best_alignment1.append(seq1[-1])
            best_alignment2.append("-")
            seq1 = seq1[:-1]


def score(a1, a2): #calculates score of alignment of two letters
    if a1 == a2 == "A":
        return 3
    elif a1 == a2 == "C":
        return 2
    elif a1 == a2 == "G":
        return 1
    elif a1 == a2 == "T":
        return 2
    elif a1 == "-" or a2 == "-":
        return -4
    else:
        return -3


def displayAlignment(alignment):
    string1 = alignment[0]
    string2 = alignment[1]
    string3 = ''
    for i in range(min(len(string1),len(string2))):
        if string1[i]==string2[i]:
            string3=string3+"|"
        else:
            string3=string3+" "
    print('Alignment ')
    print('String1: '+string1)
    print('         '+string3)
    print('String2: '+string2+'\n\n')


file1 = open(sys.argv[1], 'r')
seq1=file1.read()
file1.close()
file2 = open(sys.argv[2], 'r')
seq2=file2.read()
file2.close()
start = time.time()


best_alignment1 = []
best_alignment2 = []
best_score = []
matrices(seq1, seq2)
best_alignment1 = "".join(best_alignment1[::-1])
best_alignment2 = "".join(best_alignment2[::-1])
best_alignment = [best_alignment1, best_alignment2]


stop = time.time()
time_taken = stop - start


print('Time taken: ' + str(time_taken))
print('Best (score '+str(best_score)+'):')
displayAlignment(best_alignment)