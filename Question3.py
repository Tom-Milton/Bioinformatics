import numpy as np
from tabulate import tabulate

# cd \Users\Tom\Documents\GitHub\Bioinformatics-Coursework
# python
# from Question3 import NJ
# NJ("Matrix1")

def NJ(file_name): #opens and stores txt file as arrays
    with open(file_name+".txt") as f:
        lines = [l.rstrip() for l in f] #removes blank rows
        matrix = np.array([list(l.split()) for l in lines if l]) #removes blank spaces and stores in multidimensional array
        heading = matrix[:1, 1:].astype(np.str).flatten().tolist() #extracts letters and stores in array
        matrix = matrix[1:, 1:].astype(np.float) #extracts numbers and stores in array
    f.close()
    return qMatrix(heading, matrix)


def qMatrix(heading, matrix): #calculates qscore of numbers matrix
    sums = [sum(i) for i in matrix] #calculates row sums
    q_matrix = matrix.copy()
    r = np.size(matrix, 0) #calculates number of rows
    for x, v in np.ndenumerate(matrix): #calculates qscores except down leading diagonal
        i, j = x[0], x[1]
        if i == j:
            continue
        else:
            q_matrix[i, j] = (r-1)*v-sums[i]-sums[j]
    printing(heading, matrix, q_matrix, sums)
    if r == 2: #checks for recursive base case
        return
    else:
        return group(heading, matrix, q_matrix)


def group(heading, matrix, q_matrix):
    i, j = np.unravel_index(q_matrix.argmin(), q_matrix.shape) #finds minimum in numbers matrix
    for x, v in np.ndenumerate(matrix[i]): #calculates new scores for new combined species except down leading diagonal
        x = x[0]
        if x == i or x == j:
            continue
        else:
            matrix[i, x] = (matrix[x, i] + matrix[x, j] - matrix[i, j])/ 2
    matrix[..., i] = matrix[i] #copies new scores from column to row
    matrix = np.delete(np.delete(matrix, [j], 0), [j], 1) #removes extra column
    new_heading = heading.copy()
    new_heading[i] = heading[i] + heading[j] #combines letters
    del new_heading[j] #deletes extra column
    return qMatrix(new_heading, matrix)


def printing(heading, matrix, q_matrix, sums):
    heading.append("Row Sums")
    print("distance matrix:\n" + tabulate(np.column_stack((matrix, sums)), heading, tablefmt="fancy_grid"))
    print("Q matrix:\n" + tabulate(q_matrix, heading, tablefmt="fancy_grid"))
