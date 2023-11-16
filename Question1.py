import time
import sys

# cd \Users\Tom\Documents\GitHub\Bioinformatics-Coursework
# py Question1.py length9_A.txt length9_B.txt

def alignments(seq1, seq2, a1, a2, s):
    global best_score
    global best_alignment
    global num_alignments
    if len(a1) != 0:
        s += score(a1[-1], a2[-1]) #calculates score during alignment generation
    l1 = len(seq1)
    l2 = len(seq2)
    if l1 == 0 and l2 == 0: #checks for base case, if both sequences are empty
        num_alignments += 1
        if s > best_score: #updates best score and alignment if a new best score is found
            best_score = s
            best_alignment = [a1, a2]
        return
    else:
        if l1 == 0: #if seq1 is empty, can only receive "-" from seq1
            option2 = seq2[0]
            all = alignments(seq1, seq2[1:], a1 + "-", a2 + option2, s)
        elif l2 == 0: #if seq2 is empty, can only receive "-" from seq2
            option1 = seq1[0]
            all = alignments(seq1[1:], seq2, a1 + option1, a2 + "-", s)
        else:
            option1 = seq1[0]
            option2 = seq2[0]
            all = alignments(seq1[1:], seq2[1:], a1 + option1, a2 + option2, s) #receive first letter from both sequences
            all = alignments(seq1[1:], seq2, a1 + option1, a2 + "-", s) #receive first letter from seq1 and "-" from seq2
            all = alignments(seq1, seq2[1:], a1 + "-", a2 + option2, s) #receive first letter from seq2 and "-" from seq1
        return all


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


num_alignments = 0
best_score = -999999
best_alignment = ""
alignments(seq1, seq2, "", "", s=0)


stop = time.time()
time_taken=stop-start


print('Alignments generated: '+str(num_alignments))
print('Time taken: '+str(time_taken))
print('Best (score '+str(best_score)+'):')
displayAlignment(best_alignment)