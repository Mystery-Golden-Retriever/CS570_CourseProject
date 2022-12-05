# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import numpy as np
import time
import tracemalloc
import sys
from resource import *
import time
import psutil


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press ⌘F8 to toggle the breakpoint.

def conv(c):
    if c == 'A':
        return 0
    elif c == 'C':
        return 1
    elif c == 'G':
        return 2
    else:
        return 3


def dp(A, B, simMatrix, gapPenalty):
    # The Needleman-Wunsch algorithm

    # Stage 1: Create a zero matrix and fills it via algorithm
    n, m = len(A), len(B)
    mat = []
    for i in range(n + 1):
        mat.append([0] * (m + 1))
    for j in range(m + 1):
        mat[0][j] = gapPenalty * j
    for i in range(n + 1):
        mat[i][0] = gapPenalty * i
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            mat[i][j] = min(mat[i - 1][j - 1] + simMatrix[conv(A[i - 1])][conv(B[j - 1])],
                            mat[i][j - 1] + gapPenalty, mat[i - 1][j] + gapPenalty)

    # Stage 2: Computes the final alignment, by backtracking through matrix
    alignmentA = ""
    alignmentB = ""
    i, j = n, m
    while i and j:
        score, scoreDiag, scoreUp, scoreLeft = mat[i][j], mat[i - 1][j - 1], mat[i - 1][j], mat[i][j - 1]
        if score == scoreDiag + simMatrix[conv(A[i - 1])][conv(B[j - 1])]:
            alignmentA = A[i - 1] + alignmentA
            alignmentB = B[j - 1] + alignmentB
            i -= 1
            j -= 1
        elif score == scoreUp + gapPenalty:
            alignmentA = A[i - 1] + alignmentA
            alignmentB = '_' + alignmentB
            i -= 1
        elif score == scoreLeft + gapPenalty:
            alignmentA = '_' + alignmentA
            alignmentB = B[j - 1] + alignmentB
            j -= 1
    while i:
        alignmentA = A[i - 1] + alignmentA
        alignmentB = '_' + alignmentB
        i -= 1
    while j:
        alignmentA = '_' + alignmentA
        alignmentB = B[j - 1] + alignmentB
        j -= 1
    # Now return result in format: [1st alignment, 2nd alignment, similarity]
    return [alignmentA, alignmentB, mat[n][m]]
# Press the green button in the gutter to run the script.
def forwards(x, y, simMatrix, gapPenalty):
    # This is the forwards subroutine.
    n, m = len(x), len(y)
    mat = []
    for i in range(n+1):
        mat.append([0]*(m+1))
    for j in range(m+1):
        mat[0][j] = gapPenalty*j
    for i in range(1, n+1):
        mat[i][0] = mat[i-1][0] + gapPenalty
        for j in range(1, m+1):
            mat[i][j] = min(mat[i-1][j-1] + simMatrix[conv(x[i-1])][conv(y[j-1])],
                            mat[i-1][j] + gapPenalty,
                            mat[i][j-1] + gapPenalty)
        # Now clear row from memory.
        mat[i-1] = []
    return mat[n]
def backwards(x, y, simMatrix, gapPenalty):
    # This is the backwards subroutine.
    n, m = len(x), len(y)
    mat = []
    for i in range(n+1):
        mat.append([0]*(m+1))
    for j in range(m+1):
        mat[0][j] = gapPenalty*j
    for i in range(1, n+1):
        mat[i][0] = mat[i-1][0] + gapPenalty
        for j in range(1, m+1):
            mat[i][j] = min(mat[i-1][j-1] + simMatrix[conv(x[n-i])][conv(y[m-j])],
                            mat[i-1][j] + gapPenalty,
                            mat[i][j-1] + gapPenalty)
        # Now clear row from memory.
        mat[i-1] = []
    return mat[n]
def hirschberg(x, y, simMatrix, gapPenalty):
    # This is the main Hirschberg routine.
    n, m = len(x), len(y)
    if n<2 or m<2:
        # In this case we just use the N-W algorithm.
        return dp(x, y, simMatrix, gapPenalty)
    else:
        # Make partitions, call subroutines.
        F, B = forwards(x[:int(n/2)], y, simMatrix, gapPenalty), backwards(x[int(n / 2):], y, simMatrix, gapPenalty)
        partition = [F[j] + B[m-j] for j in range(m+1)]
        cut = partition.index(min(partition))
        # Clear all memory now, so that we don't store data during recursive calls.
        F, B, partition = [], [], []
        # Now make recursive calls.
        callLeft = hirschberg(x[:int(n/2)], y[:cut], simMatrix, gapPenalty)
        callRight = hirschberg(x[int(n/2):], y[cut:], simMatrix, gapPenalty)
        # Now return result in format: [1st alignment, 2nd alignment, similarity]
        return [callLeft[r] + callRight[r] for r in range(3)]

if __name__ == '__main__':
    # start = time.perf_counter()
    start_time = time.time()
    process = psutil.Process()
    # tracemalloc.clear_traces()
    # tracemalloc.start(25)
    f = open("./input111.txt")
    line = f.readline()
    base1 = line.strip('\n')
    line = f.readline()
    while line:
        x = line.strip('\n')
        if x.isdigit():
            base1 = base1[:int(x)+1] + base1 + base1[int(x)+1:]
        else:
            break
        line = f.readline()
    base2 = line.strip('\n')
    line = f.readline()
    while line:
        x = line.strip('\n')
        if x.isdigit():
            base2 = base2[:int(x)+1] + base2 + base2[int(x)+1:]
        else:
            break
        line = f.readline()
    f.close()

    # print(base1)
    # print(base2)
    skip_cost = 30
    mis_match_cost = np.array([[0, 110, 48, 94], [110, 0, 118, 48], [48, 118, 0, 110], [94, 48, 110, 0]])
    # m = len(base1)
    # n = len(base2)
    # dp = np.zeros((m+1, n+1), dtype=int)

    # 1:diagonal 2:up 3:left
    # alignment = np.zeros((m+1, n+1), dtype=int)

    # for q in range(1, m+1):
    #     dp[q][0] = q * skip_cost
    #     # alignment[q][0] = 2
    # for p in range(1, n+1):
    #     dp[0][p] = p * skip_cost
    #     # alignment[0][p] = 3
    # for i in range(1, m+1):
    #     for j in range(1, n+1):
    #         if base1[i-1] == base2[j-1]:
    #             dp[i][j] = dp[i-1][j-1]
    #         else:
    #             dp[i][j] = min(dp[i-1][j-1] + mis_match_cost[conv(base1[i-1])][conv(base2[j-1])],
    #                            dp[i-1][j] + skip_cost,
    #                            dp[i][j-1] + skip_cost
    #                            )
    # print(dp[m][n])
    # res = dp(base1, base2, mis_match_cost, skip_cost)
    # align1 = res[0]
    # align2 = res[1]
    # align3 = res[2]
    # print("Align A:",align1)
    # print("Align B:",align2)
    # print("cost:",align3)
    # #TEST
    # sum = 0
    # for i in range(0, len(align1)):
    #     if(align1[i] == '_' or align2[i] == '_'):
    #         sum += skip_cost
    #     elif(align1[i] != align2[i]):
    #         sum+=mis_match_cost[conv(align1[i])][conv(align2[i])]
    # print(sum)
    #print("First sequence:", base1)
    #print("Second sequence:", base2)
    #print("Calculating alignment distance by Hirschberg method...")
    z = hirschberg(base1, base2, mis_match_cost, skip_cost)
    #print("Alignment of A: ", z[0])
    #print("Alignment of B: ", z[1])
    #print("Alignment Cost: ", z[2], '\n')
    # end = time.perf_counter()
    end_time = time.time()
    memory_info = process.memory_info()
    # size, peak = tracemalloc.get_traced_memory()
    #print('memory: ', peak / 1024)
    #print("running time: %s seconds"%(end - start))
    #TEST
    # sum = 0
    # for i in range(0, len(z[0])):
    #     if(z[0][i] == '_' or z[1][i] == '_'):
    #         sum += skip_cost
    #     elif(z[0][i] != z[1][i]):
    #         sum+=mis_match_cost[conv(z[0][i])][conv(z[1][i])]
    # print(sum)
    g = open("./alignments.txt", 'w')
    #g.write(str(z[2]) + "\n")
    g.write(z[0][:50] + " " + z[0][-50:] + "\n")
    g.write(z[1][:50] + " " + z[1][-50:] + "\n")
    # g.write(str(end - start) + "\n")
    time_taken = (end_time - start_time) * 1000
    g.write(str(time_taken) + "\n")
    memory_consumed = int(memory_info.rss / 1024)
    g.write(str(memory_consumed) + '\n')
    # g.write(str(peak/1024)+'\n')
    g.write("\n")

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
