import os
import sys
import time
import psutil

COST_TABLE = {
    'delta' : 30,
    'alpha' : {
        'A' : {'A' : 0, 'C' : 110, 'G' : 48, 'T' : 94},
        'C' : {'A' : 110, 'C' : 0, 'G' : 118, 'T' : 48},
        'G' : {'A' : 48, 'C' : 118, 'G' : 0, 'T' : 110},
        'T' : {'A' : 94, 'C' : 48, 'G' : 110, 'T' : 0},
    },
}

total_cost = 0
cost = 0

def seq_generator(s_init, gen_list):
    '''
    Take a root seq and generate the long seq
    '''
    new_seq = s_init
    for idx in gen_list:
        new_seq = new_seq[:idx+1] + new_seq + new_seq[idx+1:]
    return new_seq


def input_parser(input_file):
    '''
    Parse input file to a sequence of string
    '''
    X,Y = '',''
    X_init, Y_init = '',''
    X_gen_list, Y_gen_list = [],[]
    with open(input_file) as fp:
        lines = fp.read().splitlines()
        cnt = 0
        for line in lines:
            if not line.isdigit():
                if cnt == 0:
                    X_init = line
                elif cnt == 1:
                    Y_init = line
                cnt += 1
            else:
                if cnt == 1:
                    X_gen_list.append(int(line))
                elif cnt == 2:
                    Y_gen_list.append(int(line))
    X = seq_generator(X_init, X_gen_list)
    Y = seq_generator(Y_init, Y_gen_list)
    return X,Y


def basicAlign(str_1, str_2):
    N = len(str_1)
    M = len(str_2)
    A = [[0 for i in range(M+1)] for i in range(N+1)]

    for i in range(M+1):
        A[0][i] = COST_TABLE['delta'] * i
    
    for j in range(N+1):
        A[j][0] = COST_TABLE['delta'] * j

    for i in range(1, N+1):
        for j in range(1, M+1):
            A[i][j] = min(
                    A[i-1][j-1] + COST_TABLE['alpha'][str_1[i-1]][str_2[j-1]],
                    A[i-1][j] + COST_TABLE['delta'], 
                    A[i][j-1] + COST_TABLE['delta']
                )

    res_1 = ''
    res_2 = ''
    j = M
    i = N
    while(i >= 1 and j >= 1):
        x_i = str_1[i-1]
        y_j = str_2[j-1]
        if A[i][j] == A[i-1][j-1] + COST_TABLE['alpha'][x_i][y_j]:
            res_1 = x_i + res_1
            res_2 = y_j + res_2
            i -= 1
            j -= 1
        elif A[i][j] == A[i][j-1] + COST_TABLE['delta']:
            res_1 = '_' + res_1
            res_2 = y_j + res_2
            j -= 1
        else:
            res_1 = x_i + res_1
            res_2 = '_' + res_2
            i -= 1
    while(i >= 1):
        res_1 = str_1[i - 1] + res_1
        res_2 = "_" + res_2
        i -= 1
    while(j >= 1):
        res_1 = "_" + res_1
        res_2 = str_2[j - 1] + res_2
        j -= 1
    global cost
    cost = A[N][M]
    return [res_1, res_2]


def dcAlign(str_1, str_2):
    M = len(str_1)
    N = len(str_2)

    if (M <= 2 or N <= 2):
        tmp = basicAlign(str_1, str_2)
        global total_cost
        global cost
        total_cost += cost
        return tmp
    
    f = [[0 for i in range(2)] for j in range(M+1)]
    g = [[0 for i in range(2)] for j in range(M+1)]
    mid = int(N / 2)

    for i in range(M+1):
        f[i][0] = COST_TABLE['delta'] * i
    
    for j in range(1, mid+1):
        f[0][1] = COST_TABLE['delta'] * j
        for i in range(1, M+1):
            f[i][1] = min(
                COST_TABLE['alpha'][str_1[i-1]][str_2[j-1]] + f[i-1][0],
                COST_TABLE['delta'] + f[i-1][1],
                COST_TABLE['delta'] + f[i][0]
            )
        for i in range(M+1):
            f[i][0] = f[i][1]
    
    str_1 = str_1[::-1]
    str_2 = str_2[::-1]
    for i in range(M+1):
        g[i][0] = COST_TABLE['delta'] * i

    for j in range(mid+1, N+1):
        g[0][1] = COST_TABLE['delta'] * (j - mid)
        for i in range(1, M+1):
            x_i = str_1[i-1]
            y_j = str_2[j - (mid + 1)]
            g[i][1] = min(
                    COST_TABLE['alpha'][x_i][y_j] + g[i-1][0],
                    COST_TABLE['delta'] + g[i-1][1], 
                    COST_TABLE['delta'] + g[i][0]
                )
        for i in range(M+1):
            g[i][0] = g[i][1]

    str_1 = str_1[::-1]
    str_2 = str_2[::-1]

    q = -1
    sum = sys.maxsize
    for i in range(M+1):
        if f[i][1] + g[M-i][1] < sum:
            q = i
            sum = f[i][1] + g[M-i][1]
    left = dcAlign(str_1[0: q], str_2[0: mid])
    right = dcAlign(str_1[q: M], str_2[mid: N])
    return [left[0]+right[0], left[1]+right[1]]
  
            
def main(input_file):
    input_X, input_Y = input_parser(input_file)
    res = dcAlign(input_X, input_Y)
    return res[0], res[1]


if __name__ == "__main__":
    # parse input augments
    if len(sys.argv) != 3:
        raise Exception("Must have 3 Augments!")

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # process info. of memory
    process = psutil.Process()
    
    # main function
    tik = time.time()
    alignment_1, alignment_2 = main(input_file)
    tok = time.time()

    # memory info.
    memory_info = process.memory_info()
    memory_consumed = int(memory_info.rss/1024)

    with open(output_file, 'w') as fp:
        fp.write(f"The optimal Cost is: {total_cost}\n")
        fp.write(f"The optimal alignment for X is: {alignment_1}\n")
        fp.write(f"The optimal alignment for Y is: {alignment_2}\n")
        fp.write(f"Total time cost: {(tok-tik) * 1000}ms\n")
        fp.write(f"Total memory cost: {memory_consumed}KB\n")