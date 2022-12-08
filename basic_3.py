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

def init_DPTable(DP_table, len_x, len_y):
    for i in range(len_y+1):
        DP_table[0][i] = i * COST_TABLE['delta']
    for i in range(len_x+1):
        DP_table[i][0] = i * COST_TABLE['delta']
    return

def comp_DPTable(DP_table, len_x, len_y, X, Y):
    for i in range(1,len_x+1):
        for j in range(1,len_y+1):
            DP_table[i][j] = min(
                COST_TABLE['alpha'][X[i-1]][Y[j-1]] + DP_table[i-1][j-1],
                COST_TABLE['delta'] + DP_table[i-1][j],
                COST_TABLE['delta'] + DP_table[i][j-1]
            )
    return DP_table[-1][-1]

def retrieve_alignment(DP_table, len_x, len_y, X, Y):
    align_X_r, align_Y_r = '',''
    curr_OPT = DP_table[-1][-1]
    curr_idx_X, curr_idx_Y = len_x, len_y
    
    while (curr_idx_X != 0 or curr_idx_Y != 0):
        # both x_i and y_j in the opt. alignment:
        if curr_OPT == COST_TABLE['alpha'][X[curr_idx_X-1]][Y[curr_idx_Y-1]]\
                        + DP_table[curr_idx_X-1][curr_idx_Y-1]:
            align_X_r += X[curr_idx_X-1]
            align_Y_r += Y[curr_idx_Y-1]
            curr_idx_X -= 1
            curr_idx_Y -= 1
            curr_OPT = DP_table[curr_idx_X][curr_idx_Y]
        # x_i mismatch, put x_i in the X alignment and _ in the Y alignment
        elif curr_OPT == COST_TABLE['delta'] + DP_table[curr_idx_X-1][curr_idx_Y]:
            align_X_r += X[curr_idx_X-1]
            align_Y_r += '_'
            curr_idx_X -= 1
            curr_OPT = DP_table[curr_idx_X][curr_idx_Y]
        # y_j mismatch, put y_j in the Y alignment and _ in the X alignment
        elif curr_OPT == COST_TABLE['delta'] + DP_table[curr_idx_X][curr_idx_Y-1]:
            align_Y_r += Y[curr_idx_Y-1]
            align_X_r += '_'
            curr_idx_Y -= 1
            curr_OPT = DP_table[curr_idx_X][curr_idx_Y]
    
    return align_X_r[::-1], align_Y_r[::-1]
            


def main(input_file, output_file):
    input_X, input_Y = input_parser(input_file)
    len_x = len(input_X)
    len_y = len(input_Y)

    # initialize DP table of size [len_x+1][len_y+1]
    # stored in the order of table[xi][yj]
    OPT_table = [[0 for _ in range(len_y+1)] for _ in range(len_x+1)]
    init_DPTable(OPT_table, len_x, len_y)

    # compute the DP Table and the optimal value
    OPT_value = comp_DPTable(OPT_table, len_x, len_y, input_X, input_Y)
    # print(f"The optimal Cost is: {OPT_value}\n")

    # retrieve the alignement from the DP table
    alignment_X, alignment_Y = retrieve_alignment(OPT_table, len_x, len_y, input_X, input_Y)
    # print(f"The optimal alignment for X is: {alignment_X}\n")
    # print(f"The optimal alignment for Y is: {alignment_Y}\n")

    return OPT_value, (alignment_X, alignment_Y)

if __name__ == "__main__":
    # # parse input augments
    if len(sys.argv) != 3:
        raise Exception("Must have 3 Augments!")

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # process info. of memory
    process = psutil.Process()

    # main function
    tik = time.time()
    OPT_value, alignments = main(input_file, output_file)
    tok = time.time()

    # memory info.
    memory_info = process.memory_info()
    memory_consumed = int(memory_info.rss/1024)

    with open(output_file, 'w') as fp:
        fp.write(f"The optimal Cost is: {OPT_value}\n")
        fp.write(f"The optimal alignment for X is: {alignments[0]}\n")
        fp.write(f"The optimal alignment for Y is: {alignments[1]}\n")
        fp.write(f"Total time cost: {(tok-tik) * 1000}ms\n")
        fp.write(f"Total memory cost: {memory_consumed}KB\n")