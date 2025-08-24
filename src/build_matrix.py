# This is the file to build the binary matrix

from utils import lowest_set_bit

# Create the sparse matrix:
# Consider the dense binary matrix A[i][j]
# The sparse matrix M is defined by M[i] contains j if and only if A[i][j] = 1
# For large numbers to factor, this reduces that memory requirement for the linear algebra part
# Each prime in the factor base corresponds to one line
# Each relation corresponds to one column
def build_sparse_matrix(relations,primes):
    bin_matrix = []
    line = []
    for i in range(len(relations)):
        if relations[i] < 0: line.append(i)
    bin_matrix.append(line)
    for i in range(len(primes)):
        line = []
        for j in range(len(relations)):
            tmp = 0
            tmp2 = primes[i]
            while not relations[j]%tmp2:
                tmp ^= 1
                tmp2 *= primes[i]
            if tmp: line.append(j)
        bin_matrix.append(line)
    while bin_matrix[-1] == []: del bin_matrix[-1]
    return bin_matrix
    
# Reduce the sparse matrix, according to various rules:
# 1. If a line contains only one non-zero value, then it is deleted as well as the column containing the corresponding non-zero value
# 2. If a line contrains no non-zero value, then it is deleted
# 3. We only keep 10 more columns than lines, to ensure we still have solutions while reducing the matrix size
def reduce_sparse_matrix(matrix,relations,smooth):
    flag = True
    while flag:
        flag = False
        i = 0
        while i < len(matrix):
            if len(matrix[i]) == 1:
                flag = True
                coeff = matrix[i][0]
                for j in range(len(matrix)):
                    if j != i:
                        stored = -1
                        for z in range(len(matrix[j])):
                            if matrix[j][z] == coeff: stored = z
                            elif matrix[j][z] > coeff: matrix[j][z] -= 1
                        if stored > -1: del matrix[j][stored]
                del relations[coeff]
                del smooth[coeff]
                del matrix[i]
            else: i += 1
            
        i = 0
        while i < len(matrix):
            if matrix[i] == []:
                del matrix[i]
                flag = True
            else: i += 1
        length = len(matrix)+10
        for i in range(len(matrix)):
            if matrix[i][0] >= length:
                flag = True
                matrix[i] = []
            else:
                for j in range(len(matrix[i])):
                    if matrix[i][-j-1] < length:
                        if j > 0:
                            flag = True
                            matrix[i] = matrix[i][:len(matrix[i])-j]
                        break
                        
# Build the dense binary matrix by computing explicitely every A[i][j]
def build_dense_matrix(relations,primes):
    bin_matrix = []
    for rel in relations:
        line = [0]*(len(primes)+1)
        
        line[0] = 1*(rel < 0)
            
        for i in range(len(primes)):
            tmp = primes[i]
            while not rel%tmp:
                tmp *= primes[i]
                line[i+1] ^= 1
                
        bin_matrix.append(line)
        
    return bin_matrix
    
def siqs_build_matrix_opt(M):
    """Convert the given matrix M of 0s and 1s into a list of numbers m
    that correspond to the columns of the matrix.
    The j-th number encodes the j-th column of matrix M in binary:
    The i-th bit of m[i] is equal to M[i][j].
    """
    m = len(M[0])
    cols_binary = [""] * m
    for mi in M:
        for j, mij in enumerate(mi):
            cols_binary[j] += "1" if mij else "0"
    return [int(cols_bin[::-1], 2) for cols_bin in cols_binary], len(M), m