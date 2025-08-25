# This is the file to build the binary matrix

from utils import lowest_set_bit

# Create the sparse matrix:
# Consider the dense binary matrix A[i][j]
# The sparse matrix M is defined by M[i] contains j if and only if A[i][j] = 1
# For large numbers to factor, this reduces that memory requirement for the linear algebra part
# Each prime in the factor base corresponds to one line
# Each relation corresponds to one column
def build_sparse_matrix(relations, primes):
    bin_matrix = []
    line = []
    for i in range(len(relations)):
        if relations[i] < 0: line.append(i)
    bin_matrix.append(set(line))
    for i in range(len(primes)):
        line = []
        for j in range(len(relations)):
            tmp = 0
            tmp2 = primes[i]
            while not relations[j]%tmp2:
                tmp ^= 1
                tmp2 *= primes[i]
            if tmp: line.append(j)
        bin_matrix.append(set(line))
    while bin_matrix[-1] == []: del bin_matrix[-1]
    return bin_matrix
    
# Reduce the sparse matrix, according to the following rules:
# 1. Delete empty rows, i.e. primes that do not impact square formation
# 2. Delete singletons, i.e. row where only one relations has non zero value
# 3. If there it is possible to delete some columns, delete columns in rows with as few elements as possible
#   This maximises the chance to remove a row and thus reduce the computation cost of finding kernel vectors
def reduce_sparse_matrix(matrix, relations, smooth):
    flag = True
    while flag:
        flag = False

        singleton_queue = [index for index, row in enumerate(matrix) if len(row) == 1]
        active_cols = set(range(len(relations)))

        flag = len(singleton_queue)

        while singleton_queue:
            index = singleton_queue.pop()

            if len(matrix[index]) != 1: continue
            coeff = next(iter(matrix[index]))

            for j, row in enumerate(matrix):
                if j != index and coeff in row:
                        row.remove(coeff)
                        if len(row) == 1: singleton_queue.append(j)

            matrix[index].clear()
            active_cols.discard(coeff)

        matrix = [row for row in matrix if row]

        relations = [relations[index] for index in sorted(active_cols)]
        smooth = [smooth[index] for index in sorted(active_cols)]

        mapping = {old: new for new, old in enumerate(sorted(active_cols))}
        matrix = [{mapping[c] for c in row if c in mapping} for row in matrix]

        active_cols = set(range(len(relations)))
        target_n_cols = len(matrix)+10
        to_delete = len(relations) - target_n_cols
        current_len_row = 2

        while to_delete >= current_len_row-1:

            index = -1
            for i in range(len(matrix)):
                if len(matrix[i]) == current_len_row:
                    flag = True
                    index = i
                    break

            if index != -1:
                for coeff in matrix[index]:
                    active_cols.discard(coeff)

                    for j, row in enumerate(matrix):
                        if j != index and coeff in row:
                            row.discard(coeff)

                to_delete -= current_len_row-1

                matrix[index].clear()

            else:
                current_len_row += 1

        matrix = [row for row in matrix if row]

        relations = [relations[index] for index in sorted(active_cols)]
        smooth = [smooth[index] for index in sorted(active_cols)]

        mapping = {old: new for new, old in enumerate(sorted(active_cols))}
        matrix = [{mapping[c] for c in row if c in mapping} for row in matrix]

    return matrix, relations, smooth
                        
# Build the dense binary matrix by computing explicitely every A[i][j]
def build_dense_matrix(relations, primes):
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