# This is the file conatining functions to generate the solutions

# Used when the gaussian pivot is used, convert sparse vector of indices into dense binary vector
def convert_to_binary(z, smooth_number):
    res = [0]*len(smooth_number)
    for index in z: res[index] ^= 1
    return res
    
def convert_to_binary_lanczos(z, smooth_number):
    res = [0]*len(smooth_number)
    for i in range(len(smooth_number)):
        if (z >> len(smooth_number) - i - 1)&1: res[i] = 1
    return res
    
# Compute the x and y from the null space vector
def compute_solution(relation_set, smooth, z, n, primes):
    vector_y = [0]*(len(primes)+1)
    x = 1
    y = 1
    for e in range(len(z)):
        if z[e]:
            x *= smooth[e]
            x %= n
            if relation_set[e] < 0:
                vector_y[0] += 1
            for i in range(len(primes)):
                tmp2 = primes[i]
                while not relation_set[e]%tmp2:
                    vector_y[i+1] += 1
                    tmp2 *= primes[i]
    if (vector_y[0]>>1)&1:
        y = -1
    for k in range(1, len(vector_y)):
        y *= pow(primes[k-1], vector_y[k]>>1, n)
        y %= n
    return x,y