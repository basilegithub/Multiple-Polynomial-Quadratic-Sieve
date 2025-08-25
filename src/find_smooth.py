# This is the file to find smooth numbers

import math, random
from utils import fermat_primality

# batch smooth test using a product tree
# more efficient than testing for smoothness every candidate : should be prefered
# Note that in this function, we test for smoothness abs(candidates), see QS.sieve_and_batch_smooth function
# Source: https://cr.yp.to/factorization/smoothparts-20040510.pdf (Bernstein)
def batch_smooth_test(candidates, prod_primes, cst_1, cst_2):
    # compute the smallest e such that 2**(2**e) >= max(candidates), ie e = int(log2(log2(max(candidates)))
    e = 0
    tmp = 2
    tmp_max = max(candidates)
    while tmp < tmp_max:
        tmp = tmp*tmp
        e += 1
        
    # build the product tree
    tmp = candidates.copy()
    prod_tree = [tmp]
    
    while len(prod_tree[-1]) > 1:
        line = [0]*(len(prod_tree[-1])>>1)

        for i in range(1,len(prod_tree[-1]),2):
            tmp = prod_tree[-1][i]*prod_tree[-1][i-1]

            if tmp <= prod_primes: line[i>>1] = tmp
            else: line[i>>1] = prod_primes+1

        if len(prod_tree[-1]) & 1: line.append(prod_tree[-1][-1])
        prod_tree.append(line)
    
    # update the product tree by computing the adequate remainders

    prod_tree[-1][0] = prod_primes%prod_tree[-1][0]
    for i in range(len(prod_tree)-1):

        for j in range(len(prod_tree[-i-1])-1):
            prod_tree[-i-2][(j<<1)] = prod_tree[-i-1][j]%prod_tree[-i-2][j<<1]
            prod_tree[-i-2][(j<<1)+1] = prod_tree[-i-1][j]%prod_tree[-i-2][(j<<1)+1]

        tmp = len(prod_tree[-i-1])-1
        prod_tree[-i-2][tmp<<1] = prod_tree[-i-1][tmp]%prod_tree[-i-2][tmp<<1]

        if (tmp<<1)+1 < len(prod_tree[-i-2]):
            prod_tree[-i-2][(tmp<<1)+1] = prod_tree[-i-1][tmp]%prod_tree[-i-2][(tmp<<1)+1]

    # test for smoothness
    smooth = [0]*len(candidates)
    for i in range(len(prod_tree[0])):
        tmp = 0
        j = prod_tree[0][i]
        while j and tmp < e:
            j = (j*j)%candidates[i]
            tmp += 1
        if not j: smooth[i] = [True]
        else:
            tmp = candidates[i]//math.gcd(candidates[i],j)
            if tmp < cst_1:
                smooth[i] = ["large", 1, tmp]
            elif tmp < cst_2 and fermat_primality(tmp):
                tmp = pollard_rho(tmp)
                smooth[i] = ["large", min(tmp), max(tmp)]
            else: smooth[i] = [False]
            
    return smooth

# Naive smoothness test: trial division by every prime in the base  
# Note that in this function, we test for smoothness (candidate) and not abs(candidate), see QS.sieve_and_smooth
def smooth_test(candidate, prime_base, const_1, const_2):
    if candidate < 0:
        candidate = -candidate
        
    for p in prime_base:
        while not candidate%p:
            candidate //= p
        if candidate==1: return [True]
        
    # candidate with one large prime
    if candidate < const_1:
        return ["large", 1, candidate]
        
    # candidate with two large primes
    elif candidate < const_2 and fermat_primality(candidate):
        tmp = pollard_rho(candidate)
        return ["large", min(tmp), max(tmp)]
        
    return [False]
    
    
# Pollard rho algorithm to factor small primes
# This is used to find the two large primes when needed
def pollard_rho(n):
    a = int(random.randint(1, n-3))
    s = int(random.randint(0, n-1))
    x = s
    y = s
    d = 1
    while d == 1:
        e = 1
        X = x
        Y = y
        for k in range(0, 100):
            x = (x*x+a)%n
            y = (y*y+a)%n
            y = (y*y+a)%n
            e = e*(x-y)%n
        d = math.gcd(e, n)
    if d == n:
        x = X
        y = Y
        d = 1
        while d == 1:
            x = (x*x+a)%n
            y = (y*y+a)%n
            y = (y*y+a)%n
            d = math.gcd(x-y, n)
    if d == n:
        return pollard_rho(n)
    return d, n//d