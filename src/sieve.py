# This is the file that executes the sieve algorithm

from utils import invmod

# we do not sieve with small primes. We need to compute the average contribution of these skipped primes
def compute_skipped(skipped, logs, primes, prime_start):
    for i in range(prime_start): skipped += 2*logs[i]/primes[i]
    skipped = round(skipped)
    
    return skipped

# Sieve and return candidates for rigorous smoothness test
def sieve(L,poly_used,logs,primes,a,n,param,sieve_len,skipped,prime_start,tmp1,tmp2,tmp3):
    # Sieve
    sieve = [0]*sieve_len
    for i in range (prime_start, len(primes)):
        p = primes[i]
        k = logs[i]
        if param[0][i]:
            z = ((-poly_used[1]+a[i])*param[0][i]+L)%p

            for j in range(z, sieve_len,p):
                sieve[j] += k
            z = (z+param[1][i])%p
        else:
            z = (-poly_used[2]*invmod(tmp2, p)+L)%p
        
        for j in range(z, sieve_len, p): sieve[j] += k
         
    # identify candidates
    poly_eval = tmp3-tmp2*L+poly_used[2]
    B_smooth = []
    for x in range(L):
        k = x-L
        log_eval = poly_eval.bit_length()-1-skipped
        if sieve[x] >= log_eval : B_smooth.append(k)
        if sieve[-1-x] >= log_eval : B_smooth.append(-k)
        poly_eval += tmp1*k+poly_used[0]+tmp2
        
    log_eval = poly_eval.bit_length()-1-skipped
    if sieve[L] >= log_eval : B_smooth.append(0)
    
    return B_smooth