# This is the file containing functions handling the polynomials generation

import random, math
from utils import invmod, CRT

# We use multi-polynomail QS, these functions generate efficiently the required polynomials


# sample the primes to build the a coefficient
def create_polynomial(n,prime_base,partial_solutions,bounds,target_log):
    best_bound, best_selected, best_sol,best_where, best_a = None, [], [], [], 1
    
    # we sample 7 candidates for the a coefficient
    for _ in range(7):
        where, bound =  [], target_log
        # sample primes of appropriate size so that the a coefficient has the appropriate size
        while bound > 1:
            if bound > 3:
                location = random.randint(bounds[2],bounds[3]-1)
            elif bound > 2:
                location = random.randint(bounds[1],bounds[2]-1)
            else:
                location = random.randint(bounds[0],bounds[1]-1)
            if location in where : continue # do not sample the same prime twice
            where.append(location)
            bound -=  math.log10(prime_base[location])
        if best_bound is None or abs(bound) < best_bound:
            best_bound = abs(bound)
            best_where = where
            
    for i in best_where :
        best_a *= prime_base[i]
        best_selected.append(prime_base[i])
        best_sol.append(partial_solutions[i])
        
    second_part = [best_a//i*invmod(best_a//i,i)%best_a for i in best_selected]
    return (best_a,best_sol,second_part,best_where,best_selected)
    
# compute relevant quantities for efficiently initialize sieving
def initialization(a,primes,residues,locations):
    inverse_a, way_to_root = [], []
    for i in range(len(primes)):
        if i in locations:
            inverse_a.append(0)
            way_to_root.append(0)
        else:
            tmp = invmod(a,primes[i])
            inverse_a.append(tmp)
            way_to_root.append(-2*residues[i]*tmp%primes[i])
    return inverse_a, way_to_root
    
# master function calling create_polynomial and initialization
# finds a suitable polynomial, compute its coefficients and return all the required parameters
def find_poly(n,primes,a,bounds,target):
    polynomials = create_polynomial(n,primes,a,bounds,target)
    poly_selected = [0,0,0,0]
    poly_selected[0],poly_selected[3] = polynomials[0], polynomials[3]
    param = initialization(poly_selected[0],primes,a,poly_selected[3])
    poly_index = 0
    threshold = (1<<(len(polynomials[3])-1))-1
    moduli = [polynomials[1][i] for i in range(len(poly_selected[3]))]

    tmp = CRT(moduli,polynomials[0],polynomials[2])
    if tmp >= poly_selected[0]//2: tmp = polynomials[0]-tmp
    poly_selected[1],poly_selected[2] = tmp,(tmp*tmp-n)//poly_selected[0]
    return poly_selected, param, moduli, polynomials[2], threshold
    
# finds the next polynomial in the list of polynomials defined by the coefficient a
def find_next_poly(poly_index,moduli,n,poly_selected,needed):
    # poly_index bits encode the positive/negative moduli
    tmp = poly_index^(poly_index+1) # encodes the moduli that change sign
    ind = 0
    while tmp:
        moduli[-1-ind] = -moduli[-1-ind]
        tmp >>= 1
        ind += 1
    
    tmp = CRT(moduli,poly_selected[0],needed)
    if tmp >= poly_selected[0]//2: tmp = poly_selected[0]-tmp
    poly_selected[1],poly_selected[2] = tmp,(tmp*tmp-n)//poly_selected[0] # compute the updated polynomial coefficients, coefficient a is not changed
    return poly_selected[1],poly_selected[2]