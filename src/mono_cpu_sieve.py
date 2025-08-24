# This file contains the functions that realize the sieve step in the case where only one cpu is used

from datetime import datetime
import sys
import math
import log
import sieve
import find_smooth
from sieve_polynomials_functions import *
from relations import *
from utils import format_duration

def sieve_and_batch_smooth(relations, smooth_number, full_found, partial_relations, partial_found, possible_smooth, graph, a, b, poly_selected, coeff, logs, primes, param, sieve_len, skipped, prime_start, prod_primes, const_1, const_2, tmp1, tmp2, tmp3, block, to_batch, size_partials, parent, cycle_len, n):
    
    tmp_block = sieve.sieve(b,poly_selected,logs,primes,a,n,param,sieve_len,skipped,prime_start,tmp1,tmp2,tmp3)
    to_batch += [abs(poly_selected[0]*i**2+tmp2*i+poly_selected[2]) for i in tmp_block]
    block += tmp_block
    coeff.append([poly_selected[0],poly_selected[1],len(tmp_block)])
    
    if len(block) > 512:
        smooth = find_smooth.batch_smooth_test(to_batch,prod_primes,const_1,const_2)
        tmp_index = 0
        for z in range(len(coeff)):
            coeff1 = coeff[z][0]
            coeff2 = coeff[z][1]
            for i in range(coeff[z][2]):
                tmp_smooth = smooth[tmp_index+i]
                value = coeff1*block[tmp_index+i]+coeff2
                
                relations, smooth_number, partial_relations, possible_smooth, full_found, partial_found, graph, size_partials, parent = handle_possible_smooth(value,tmp_smooth,full_found,partial_found,relations,smooth_number,partial_relations,possible_smooth,graph,size_partials,parent,cycle_len,n)
                
            tmp_index += coeff[z][2]
            
        block = []
        coeff = []
        to_batch = []
        
    return relations, smooth_number, full_found, partial_relations, partial_found, possible_smooth, block, coeff, to_batch, graph, size_partials, parent

   
def sieve_and_smooth(relations, smooth_number, full_found, partial_relations, partial_found, possible_smooth, graph, a, b, poly_selected, logs, primes, param, sieve_len, skipped, prime_start, const_1, const_2, tmp1, tmp2, tmp3, size_partials, parent, cycle_len, n):
    smooth = sieve.sieve(b,poly_selected,logs,primes,a,n,param,sieve_len,skipped,prime_start,tmp1,tmp2,tmp3)
    
    coeff1 = poly_selected[0]
    coeff2 = poly_selected[1]
    coeff3 = poly_selected[2]
  
    for i in smooth:
        tmp_smooth = find_smooth.smooth_test(coeff1*i**2+tmp2*i+coeff3,primes,const_1,const_2)
        value = coeff1*i+coeff2
        
        relations, smooth_number, partial_relations, possible_smooth, full_found, partial_found, graph, size_partials, parent = handle_possible_smooth(value,tmp_smooth,full_found,partial_found,relations,smooth_number,partial_relations,possible_smooth,graph,size_partials,parent,cycle_len,n)
        
    return relations, smooth_number, full_found, partial_relations, partial_found, possible_smooth, graph, size_partials, parent

def find_relations(primes, const, prod_primes, bounds, target, logs, a, b, flag_use_batch_smooth_test, n, LOG_PATH):
    log.write_log(LOG_PATH, "sieving...")
    log.write_log(LOG_PATH, "need to find at least "+str(len(primes)+10)+" relations")
    
    poly_selected, threshold, poly_index = [0,0,0,0], 0, 0
    relations, smooth_number, partial_relations, possible_smooth, graph, cycle_len = [], [], {}, {}, {}, [0]*10
    size_partials = 0
    parent = {}
    
    partial_found, full_found, skipped, last, sieve_len = 0, 0, 0, primes[-1], (b<<1)+1
    const_1 = const*primes[-1]
    const_2 = const*primes[-1]**2
    skipped += int(math.log2(const_2))
    prime_start = 30
    skipped = sieve.compute_skipped(skipped, logs, primes, prime_start)
    block = []
    coeff = []
    to_batch = []
    
    time_1 = datetime.now()
        
    sys.stdout.write('\r'+"0/("+str(len(primes)+1)+"+10) relations found")
    while len(relations) <= len(primes)+10:
        if poly_index >= threshold:
            poly_selected, param, moduli, needed, threshold = find_poly(n,primes,a,bounds,target)
            tmp1, tmp3 = poly_selected[0]<<1, poly_selected[0]*b**2
            poly_index = 0
        else:
            poly_selected[1],poly_selected[2] = find_next_poly(poly_index,moduli,n,poly_selected,needed)
            poly_index += 1
        tmp2 = poly_selected[1]<<1
        
        if flag_use_batch_smooth_test:
            relations, smooth_number, full_found, partial_relations, partial_found, possible_smooth, block, coeff, to_batch, graph, size_partials, parent = sieve_and_batch_smooth(relations,smooth_number,full_found,partial_relations,partial_found,possible_smooth,graph,a,b,poly_selected,coeff,logs,primes,param,sieve_len,skipped,prime_start,prod_primes,const_1,const_2,tmp1,tmp2,tmp3,block,to_batch,size_partials,parent,cycle_len,n)
 
        else:
            relations, smooth_number, full_found, partial_relations, partial_found, possible_smooth, graph, size_partials, parent = sieve_and_smooth(relations,smooth_number,full_found,partial_relations,partial_found,possible_smooth,graph,a,b,poly_selected,logs,primes,param,sieve_len,skipped,prime_start,const_1,const_2,tmp1,tmp2,tmp3,size_partials,parent,cycle_len,n)
                
        sys.stdout.write('\r'+str(len(smooth_number))+"/("+str(len(primes)+1)+"+10) relations found : full = "+str(full_found)+" ; partial = "+str(partial_found)+ " ("+str(size_partials)+")")
    print("\n")
    
    time_2 = datetime.now()
    
    log.write_log(LOG_PATH, "sieving done in "+format_duration(time_2-time_1)+".\n")
    log.write_log(LOG_PATH, str(len(smooth_number)) + " relations found")
    log.write_log(LOG_PATH, str(full_found)+" full relations found")
    log.write_log(LOG_PATH, str(partial_found)+" partial relations found")
    log.write_log(LOG_PATH, str(len(possible_smooth))+" smooth with large prime(s) found")
    log.write_log(LOG_PATH, "Distribution of cycle length:")
    for i in range(9):
        log.write_log(LOG_PATH, str(i+2)+"-cycle: "+str(cycle_len[i]))
    log.write_log(LOG_PATH, "11+-cycle: "+str(cycle_len[-1])+"\n")
                
    return relations, smooth_number 