# This file contains the functions that realize the sieve step when multiprocessing is on (NB_CPU > 1)

from datetime import datetime
import time
import sys
import math
import log
import sieve
import find_smooth
from sieve_polynomials_functions import *
from relations import *
from utils import format_duration
import multiprocessing
    
def siever_batch(polynomials, tmp_rels, b, n, logs, primes, a, sieve_len, skipped, prime_start, const_1, const_2, prod_primes):
    to_batch, values, coeff, block = [], [], [], []
    while True:

        try:
            poly_selected, param, moduli, needed, threshold = polynomials.get(timeout=0)
            poly_index = 0
            tmp1, tmp2, tmp3 = poly_selected[0]<<1, poly_selected[1]<<1, poly_selected[0]*b**2
            
            tmp_block = sieve.sieve(b, poly_selected, logs, primes, a, n, param, sieve_len, skipped, prime_start, tmp1, tmp2, tmp3)
        
            to_batch += [abs(poly_selected[0]*i**2+tmp2*i+poly_selected[2]) for i in tmp_block]
            values += [poly_selected[0]*i+poly_selected[1] for i in tmp_block]
            coeff.append([poly_selected[0],poly_selected[1],len(tmp_block)])
            
            block += tmp_block
            
            while poly_index < threshold:
                poly_selected[1], poly_selected[2] = find_next_poly(poly_index, moduli, n, poly_selected, needed)
                tmp2 = poly_selected[1]<<1
                poly_index += 1
                
                tmp_block = sieve.sieve(b, poly_selected, logs, primes, a, n, param, sieve_len, skipped, prime_start, tmp1, tmp2, tmp3)
            
                to_batch += [abs(poly_selected[0]*i**2+tmp2*i+poly_selected[2]) for i in tmp_block]
                values += [poly_selected[0]*i+poly_selected[1] for i in tmp_block]
                coeff.append([poly_selected[0],poly_selected[1],len(tmp_block)])
                
                block += tmp_block
                
                if len(block) > 512:
                    smooth = find_smooth.batch_smooth_test(to_batch, prod_primes, const_1, const_2)
                    
                    rels = []
                    
                    tmp_index = 0
                    for z in range(len(coeff)):
                        for i in range(coeff[z][2]):
                            tmp_smooth = smooth[tmp_index+i]
                            value = values[tmp_index+i]

                            if tmp_smooth[0] == True:
                                rels.append([value, value*value-n, True])
                            elif tmp_smooth[0] == "large":
                                rels.append([value, value*value-n, False, tmp_smooth[1], tmp_smooth[2]])
                                
                        tmp_index += coeff[z][2]

                    tmp_rels.put(rels)
                        
                    block = []
                    coeff = []
                    to_batch = []
                    values = []

        except:
            pass

def siever(polynomials, tmp_rels, b, n, logs, primes, a, sieve_len, skipped, prime_start, const_1, const_2):
    while True:
        try:
            poly_selected, param, moduli, needed, threshold = polynomials.get(timeout=0)
            poly_index = 0
            tmp1, tmp2, tmp3 = poly_selected[0]<<1, poly_selected[1]<<1, poly_selected[0]*b**2
            
            coeff1 = poly_selected[0]
            coeff2 = poly_selected[1]
            coeff3 = poly_selected[2]
            tmp = coeff2<<1
            
            smooth = sieve.sieve(b, poly_selected, logs, primes, a, n, param, sieve_len, skipped, prime_start, tmp1, tmp2, tmp3)

            rels = []

            for i in smooth:
                tmp_smooth = find_smooth.smooth_test(coeff1*i*i+tmp*i+coeff3, primes, const_1, const_2)
                value = coeff1*i+coeff2
        
                if tmp_smooth[0] == True:
                    rels.append([value, value*value-n, True])
                elif tmp_smooth[0] == "large":
                    rels.append([value, value*value-n, False, tmp_smooth[1], tmp_smooth[2]])
            
            tmp_rels.put(rels)
            
            while poly_index < threshold:
                poly_selected[1],poly_selected[2] = find_next_poly(poly_index, moduli, n, poly_selected, needed)
                tmp2 = poly_selected[1]<<1
                poly_index += 1
                
                coeff1 = poly_selected[0]
                coeff2 = poly_selected[1]
                coeff3 = poly_selected[2]
                tmp = coeff2<<1
                
                smooth = sieve.sieve(b, poly_selected, logs, primes, a, n, param, sieve_len, skipped, prime_start, tmp1, tmp2, tmp3)
            
                for i in smooth:
                    tmp_smooth = find_smooth.smooth_test(coeff1*i*i+tmp*i+coeff3, primes, const_1, const_2)
                    value = coeff1*i+coeff2
            
                    if tmp_smooth[0] == True:
                        tmp_rels.put([value, value*value-n, True])
                    elif tmp_smooth[0] == "large":
                        tmp_rels.put([value, value*value-n, False, tmp_smooth[1], tmp_smooth[2]])
        except:
            pass

def find_relations(primes, CONST_LARGE_PRIME, prod_primes, bounds, target, logs, a, b, FLAG_USE_BATCH_SMOOTH_TEST, n, LOG_PATH, NB_CPU):
    log.write_log(LOG_PATH, "sieving...")
    log.write_log(LOG_PATH, "need to find at least "+str(len(primes)+10)+" relations")
    
    relations, smooth_number, partial_relations, possible_smooth, graph, cycle_len = [], [], {}, {}, {}, [0]*10
    size_partials = 0
    parent = {}
    
    partial_found, full_found, skipped, sieve_len = 0, 0, 0, (b<<1)+1
    const_1 = CONST_LARGE_PRIME*primes[-1]
    const_2 = CONST_LARGE_PRIME*primes[-1]**2
    skipped += int(math.log2(const_2))
    prime_start = 30
    skipped = sieve.compute_skipped(skipped, logs, primes, prime_start)
    
    cpu, sievers = min(NB_CPU, multiprocessing.cpu_count()), []

    polynomials = multiprocessing.Queue()
    tmp_rels = [multiprocessing.Queue() for _ in range(cpu-1)]

    for _ in range(2*(cpu-1)):
        polynomials.put(find_poly(n, primes, a, bounds, target))
    
    for k in range(cpu-1):
        if FLAG_USE_BATCH_SMOOTH_TEST:
            pro = multiprocessing.Process(target=siever_batch,args=(polynomials, tmp_rels[k], b, n, logs, primes, a, sieve_len,
                                                                    skipped, prime_start, const_1, const_2, prod_primes))
        else:
            pro = multiprocessing.Process(target=siever,args=(polynomials, tmp_rels[k], b, n, logs, primes, a, sieve_len,
                                                              skipped, prime_start, const_1, const_2))
        sievers.append(pro)
        pro.start()
    
    time_1 = datetime.now()
        
    sys.stdout.write('\r'+"0/("+str(len(primes)+1)+"+10) relations found")
    
    while len(relations) <= len(primes)+10:
        while polynomials.qsize() < 2*(cpu-1): polynomials.put(find_poly(n, primes, a, bounds, target))
        
        for k in range(cpu-1):
            try:
                rels = tmp_rels[k].get(timeout=0)

                for rel in rels:
                    value = rel[0]
                    if rel[2]:
                        tmp_smooth = [True]
                    else:
                        tmp_smooth = ["large", rel[3], rel[4]]
                        
                    relations, smooth_number, partial_relations, possible_smooth, full_found, partial_found, graph, size_partials, parent = handle_possible_smooth(value, tmp_smooth, full_found, partial_found, relations, smooth_number,
                                                                                                                                                                partial_relations, possible_smooth, graph, size_partials, parent,
                                                                                                                                                                cycle_len, n)
                    if len(relations) > len(primes)+10: break
            except:
                pass
                
        sys.stdout.write('\r'+str(len(smooth_number))+"/("+str(len(primes)+1)+"+10) relations found : full = "+str(full_found)+" ; partial = "+str(partial_found)+ " ("+str(len(possible_smooth))+")")
    print("\n")
    
    time_2 = datetime.now()

    while len(tmp_rels):
        tmp_rels[0].close()
        del tmp_rels[0]
    
    for p in sievers:
        p.terminate()
        p.join()

    polynomials.close()
    
    log.write_log(LOG_PATH, "sieving done in "+format_duration(time_2-time_1)+".\n")
    log.write_log(LOG_PATH, str(len(smooth_number)) + " relations found")
    log.write_log(LOG_PATH, str(full_found)+" full relations found")
    log.write_log(LOG_PATH, str(partial_found)+" partial relations found")
    log.write_log(LOG_PATH, str(len(possible_smooth))+" smooth with large prime(s) found")
    log.write_log(LOG_PATH, "Distribution of cycle length:")
    for i in range(9):
        log.write_log(LOG_PATH, str(i+2)+"-cycle: "+str(cycle_len[i]))
    log.write_log(LOG_PATH, "11+-cycle: "+str(cycle_len[-1])+"\n")
    
    time.sleep(1e-1)
                
    return relations, smooth_number 