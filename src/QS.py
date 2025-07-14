# This is the file that execute the Quadratic Sieve algorithm

## import libraries

import math, time
from datetime import datetime
import log
import parse_config
from utils import *
import mono_cpu_sieve
import multi_cpu_sieve
import generate_primes
import sieve
from sieve_polynomials_functions import *
import find_smooth
from relations import *
from build_matrix import *
from gaussian_elimination import *
from wiedemann import *
from block_lanczos import *
import compute_solutions

## Set path to find config file

CONFIG_PATH = "C:\\Users\\basil\\OneDrive\\Bureau\\Git\\Quadratic Sieve Python\\config\\config.ini"

## display functions
    
def print_final_message(x, y, n, time_1, LOG_PATH):
    time_2 = datetime.now()
    log.skip_line(LOG_PATH)
    log.write_log(LOG_PATH, "found factor : "+str(math.gcd(x-y,n)))
    log.write_log(LOG_PATH, "found factor : "+str(math.gcd(x+y,n)))
    log.skip_line(LOG_PATH)
    log.write_log(LOG_PATH, "null space found in "+format_duration(time_2-time_1)+".\n")
    
## QS functions
    

def set_bounds(primes, n):        
    low_bound_10,low_bound_100,low_bound_1000,up_bound_3000 = -1,-1,-1,-1
    for i in range(len(primes)):
        if primes[i] > 10 and low_bound_10 == -1: low_bound_10 = i
        elif primes[i] > 100 and low_bound_100 == -1: low_bound_100 = i
        elif primes[i] > 1000 and low_bound_1000 == -1: low_bound_1000 = i
        elif primes[i] > 3000 and up_bound_3000 == -1:
            up_bound_3000 = i
            break
    if low_bound_1000 != -1 and up_bound_3000 == -1 and len(primes)-1-low_bound_1000 > 20:
        up_bound_3000 = len(primes)-1
    elif up_bound_3000 == -1: return "try normal QS"
    return [low_bound_10,low_bound_100,low_bound_1000,up_bound_3000]
    
def initialize(n):
    
    b = int(math.exp(math.sqrt(math.log(n)*math.log(math.log(n)))/2)//3.2)
    
    prime_candidates = generate_primes.create_smooth_primes_base(b)
    
    primes = []
    a = []
    logs = []
    prod_primes = 1
    for p in prime_candidates:
        if compute_legendre_character(n,p) == 1:
            primes.append(p)
            prod_primes *= p
            a.append(compute_sqrt_mod_p(n,p))
            logs.append(round(math.log2(p)))
    target = math.log10(2)/2+math.log10(n)/2-math.log10(b)
    
    return b, primes, a, logs, prod_primes, target
    
def find_null_space_and_compute_factors(relations, smooth_number, primes, n, flag_gaussian_pivot, flag_lanczos, BLOCK_SIZE, LOG_PATH):
    if not flag_gaussian_pivot:
    
        bin_matrix = build_sparse_matrix(relations,primes)
        
        log.write_log(LOG_PATH, "matrix created "+str(len(primes)+1)+"x"+str(len(smooth_number))+" solving...")
    
        reduce_sparse_matrix(bin_matrix,relations,smooth_number)
        e = (len(bin_matrix)+10)
        relations = relations[0:e]
        smooth_number = smooth_number[0:e]
        log.write_log(LOG_PATH, "matrix reduced to "+str(len(bin_matrix))+"x"+str(len(relations))+" solving...\n")
        
        nb_attempts = 1
        time_1 = datetime.now()
        
        mini_poly_estim = 1
        
        while True:
            if flag_lanczos:
                null_space = block_lanczos(bin_matrix, len(primes)+1, len(smooth_number), BLOCK_SIZE, LOG_PATH)
            else:
                null_space, mini_poly_estim = wiedemann(bin_matrix,len(relations),BLOCK_SIZE,mini_poly_estim)
                reduce_null_space_vectors(null_space,len(primes)+1)
            
            log.write_log(LOG_PATH, "attempt "+str(nb_attempts)+": "+str(len(null_space))+" kernel vectors found")
            for vector in null_space:
                if flag_lanczos:
                    vector = compute_solutions.convert_to_binary_lanczos(vector, smooth_number)
                x,y = compute_solutions.compute_solution(relations,smooth_number,vector,n,primes)
                if x != y and math.gcd(x-y,n) != 1 and math.gcd(x+y,n) != 1:
                    print_final_message(x, y, n, time_1, LOG_PATH)
                    return str(time.localtime()[3])+":"+str(time.localtime()[4])+":"+str(time.localtime()[5]), math.gcd(x-y,n), math.gcd(x+y,n)
            nb_attempts += 1

    else:
        bin_matrix = build_dense_matrix(relations, primes)
        bin_opt, bin_n,bin_m = siqs_build_matrix_opt(bin_matrix)
        time_1 = datetime.now()
        log.write_log(LOG_PATH, "matrix created "+str(len(primes)+1)+"x"+str(len(smooth_number))+" solving...")
        null_space = siqs_solve_matrix_opt(bin_opt,bin_n,bin_m)
        log.write_log(LOG_PATH, "matrix created "+str(len(primes)+1)+"x"+str(len(smooth_number))+" null space found !")
        
        for z in null_space:
            bin_encoding = compute_solutions.convert_to_binary(z, smooth_number)

            x,y = compute_solutions.compute_solution(relations,smooth_number,bin_encoding,n,primes)

            if x != y and math.gcd(x-y,n) != 1 and math.gcd(x+y,n) != 1:
                print_final_message(x, y, n, time_1, LOG_PATH)
                return str(time.localtime()[3])+":"+str(time.localtime()[4])+":"+str(time.localtime()[5]), math.gcd(x-y,n), math.gcd(x+y,n)
    
def QS(n):
    now = datetime.now()
    
    LOG_PATH = "C:\\Users\\basil\\OneDrive\\Bureau\\Git\\Quadratic Sieve Python\\logs\\log_"+str(now.year)+str(now.month)+str(now.day)+"_"+str(now.hour)+str(now.minute)+str(now.second)+".txt"
    
    parameters = parse_config.parse_config(CONFIG_PATH)
    
    flag_use_batch_smooth_test = parameters[0].lower() in ["true"]
    flag_gaussian_pivot = parameters[1].lower() in ["true"]
    flag_lanczos = parameters[2].lower() in ["true"]
    const = int(parameters[3])
    BLOCK_SIZE = int(parameters[4])
    NB_CPU = int(parameters[5])
    
    n = int(n)
    
    b, primes, a, logs, prod_primes, target = initialize(n)
    
    bounds = set_bounds(primes, n)
    if bounds == "try normal QS":
        log.write_log(LOG_PATH, "Try normal QS")
        return 0
    
    log.write_log(LOG_PATH, "factor base created : "+str(len(primes))+" primes\n")
    
    if NB_CPU == 1:
        relations, smooth_number = mono_cpu_sieve.find_relations(primes, const, prod_primes, bounds, target, logs, a, b, flag_use_batch_smooth_test, n, LOG_PATH)
    
    else:
        relations, smooth_number = multi_cpu_sieve.find_relations(primes, const, prod_primes, bounds, target, logs, a, b, flag_use_batch_smooth_test, n, LOG_PATH, NB_CPU)
    
    e = (len(primes)+11)
    relations = relations[0:e]
    smooth_number = smooth_number[0:e]
    
    return find_null_space_and_compute_factors(relations,smooth_number,primes,n,flag_gaussian_pivot, flag_lanczos,BLOCK_SIZE,LOG_PATH)