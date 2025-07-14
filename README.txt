Thanks for looking at this project.

This is a project I started a few years ago, during the first half of 2019.
I re-worked it recently to merge some versions, make it more readable, add comments, and get all the sources I used.


##### Introduction ######

This project is a Python implementation of the Quadratic Sieve to factor integers.

Here are some detailed characteristics :
- Multiple polynomials used for sieving
- Double large prime variation used to collect relations
- Batch smoothness test and naive smoothness test available
- Parallel sieving allowed
- gaussian elimination, block Lanczos, Wiedemann algorithms available for the linear algebra step (no parallelization available for now)

##### Sources #####

Overall algorithm:
- "Prime numbers, a computational perspective" by Richard Crandall and Carl Pomerance (really good): https://link.springer.com/book/10.1007/0-387-28979-8
- "The Quadratic Sieve Factoring Algorithm" by Eric Landquist: https://www.cs.virginia.edu/crab/QFS_Simple.pdf

Multiple polynomials and double large primes
- "Factoring Integers with Large-Prime Variations of the Quadratic Sieve" by Henk Boender and Herman J. J. te Riele: https://projecteuclid.org/journals/experimental-mathematics/volume-5/issue-4/Factoring-integers-with-large-prime-variations-of-the-quadratic-sieve/em/1047565445.pdf

Batch smoothness test:
- "HOW TO FIND SMOOTH PARTS OF INTEGERS" by DANIEL J. BERNSTEIN: https://cr.yp.to/factorization/smoothparts-20040510.pdf

Gaussian elimination:
- Took the gaussian elimination code from https://github.com/skollmann/PyFactorise/blob/master/factorise.py#L39

Block Lanczos:
- "A Block Lanczos Algorithm for Finding Dependencies over GF(2)" by Peter L. Montgomery: https://scispace.com/pdf/a-block-lanczos-algorithm-for-finding-dependencies-over-gf-2-ezdu2qt0pp.pdf
- "A modified block Lanczos algorithm with fewer vectors" by Emmanuel Thomé (really good): https://eprint.iacr.org/2016/329.pdf (really good)

Wiedemann algorithm:
- "SOLVING HOMOGENEOUS LINEAR EQUATIONSOVER GF(2) VIA BLOCK WIEDEMANN ALGORITHM" by Don Coppersmith: https://www.ams.org/journals/mcom/1994-62-205/S0025-5718-1994-1192970-7/S0025-5718-1994-1192970-7.pdf

##### running the algortihm #####

When running the main.py file, the required argument is --n the number you want to factor.

A suitable command looks like "python main.py --n [your_number]"

##### config parameters #####

I have defined a few parameters directly you can edit in the config file.

- batch_smooth_test: If True, the batch smoothness test from Bernstein will be used. Otherwise, the naive (trial division from the prime factor base) approach will be used. Batch smoothness test is generally faster
- gaussian_pivot: If True, gaussian elimination will ALWAYS be used in the linear algebra step. If you want to use block Lanczos of Wiedemann algorithms, this has to be set to False.
- lanczos: If True, the block Lanczos algorithm will be used. If False, the Wiedemann algorithm will be used. No matter its value, if gaussian_pivot is True then gaussian elimination will be performed.
- large_primes_constant: Define the constant that will be multiplied by the last prime in the factor base to obtain the bound for the single large primes. For the double large primes, the bound is this constant multiplied by the last prime in the factor base squared.
- block_size: If Block Lanczos or Wiedemann algorithm are used, this sets the block size. In the Wiedemann algorithm, this is only used to compute many matrix-vector products at once, and not to compute the matrix generator of the matrix sequence you obtain.
- NB_CPU: Sets the number of cpu used for sieving. One cpu is always kept as a "leader" that collects relations from the sievers, and the rest of the cpus sieve, test for smoothness, and send their candidate relations to the leader.

##### General discussion #####

I tried overall to use as few exernal libraries as possible. While this makes the code way harder to write for me, and understand for you,
I viewed it as a good challenge to write good code, and get to know hard and foreign concepts to me. It took a lot of work, but I got a bug free code,
and got much more knowledgable on the subroutines thanks to this.

Here are some points I consider working on at some point:
- The data architecture for the large primes relations graph may not be optimal.
	Right now it is an array of arrays. The first index always indicates the smallest large prime of the two.
	The global array is always sorted according to the first value of its arrays.
	This allows for efficient cycle search.
	However, this leads to a very heavy code.
	I consider giving a try to representing the graph in a dict of array, which would allow finding a large prime in the graph very fast.
	However, I have to check the complecity of checking IF a large prime exists in the graph.
	If this is as good as binary search, then I will leave it as it is now.

- No matter how hard I tried, I am stuck on understanding the block Wiedemann algorithm. For now, the best I can do is the scalar one, with some optimizations.
	Namely, I use binary encoding of the blocks of vectors to compute matrix-vector product very efficiently. Then, for each scalar sequence, I compute
	its generator using scalar Berlekamp-Massey algorithm, and use it to update the minimal polynomial of the matrix. I consider working some way to
	compute the matrix generator of the matrix sequence through some algorithm (matrix Berlekamp-Massey or some other).

Here are the next steps:

- Some optimizations may be possible here and there, but I have already a big amount of work to get fast code. The next step on this topic is switching to C.
	Python is intrisically limited on the topic of speed, as it is compiled on the fly. Just like this project, I have already done part of the C project
	a while ago. I have to rework it to get something on the same level of quality as this Python version.

- I will work on the General Number Field Sieve. Just like this project, I have already done part of the GNFS project a while ago. I have to rework it to get
	something of the same level of quality as this QS algorithm.