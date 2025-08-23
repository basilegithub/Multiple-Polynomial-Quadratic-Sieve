Thanks for looking at this project.

This is a project I started a few years ago, during the first half of 2019.
I re-worked it recently to merge some versions, make it more readable, add comments, and get all the sources I used.


##### Introduction ######

This project is a Python implementation of the Quadratic Sieve to factor integers.

Here are some detailed characteristics :
- Multiple polynomials used for sieving
- Double large prime variation used to collect relations
- Union-find algorithm to detect cycle in the partial relations graph
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
- "A modified block Lanczos algorithm with fewer vectors" by Emmanuel Thomé (really good): https://eprint.iacr.org/2016/329.pdf

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

I tried overall to use as few exernal libraries as possible. While this makes the code 
	way harder to write for me, and understand for you, I viewed it as a good challenge
	to write good code, and get to know hard and foreign concepts to me. It took a lot
	of work, but I got a bug free code, and got much more knowledgable on the subroutines
	thanks to this.

I have a Python implementation of the GNFS algorithm on my Github.

Here are some points I consider working on at some point:
- The implementation of the Union-find algorithm for might not be optimal.
	Joining two connected componentes takes O(n) times.
	I have seen how to modify it, I just have to take the time to implement it.

- No matter how hard I tried, I am stuck on understanding the block Wiedemann algorithm.
	For now, the best I can do is the scalar one, with some optimizations. Namely, I
	use binary encoding of the blocks of vectors to compute matrix-vector product very
	efficiently. Then, for each scalar sequence, I compute its generator using scalar
	Berlekamp-Massey algorithm, and use it to update the minimal polynomial of the matrix.
	I consider working some way to compute the matrix generator of the matrix sequence
	through some algorithm (matrix Berlekamp-Massey or some other).

- I have to make the linear algebra step more parallel. Even though it is the fastest of the
	two main steps of the algorithm, no one likes to wait if it is possible to go faster.

- I had seen that there was some work on how to chose k such that k*n is much easier to factor.
	If I remember correctly, you test multiple small values of k, and see for which
	value of k, k*n has the most quadratic residues. You weight the value of the quadratic
	residues by some function of p. I have yet to choose if I do it or not, as it implies
	many modifications for a speedup that I don't have a good estimate for. This is akin
	to the polynomial selection step of the GNFS.

Here are the next steps:

- Some optimizations may be possible here and there, but I have already a big amount of work
	to get fast code. The next step on this topic is switching to C. Python is intrisically
	limited on the topic of speed, as it is compiled on the fly. Just like this project, I
	have already done part of the C project a while ago. I have to rework it to get
	something on the same level of quality as this Python version.

- Likewise, implement the GNFS algorithm in C.