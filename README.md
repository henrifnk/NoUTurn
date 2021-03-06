# NoUTurn
No-U-Turn-Sampler introduced by Hoffman and Gelman in 2014

## Library used:

The No-U-Turn Sampler:  Adaptively Setting Path Lengths in Hamiltonian Monte Carlo (Hofmann 2014)
	Topics: Main Paper introducing NUTS
	Link:	http://www.stat.columbia.edu/~gelman/research/published/nuts.pdf
	
NIPS Learning (2011 Hofmann) [Youtube Video]
	Topics: Supportiv Presentation of NUTS (especially tree building)
	Link:	

Baysian Data Analysis (Gelman 3rd Edition), Chapter 11, 12
	Topics: MCMC, Metropolis(-Hastings), Hamiltonian MC
	Link:	https://statisticalsupportandresearch.files.wordpress.com/2017/11/bayesian_data_analysis.pdf
	
Markov Chains: Why Walk When You Can Flow? (Richards 2017)
	Topics: Vizualization Problems of MH and HMC, showing how NUTS can solve them
	Link:	https://elevanth.org/blog/2017/11/28/build-a-better-markov-chain/

Hamiltonian Monte Carlo explained
  Topics: From MCMC to hamiltonian Monte Carlo with nice vizualizations
  Link: https://arogozhnikov.github.io/2016/12/19/markov_chain_monte_carlo.html
	
A  Conceptual  Introduction  to Hamiltonian  Monte  Carlo (Betancourt 2018), Acknowledges A
	Topics: Technical implementation of trajectory and transition
	Link:	https://arxiv.org/pdf/1701.02434.pdf

Python: Hamiltonian Monte Carlo from scratch (Moore, 2020)
	Topics: Code implementation of HMC algorithm in python
	Link:	https://towardsdatascience.com/python-hamiltonian-monte-carlo-from-scratch-955dba96a42d

MCMC from scratch (Moore, 2020)  [Colab Reasearch Doc]
	Topics: Code implementation and Regression for MCMC/ Good Intuition for HMC energy analogy
	Link:	https://colab.research.google.com/drive/1YQBSfS1Nb8a9TAMsV1RjWsiErWqXLbrj#scrollTo=HUoQiEtDX_xL
	
Slice Sampling (Neal 2003)
  Topics: Slice Sampler for univariate and multivariate Problems
  Link: https://projecteuclid.org/download/pdf_1/euclid.aos/1056562461
  
Stan:A Probabilistic Programmin Language
  Topics: Stanis a probabilistic programming language for specifying statistical models
  Link: http://www.stat.columbia.edu/~gelman/research/published/stan-paper-revision-feb2015.pdf
  
Regression (Fahrmeir, 2013)
  Topics: Bayesian Generalized Regression
  
MCMC Interactive Gallery (Chi Feng, 2016)
  Topics: Vizualization of different MCMC algorithms
  Link: https://chi-feng.github.io/mcmc-demo/app.html?algorithm=MALA&target=banana



## Proposed Structure
	
	1) Brief Intro to MH-Algorithms random walk vs. HMC Idea of supressing local random walk behavior
	2) Brief Intro to Hamiltonian HMC with focus on tuning parameters L & \epsilon 
		a) Leapfrog algorithm and physical analogy (frictionless bowl as intuition)
		b) Influence and setting heurictics for the tuning parameters stepsize and Leapfrog-steps (locally adaptive HMC)
		   Supply typical problems for wrongly setting L to high/low (like trajectory bites own tail)
	3) Introduce the "No-U-Turn" Idea
	4) Provide Code for 
		a) naive NUTS algorithm
		b) effective NUTS algorithm
	5) Introduce "on the fly" tuning of stepsize
	  a) find reasonable epsilon
	  b) dual averaging
	*Idea 4) and 5)* show developement of algorithm by diffs from Metropolis Hastings to efficient NUTS with dual averaging
	6) Synthetic simple examples for Regression and Classification
		a) Benchmarking algorthims (traceplots, convergence diagnostic etc) => Pros and Cons
	7) RealWorld example/ Maybe little outlook on apllications like Stan!?
