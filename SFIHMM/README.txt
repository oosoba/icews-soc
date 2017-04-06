SFIHMM 
http://bit.ly/sfihmm

high-speed C code for the estimation of Hidden Markov Models (finite state machines) on arbitrary time-series, for Viterbi Path Reconstruction, PCCA+ (Perron-Cluster Cluster Analysis), and for the generation of simulated data from HMMs. Based on the code used in "Conflict and Computation on Wikipedia: a Finite-State Machine Analysis of Editor Interactions”, where these methods were used to detect critical transitions in symbolic time series.

Now includes glue code for parallel use of multiple cores.

This is a pre-release version of the code; feedback and requests for additional features solicited; please write simon[at]santafe.edu, or ping me on Twitter at @simondedeo  

1. Installing 
2. Getting started: estimating an HMM from data
3. How many hidden states?
4. Parallel SFIHMM
5. Viterbi Path Reconstruction
6. Module Identification — Critical Transition detection 
7. Generating new sequences
8. Polishing solutions
9. Notes / Citing / Acknowledgements

Installing
SFIHMM works on Mac OS X and UNIX systems such as GNU/Linux. You’ll need the gsl libraries. On the Mac, you can load these in using tools such as macports (sudo port install gsl) or homebrew. You’ll know it’s working if you can run the command gsl-config. We know it works on both Ubuntu and Gentoo Linux, and CentOS, as well as Mac OS X; for Linux boxes you may need to install (via apt-get) gsl-bin, libgsl0-dev, and libgsl0ldb.

Once you’ve installed gsl, download SFIHMM.zip, and uncompress by typing "unzip SFIHMM.zip". Type "cd SFIHMM" (i.e., go to the SFIHMM directory), and type "make". You should now be able to run ./hmm when inside the SFIHMM directory. You’re off and running! I’ve included a test file, test_data, so you can experiment with the different options available for the code. 
 
Estimating a Hidden Markov Model (pre-specified number of states)
The most basic task: given a sequence, find the best-fitting Hidden Markov Model that describes it.
ayerie:~simon$ ./hmm -f [N_iterations] [N_states] [name of data file] 
 
Input
data file. An ASCII text file with two lines; I’ve given you a sample file, test_data, to play with, but you’ll want to use your own eventually!
Line One: the total number of observations in the sequence. 
Line Two: the sequence itself. The code reads in all the symbols, and orders them alphabetically (by ASCII order). Don’t use multi-letter symbols!—you may need to pre-process your data to create a proper input file to SFIHMM. For example, you can use all of the letters from A to Z, and the all of the lowercase letters from a to z, to get up to 52 different symbols. If a symbol doesn’t appear in the input stream, SFIHMM can’t know it exists, so (e.g.) if you have A, B, and D in the data, then SFIHMM will call A symbol one, B symbol two, and D symbol three. Beware!

Command line options
[N_iterations]: the number of times to run SFIHMM’s EM algorithm from random initial conditions. Landscapes are "rugged", meaning that there can be multiple solutions, some of which are sub-optimal. By increasing the number of times you run SFIHMM, you increase your chances of finding the optimal machine—at the cost of running the algorithm for longer, of course!
[N_states]: the number of states in the final output. 

Example:
ayerie:~simon$ ./hmm -f 100 5 test_data
Fit a five-state machine; run the algorithm 100 times; return the highest-likelihood machine. 
 
Output
What you want is sent to the file [name of file]_OUT_[N_states]states ; for the example above, that would be test_data_OUT_5states
It is an ASCII text file listing the log-likelihood, the number of states, the number of output symbols, the state-to-state transition probabilities, the symbol emission probabilities for each state, the list of symbols in order, the stationary distribution over the hidden states, and λ2, the second eigenvalue (see Equation 1 of the Wikipedia paper).

Estimating a Hidden Markov Model (without specifying the number of states)—recommended!
What happens when you don’t have a good idea how many states you should use? Too few, and you miss important patterns; too many, and you are in danger of overfitting to noise. The general problem of how to select the correct number of states is called "model selection". One model selection solution (far from the only one) is to use AIC (the "Akaike Information Criterion"). SFIHMM will do this for you—all you have to do is leave off the [N_states] option:
ayerie:~simon$ ./hmm -f [N_iterations] [name of data file]
For example:
ayerie:~simon$ ./hmm -f 100 test_data
I’d recommend using the AIC model selection tool rather than pre-specifying the number of states. It is hard to get an intuition about the "right" number of states to use; tests in our paper suggest that AIC does very well in approximating the true number of states.

Parallel SFIHMM: fitting HMMs with multiple cores
SFIHMM now has a little ruby glue code to run the -f option in parallel. This is in test mode! It should be tidy and clean up after itself.
ayerie:~simon$ ./parallel_glue.rb [number of cores] [number of iterations per core] [filename]
For example:
ayerie:~simon$ ./parallel_glue.rb 8 10 test_data 
To state the obvious: you can either run on multiple cores with the same number of iterations (takes the same amount of time, but better fits); or with a smaller number of iterations (if the number of cores times the number of iterations is the same, the fit quality won’t decline, but it will run ncores times faster.)

Viterbi path reconstruction
Say you’ve found a good model to describe your data. As the timeseries unfolds, what internal state is the system most likely to be in over time? Viterbi Path Reconstruction finds the most likely path through the system, given a set of data, and an accompanying model. Note that the data you feed it doesn’t have to be the same data you used to build the model! (But it can’t have any new symbols the original model didn’t see.) Use the -v ("Viterbi") option. 
ayerie:~simon$ ./hmm -v [name of data file] [name of HMM file]
For example:
ayerie:~simon$ ./hmm -v test_data test_data_OUT_4states
The output is a list of numbers, specifying which hidden state the system is most likely to be in at each step in the time-series.

Module identification and path reconstruction
Depending on your system, the relaxation time may be very long, and the system may be trapped in a subspace of the available states. Using the methods described in the Wikipedia paper (see Sec. 1.4), we can use the second eigenvector to find the two main system modules. Using the -m ("Module") option, SFIHMM will identify these modules; it will name one H ("Higher probability of emitting symbol 1") and the other L, and output the path the system takes between them
ayerie:~simon$ ./hmm -m [name of data file] [name of HMM file]
For example:
ayerie:~simon$ ./hmm -m test_data test_data_OUT_4states
The output is a list of letters (H or L), specifying which module the system is most likely to be in at each step in the time-series. It will also tell you the state membership in the two modules, the probability of emitting the first symbol conditional upon being within the module (using the stationary distribution), and the empirical trapping times in the two states.

Generating new sequences
Say you’ve found a good model to describe your data, and you’d like to generate more of it. SFIHMM is happy to help; with the -g ("generate") option, tell SFIHMM the number of letters you’d like, and feed it the HMM you desire.  
ayerie:~simon$ ./hmm -g [N_letters] [name of HMM file]
For example:
ayerie:~simon$ ./hmm -g 100 test_data_OUT_4states

Polishing solutions
The standard SFIHMM algorithm terminates when it can not improve the solution beyond a tolerance of 10-3 in log-likelihood. After finding a solution, however, you can polish it using the -p option:
ayerie:~simon$ ./hmm -p [name of data file] [name of HMM file] [tolerance]
For example:
ayerie:~simon$ ./hmm -p test_data test_data_OUT_4states 1e-6

Notes
I’ve written SFIHMM with the goal of maximal simplicity. For example, you don’t specify the symbol space; SFIHMM figures it out from the data itself. And its output HMM descriptions are also valid file formats for input.

Citation
SFIHMM was written for the calculations first described and presented in this paper, including the use of Viterbi path reconstruction, and λ2 / trapping time characterization.

BibTeX
@Article{dedeo16,
AUTHOR = {DeDeo, Simon},
TITLE = {Conflict and Computation on {W}ikipedia: A Finite-State Machine Analysis of Editor Interactions},
JOURNAL = {Future Internet},
VOLUME = {8},
YEAR = {2016},
NUMBER = {3},
PAGES = {31},
URL = {http://www.mdpi.com/1999-5903/8/3/31},
ISSN = {1999-5903},
DOI = {10.3390/fi8030031}
}

Acknowledgements
The development of SFIHMM was supported in part by the National Science Foundation under Grant No. EF-1137929, "the small-number limit of biological information-processing", and by the Santa Fe Institute through an Omidyar Fellowship award to Simon DeDeo. Any opinions, findings and conclusions or recommendations expressed by SFIHMM are those of SFIHMM and do not necessarily reflect the views of the National Science Foundation (NSF) or any other funders. I thank Anne Sallaska, David Slater, and Haven Liu (MITRE Corp) for beta-testing this code.
