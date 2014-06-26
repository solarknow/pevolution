-------------------------------------------------------------------------------------------------
Basic usage: 
 - Sequential version: 
        java -jar prottest-2.1.jar -i alignm_file [OPTIONS]
 - Parallel version: 
        mpjrun.sh -wdir $PWD/ -np [NUM_PROCS] -jar ModelTest-2.1.jar -i alignm_file [OPTIONS]
OPTIONS:
 -i alignment_filename
            Alignment input file (required)
 -t tree_filename
            Tree file       (optional) [default: NJ tree]
 -o output_filename
            Output file     (optional) [default: standard output]
 -log enabled/disabled
            Enables / Disables PhyML logging into log directory (see prottest.properties)
 -[matrix]
            Include matrix (Amino-acid) = JTT LG DCMut MtREV MtMam MtArt Dayhoff WAG 
	                                  RtREV CpREV Blosum62 VT HIVb HIVw FLU 
                If you don't specify any matrix, all matrices displayed above will
                be included.
 -I
            Include models with a proportion of invariable sites
 -G
            Include models with rate variation among sites and number of categories
 -IG
            include models with both +I and +G properties
 -all-distributions
            Include models with rate variation among sites, number of categories and both
 -ncat number_of_categories
            Define number of categories for +G and +I+G models [default: 4]
 -F
            Include models with empirical frequency estimation 
 -AIC
            Display models sorted by Akaike Information Criterion (AIC)
 -BIC
            Display models sorted by Bayesian Information Criterion (BIC)
 -AICC
            Display models sorted by Corrected Akaike Information Criterion (AICc)
 -DT
            Display models sorted by Decision Theory Criterion
 -all
            Displays a 7-framework comparison table
 -S optimization_strategy
            Optimization strategy mode: [default: 0]
             		0: Fixed BIONJ JTT
             		1: BIONJ Tree
             		2: Maximum Likelihood tree
             		3: User defined topology
 -s moves
            Tree search operation for ML search: 
            NNI (fastest), SPR (slowest), BEST (best of NNI and SPR) [default: NNI]
 -t1      				
            Display best-model's newick tree [default: false]
 -t2      				
            Display best-model's ASCII tree  [default: false]
 -tc consensus_threshold 
            Display consensus tree with the specified threshold, between 0.5 and 1.0
            [0.5 = majority rule consensus ; 1.0 = strict consensus]
 -threads number_or_threads			
            Number of threads requested to compute (only if MPJ is not used) [default: 1]
 -verbose
            Verbose mode [default: false]
-------------------------------------------------------------------------------------------------
Example: 
- Sequential version:
    java -jar ModelTest-2.1.jar -i alignm_file -t tree_file -S 0 -all-distributions -F -AIC -BIC -tc 0.5 > output
- Parallel version:
    mpjrun.sh -wdir $PWD/ -np 2 -jar ModelTest-2.1.jar -i alignm_file -t tree_file -S 0 -all-distributions -F -AIC -BIC -tc 0.5
