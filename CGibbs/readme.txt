vector_CGibbs_min_e_cnt_d_nbEC_noG

non-zero min_e_ph,.. context descriptors

output: ++ for all positioned ecs at each perturbation 

co_occurence fixed: do not accumulate difference for all genes. No matter how many genes sit on nodes, just consider the occupied ec connections.
median_E co-occurence: unpenalize EC-EC co-occurence for ECs that do not have any possible EC connections. (EC_EC_link[ec].size() == 1: itself) If there are possible multiple nb ECs but no genes are on those nb ECs, that case is not the case.

context descriptor: Whenever a gene is in the network (in-node), the minimun score is min_e_X, not 0 (no matter EC_EC_link[ec].size() == 1 or > 1, if the gene doesn't have nb genes to calculate gene-gene context, give min_e_X to the gene). 

log p mapping function: becareful. b parameters are set for 10-base calculattion. Not just for 3 context descriptors and homology, the orthology & EC co-occurence b values are also for 10 to the power of X.


