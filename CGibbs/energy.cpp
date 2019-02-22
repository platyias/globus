//
//  energy.cpp
//
//  Created by Jeewoen Sin on 5/3/17.
//

#include "energy.h"

energy::energy(){
}

void energy::initialize(const double &min_ph, const double &min_gc, const double &min_ce){

    // TODO //
    min_e_ph = min_ph;
    min_e_gc = min_gc;
    min_e_ce = min_ce;

    homology = 0.0; // 0
    orthology = 0; // 8

    //co_occur = 0.0;  // the average correlation between the EC activity of the assigned gene and the EC activities for all its network neighbors

    phylo_corr = min_e_ph; //0.0; //min_e_ph; // 2
    max_phylo_corr_gene_index = -1; // -1 == not assigned yet.

    gene_cluster = min_e_gc; //0.0; //min_e_gc; // 3
    max_gene_cluster_gene_index = -1; // -1 == not assigned yet.

    co_exprsn = min_e_ce; //0.0; //min_e_ce; // 4
    max_co_exprsn_gene_index = -1; // -1 == not assigned yet.

    no_EC_penalty = 0.0; // either 1 or 0

    // neibhors.. Do I need the following?
    // 5,6,7 = ph, gc, co_exp
}

