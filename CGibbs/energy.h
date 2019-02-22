//
//  energy.h
//
//  Created by Jeewoen Sin on 5/3/17.
//

#ifndef energy_h
#define energy_h

class energy{

public:
	// variables:
    double homology;
    int orthology;

    //double co_occur;

    double phylo_corr;
    int max_phylo_corr_gene_index;

    double gene_cluster;
    int max_gene_cluster_gene_index;

    double co_exprsn;
    int max_co_exprsn_gene_index;

    double no_EC_penalty;

    // parameters:
    double min_e_ph;
    double min_e_gc;
    double min_e_ce;

    // functions:
    energy();

    void initialize(const double &min_ph, const double &min_gc, const double &min_ce);

};
#endif

