//
//  fitness_calculate.h
//
//  Created by Jeewoen Shin on 5/3/17.
//

#ifndef fitness_calculate_h
#define fitness_calculate_h

#include "energy.h"

#include <vector>
using std::vector;

#include <vector>
using std::vector;

#include "tnt_126/tnt_vec.h"
using TNT::Vector;

#include "tnt_126/tnt_cmat.h"
using TNT::Matrix;



void initial_gene_ec_assignment(
        const int &num_genes,
        const int &num_ECs,
        const vector< vector <int> > &EC_EC_link,
        const vector< vector <int> > &gene_EC_link,
        const vector< vector <double> > &EC_EC_co_occur,
        const double &median_ec_co_occur,
        const Matrix<double> &gene_gene_phylo_corr,
        const Matrix<double> &gene_gene_cluster_corr,
        const Matrix<double> &gene_gene_coexp_corr,
        const vector< vector <double> > &gene_EC_homol,
        const vector< vector <int> > &gene_EC_orth,
        vector< vector<int> > &genes_on_EC,
        vector<int> &assigned_gene2ec,
        vector< energy > &gene_energy,
        vector<double> &EC_co_occur_vec,
        vector<int> &EC_cnt_nbEC_vec
    );

void pull_out_gene(
        const double &min_e_ph, const double &min_e_gc,  const double &min_e_ce,
        const int &num_ECs,
        const int &gene, const int &org_ec_loc,
        vector<energy> &gene_energy,
        const vector< vector <int> > &EC_EC_link,
        const vector< vector <double> > &EC_EC_co_occur,
        const double &median_ec_co_occur,
        vector< vector<int> > &genes_on_EC,
        const Matrix<double> &gene_gene_phylo_corr,
        const Matrix<double> &gene_gene_cluster_corr,
        const Matrix<double> &gene_gene_coexp_corr,
        vector<double> &EC_co_occur_vec, vector<int> &EC_cnt_nbEC_vec,
        energy &delta_E, double &delta_EC_co_occur
    );

void put_back_new_loc(
        const int &num_ECs,
        const double &min_e_ph, const double &min_e_gc, const double &min_e_ce,
        const int &gene, const int &new_ec_loc, const vector< energy > &gene_energy,
        const vector< vector <int> > &EC_EC_link,
        const vector< vector <int> > &gene_EC_link,
        vector< vector<int> > &genes_on_EC,
        const vector< vector <double> > &gene_EC_homol,
        const vector< vector <int> > &gene_EC_orth,
        const Matrix<double> &gene_gene_phylo_corr,
        const Matrix<double> &gene_gene_cluster_corr,
        const Matrix<double> &gene_gene_coexp_corr,
        const vector< vector <double> > &EC_EC_co_occur,
        const double &median_ec_co_occur,
        vector< vector<int> > &updated_ph_index,
        vector< double > &updated_ph,
        vector< vector<int> > &updated_gc_index,
        vector< double > &updated_gc,
        vector< vector<int> > &updated_ce_index,
        vector< double > &updated_ce,
        vector<int> &updated_co_occur_EC_vec,
        vector< double > &updated_co_occur,
        vector<int> &updated_co_occur_nb_cnt,
        const vector<double> &EC_co_occur_vec, const vector<int> &EC_cnt_nbEC_vec,
        energy &delta_E, double &delta_EC_co_occur
    );

 void switch_loc_genes(
        const int &num_ECs,
        const double &min_e_ph, const double &min_e_gc, const double &min_e_ce,
        const double &coeff_hm, const double &coeff_ort,
        const double &coeff_ph, const double &coeff_gc, const double &coeff_ce,
        const double &p_ec_co_occur,
        const double &p_no_EC,
        const int &rand_gene_index,
        const vector< vector <int> > &EC_EC_link,
        const vector< vector <int> > &gene_EC_link,
        const vector< vector <double> > &EC_EC_co_occur,
        const double &median_ec_co_occur,
        const vector< vector <double> > &gene_EC_homol,
        const vector< vector <int> > &gene_EC_orth,
        const Matrix<double> &gene_gene_phylo_corr,
        const Matrix<double> &gene_gene_cluster_corr,
        const Matrix<double> &gene_gene_coexp_corr,
        vector< vector<int> > &genes_on_EC,
        int &new_ec_loc_index,
        vector<int> &assigned_gene2ec,
        vector< energy > &gene_energy,
        vector<double> &EC_co_occur_vec, vector<int> &EC_cnt_nbEC_vec
    );
/*
void calculate_ec_co_occur(
        const vector< vector<int> > &genes_on_EC,
        const vector <int> &new_ec_nb_ecs,
        const vector <double> &new_ec_co_occur,
        double &delta_co_occur
    );
*/
//void remove_gene_on_EC(vector< vector<int> > &genes_on_EC, const int &ec_loc, const int &gene);

void remove_gene_on_EC(vector<int> &genes_on_EC_vec, const int &ec_loc, const int &gene);

void find_index_of_id(const vector<int> &id_vec, const int &id, int &index);

void select_new_loc(
        const double &coeff_hm,
        const double &coeff_ort,
        const double &coeff_ph,
        const double &coeff_gc,
        const double &coeff_ce,
        const double &p_ec_co_occur,
        const double &p_no_EC,
        const vector<energy> &delta_E_vec,
        const vector<double> &delta_co_occur_vec,
        int &new_ec_loc_index
    );

//int locate_gene(const vector<double> &vec_prob);
void locate_gene(const vector<double> &vec_prob, int &index);

void unique_EC_occupancy(const vector<int> &assigned_gene2ec, Vector<int> &count_ec_on_init);

#endif
