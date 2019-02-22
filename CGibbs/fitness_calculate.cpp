//
//  fitness_calculate.cpp
//
//  Created by Jeewoen Shin on 5/3/17.
//

#include "parameters.h"
#include "fitness_calculate.h"
#include "random.h"

#include <cmath> // exp
#include <algorithm> // sort, unique
#include <cmath>       /* isfinite, sqrt */
#include <iostream>
using std::endl;
using std::cerr;


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
        ){
    // gene_energy[i]'s are all initialized in main.cpp
    for (int i=0;i<num_genes;i++) {
        int index;

        // gene_EC_link is always >0 (1 if only outside node)
        if(gene_EC_link[i].size() > 1) {
            // if this gene has possible EC assignments within the network
            // initially assign all genes to node in network. not to the out-network-node ??
            // the order of EC index is the same for gene_EC_link, gene_EC_homol, gene_EC_orth
            index = int_rand( gene_EC_link[i].size()-1 ); // index = 0, .., gene_EC_link[i].size()-1-1
        }
        else if(gene_EC_link[i].size() == 1){ // only can be at out-network-node
            index = 0; // where gene_EC_link[i][0] is num_ECs-1
        }

        int ec_index = gene_EC_link[i][index];
        // update gene-ec locations:
        // all genes are assigned to nodes within the network initially
        genes_on_EC[ec_index].push_back(i); // genes_on_EC.size() = num_ECs
        assigned_gene2ec[i] = ec_index;

    // calculate sequence based E:
        gene_energy[i].homology = gene_EC_homol[i][index];
        gene_energy[i].orthology = gene_EC_orth[i][index];
    }

    // After assignment, calculate Energy:

    //EC_co_occur_vec length = num_ECs-1
    for (int i=0;i<num_ECs-1;i++){
        int cnt_ecs = 0;
        double co_occur_val = 0.0;
        //EC_EC_link; // linked ECs for each EC including itself (from global network file) (ECnb)
	    for(int j=0;j<(EC_EC_link[i]).size();j++){ // for EC neighbors (0-th = itself) to the assigned EC
            int nb_ec_index = EC_EC_link[i][j];
            if(j>0){ // for EC neighbors (0-th = itself)
                if(genes_on_EC[nb_ec_index].size()>0){ //if any gene is on nb_ec_index, then add EC_EC_co_occur
                    co_occur_val += EC_EC_co_occur[i][j];
                    cnt_ecs ++;
                }
            }
	    }
	    EC_cnt_nbEC_vec[i] = cnt_ecs;
	    if(EC_EC_link[i].size()==0){
            cerr << "EC_EC_link[" << i << "]).size()==0" << endl;
        }
        else if(EC_EC_link[i].size() == 1){ // if EC_EC_link[i].size() == 1, cnt_ecs can not be >0
            EC_co_occur_vec[i] = median_ec_co_occur;
        }
        else{
            if(cnt_ecs>0){
                EC_co_occur_vec[i] = (co_occur_val/ (double) cnt_ecs);
            }
            else{
                EC_co_occur_vec[i] = 0.0;
            }
        }
    }

    for (int i=0;i<num_genes;i++) {
        // gene_EC_link is always >0 (1 if only outside node)

        int ec_index = assigned_gene2ec[i];

        if(ec_index < num_ECs-1){ // EC is within a network (not an outside node)
            //gene_energy[i].co_occur = 0.0;
            for(int j=0;j<(EC_EC_link[ec_index]).size();j++){ // for EC neighbors (0-th = itself) to the assigned EC
                //j==0: nb_ec_index == ec_index

                int nb_ec_index = EC_EC_link[ec_index][j];
                //
                // TODO: j: close current for loop here and open a new one?
                //
                for(int k=0;k<(genes_on_EC[nb_ec_index]).size();k++){// genes sitting on the nb_ec
                    int nb_gene_index = genes_on_EC[nb_ec_index][k];
                    if(i != nb_gene_index){ // for the neighbor genes other than itself
                        double e_ph = gene_gene_phylo_corr[i][nb_gene_index];
                        if( e_ph > gene_energy[i].phylo_corr ){
                            gene_energy[i].phylo_corr = e_ph;
                            gene_energy[i].max_phylo_corr_gene_index = nb_gene_index;
                        }

                        double e_gc = gene_gene_cluster_corr[i][nb_gene_index];
                        if( e_gc > gene_energy[i].gene_cluster ){
                            gene_energy[i].gene_cluster = e_gc;
                            gene_energy[i].max_gene_cluster_gene_index = nb_gene_index;
                        }

                        double e_ce = gene_gene_coexp_corr[i][nb_gene_index];
                        if( e_ce > gene_energy[i].co_exprsn ){
                            gene_energy[i].co_exprsn = e_ce;
                            gene_energy[i].max_co_exprsn_gene_index = nb_gene_index;
                        }
                    }
                }
            }
        }
        else{ // gene at the out-node
            // gene_energy[i]'s are all initialized in main.cpp
            gene_energy[i].phylo_corr = 0.0;//gene_energy[i].min_e_ph;
            gene_energy[i].max_phylo_corr_gene_index = -1;
            gene_energy[i].gene_cluster = 0.0; //gene_energy[i].min_e_gc;
            gene_energy[i].max_gene_cluster_gene_index = -1;
            gene_energy[i].co_exprsn = 0.0; //gene_energy[i].min_e_ce;
            gene_energy[i].max_co_exprsn_gene_index = -1;
            //gene_energy[i].co_occur = 0.0;
            gene_energy[i].no_EC_penalty = 1.0;
        }
    } // done all genes


}

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
        energy &delta_E, double &delta_EC_co_occur){
        // Update gene_energy of all genes to state of no gene at the org_ec_loc.

    delta_E.initialize(0.0,0.0,0.0); //(min_e_ph, min_e_gc, min_e_ce); // set scores to zero
    delta_EC_co_occur = 0.0;
    // context : E without the gene - E with the gene
    // co_occur : score for having this gene at org_ec_loc ( E with the gene-E without the gene)
    // homol, orth = gene's score

    remove_gene_on_EC(genes_on_EC[org_ec_loc], org_ec_loc, gene);

    bool change_co_occur = false;

    // calculate energy difference
    if(org_ec_loc < num_ECs-1){ // org EC was within a network (not an outside node)

        // First layer:
        // pull_out_delta_E's context = Delta_E_R = E_R - E_0
        if(genes_on_EC[org_ec_loc].size() == 0){ // genes_on_EC is reduced above (i.e., if it still have at least one gene on the org_ec_loc, the network does not lose its connection. co-occurence change only if genes_on_EC[org_ec_loc] == 0)
            change_co_occur = true;
            EC_cnt_nbEC_vec[org_ec_loc] = 0;
            delta_EC_co_occur += EC_co_occur_vec[org_ec_loc];
            EC_co_occur_vec[org_ec_loc] = 0.0; // the gene is not assigned to any EC-> no score
        }

        delta_E.homology = gene_energy[gene].homology; // score for having gene at org_ec_loc
        delta_E.orthology = gene_energy[gene].orthology;

        gene_energy[gene].homology = 0.0;
        gene_energy[gene].orthology = 0.0;

         // First layer:
        delta_E.phylo_corr += (min_e_ph - gene_energy[gene].phylo_corr);
        delta_E.gene_cluster += (min_e_gc - gene_energy[gene].gene_cluster);
        delta_E.co_exprsn += (min_e_ce - gene_energy[gene].co_exprsn);

        // Second layer:
        vector< vector<int> > updated_ph_index; // 1st = old gene, 2nd = updated gene
        vector< double > updated_ph; // new phylo corr val with the updated gene

        vector< vector<int> > updated_gc_index;
        vector< double > updated_gc;

        vector< vector<int> > updated_ce_index;
        vector< double > updated_ce;

        // As the gene is removed, the neighbor loses its neighbor.
        // for each gene sitting on the neighboring EC, measure gene-gene context
        // the pulled out gene was already removed from genes_on_EC[index_nb_ec] above.

        // for each neighbor EC position:
        for(int i=0;i<EC_EC_link[org_ec_loc].size();i++){ //i==0: index_nb_ec==org_ec_loc

            int index_nb_ec = EC_EC_link[org_ec_loc][i]; // a neighbor (can be itself)

            double new_co_occur;
            // EC-EC co_occur: (Second layer)
            if(i>0 && change_co_occur == true){

                if(genes_on_EC[index_nb_ec].size()>0){ // it loses one nb

                    EC_cnt_nbEC_vec[index_nb_ec] --; //reduce the number of neighboring ECs that are occupied.

                    if(EC_cnt_nbEC_vec[index_nb_ec] > 0){ // even after pulling out a nb gene, it still has nbs.
                        double org_co_occur = EC_co_occur_vec[index_nb_ec] * (EC_cnt_nbEC_vec[index_nb_ec]+1.0); // +1 for the pulled out gene

                        new_co_occur = org_co_occur - EC_EC_co_occur[org_ec_loc][i];
                        new_co_occur /= (double) EC_cnt_nbEC_vec[index_nb_ec]; // updated co-occur for gene_on_nbEC_index
                    }
                    else{ // that pulled out gene was the only nb gene for gene_on_nbEC_index
                        new_co_occur = 0.0; //median_ec_co_occur;
                    }

                    delta_EC_co_occur += (EC_co_occur_vec[index_nb_ec] - new_co_occur);
                    EC_co_occur_vec[index_nb_ec] = new_co_occur;
                }
            }

            for(int j=0;j<genes_on_EC[index_nb_ec].size();j++){

                int gene_on_nbEC_index = genes_on_EC[index_nb_ec][j]; // one of the neighbor gene

                // if gene-gene_on_nbEC_index is the highest score
                if( gene_on_nbEC_index != gene && gene == gene_energy[gene_on_nbEC_index].max_phylo_corr_gene_index ){

                    delta_E.phylo_corr -= gene_energy[gene_on_nbEC_index].phylo_corr;

                    // search this neighbor gene's neighbors for the second best gene pair
                    // After take out the gene, a nb gene may not have any nb genes
                    double second_max_ph = min_e_ph; // -1.0;
                    int second_max_ph_index = -1; // means isolated gene, i.e., no nb

                    //for(int u=1;u<EC_EC_link[index_nb_ec].size();u++) // nb EC's nb ECs
                    for(int u=0;u<EC_EC_link[index_nb_ec].size();u++){ // nb EC's nb ECs

                        int index_nb_nb_ec = EC_EC_link[index_nb_ec][u];

                        for(int v=0;v<genes_on_EC[index_nb_nb_ec].size();v++){ // genes on the nb EC's nb ECs
                            // After take out the gene, a nb gene may not have any nb genes -> this loop is not run
                            int nb_of_nb_gene_index = genes_on_EC[index_nb_nb_ec][v];

                            if( gene != nb_of_nb_gene_index ){
                                if(second_max_ph < gene_gene_phylo_corr[gene_on_nbEC_index][nb_of_nb_gene_index]){
                                    second_max_ph = gene_gene_phylo_corr[gene_on_nbEC_index][nb_of_nb_gene_index];
                                    second_max_ph_index = nb_of_nb_gene_index;
                                }
                            }
                        }
                    }

                    delta_E.phylo_corr += second_max_ph;

                    // updated_ph_index =  1st = old gene, 2nd = updated neighbor gene
                    // updated_ph = new phylo corr val with the updated gene
                    updated_ph.push_back(second_max_ph);
                    int tmp_arr[] = {gene_on_nbEC_index, second_max_ph_index};
                    vector<int> tmp(tmp_arr, tmp_arr+sizeof(tmp_arr)/sizeof(tmp_arr[0]));
                    updated_ph_index.push_back(tmp);

                } // phylo_corr

                // if gene-gene_on_nbEC_index is the highest score
                if( gene_on_nbEC_index != gene && gene == gene_energy[gene_on_nbEC_index].max_gene_cluster_gene_index ){

                    delta_E.gene_cluster -= gene_energy[gene_on_nbEC_index].gene_cluster;

                    // search this neighbor gene's neighbors for the second best gene pair
                    double second_max_gc = min_e_gc; //-1.0;
                    int second_max_gc_index = -1;

                    //for(int u=1;u<EC_EC_link[index_nb_ec].size();u++) // nb EC's nb ECs
                    for(int u=0;u<EC_EC_link[index_nb_ec].size();u++){ // nb EC's nb ECs

                        int index_nb_nb_ec = EC_EC_link[index_nb_ec][u];

                        for(int v=0;v<genes_on_EC[index_nb_nb_ec].size();v++){ // genes on the nb EC's nb ECs

                            int nb_of_nb_gene_index = genes_on_EC[index_nb_nb_ec][v];

                            if( gene != nb_of_nb_gene_index ){
                                if(second_max_gc < gene_gene_cluster_corr[gene_on_nbEC_index][nb_of_nb_gene_index]){
                                    second_max_gc = gene_gene_cluster_corr[gene_on_nbEC_index][nb_of_nb_gene_index];
                                    second_max_gc_index = nb_of_nb_gene_index;
                                }
                            }
                        }
                    }

                    delta_E.gene_cluster += second_max_gc;

                    updated_gc.push_back(second_max_gc);
                    int tmp_arr[] = {gene_on_nbEC_index, second_max_gc_index};
                    vector<int> tmp(tmp_arr, tmp_arr+sizeof(tmp_arr)/sizeof(tmp_arr[0]));
                    updated_gc_index.push_back(tmp);

                } // gene cluster

                // if gene-gene_on_nbEC_index is the highest score
                if( gene_on_nbEC_index != gene && gene == gene_energy[gene_on_nbEC_index].max_co_exprsn_gene_index ){

                    delta_E.co_exprsn -= gene_energy[gene_on_nbEC_index].co_exprsn;

                    // search this neighbor gene's neighbors for the second best gene pair
                    double second_max_ce = min_e_ce; //-1.0;
                    int second_max_ce_index = -1;

                    //for(int u=1;u<EC_EC_link[index_nb_ec].size();u++) // nb EC's nb ECs
                    for(int u=0;u<EC_EC_link[index_nb_ec].size();u++){ // nb EC's nb ECs

                        int index_nb_nb_ec = EC_EC_link[index_nb_ec][u];

                        for(int v=0;v<genes_on_EC[index_nb_nb_ec].size();v++){ // genes on the nb EC's nb ECs

                            int nb_of_nb_gene_index = genes_on_EC[index_nb_nb_ec][v];

                            if( gene != nb_of_nb_gene_index ){
                                if(second_max_ce < gene_gene_coexp_corr[gene_on_nbEC_index][nb_of_nb_gene_index]){
                                    second_max_ce = gene_gene_coexp_corr[gene_on_nbEC_index][nb_of_nb_gene_index];
                                    second_max_ce_index = nb_of_nb_gene_index;
                                }
                            }
                        }
                    }

                    delta_E.co_exprsn += second_max_ce;

                    updated_ce.push_back(second_max_ce);
                    int tmp_arr[] = {gene_on_nbEC_index, second_max_ce_index};
                    vector<int> tmp(tmp_arr, tmp_arr+sizeof(tmp_arr)/sizeof(tmp_arr[0]));
                    updated_ce_index.push_back(tmp);

                } // co_exprsn

            } // the gene's neighboring genes

        } // the gene's neighboring EC locations


        // EC-EC co-occur Second layer
        if(change_co_occur == true){
            if(EC_EC_link[org_ec_loc].size() == 1){
                delta_EC_co_occur = median_ec_co_occur;
            }
        }
        else{
            delta_EC_co_occur = 0.0;
        }

        // update the gene's context to second best become gene's context score
        for(int i=0;i<updated_ph.size();i++){
            gene_energy[updated_ph_index[i][0]].phylo_corr = updated_ph[i];
            gene_energy[updated_ph_index[i][0]].max_phylo_corr_gene_index = updated_ph_index[i][1];
        }
        for(int i=0;i<updated_gc.size();i++){
            gene_energy[updated_gc_index[i][0]].gene_cluster = updated_gc[i];
            gene_energy[updated_gc_index[i][0]].max_gene_cluster_gene_index = updated_gc_index[i][1];
        }
        for(int i=0;i<updated_ce.size();i++){
            gene_energy[updated_ce_index[i][0]].co_exprsn = updated_ce[i];
            gene_energy[updated_ce_index[i][0]].max_co_exprsn_gene_index = updated_ce_index[i][1];
        }

    }
    else{ // originally outside of the EC network
        delta_E.no_EC_penalty = -1.0;
        // delta_E.initialize(0.0,0.0,0.0);
        // pull_out_delta_E's context = Delta_E_R = E_R - E_0
        // E_R = min_e_X, E_0 = 0
        delta_E.phylo_corr = min_e_ph;
        delta_E.gene_cluster = min_e_gc;
        delta_E.co_exprsn = min_e_ce;
    }

    // no matter pulling out in-node or out-node, at E_R state:
    gene_energy[gene].phylo_corr = min_e_ph;
    gene_energy[gene].max_phylo_corr_gene_index = -1;
    gene_energy[gene].gene_cluster = min_e_gc;
    gene_energy[gene].max_gene_cluster_gene_index = -1;
    gene_energy[gene].co_exprsn = min_e_ce;
    gene_energy[gene].max_co_exprsn_gene_index = -1;

}

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
        energy &delta_E, double &delta_EC_co_occur){

    delta_E.initialize(0.0,0.0,0.0); // (min_e_ph, min_e_gc, min_e_ce);
    delta_EC_co_occur = 0.0;

    updated_ph_index.clear();
    updated_ph.clear();
    updated_gc_index.clear();
    updated_gc.clear();
    updated_ce_index.clear();
    updated_ce.clear();
    updated_co_occur_EC_vec.clear();
    updated_co_occur.clear();
    updated_co_occur_nb_cnt.clear();

    // new_ec_loc can be out-node
    genes_on_EC[new_ec_loc].push_back(gene); // undo is required at the end

    int index;
    find_index_of_id(gene_EC_link[gene], new_ec_loc, index); // index of new_ec_loc in gene_EC_link[gene]
    if(index == -1){ cerr << "no " << new_ec_loc << " in gene_EC_link[" << gene << "]" << endl;
    }

    delta_E.homology = gene_EC_homol[gene][index]; // out-node column = 0.0
    delta_E.orthology = gene_EC_orth[gene][index]; // out-node column = 0

    bool change_co_occur = false;
    double new_ec_co_occur = 0.0;
    int new_ec_nbs = 0; // will be added to EC_cnt_nbEC_vec[new_ec_loc] outside the for loop

    // std::vector<>.push_back() copies element.
    //int tmp_arr[] = {gene_on_nbEC_index, gene};
    int tmp_arr[] = {0, gene};
    vector<int> tmp(tmp_arr, tmp_arr+sizeof(tmp_arr)/sizeof(tmp_arr[0]));
    vector<int> tmp2(tmp_arr, tmp_arr+sizeof(tmp_arr)/sizeof(tmp_arr[0]));

    if(new_ec_loc < num_ECs-1 ){ // the new location is within the network

        double max_ph = min_e_ph;
        int max_ph_index = -1; // the gene index the new_ec_loc may be an isolate location. i.e., no neighboring genes at all
        double max_gc = min_e_gc;
        int max_gc_index = -1; // the gene index
        double max_ce = min_e_ce;
        int max_ce_index = -1; // the gene index

        //No matter genes were already there, the gene will get its own co_occur score. Only the nb genes will not update their co_occurs.
        updated_co_occur_EC_vec.push_back(new_ec_loc);
        updated_co_occur.push_back(0.0); //median_ec_co_occur); // give the true value at the end. If it doesn't have nbs, its val= median_ec_co_occur
        updated_co_occur_nb_cnt.push_back(0);

        if(genes_on_EC[new_ec_loc].size() == 1){ // the new gene is added above. i.e., no genes were at new_ec_loc
            change_co_occur = true;
        }

        // for each neighbor EC position:
        for(int i=0;i<EC_EC_link[new_ec_loc].size();i++){ // i==0: index_nb_ec==new_ec_loc

            int index_nb_ec = EC_EC_link[new_ec_loc][i];
            double delta_co_occur;
            double new_co_occur;

            if(i>0){
                if(change_co_occur == true){

                    if(genes_on_EC[index_nb_ec].size() > 0){
                        // First layer EC co-occur
                        new_ec_co_occur += EC_EC_co_occur[new_ec_loc][i];
                        new_ec_nbs ++;

                        // Second layer EC co-occur (nb gene acquires one EC-EC with new_ec_loc)
                        int nb_ecs = EC_cnt_nbEC_vec[index_nb_ec] + 1; // after adding the gene, the number of nb genes of gene_on_nbEC_index
                        if(nb_ecs > 1){
                            double org_co_occur = EC_co_occur_vec[index_nb_ec] * (nb_ecs-1); // -1 for the adding a gene

                            new_co_occur = org_co_occur + EC_EC_co_occur[new_ec_loc][i];
                            new_co_occur /= (double) nb_ecs; // updated co-occur for gene_on_nbEC_index
                            delta_co_occur = new_co_occur - EC_co_occur_vec[index_nb_ec];
                        }
                        else{ // the nb gene gene_on_nbEC_index didn't have nbs.
                            new_co_occur = EC_EC_co_occur[new_ec_loc][i];
                            delta_co_occur = new_co_occur; // new_co_occur -median_ec_co_occur;

                            if(nb_ecs == 0){
                                cerr << "No EC-EC connection even after adding a gene" << endl;
                            }
                        }
                        delta_EC_co_occur += delta_co_occur;

                        updated_co_occur_EC_vec.push_back(index_nb_ec);
                        updated_co_occur.push_back(new_co_occur);
                        updated_co_occur_nb_cnt.push_back(nb_ecs);
                    }
                }
            }

            for(int j=0;j<genes_on_EC[index_nb_ec].size();j++){

                int gene_on_nbEC_index = genes_on_EC[index_nb_ec][j]; // one of the neighbor gene
                tmp[0] = gene_on_nbEC_index;

                if(gene_on_nbEC_index != gene){
                    // First layer:
                    if(max_ph<gene_gene_phylo_corr[gene][gene_on_nbEC_index]){
                        max_ph = gene_gene_phylo_corr[gene][gene_on_nbEC_index];
                        max_ph_index = gene_on_nbEC_index;
                    }
                    if(max_gc<gene_gene_cluster_corr[gene][gene_on_nbEC_index]){
                        max_gc = gene_gene_cluster_corr[gene][gene_on_nbEC_index];
                        max_gc_index = gene_on_nbEC_index;
                    }
                    if(max_ce<gene_gene_coexp_corr[gene][gene_on_nbEC_index]){
                        max_ce = gene_gene_coexp_corr[gene][gene_on_nbEC_index];
                        max_ce_index = gene_on_nbEC_index;
                    }

                    // Second layer:
                    // if cotext E with the newly added gene is higher than what it had before:
                    if(gene_energy[gene_on_nbEC_index].phylo_corr < gene_gene_phylo_corr[gene_on_nbEC_index][gene]){
                        updated_ph.push_back(gene_gene_phylo_corr[gene_on_nbEC_index][gene]);
                        updated_ph_index.push_back(tmp);
                    }
                    if(gene_energy[gene_on_nbEC_index].gene_cluster < gene_gene_cluster_corr[gene_on_nbEC_index][gene]){
                        updated_gc.push_back(gene_gene_cluster_corr[gene_on_nbEC_index][gene]);
                        updated_gc_index.push_back(tmp);
                    }
                    if(gene_energy[gene_on_nbEC_index].co_exprsn < gene_gene_coexp_corr[gene_on_nbEC_index][gene]){
                        updated_ce.push_back(gene_gene_coexp_corr[gene_on_nbEC_index][gene]);
                        updated_ce_index.push_back(tmp);
                    }
                }
            } // for each gene on neighbor EC

        } // neighbor EC

        if(change_co_occur == true){
            if(EC_EC_link[new_ec_loc].size()==1){  // if EC_EC_link[ec_index].size() == 1, new_ec_nbs can not be >0
                delta_EC_co_occur = median_ec_co_occur;
                updated_co_occur[0] = median_ec_co_occur;
            }
            else{
                if(new_ec_nbs > 0){ // the added gene has neighboring genes
                    updated_co_occur[0] = new_ec_co_occur/ (double)new_ec_nbs;
                    delta_EC_co_occur += updated_co_occur[0];
                    updated_co_occur_nb_cnt[0] = new_ec_nbs;
                }
            }
        }
        else{
            updated_co_occur[0] = EC_co_occur_vec[new_ec_loc];
            updated_co_occur_nb_cnt[0] = EC_cnt_nbEC_vec[new_ec_loc];
        }

        updated_ph.push_back(max_ph);
        //int tmp_arr[] = {gene, max_ph_index};
        //vector<int> tmp(tmp_arr, tmp_arr+sizeof(tmp_arr)/sizeof(tmp_arr[0]));
        tmp[0] = gene;
        tmp[1] = max_ph_index;
        updated_ph_index.push_back(tmp);

        updated_gc.push_back(max_gc);
        //vector<int> tmp {gene, max_gc_index};
        tmp[1] = max_gc_index;
        updated_gc_index.push_back(tmp);

        updated_ce.push_back(max_ce);
        //vector<int> tmp {gene, max_ce_index};
        tmp[1] = max_ce_index;
        updated_ce_index.push_back(tmp);

        // if new_ec_loc is isolated one, add the minimum.
        delta_E.phylo_corr += (max_ph-min_e_ph);
        delta_E.gene_cluster += (max_gc-min_e_gc);
        delta_E.co_exprsn += (max_ce-min_e_ce);

        for(int i=0;i<updated_ph.size()-1;i++){ // i=0: the new gene its own E: already added to delta_E above
            delta_E.phylo_corr -= gene_energy[updated_ph_index[i][0]].phylo_corr;
            delta_E.phylo_corr += updated_ph[i];
        }
        for(int i=0;i<updated_gc.size()-1;i++){ // i=0: the new gene its own E: already added to delta_E above
            delta_E.gene_cluster -= gene_energy[updated_gc_index[i][0]].gene_cluster;
            delta_E.gene_cluster += updated_gc[i];
        }
        for(int i=0;i<updated_ce.size()-1;i++){ // i=0: the new gene its own E: already added to delta_E above
            delta_E.co_exprsn -= gene_energy[updated_ce_index[i][0]].co_exprsn;
            delta_E.co_exprsn += updated_ce[i];
        }

    }
    else{ // put back to out-node:
        delta_E.no_EC_penalty = 1.0;
        // delta_E_vec = E_A -E_R
        // E_R = min_e_X, E_A = 0 (no gene-gene context score for outside node)
        delta_E.phylo_corr = -min_e_ph;
        delta_E.gene_cluster = -min_e_gc;
        delta_E.co_exprsn = -min_e_ce;

        //int tmp_arr[] = {gene, -1};

        updated_ph.push_back(0.0);
        //vector<int> tmp(tmp_arr, tmp_arr+sizeof(tmp_arr)/sizeof(tmp_arr[0]));
        tmp[0] = gene;
        tmp[1] = -1;
        updated_ph_index.push_back(tmp);

        updated_gc.push_back(0.0);
        //int tmp_arr[] = {gene, -1};
        //vector<int> tmp(tmp_arr, tmp_arr+sizeof(tmp_arr)/sizeof(tmp_arr[0]));
        updated_gc_index.push_back(tmp);

        updated_ce.push_back(0.0);
        //int tmp_arr[] = {gene, -1};
        //vector<int> tmp(tmp_arr, tmp_arr+sizeof(tmp_arr)/sizeof(tmp_arr[0]));
        updated_ce_index.push_back(tmp);
    }

    // undo putting this gene to the candidate EC loc
    remove_gene_on_EC(genes_on_EC[new_ec_loc], new_ec_loc, gene);

}

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
        ){

    // switch a gene's EC location to all possible locations and re-calculate the energy differences -> choose new_ec_loc_index

    int org_ec_loc = assigned_gene2ec[rand_gene_index]; // original loc

    energy pull_out_delta_E;
    double pull_out_delta_co_occur;

    vector<energy> delta_E_vec(gene_EC_link[rand_gene_index].size());
    for(int i=0;i<gene_EC_link[rand_gene_index].size();++i){
        delta_E_vec[i].initialize(min_e_ph, min_e_gc, min_e_ce);
    }
    vector<double> delta_co_occur_vec(gene_EC_link[rand_gene_index].size());

    // update gene_energy and genes_on_EC (remove a gene-EC pair from genes_on_EC)

    // As the gene is removed, the energy of the network changes:
    vector< vector<int> > updated_ph_index;
    vector< double > updated_ph;
    vector< vector<int> > updated_gc_index;
    vector< double > updated_gc;
    vector< vector<int> > updated_ce_index;
    vector< double > updated_ce;
    vector<int> updated_co_occur_EC_vec;
    vector< double > updated_co_occur;
    vector< int > updated_co_occur_nb_cnt;

    int new_ec_loc;

    // Update gene_energy of all genes to state of no gene at the org_ec_loc.
    // out-node vs. in-node / out-node: pull_out_delta_E: all zeros except penalty term
    // Every gene has to be updated to E_R state (update genes_on_EC and gene_energy of every gene if necessary)

    // pull_out_delta_E:

        // context, no_EC_p : E without the gene - E with the gene (== Delta_E_R)
        // co_occur : score for having this gene at org_ec_loc ( E with the gene-E without the gene) ( == -Delta_E_R )
        // homol, orth = gene's score ( == E_0  == -Delta_E_R since E_R==0)
    pull_out_gene(
        min_e_ph, min_e_gc, min_e_ce,
        num_ECs,
        rand_gene_index, org_ec_loc, gene_energy,
        EC_EC_link, EC_EC_co_occur, median_ec_co_occur, genes_on_EC,
        gene_gene_phylo_corr, gene_gene_cluster_corr, gene_gene_coexp_corr,
        EC_co_occur_vec, EC_cnt_nbEC_vec,
        pull_out_delta_E, pull_out_delta_co_occur); // gene_energy and genes_on_EC are updated

    // for all possible new EC location (including out of the network node):
    for(int i=0;i<gene_EC_link[rand_gene_index].size();++i){

        new_ec_loc = gene_EC_link[rand_gene_index][i];

        if(org_ec_loc == new_ec_loc){ //delta_E_vec[i] <= -Delta_E_R (== score for having it at this loc)
            delta_E_vec[i].homology = pull_out_delta_E.homology;
            delta_E_vec[i].orthology = pull_out_delta_E.orthology;

            delta_E_vec[i].phylo_corr = -1.0*pull_out_delta_E.phylo_corr; // -1*pull_out_delta_E == score for having gene at rand_gene_index relative to the state without having gene within the network (== -Delta_E_R)
            delta_E_vec[i].gene_cluster = -1.0*pull_out_delta_E.gene_cluster; // -Delta_E_R
            delta_E_vec[i].co_exprsn = -1.0*pull_out_delta_E.co_exprsn;
            //delta_E_vec[i].co_occur = pull_out_delta_E.co_occur; // -Delta_E_R
            delta_co_occur_vec[i] = pull_out_delta_co_occur; // -Delta_E_R
            delta_E_vec[i].no_EC_penalty = -1.0*pull_out_delta_E.no_EC_penalty;
        }
        else{
            // do not update assigned_gene2ec, genes_on_EC and gene_energy. Do it after select a new EC location
            // delta_E_vec[i] = delta E with a gene at a new loc relative to the state without the gene assigned in the network. (== E_A - E_R == Delta_E_A)
            put_back_new_loc(num_ECs, min_e_ph, min_e_gc, min_e_ce,
                rand_gene_index, new_ec_loc, gene_energy,
                EC_EC_link, gene_EC_link, genes_on_EC,
                gene_EC_homol, gene_EC_orth,
                gene_gene_phylo_corr, gene_gene_cluster_corr, gene_gene_coexp_corr,
                EC_EC_co_occur, median_ec_co_occur,
                updated_ph_index, updated_ph, updated_gc_index, updated_gc, updated_ce_index, updated_ce,
                updated_co_occur_EC_vec, updated_co_occur, updated_co_occur_nb_cnt,
                EC_co_occur_vec, EC_cnt_nbEC_vec,
                delta_E_vec[i], delta_co_occur_vec[i]);
        }
    }

    select_new_loc(
        coeff_hm, coeff_ort, coeff_ph, coeff_gc, coeff_ce, p_ec_co_occur, p_no_EC,
        delta_E_vec, delta_co_occur_vec, new_ec_loc_index);

    // update the gene's new location to new_ec_loc:
    new_ec_loc = gene_EC_link[rand_gene_index][new_ec_loc_index];

    // delta_E_vec[new_ec_loc_index].co_occur == Delta_E_A
    // update not only rand_gene_index but also neighboring genes using updated_....
    energy tmp_delta_E;
    double tmp_delta_co_occur;
    put_back_new_loc(num_ECs, min_e_ph, min_e_gc, min_e_ce,
        rand_gene_index, new_ec_loc, gene_energy,
        EC_EC_link, gene_EC_link, genes_on_EC,
        gene_EC_homol, gene_EC_orth,
        gene_gene_phylo_corr, gene_gene_cluster_corr, gene_gene_coexp_corr,
        EC_EC_co_occur, median_ec_co_occur,
        updated_ph_index, updated_ph, updated_gc_index, updated_gc, updated_ce_index, updated_ce,
        updated_co_occur_EC_vec, updated_co_occur, updated_co_occur_nb_cnt,
        EC_co_occur_vec, EC_cnt_nbEC_vec,
        tmp_delta_E, tmp_delta_co_occur); // redo put_back_new_loc just to get updated_xxx at the selected new_ec_loc

    assigned_gene2ec[rand_gene_index] = new_ec_loc;

    genes_on_EC[new_ec_loc].push_back(rand_gene_index);

    gene_energy[rand_gene_index].homology = tmp_delta_E.homology;
    gene_energy[rand_gene_index].orthology = tmp_delta_E.orthology;

    int old_at; int new_at; double new_val;
    for(int i=0;i<updated_ph.size();++i){
        old_at = updated_ph_index[i][0];
        new_at = updated_ph_index[i][1];
        new_val = updated_ph[i];

        gene_energy[old_at].phylo_corr = new_val;
        gene_energy[old_at].max_phylo_corr_gene_index = new_at;
    }
    for(int i=0;i<updated_gc.size();++i){
        old_at = updated_gc_index[i][0];
        new_at = updated_gc_index[i][1];
        new_val = updated_gc[i];

        gene_energy[old_at].gene_cluster = new_val;
        gene_energy[old_at].max_gene_cluster_gene_index = new_at;
    }
    for(int i=0;i<updated_ce.size();++i){
        old_at = updated_ce_index[i][0];
        new_at = updated_ce_index[i][1];
        new_val = updated_ce[i];

        gene_energy[old_at].co_exprsn = new_val;
        gene_energy[old_at].max_co_exprsn_gene_index = new_at;
    }
    for(int i=0;i<updated_co_occur_EC_vec.size();++i){ // include all ECs including itself
        EC_co_occur_vec[updated_co_occur_EC_vec[i]] = updated_co_occur[i];
        EC_cnt_nbEC_vec[updated_co_occur_EC_vec[i]] = updated_co_occur_nb_cnt[i];
    }
}

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
        int &new_ec_loc_index){
    // delta_E_vec:
        // context : score for having this gene at the new loc ( == Delta_E_A)
        // co_occur, no_EC_p: gene's E ( == Delta_E_A )
        // homol, orth = gene's score ( == E_A  == Delta_E_A since E_R==0)

    // Sum all score changes -> total score of the network -> fitness
    int num_candi_ec = delta_E_vec.size();
    vector<double> total_score_vec(num_candi_ec);

    for(int i=0;i<num_candi_ec;++i){
        total_score_vec[i] = 0.0;
        total_score_vec[i] += (coeff_hm * delta_E_vec[i].homology);
        total_score_vec[i] += (coeff_ort * delta_E_vec[i].orthology);
        total_score_vec[i] += (coeff_ph * delta_E_vec[i].phylo_corr);
        total_score_vec[i] += (coeff_gc * delta_E_vec[i].gene_cluster);
        total_score_vec[i] += (coeff_ce * delta_E_vec[i].co_exprsn);
        //total_score_vec[i] += (p_ec_co_occur * delta_E_vec[i].co_occur);
        total_score_vec[i] += (p_ec_co_occur * delta_co_occur_vec[i]);
        total_score_vec[i] += (p_no_EC*delta_E_vec[i].no_EC_penalty); // p_no_EC (-)
        if(!std::isfinite(total_score_vec[i])){
            cerr << "843: total_score_vec NAN" << endl;
        }
    }

    locate_gene(total_score_vec, new_ec_loc_index);

}

#include <iomanip> // setprecision
void locate_gene(const vector<double> &vec_prob, int &index){

    index = -1;
    int num_loc = vec_prob.size();
    vector<double> normalized_vec_prob = vec_prob;

    double max_prob = *(std::max_element(normalized_vec_prob.begin(), normalized_vec_prob.end()));
    // replaced
    //double max_prob = normalized_vec_prob[0];
    //for(int i=1;i<num_loc;++i){
    //    if(normalized_vec_prob[i]>max_prob){
    //        max_prob = normalized_vec_prob[i];
    //    }
    //}

    double sum_prob = 0.0;
    for(int i=0;i<num_loc;++i){
        normalized_vec_prob[i] = pow(10,normalized_vec_prob[i]-max_prob); // fitness scores are shifted so that the maximum becomes 1. otherwise, a fitness score can be an INF value
        //normalized_vec_prob[i] = exp(normalized_vec_prob[i]-max_prob); // fitness scores are shifted so that the maximum becomes 1. otherwise, a fitness score can be an INF value
        sum_prob += normalized_vec_prob[i];
    }

    double posi; // (0,1]
    do{
        posi = uniform_rand();
    }while(posi>1 || posi==0.0);


    for(vector<double>::iterator itr=normalized_vec_prob.begin(), itr_end=normalized_vec_prob.end(); itr != itr_end; ++itr){
        *itr = (*itr) / sum_prob;
    }
    //for(int i=0;i<num_loc;++i){
    //    normalized_vec_prob[i] = (normalized_vec_prob[i])/sum_prob;
    //}

    double tmp_sum = 0.0;
    for(int i=0, end_i=num_loc-1; i<end_i; ++i){
        tmp_sum += normalized_vec_prob[i];
        if( posi <= tmp_sum ){
            index = i;
            break;
        }
    }
    if(index < 0){
        index = num_loc - 1;
    }
    //double tmp_sum = 0.0;
    //double tmp_sum_next = 0.0;
    //for(int i=0;i<num_loc;++i){ // where does posi fall in? posi: (0,sum_prob]
    //    if(i<num_loc-1){
    //        tmp_sum_next = tmp_sum + normalized_vec_prob[i];
    //    }
    //    else{
    //        tmp_sum_next = 1;
    //    }

    //    if( (tmp_sum < posi) && (posi <= tmp_sum_next) ){
    //        index = i;
    //        break;
    //    }
    //    else{
    //        tmp_sum = tmp_sum_next;
    //    }
    //}
    //
    //if(index < 0){
    //    index = num_loc-1;
    //    cerr << "locate_gene: wrong new_ec_loc_index was selected." << endl;
    //}

}

void remove_gene_on_EC(vector<int> &genes_on_EC_vec, const int &ec_loc, const int &gene){
    for(vector<int>::iterator itr=genes_on_EC_vec.begin(), itr_end=genes_on_EC_vec.end(); itr != itr_end; ++itr){
        if(*itr == gene){
            genes_on_EC_vec.erase(itr);
            break;
        }
    }
    //for(int i=0;i<genes_on_EC_vec.size();++i){
    //    if(genes_on_EC_vec[i] == gene){
    //        genes_on_EC_vec.erase(genes_on_EC_vec.begin()+i);
    //        break;
    //    }
    //}
}

void find_index_of_id(const vector<int> &id_vec, const int &id, int &index){
    vector<int>::const_iterator itr = std::find(id_vec.begin(), id_vec.end(), id);
    if(itr == id_vec.end()){
        index = -1;
    }
    else{
        index = std::distance(id_vec.begin(), itr);
    }
    //
    //bool check = false; // True if id exists in the id_vec
    //for(int i=0;i<id_vec.size();++i){
    //    if(id_vec[i] == id){
    //        index = i;
    //        check = true;
    //        break;
    //    }
    //}
    //if(check == false){ index = -1; }

}

/////////////

void unique_EC_occupancy(const vector<int> &assigned_gene2ec, Vector<int> &count_ec_on_init){
    vector<int> tmp_assigned_gene2ec = assigned_gene2ec;

    std::sort(tmp_assigned_gene2ec.begin(), tmp_assigned_gene2ec.end());

    vector<int>::iterator it = std::unique (tmp_assigned_gene2ec.begin(), tmp_assigned_gene2ec.end());

    vector<int>::iterator it_p;
    for(it_p=tmp_assigned_gene2ec.begin(); it_p != it; ++it_p){
        count_ec_on_init[*it_p] ++;
    }
}


