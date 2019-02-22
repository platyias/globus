//
//  main.cpp
//
//  Created by Jeewoen Shin on 4/14/17.
//

#include "parameters.h"
#include "data_table_io.h"
#include "random.h"
#include "fitness_calculate.h"

#include <cmath>
#include <cstdlib> // strtod
#include <ctime>
//#include <chrono>

#include <fstream>
#include <string>
#include <algorithm>


int main(int args, char** argv){
//int main(){

    //std::chrono::time_point<std::chrono::system_clock> start, end;

	srand(static_cast<unsigned int>( time(0) + atoi(argv[1]) )); // number of independent runs with different initial gene-ec assignments.
	//srand(static_cast<unsigned int>(time(0)));

	/*
    if(args < arg_num)
	{
		std::cerr<<"usage: "<<argv[0]<<" $SGE_TASK_ID"<<endl;
		return -1;
	}
    */
    ifstream fin;

// input parameters
    std::string in_txt_filename = argv[2];

    string genome_name; // arg 2 = $1
    string genome_dir; // arg 3 = .../results/$1/
    string global_net_file; // arg 4 = .../database/GRCENEECnet.txt

    double coeff_hm; // coeff for homology (sequence identity score) (+)
    double coeff_ph; // coeff for phylogenetic corr (+)
    double coeff_gc; // coeff for gene clustering (+)
    double coeff_ce; // coeff for co-expression (+)
    double coeff_ort; //coeff for othology (+)
    double p_ec_co_occur; // coeff for EC-EC co-occur (+)
    double p_no_EC; // coeff for not assigning EC to a gene (penalty) (-)
    double median_ec_co_occur; // coeff of median ec-ec occurence (coff_median_ec_occur = median of all BN ec-ec occurence.

    fin.open( in_txt_filename.c_str() );
	if( ! fin.is_open() )
	{
		cerr<<"Failed to open "<< in_txt_filename <<" file."<<endl;
		return -1;
	}
    read_input_param(fin, genome_name, genome_dir, global_net_file,
                     coeff_hm, coeff_ph, coeff_gc, coeff_ce, coeff_ort, p_ec_co_occur, p_no_EC, median_ec_co_occur);

    fin.close();

// Output file preparation here //

    std::string postfix1(".GLOBUS.Gene_Prob.txt");
    std::string postfix2(".GLOBUS.EC_Prob.txt");
    std::string swissprot = argv[3];

    std::string out_filename1;
    out_filename1.append(genome_dir);
    out_filename1.append(genome_name);
    out_filename1.append(".r");
    out_filename1.append(argv[1]);
    if (swissprot.compare("1") == 0) { // if argv[3] == 1
        out_filename1.append(".old");
    }
    else{
        out_filename1.append(".new");
    }
    out_filename1.append(postfix1);

    std::string out_filename2;
    out_filename2.append(genome_dir);
    out_filename2.append(genome_name);
    out_filename2.append(".r");
    out_filename2.append(argv[1]);
    if (swissprot.compare("1") == 0) { // if argv[3] == 1
        out_filename2.append(".old");
    }
    else{
        out_filename2.append(".new");
    }
    out_filename2.append(postfix2);
    //cerr << "myfile2" << out_filename2 << endl;

    std::string out_filename3;
    out_filename3.append(genome_dir);
    out_filename3.append(genome_name);
    out_filename3.append(".all");
    out_filename3.append(".r");
    out_filename3.append(argv[1]);
    if (swissprot.compare("1") == 0) { // if argv[3] == 1
        out_filename3.append(".old");
    }
    else{
        out_filename3.append(".new");
    }
    out_filename3.append(postfix1);
    //cerr << "myfile3" << out_filename3 << endl;

    std::ofstream myfile1; // acc_count_gene_on_ec_act
	std::ofstream myfile2; // count_ec_on_act
    std::ofstream myfile3; // A row: A iteragtion / At a row, num_genes pairs of gene-EC that were assigened at the iteration. (iterations include both initial and actual runs.)

    myfile1.open(out_filename1.c_str());
    myfile2.open(out_filename2.c_str());
    myfile3.open(out_filename3.c_str());
    
    std::string postfix3 = "";
    std::string postfix4 = "";
    std::string postfix5 = "";
    std::string postfix6 = "";
    std::string postfix7 = "";
// Input file preparation here //
    if (swissprot.compare("1") == 0) { // if argv[3] == 1
        postfix3 = ".gene_sw1";
        postfix4 = ".ort_sw1";
        postfix5 = ".pc.Z_sw1";
        postfix6 = ".gc.Z_sw1";
        postfix7 = ".pc.Z_sw1";
    }
    else{
        postfix3 = ".gene";
        postfix4 = ".ort";
        postfix5 = ".pc.Z";
        postfix6 = ".gc.Z";
        postfix7 = ".pc.Z";
    }
    std::string in_filename1;
    std::string in_filename2;
    std::string in_filename3;
    std::string in_filename4;
	std::string in_filename5;
    std::string in_filename6;

    in_filename1 = global_net_file; // database/GRCENEECnet.txt

    // results/$1/$1.gene (gene list & gene-ec homology)
    in_filename2.append(genome_dir);
    in_filename2.append(genome_name);
    in_filename2.append("/");
    in_filename2.append(genome_name);
    in_filename2.append(postfix3);

    // results/$1/$1.ort (gene-ec orthology)
    in_filename3.append(genome_dir);
    in_filename3.append(genome_name);
    in_filename3.append("/");
    in_filename3.append(genome_name);
    in_filename3.append(postfix4);

    // results/$1/$1.pc.Z (gene-gene phylogenetic corr)
    in_filename4.append(genome_dir);
    in_filename4.append(genome_name);
    in_filename4.append("/");
    in_filename4.append(genome_name);
    in_filename4.append(postfix5);

    // results/$1/$1.gc.Z (gene-gene gene clustering corr)
    in_filename5.append(genome_dir);
    in_filename5.append(genome_name);
    in_filename5.append("/");
    in_filename5.append(genome_name);
    in_filename5.append(postfix6);

    // results/$1/$1.pc.Z (gene-gene co-expression corr)
    in_filename6.append(genome_dir);
    in_filename6.append(genome_name);
    in_filename6.append("/");
    in_filename6.append(genome_name);
    in_filename6.append(postfix7);

// gene clustering z score to prob. of two being neighbor
    map<int, double> gcZ_Pnb;
    map_gcZ_Pnb(gcZ_Pnb);

// phylogenetic correlation z score to prob. of two being neighbor
    map<int, double> phZ_Pnb;
    map_phZ_Pnb(phZ_Pnb);

// gene co-expression z score to prob. of two being neighbor
    map<int, double> ceZ_Pnb;
    map_ceZ_Pnb(ceZ_Pnb);
    //cerr << "map_ceZ_Pnb done" << endl;

// sequence identity score (0~100) to prob. of being correct assignment
    map<int, double> hmZ_Ptr;
    map_hmZ_Ptr(hmZ_Ptr);

    double min_e_ph = phZ_Pnb[1];
    double min_e_gc = gcZ_Pnb[0];
    double min_e_ce = ceZ_Pnb[0];

// Read input files and initialize Variables here //
    //ifstream fin;

// global network file
	fin.open( in_filename1.c_str() );
	if( ! fin.is_open() )
	{
		cerr<<"Failed to open "<< in_filename1 <<" file."<<endl;
		return -1;
	}

    int num_ECs; // from global EC network
    vector<string> EC_names;
	vector< vector <int> > EC_EC_link; // linked ECs for each EC including itself (from global network file) (ECnb)
	vector< vector <double> > EC_EC_co_occur; // for an EC,EC co-occurence with neighbor ECs (including itself) from global network file
    //between each pair of EC numbers including self-correlation, phylogenetic correlation coeff (from global network file) -99 if nan

    read_global_network(fin, num_ECs, EC_names, EC_EC_link, EC_EC_co_occur);
    // self co-occurence is not used.

    fin.close();

    num_ECs += 1; // out-network-node

// genome file
	fin.open( in_filename2.c_str() );
	if( ! fin.is_open() )
	{
		cerr<<"Failed to open "<< in_filename2 <<" file."<<endl;
		return -1;
	}

    int num_genes = 0;
    vector<string> gene_names;
	vector< vector <int> > gene_EC_link; // for a gene, possible EC number indices (G2EC) <- This includes num_ECs-1 which is for the out-network-node
	vector< vector <double> > gene_EC_homol; // for each gene, sequence identity value for those possible EC numbers (log(prob))
	vector< vector <int> > known_answer; // for each gene, if there are known EC number answers push_back those EC number indices

    read_genome(fin, hmZ_Ptr, num_ECs, num_genes, gene_names, gene_EC_link, gene_EC_homol, known_answer); // <- This includes num_ECs-1 which is for the out-network-node
    fin.close();

    vector<int> gene_index_vec(num_genes);
    for(int i=0;i<num_genes;++i){
        gene_index_vec[i] = i;
    }

// orthology file
    fin.open( in_filename3.c_str() );
	if( ! fin.is_open() )
	{
		cerr<<"Failed to open "<< in_filename3 <<" file."<<endl;
		return -1;
	}

    vector< vector <int> > gene_EC_orth; // 0 or 1

    //read_orthology(fin, gene_EC_orth); // TODO (FIX data_table_io.h, gene_names and known_answer are also read in read_genome) read_orthology can be added to read_genome  //
    read_orthology(fin, num_genes, gene_names, gene_EC_orth, known_answer);

    fin.close();

// gene-gene phylogenetic correlation file
    fin.open( in_filename4.c_str() ); // results/$1/$1.pc.Z

	if( ! fin.is_open() )
	{
		cerr<<"Failed to open "<< in_filename4 <<" file."<<endl;
		return -1;
	}

	Matrix<double> gene_gene_phylo_corr(num_genes,num_genes,0.0); // log(prob)

    read_context_descriptor(fin, num_genes, phZ_Pnb, phZ_Pnb[1], gene_gene_phylo_corr);

    fin.close();

// gene-gene gene clustering file
    fin.open( in_filename5.c_str() ); // results/$1/$1.gc.Z

	if( ! fin.is_open() )
	{
		cerr<<"Failed to open "<< in_filename5 <<" file."<<endl;
		return -1;
	}

    Matrix<double> gene_gene_cluster_corr(num_genes,num_genes,0.0); // log(prob)
	read_context_descriptor(fin, num_genes, gcZ_Pnb, gcZ_Pnb[0], gene_gene_cluster_corr);

    fin.close();

// gene-gene co-expression file  // TODO copied from phylo_corr $1.pc.Z
    fin.open( in_filename6.c_str() ); // results/$1/$1.pc.Z

	if( ! fin.is_open() )
	{
		cerr<<"Failed to open "<< in_filename6 <<" file."<<endl;
		return -1;
	}

    Matrix<double> gene_gene_coexp_corr(num_genes,num_genes,0.0);
    read_context_descriptor(fin, num_genes, ceZ_Pnb, ceZ_Pnb[0], gene_gene_coexp_corr); // log(prob)

    fin.close();

// Done reading


// Updated data set
	vector< vector<int> > genes_on_EC(num_ECs); // for each EC, which genes are sitting on it. (CGenes)

    vector<int> assigned_gene2ec(num_genes); // for each gene, the index of the assigned EC number (CEC) (-1 if gene is not assigned)

// Output data
    Vector<int> count_ec_on_init(num_ECs, 0); // number of times that an EC isoccupied by "any" genes (max = num_iter)

    // acc_count_gene_on_ec_init[i][j]: # of times that gene i is on EC j
    vector<vector <int> > acc_count_gene_on_ec_init(num_genes);
    for(int i=0;i<num_genes;++i){
        for(int j=0;j<gene_EC_link[i].size();++j){
        acc_count_gene_on_ec_init[i].push_back(0);
        }
    }

    Vector<int> count_ec_on_act(num_ECs, 0); // number of times that an EC isoccupied by "any" genes (max = num_iter)

    vector<vector <int> > acc_count_gene_on_ec_act(num_genes);
    for(int i=0;i<num_genes;++i){
        for(int j=0;j<gene_EC_link[i].size();++j){
        acc_count_gene_on_ec_act[i].push_back(0);
        }
    }

    vector<double> EC_co_occur_vec(num_ECs-1); // not include out-node
    vector<int> EC_cnt_nbEC_vec(num_ECs-1); // not include out-node (number of neighboring ECs that genes are on.)

    energy tot_energy;
    vector< energy > gene_energy(num_genes); // for each gene, a vector of energy (GNE)

    tot_energy.initialize(min_e_ph, min_e_gc, min_e_ce);
    for(int i=0;i<num_genes;++i){
        gene_energy[i] = tot_energy;
    }

// Initial gene location assignment
    // initially assign all genes to node in network. not to the out-network-node ??
    initial_gene_ec_assignment(
        num_genes, num_ECs,
        EC_EC_link, gene_EC_link, EC_EC_co_occur, median_ec_co_occur,
        gene_gene_phylo_corr, gene_gene_cluster_corr, gene_gene_coexp_corr,
        gene_EC_homol, gene_EC_orth, genes_on_EC,
        assigned_gene2ec, gene_energy, EC_co_occur_vec, EC_cnt_nbEC_vec
    );

cerr << "done reading data" << endl;


// Exercise (is not recorded)
    //start = std::chrono::system_clock::now();
    std::clock_t c_start = std::clock();

    for(int t=0;t<num_init_run;++t){

        vector<int> tmp_gene_index_vec = gene_index_vec;

        std::random_shuffle( tmp_gene_index_vec.begin(), tmp_gene_index_vec.end() );

        for(int g=0;g<num_genes;++g){
            if(gene_EC_link[tmp_gene_index_vec[g]].size()<2){ // if possible ECs for the gene is 0 or 1, no new locations are available to move.
                // gene_EC_link is always >=1 (1 if only outside node)
                cerr << "gene_EC_link[" << tmp_gene_index_vec[g] << "].size()<2" << endl;
                int rand_gene_index;
                do{
                    rand_gene_index = int_rand(num_genes);
                }while(gene_EC_link[rand_gene_index].size()<2); // if possible ECs for the gene is 0 or 1, no new locations are available to move.
                tmp_gene_index_vec[g] = rand_gene_index;
            }
        }

        for(int g=0;g<num_genes;++g){
            // choose a random gene and switch its EC location:

            int rand_gene_index = tmp_gene_index_vec[g];

            //int new_ec_loc = assigned_gene2ec[rand_gene_index]; // original loc
            int new_ec_loc_index; // index of new_ec_loc in gene_EC_link[rand_gene_index]

            switch_loc_genes(
                num_ECs, min_e_ph, min_e_gc, min_e_ce,
                coeff_hm, coeff_ort, coeff_ph, coeff_gc, coeff_ce, p_ec_co_occur, p_no_EC,
                rand_gene_index,
                EC_EC_link, gene_EC_link, EC_EC_co_occur, median_ec_co_occur,
                gene_EC_homol, gene_EC_orth,
                gene_gene_phylo_corr, gene_gene_cluster_corr, gene_gene_coexp_corr,
                genes_on_EC, new_ec_loc_index, assigned_gene2ec, gene_energy,
                EC_co_occur_vec, EC_cnt_nbEC_vec
            );

            if(gene_EC_link[rand_gene_index][new_ec_loc_index] != assigned_gene2ec[rand_gene_index]){
                cerr << "new_ec_loc_index does not match in assigned_gene2ec" << endl;
            }

            acc_count_gene_on_ec_init[rand_gene_index][new_ec_loc_index] ++;
            myfile3 << rand_gene_index << " " << new_ec_loc_index << ", ";

        }

        // one iteration is done
        myfile3 << endl;

        // check occupancy of ECs in the network.
        // unique ECs in assigned_gene2ec vector = occupied ECs
        unique_EC_occupancy(assigned_gene2ec, count_ec_on_init);
    }

// Actual runs for recording
    for(int t=0;t<num_actual_run;++t){

        vector<int> tmp_gene_index_vec = gene_index_vec;
        std::random_shuffle( tmp_gene_index_vec.begin(), tmp_gene_index_vec.end() );
        for(int g=0;g<num_genes;++g){
            if(gene_EC_link[tmp_gene_index_vec[g]].size()<2){ // if possible ECs for the gene is 0 or 1, no new locations are available to move.
                // gene_EC_link is always >=1 (1 if only outside node)
                cerr << "gene_EC_link[" << tmp_gene_index_vec[g] << "].size()<2" << endl;
                int rand_gene_index;
                do{
                    rand_gene_index = int_rand(num_genes);
                }while(gene_EC_link[rand_gene_index].size()<2); // if possible ECs for the gene is 0 or 1, no new locations are available to move.
                tmp_gene_index_vec[g] = rand_gene_index;
            }
        }

        for(int g=0;g<num_genes;++g){
            // choose a random gene and switch its EC location:
            int rand_gene_index = tmp_gene_index_vec[g];

            //int new_ec_loc = assigned_gene2ec[rand_gene_index]; // original loc
            int new_ec_loc_index; // index of new_ec_loc in gene_EC_link[rand_gene_index]

            switch_loc_genes(
                num_ECs, min_e_ph, min_e_gc, min_e_ce,
                coeff_hm, coeff_ort, coeff_ph, coeff_gc, coeff_ce, p_ec_co_occur, p_no_EC,
                rand_gene_index,
                EC_EC_link, gene_EC_link, EC_EC_co_occur, median_ec_co_occur,
                gene_EC_homol, gene_EC_orth,
                gene_gene_phylo_corr, gene_gene_cluster_corr, gene_gene_coexp_corr,
                genes_on_EC, new_ec_loc_index, assigned_gene2ec, gene_energy,
                EC_co_occur_vec, EC_cnt_nbEC_vec
            );

            if(gene_EC_link[rand_gene_index][new_ec_loc_index] != assigned_gene2ec[rand_gene_index]){
                cerr << "new_ec_loc_index does not match in assigned_gene2ec" << endl;
            }

            //acc_count_gene_on_ec_act[rand_gene_index][new_ec_loc_index] ++;
            for(int i=0;i<num_genes;++i){ // NEW
                int index = -1;
                find_index_of_id(gene_EC_link[i], assigned_gene2ec[i], index);
                acc_count_gene_on_ec_act[i][index] ++;
            } // NEW END
            myfile3 << rand_gene_index << " " << new_ec_loc_index << ", ";

            unique_EC_occupancy(assigned_gene2ec, count_ec_on_act); // NEW
        }

        // one iteration is done
        myfile3 << endl;

        // check occupancy of ECs in the network.
        // unique ECs in assigned_gene2ec vector = occupied ECs
        //unique_EC_occupancy(assigned_gene2ec, count_ec_on_act);

    }

    //end = std::chrono::system_clock::now();
    //std::chrono::duration<double> elapsed_seconds = end-start;
    //std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

    std::clock_t c_end = std::clock();

    long double time_elapsed = (c_end-c_start) / CLOCKS_PER_SEC;
    std::cout << "CPU time used: " << time_elapsed << " s\n";

    myfile1 << num_genes << "\t" << num_actual_run << endl;
    for(int i=0;i<num_genes;++i){
        myfile1 << i << "\t" << gene_EC_link[i].size()-1 << "\t";
        for(int j=0;j<gene_EC_link[i].size();++j){
            myfile1 << acc_count_gene_on_ec_act[i][j] << "\t";
        } myfile1 << endl;
    }
    myfile1.close();


    myfile2 << num_ECs << "\t" << num_actual_run << endl;
    for(int i=0;i<num_ECs;++i){
        myfile2 << i << "\t" << count_ec_on_act[i] << "\t" << num_actual_run-count_ec_on_act[i] << endl;
    }
    myfile2.close();
    myfile3.close();

    return 0;
}
