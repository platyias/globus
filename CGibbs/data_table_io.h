//
//  data_table_io.h
//
//  Created by Jeewoen Shin on 4/13/17.
//

#ifndef data_table_IO
#define data_table_IO

#include "tnt_126/tnt_cmat.h"
using TNT::Matrix;

#include "tnt_126/tnt_vec.h"
using TNT::Vector;

#include <vector>
using std::vector;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <fstream>
using std::ofstream;
using std::ifstream;

#include <sstream>
using std::stringstream;
using std::istringstream;

#include <iomanip> // setprecision

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <cstdlib> // strtod

//#include "z_score_map.h"
#include "z_score_map_log10.h"

void naming( string& base, int num, int maxnum )
{
  int sub = maxnum;
  stringstream filename;
  filename<<base<<"_";
  int count = 0;
  while( sub > 9 )
  {
    ++count;
    sub /= 10;
  }
  sub = num;
  while( sub > 9 )
  {
    --count;
    sub /= 10;
  }
  for( int i = 0; i < count; ++i )
  {
    filename<<0;
  }
  filename<<num;

  base = filename.str();
}

// READ
void read_input_param(
        ifstream& fp,
        string &genome_name,
        string &genome_dir,
        string &global_net_file,
        double &coeff_hm,
        double &coeff_ph, double &coeff_gc, double &coeff_ce,
        double &coeff_ort,
        double &p_ec_co_occur,  double &p_no_EC,
        double &median_ec_co_occur){

    string line;

    if( getline(fp, line) ){ // 1
        genome_name = line;
    }
    if( getline(fp, line) ){ // 2
        genome_dir = line;
    }
    if( getline(fp, line) ){ // 3
        global_net_file = line;
        //cout << global_net_file << endl;

        global_net_file.erase(global_net_file.find_last_not_of(" \n\r\t")+1);
        //cerr << global_net_file << endl;
    }
    if( getline(fp, line) ){ // 4
        coeff_hm = strtod(line.c_str(), NULL);
    }
    if( getline(fp, line) ){ // 5
        coeff_ph = strtod(line.c_str(), NULL);
    }
    if( getline(fp, line) ){ // 6
        coeff_gc = strtod(line.c_str(), NULL);
    }
    if( getline(fp, line) ){ // 7
        coeff_ce = strtod(line.c_str(), NULL);
    }
    if( getline(fp, line) ){ // 8
        coeff_ort = strtod(line.c_str(), NULL);
    }
    if( getline(fp, line) ){ // 9
        p_ec_co_occur = strtod(line.c_str(), NULL);
    }
    if( getline(fp, line) ){ // 10
        p_no_EC = strtod(line.c_str(), NULL);
    }
    if( getline(fp, line) ){ // 10
        median_ec_co_occur = strtod(line.c_str(), NULL);
    }
}

//T1 = int
//T2 = double
//T3 = string
//template< typename T1 >
//template< typename T1 , typename T2, typename T3>
void read_global_network(
        ifstream& fp,
        int &num_ECs,
        vector<string> &EC_names,
        vector< vector <int> > &EC_EC_link,
        vector< vector <double> > &EC_EC_co_occur){

    vector<int> tmp_EC_index;
    vector<string> tmp_EC_names;
    vector< vector <int> > tmp_EC_EC_link;
    vector< vector <double> > tmp_EC_EC_co_occur;

    string delim1 = "///";
    string delim2 = "\t";
    num_ECs = 0;
    // Each line is an EC number
    while(fp){
        string line;
        vector<string> section; // vector length = 3
        int posi = 0;

        if( !getline(fp, line) ){
            break;
        }
        num_ECs ++;

        while( (posi = line.find(delim1)) != std::string::npos ){
            section.push_back(line.substr(0,posi));
            line.erase( 0, posi+delim1.length() );
        }

        string word;
        vector<string> word_list;

        // section[0] (neighbor EC indices, including itself)
        while( (posi = section[0].find(delim2)) != std::string::npos ) {
            word_list.push_back(section[0].substr(0, posi));
            section[0].erase( 0, posi+delim2.length() );
        }

        int my_ec_index = atoi(word_list[0].c_str());
        tmp_EC_index.push_back(my_ec_index); // index of the EC
        tmp_EC_names.push_back(word_list[1]); // EC number
        int num_ec_nb = atoi(word_list[2].c_str()); // number of neighbor ECs including itself
        //if(my_ec_index == 14){
        //    cerr << "EC 14 = " << word_list[1] << " " << "num_ec_nb = " << num_ec_nb << "word_list.size() = " << word_list.size() << endl;
        //}
        vector<int> tmp_nb_index;
        for(int i=3;i<word_list.size();++i){
            int nb_index = atoi(word_list[i].c_str());
            tmp_nb_index.push_back(nb_index);
            //EC_EC_link[my_ec_index].push_back(nb_index);
        }
        tmp_EC_EC_link.push_back(tmp_nb_index);

        // section[1] (EC-EC co-occurence, including self co-occurence)
        //vector<string> word_list; // initialize
        word_list.clear();

        while( (posi = section[1].find(delim2)) != std::string::npos ) {
            if(posi>0){
                word_list.push_back(section[1].substr(0, posi));
            }
            section[1].erase( 0, posi+delim2.length() );
        }
        vector<double> tmp_co_occur;
        for(int i=0;i<word_list.size();++i){
            if (word_list[i] != "nan"){
                double co_occur = strtod(word_list[i].c_str(), NULL);
                tmp_co_occur.push_back(co_occur);
                //EC_EC_co_occur[my_ec_index].push_back(co_occur);
            }
            else{
                tmp_co_occur.push_back(0.0); // <- NAN
            }
        }
        tmp_EC_EC_co_occur.push_back(tmp_co_occur);

        // one line ended
    }

    EC_names.resize(num_ECs);
    EC_EC_link.resize(num_ECs);
    EC_EC_co_occur.resize(num_ECs);

    for(int i=0;i<num_ECs;++i){
        EC_names[tmp_EC_index[i]] = tmp_EC_names[i];
        EC_EC_link[tmp_EC_index[i]] = tmp_EC_EC_link[i];
        //if(EC_EC_co_occur[tmp_EC_index[i]].size() > 0){
        //    EC_EC_co_occur[tmp_EC_index[i]].clear();
        //}
        EC_EC_co_occur[tmp_EC_index[i]] = vector<double> ();
        for(int j=0;j<tmp_EC_EC_co_occur[i].size();++j){ // WHY??
            EC_EC_co_occur[tmp_EC_index[i]].push_back(tmp_EC_EC_co_occur[i][j]);
            //cerr << EC_EC_co_occur[tmp_EC_index[i]][j] << " ";
        } //cerr << endl;
    }
}

//T1 = int
//T2 = double
//T3 = string
//template< typename T1 >
//template< typename T1 , typename T2, typename T3>
void read_genome(
        ifstream& fp,
        map<int, double> hmZ_Ptr,
        const int &num_ECs,
        int &num_genes,
        vector<string> &gene_names,
        vector< vector<int> > &gene_EC_link,
        vector< vector<double> > &gene_EC_homol,
        vector< vector<int> > &known_answer){

    vector<int> tmp_gene_index;
    vector<string> tmp_gene_names;
    vector< vector <int> > tmp_gene_EC_link;
    vector< vector <double> >tmp_gene_EC_homol;
    vector< vector <int> > tmp_known_answer;

    num_genes = 0;
    string delim1 = "///";
    string delim2 = "\t";

    // Each line is a gene
    while(fp){
        string line;
        vector<string> section; // vector length = 3
        int posi = 0;

        if( !getline(fp, line) ){
            break;
        }
        num_genes ++;

        while( (posi = line.find(delim1)) != std::string::npos ){
            section.push_back(line.substr(0,posi));
            line.erase( 0, posi+delim1.length() );
        }

        //cerr << section.size() << endl;

        string word;
        vector<string> word_list;

        // section[0] (For a gene, possible EC indices)
        while( (posi = section[0].find(delim2)) != std::string::npos ) {
            word_list.push_back(section[0].substr(0, posi));
            section[0].erase( 0, posi+delim2.length() );
        }

        int my_gene_index = atoi(word_list[0].c_str());
        tmp_gene_index.push_back(my_gene_index); // index of the gene
        tmp_gene_names.push_back(word_list[1]); // gene name
        int num_candi_ec = atoi(word_list[2].c_str()); // number of candidate EC numbers

        vector<int> tmp_ec_index;
        for(int i=3;i<word_list.size();i++){
            int ec_index = atoi(word_list[i].c_str());
            tmp_ec_index.push_back(ec_index);
            //gene_EC_link[my_gene_index].push_back(ec_index);
        }

        tmp_ec_index.push_back(num_ECs-1); // includes outside node

        tmp_gene_EC_link.push_back(tmp_ec_index);

        // cerr << tmp_ec_index.size() << " " << num_candi_ec << endl;
        if(tmp_ec_index.size() != num_candi_ec+1){
            cerr << "tmp_gene_EC_link.size() != num_candi_ec" << endl;
        }

        // section[1] (gene-EC homology score)
        //vector<string> word_list; // initialize
        word_list.clear();

        while( (posi = section[1].find(delim2)) != std::string::npos ) {
            if(posi>0){
                word_list.push_back(section[1].substr(0, posi));
            }
            section[1].erase( 0, posi+delim2.length() );
        }

        vector<double> tmp_seq_ident;
        //for(int i=0;i<word_list.size();i++){
        for(int i=0;i<word_list.size();i++){ //WHY??
            double seq_ident = strtod(word_list[i].c_str(), NULL);
            seq_ident = seq_ident/5;
            tmp_seq_ident.push_back( hmZ_Ptr[(int) seq_ident] );
            //tmp_seq_ident.push_back(seq_ident);
            //gene_EC_homol[my_gene_index].push_back(seq_ident);
        }
        tmp_seq_ident.push_back(0.0); // includes outside node
        tmp_gene_EC_homol.push_back(tmp_seq_ident);

        // section[2] (known yeast answer)
        //vector<string> word_list; // initialize
        word_list.clear();

        while( (posi = section[2].find(delim2)) != std::string::npos ) {
            if(posi>0){
                word_list.push_back(section[2].substr(0, posi));
            }
            section[2].erase( 0, posi+delim2.length() );
        }

        if (atoi(word_list[0].c_str()) > 0){ // the number of known ECs

            vector<int> tmp_ans_ec;
            for(int i=0;i<=atoi(word_list[0].c_str());i++){
                int ans_ec = atoi(word_list[i].c_str());
                tmp_ans_ec.push_back(ans_ec);
                //known_answer[my_gene_index].push_back(ans_ec);
            }
            tmp_known_answer.push_back(tmp_ans_ec);
        }
        else{ // No known answers
            vector<int> tmp_ans_ec;
            tmp_known_answer.push_back(tmp_ans_ec);
        }

        // one line ended
    }

    //cerr << num_genes << endl;

    gene_names.resize(num_genes);
    gene_EC_link.resize(num_genes);
    gene_EC_homol.resize(num_genes);
    known_answer.resize(num_genes);

    for(int i=0;i<num_genes;i++){
        gene_names[tmp_gene_index[i]] = tmp_gene_names[i];
        //gene_EC_link[tmp_gene_index[i]] = tmp_gene_EC_link[i];
        gene_EC_link[tmp_gene_index[i]] = vector<int> ();
        for(int j=0;j<tmp_gene_EC_link[i].size();j++){
            gene_EC_link[tmp_gene_index[i]].push_back(tmp_gene_EC_link[i][j]);
        }
        gene_EC_homol[tmp_gene_index[i]] = tmp_gene_EC_homol[i];
        known_answer[tmp_gene_index[i]] = tmp_known_answer[i];
    }
}

//T1 = int
//T3 = string
//template< typename T1, typename T3 >
void read_orthology(  // TODO // Input file can be combined with with .gene file.
        ifstream& fp,
        const int &num_genes,
        vector<string> &gene_names,
        vector< vector <int> > &gene_EC_orth){
        //vector< vector <int> > &known_answer){

    vector<int> tmp_gene_index;
    vector<string> tmp_gene_names;
    vector< vector <int> > tmp_gene_EC_orth; // binary 0: no orthologous genes in any of the organisms
    //vector< vector <int> > tmp_known_answer;

    string delim1 = "///";
    string delim2 = "\t";

    // Each line is a gene
    while(fp){
        string line;
        vector<string> section; // vector length = 3
        int posi = 0;

        if( !getline(fp, line) ){
            break;
        }

        while( (posi = line.find(delim1)) != std::string::npos ){
            section.push_back(line.substr(0,posi));
            line.erase( 0, posi+delim1.length() );
        }

        string word;
        vector<string> word_list;

        // section[0] orthology
        while( (posi = section[0].find(delim2)) != std::string::npos ) {
            word_list.push_back(section[0].substr(0, posi));
            section[0].erase( 0, posi+delim2.length() );
        }

        int my_gene_index = atoi(word_list[0].c_str());
        tmp_gene_index.push_back(my_gene_index); // index of the gene
        tmp_gene_names.push_back(word_list[1]); // gene name
        //int num_ortholog = atoi(word_list[2].c_str());
        //tmp_num_ortholog.push_back(num_ortholog); // number of orthologs

        // section[1] (orthology score for gene-ec pai)
        //vector<string> word_list; // initialize
        word_list.clear();

        while( (posi = section[1].find(delim2)) != std::string::npos ) {
            if(posi>0){
                word_list.push_back(section[1].substr(0, posi));
            }
            section[1].erase( 0, posi+delim2.length() );
        }

        vector<int> tmp_orth_index;
        //for(int i=3;i<word_list.size();i++){
        //for(int i=0;i<word_list.size();i++){
        for(int i=0;i<word_list.size();i++){ //WHY??
            int orth_index = atoi(word_list[i].c_str());
            tmp_orth_index.push_back(orth_index);
        }
        tmp_orth_index.push_back(0); // includes outside network node
        tmp_gene_EC_orth.push_back(tmp_orth_index);

        // section[2] known answer
        //vector<string> word_list; // initialize
        word_list.clear();
        /*
        while( (posi = section[2].find(delim2)) != std::string::npos ) {
            if(posi>0){
                word_list.push_back(section[2].substr(0, posi));
            }
            section[2].erase( 0, posi+delim2.length() );
        }

        if (atoi(word_list[0].c_str()) > 0){ // the number of known ECs

            vector<int> tmp_ans_ec;
            for(int i=0;i<=atoi(word_list[0].c_str());i++){
                int ans_ec = atoi(word_list[i].c_str());
                tmp_ans_ec.push_back(ans_ec);
                //known_answer[my_gene_index].push_back(ans_ec);
            }
            tmp_known_answer.push_back(tmp_ans_ec);
        }
        else{ // No known answers
            vector<int> tmp_ans_ec;
            tmp_known_answer.push_back(tmp_ans_ec);
        }
        */
        // one line ended
    }

    gene_names.resize(num_genes);
    gene_EC_orth.resize(num_genes);
    //known_answer.resize(num_genes);

    for(int i=0;i<num_genes;i++){
        gene_names[tmp_gene_index[i]] = tmp_gene_names[i];
        gene_EC_orth[tmp_gene_index[i]] = tmp_gene_EC_orth[i];
        //known_answer[tmp_gene_index[i]] = tmp_known_answer[i];
    }
}

//T1 = int
//T2 = double
//T3 = string
//template< typename T1 >
//template< typename T1 , typename T2, typename T3>
void read_context_descriptor(
    ifstream& fp,
    const int &num_genes,
    map<int, double> &Z_Pnb,
    const double &unknown_cut,
    Matrix<double> &gene_gene_mat){
    // Z_Pnb length = 41;
    // Each line is a gene

    int gene1 = 0; //row
    int gene2 = 0; //colm

    while(fp){
        string line;

        if( !getline(fp, line) ){
            break;
        }

        //for(int i=0;i<num_genes;i++){
        stringstream iss(line);
        double corr;
        int range_start_i;
        //double range_start_d;
        //string range_start_s;

        gene2 = 0;
        while(iss >> corr){ // TODO check if it reads -2 as a double
            // Z_Pnb length = 41;
            // if corr falls in [3-3.5) range, range_start = 3

            if(corr < 0){ // = -2
                gene_gene_mat[gene1][gene2] = unknown_cut;
            }
            else if(corr >= 20){
                //gene_gene_mat[gene1][gene2] = Z_Pnb["20"];
                gene_gene_mat[gene1][gene2] = Z_Pnb[40];
            }
            else{
                if ( round(corr) <= corr ){
                    range_start_i = (int) round(corr) * 2;
                    /*
                    range_start_i = (int) round(corr);
                    stringstream ss;
                    ss <<  range_start_i;
                    range_start_s = ss.str();
                    */
                }
                else{
                    double range_start_d = 2*(round(corr) - 0.5);
                    range_start_i = (int) range_start_d;
                    /*
                    range_start_d = round(corr) - 0.5;
                    stringstream ss;
                    ss << std::fixed << std::setprecision(1) << range_start_d; //X.X
                    range_start_s = ss.str();
                    */
                }
                //gene_gene_mat[gene1][gene2] = Z_Pnb[range_start_s];
                gene_gene_mat[gene1][gene2] = Z_Pnb[range_start_i];
            }

            gene2 ++;
        }
        if(gene2 != num_genes){
	        cerr<<"Error: read_context_descriptor matrix (gene2) size = " << gene2 << " " << num_genes <<endl;
        }

        gene1 ++;
    }

    if(gene1 != num_genes){
	    cerr<<"Error: read_context_descriptor matrix (gene2) size = " << gene1 << " " << num_genes <<endl;
    }

}

/*
template< typename T >
void save_single(ofstream& fp, const T& element){
	cerr<<"_single "<<element<<endl;
	fp.write( (char*)&element, sizeof(element) );
}

template< typename T >
void read_single(ifstream& fp, T& element){
	fp.read( (char*)&element, sizeof(element) );
	cerr<<"_single "<<element<<endl;
}
*/

#endif

