//
//  random.cpp
//
//  Created by Jeewoen Shin on 4/13/17.
//

#include "parameters.h"
#include "random.h"

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <ctime>
//#include <numeric>
//#include <algorithm>
//#include <random>
//#include <chrono>

//unsigned seed = std::chrono::system_clocl::now().time_since_epoch().count();
//std::default_random_engine generator(seed);
//std::default_random_engine generator;
/*
 double normal_rand(){
 //std::default_random_engine generator;
 std::normal_distribution<double> distribution(0.0,1.0);
 
 double number;
 do{
 number = distribution(generator);
 }while(number == 0.0);
 
 return number;
 }
 */
double normal_pdf(double &x){
	const double PI = 3.141592653589793;
	double p;
	p = (1/sqrt(2*PI)) * exp(-x*x/2);
		
	return p;
}

double normal_rand(){
	double u; 
    double v; 
    double r; 
    double c;
    do{
        do{
            u = ((double) rand() / (RAND_MAX)) * 2 - 1;
            v = ((double) rand() / (RAND_MAX)) * 2 - 1;
            r = u * u + v * v;
        } while (r == 0 || r > 1);
        
        c = sqrt(-2 * log(r) / r);
    }while(fabs(u*c)<tol);
    return u * c;
    
}
/*
 double uniform_rand(){
 //std::default_random_engine generator;
 std::uniform_real_distribution<double> distribution(0.0,1.0);
 
 double number = distribution(generator);
 
 return number;
 }
 */
double uniform_rand(){ // [0,1]
    double x = (double)rand()/(RAND_MAX);
    return  x;
}

/*
 int int_rand(int N){
 //std::default_random_engine generator;
 std::uniform_int_distribution<int> distribution(0,N-1);    
 int number = distribution(generator);
 
 return number;
 }
 */

int int_rand(int N) // 0, 1, 2, ..., N-1
{
    int x;
    do {
        x = (int)(rand()/(((double)RAND_MAX+1.0)/N));
    } while (x >= N);
    return x;
}

