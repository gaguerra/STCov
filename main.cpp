//
//  main.cpp
//  Covariance_Calculator
//
//  Created by Geno Guerra on 11/9/19.
//  Copyright @ Geno Guerra. All rights reserved.
//  University of California, Berkeley

#include "main.hpp"
#include <fstream>
#include <iostream>
#include<iomanip>
#include <math.h>
#include <Eigen/Dense> // Need to tell the program what the path is (aka /usr/local/include in the build settings (searchpath)
#include <vector>
#include <set>
#include <math.h>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <unistd.h>


using namespace std;
using namespace Eigen;


struct tree_info
{
    int node_num ;
    char leaf_flag;
    int child_node_nums[2] ;
    int parent_num;
    int node_size ;
    int *node_list ;
    double branch_lengths[2] ;
};

class BasicBranchNode {
public:
    
    int Offspring [2];
    int Parent;
    double Times [2];
    double Ne;
    // int Ancestors [n]; // This needs to scale up with the number of individuals.
    vector<int> Ancestors;
    double P22; // see branch node for description.
    double TP;
    double T2P;
    double Pptheta; // doesnt include the ^ 2, that is included in the composite to make Pnomut.
    Matrix<double, Dynamic, Dynamic> Ppthetas; // Pptheta for every gene, allowing for differing gene lengths.
    double  Mean; // Mean of any pair that diverge at Times[0]
    double Variance; // Variance of any pair that diverge at Times[0]
    double Pnomut; // P of no mutations of any pair that diverge at Times[0]
    
    Matrix<double,Dynamic,Dynamic> Pnomuts; // PnoMut for every gene, allowing for differing gene lengths.
    
    double TPPnomut; //expected value of t given no mutation occurred.
    double T2PPnomut; // exp val of t^2 given no mutations occurred.
    double ETpnomut; // Overall expected val of t given no mut.
    double ET2pnomut; //same as above but t^2
    double ZIpptheta; //prob of no mutation, for the zeroinflated stuff to correct so we can use the real pptheta
    double ZIPnoMut; //Same comment as above, analagous to Pnomut but zeroinflated so wrong value.
    double L; // Used in TransformToPLUS and TransformBackPlus.
};

class SpeciesTreeType{
    
public:
    // BasicBranchNode Branches [2*n-1];
    std::vector<BasicBranchNode> Branches ; //2*n-1 of them.
    //Need to resize this when you declare it.
    
};


// Topology/History data type:
class History {
public:
    int T4 [2]; // Two taxa involved in the first coalescence.
    int T3 [3]; // Two or three taxa involved in the second colalesence.
    int T2 [4]; // All four taxa are involved in the last coalesence.
    
    double ST4; //The first time T4 can occur, aka the speciation event involving the two taxa.
    double ST3; //The first time T3 can occur, the most ancient speciation event involving the 2-3 taxa.
    double ST2; //The first time T2 can occur, aka the most ancient speciation event.
    
    double Prob; // Probability of the history given the species tree.
};

//const InitializerElements IE;



// Declaring the data type for a node. Which consists of offspring, parent, times, and population sizes, along with some useful precomputed values.
class BranchNode {
public:
    //////////// Note 20 found here just needs to be n or larger. (n is the number of species/individuals) ////////
    
    int Offspring [2];  // Offspring Nodes (Note nodes are labeled 0-6 in a quartet tree.
    int Parent; // Parent node
    double Times [20]; //Times on the branch (t_i = start height of interval i, t_i+1 = end height), needs to be larger than (2*n_taxa - 1)
    double Ne [20]; // Population sizes within the times. Needs to be of same length.
    
    int Counter; // Keeps track of how many intervals there are in this branch.
    
    double FullP22; //P22 across whole branch
    double FullP33; // P33 across whole branch
    double FullP32; //P32 across the whole branch
    double P22 [20]; // Probability of 2 lineages entering an interval, and 2 leaving. P21 = 1- P22.
    double P33 [20]; // Prob 3 in 3 out
    double P32 [20]; // P 3 in 2 out, Note P31 = 1 - P33 - P32
    //double P44 [50]; // Only useful in last branch, unsure if I should calculate it here.
    //double P43 [50]; // 4 in 3 out
    //double P42 [50]; // 4 in 2 out, Note P41 = 1 - P44 - P43 - P42
    double TP [20]; // int_t_i ^t_{i+1} T*P(T|T >t_i)dT
    double TPP22 [20]; // int_t_i ^t_{i+1} T*P(T|T >t_i)*P22(N_i,(T,t_i+1))dT
    double T2P [20]; // int_t_i ^t_{i+1} T^2 *P(T|T >t_i) dT
    double TPP22N [20]; // int_t_i ^t_{i+1} T*P(T|T >t_i)*P22(N_(i),(T,t_i+1))dT when N_(i)  =/= N_i
    double PP22N [20]; // int_t_i ^t_{i+1} P(T|T >t_i)*P22(N_(i),(T,t_i+1))dT when N_(i)  =/= N_i
    double T2PP22 [20]; //int_t_i ^t_{i+1} T^2 *P(T|T >t_i)*P22(N_i,(T,t_i+1))dT
    double PP22 [20]; // int_t_i ^t_{i+1} P(T) * P22(N_i, (T, t_i+1)) dT
    double EgtTi [20]; // Expected value of T, given T is greater than t_i, the minimum of the interval. Exists
    //  as a recursion.
    double EgtT2i [20]; //expectation T^2 given greater than t.
    
    double P22A [20]; //using rate 2N
    double P22B [20]; //using rate 2/3 * N
    double P22C [20]; //using rate 1/3 * N
    double TPA [20]; //using rate 2N
    double TPB [20]; //using rate 2/3 * N
    double TPC [20]; //using rate 1/3 * N
    double PP22BA [20]; // using rate 2/3*N and rate 2N
    double PP22CB [20]; // using rate 1/3*N and rate 2/3* N
    double PP22CA [20]; // using rate 1/3 * N and Rate 2N
    double PP22CBA [20]; // using all three rates, its a mess.
    double P22H [20]; // using rate N
    double TPH [20]; // using rate N
    
    
    // New Approach, new requirements.
    double TP5 [20]; // int_{t_i}^t_{i+1} t * P(5 p coal = t |T_i)dt
    double P5NC [20]; // Prob 5 pairs, no coal in interval. e^(-5(t_i+1 - t_i)/2Ni)
    double TPP5 [20]; // int T_cd P(T_cd)e^(-5(Tcd-ti)/2Ni)dTcd
    double T2PP5 [20];// int T_cd^2 P(T_cd)e^(-5(Tcd-ti)/2Ni)dTcd
    double TPP2 [20]; // int T_cd P(T_cd)e^(-2(Tcd-ti)/2Ni)dTcd
    double T2PP2 [20]; // int T_cd^2 P(T_cd)e^(-2(Tcd-ti)/2Ni)dTcd
    double P5P2NC [20]; // int P(5p coal = t |t_i) * P(2 p no coal (t,t_i+1)dt
    double TPP2P22 [20]; //int T P(T) e^(-2(T-ti)/2Ni) * e^(-(ti+1 - T)/2Ni)dt
    double TPP5P22 [20]; //int T P(T) e^(-5(T-ti)/2Ni) * e^(-(ti+1 - T)/2Ni)dt
    double TPP5P3NC [20] ; //int T P(T) e^(-5(T-ti)/2Ni) * e^(-3(ti+1 - T)/2Ni)dt
    double P3P22 [20]; // int P(3p coal = t|ti) * e^( -(tj+1 - t)/2Ni)dt
    double TP2 [20]; // int t_i t_i+1 t * P(2 p coal = t |t i, ti+1)
    double P2NC [20]; // Prob 2 pairs, no coal in interval, e^(-2(ti+1 - ti)/2Ni)
    double TP3 [20]; // int ti ti+1 t * P(3 p coal = t | t_i)dt
    double P3NC [20]; // Prob 3 pairs, no coal in interval e^(-3(ti+1 - ti)/2Ni)
    double P6NC [20]; // Prob no coal when 4 lineages, 6 pairs.
    
    
    // Equations for shared branch length.
    double TPP5NCP3NC [20]; //t*p(t)*P(5p no coal(ti,t))P(3p no coal (t,ti+1))
};


// Declaring the data type for quartet trees. there will be four of these, with sets of ancestry that we will alter.
class Quartet {
public:
    int indiv;
    std::set <int> Ancestry;
    
};

int MatLabel(int i, int j, int s){
    int label;
    
    // s = num_individuals
    
    int m = max(i,j);
    int l = -max(-i,-j);
    
    int k = 0;
    
    for (int g =0; g <= l; g++) {
        k = k + g;
    }
    label = l * s - k + m - l - 1;
    
    return label;
}


//Function to calculate P22 across entire branch
double CalcP22 (BranchNode Branch){
    double fullP22 = 1.0;
    for (int i=0; i < Branch.Counter; i++) {
        fullP22 = fullP22 * Branch.P22[i];
    }
    return fullP22;
}

//Function to calculate P33 across entire branch
double CalcP33 (BranchNode Branch){
    double fullP33 = 1;
    for (int i=0; i <Branch.Counter; i++) {
        fullP33 = fullP33 * Branch.P33[i];
    }
    return fullP33;
}

//Function to calculate P32 across entire branch
double CalcP32 (BranchNode Branch){
    double fullP32 = 0;
    for (int i=0; i < Branch.Counter; i++) {
        double pl = 1.0; //placeholder to build on in this iteration.
        if( i > 0){
            for( int k=0 ; k < i; k++){
                pl= pl * Branch.P33[k];
            }
        }
        if ( i < Branch.Counter -1 ){
            for (int l = i+1; l < Branch.Counter; l++) {
                pl =pl * Branch.P22[l];
            }
        }
        fullP32 = fullP32 + pl * Branch.P32[i]  ; //i dont knwo if the 1/3 should be there.
    }
    
    return fullP32;
}

double GetCase1A ( BranchNode Branch , int i){
    
    double intervalSum = 0.0;
    double ProbNoCoal = 1.0;
    if ( i !=0){
        // If there are intervals to worry about before i.
        for (int j = 0; j < i; j++) {
            
            intervalSum = intervalSum + (1.0/5.0 * Branch.TP5[j]) * ProbNoCoal;
            ProbNoCoal = ProbNoCoal * Branch.P5NC[j];
        }
        
        
    }else{intervalSum = 0.0;}
    return intervalSum;
    
}


double GetCase1 ( BranchNode Branch ){
    
    double intervalSum =0.0;
    double ProbNoCoal = 1.0;
    double P5NoCoal = 1.0;
    
    for (int i =0; i < Branch.Counter; i++) {
        double Ni = Branch.Ne[i];
        double ti = Branch.Times[i];
        double C1A = GetCase1A(Branch, i);
        double C1i= Branch.TP[i] * C1A + P5NoCoal * (-2.0*Ni/25.0 * Branch.TPP5[i] - 1.0/5.0 * Branch.T2PP5[i] + 1.0/5.0 * (2.0*Ni/5.0 + ti) * Branch.TP[i] );
        
        intervalSum = intervalSum + C1i * ProbNoCoal;
        ProbNoCoal = ProbNoCoal * Branch.P22[i];
        
        P5NoCoal = P5NoCoal * Branch.P5NC[i];
    }
    
    //cout << "Case 1 = " << intervalSum << "\n";
    return intervalSum;
    
}

double GetCase1altA ( BranchNode Branch , int i){
    
    double intervalSum = 0.0;
    double ProbNoCoal = 1.0;
    if ( i !=0){
        // If there are intervals to worry about before i.
        for (int j = 0; j < i; j++) {
            
            intervalSum = intervalSum + (1.0/5.0 *(1-Branch.P5NC[j])) * ProbNoCoal;
            ProbNoCoal = ProbNoCoal * Branch.P5NC[j];
        }
        
        
    }else{intervalSum = 0.0;}
    return intervalSum;
    
}



double GetCase1alt (BranchNode Branch){
    double intervalSum =0.0;
    double ProbNoCoal = 1.0;
    double P5NoCoal = 1.0;
    
    for (int i =0; i < Branch.Counter; i++) {
        
        double C1A = GetCase1altA(Branch, i);
        double C1i= Branch.T2P[i] * C1A + P5NoCoal * (1.0/5.0 * Branch.T2P[i] - 1.0/5.0 * Branch.T2PP5[i] );
        
        intervalSum = intervalSum + C1i * ProbNoCoal;
        ProbNoCoal = ProbNoCoal * Branch.P22[i];
        
        P5NoCoal = P5NoCoal * Branch.P5NC[i];
    }
    
    // cout << "Case 1alt = " << intervalSum << "\n";
    return intervalSum;
    
    
}

double GetCase2A(BranchNode Branch, int i){
    
    double intervalSum =0.0;
    double ProbNoCoal = 1.0;
    
    if ( i !=0){
        
        for (int j = 0; j < i; j++) {
            double Nj = Branch.Ne[j];
            
            double C2ai = 0.0;
            
            C2ai = 2.0/5.0 * (-(Nj + Branch.Times[j+1]) * Branch.P5P2NC[j] + Nj * (1-Branch.P5NC[j]) + Branch.TP5[j]);
            intervalSum = intervalSum + C2ai * ProbNoCoal;
            ProbNoCoal = ProbNoCoal * Branch.P5NC[j];
            
        }
        
    }else{ intervalSum = 0.0;}
    
    
    
    
    return intervalSum;
}

double GetCase2B2 ( BranchNode Branch, int j, int i ){
    
    double intervalSum =0.0;
    //double ProbNoCoal = 1.0;
    double P2NoCoal = 1.0;
    if(j+1 < i) {
        for (int k = j+1; k < i; k++) {
            double C2b2i = 0.0;
            
            C2b2i = 1.0/2.0 * Branch.TP2[k];
            intervalSum = intervalSum + C2b2i * P2NoCoal;
            P2NoCoal = P2NoCoal * Branch.P2NC[k];
            
            
        }
    }
    intervalSum = intervalSum * Branch.TP[i] + 1.0/2.0 * P2NoCoal *( - Branch.Ne[i] * Branch.TPP2[i] - Branch.T2PP2[i] + (Branch.Ne[i] + Branch.Times[i]) * Branch.TP[i]);
    
    
    
    return intervalSum;
}

double GetCase2B( BranchNode Branch, int i){
    
    double intervalSum =0.0;
    double ProbNoCoal = 1.0;
    
    if ( i !=0){
        
        for (int j =0; j < i; j++) {
            //double Nj = Branch.Ne[j];
            //double tjp1 = Branch.Times[j+1];
            
            double C2bi = 0.0;
            double C2B2 = 0.0;
            C2B2 = GetCase2B2 ( Branch, j, i);
            
            
            C2bi = 4.0/5.0 * Branch.P5P2NC[j] * C2B2;
            intervalSum = intervalSum + C2bi * ProbNoCoal;
            ProbNoCoal = ProbNoCoal * Branch.P5NC[j];
            
        }
        
    }else{ intervalSum = 0.0;}
    
    return intervalSum;
    
    
}


double GetCase2 ( BranchNode Branch){
    
    double intervalSum =0.0;
    double ProbNoCoal = 1.0;
    double P5NoCoal = 1.0;
    
    for (int i =0; i < Branch.Counter; i++) {
        double Ni = Branch.Ne[i];
        
        double C2i = 0.0;
        double C2A = GetCase2A(Branch, i);
        double C2B = GetCase2B(Branch, i);
        
        
        C2i = P5NoCoal * ( -2.0/3.0 * Ni * Branch.TPP2[i] - 2.0/3.0 * Branch.T2PP2[i] + (8.0/75.0 * Ni)* Branch.TPP5[i] + (2.0/3.0-2.0/5.0)* Branch.T2PP5[i] + 2.0/5.0*(7.0/5.0*Ni +Branch.Times[i])*Branch.TP[i]) + Branch.TP[i] *C2A + C2B;
        
        
        intervalSum = intervalSum + C2i * ProbNoCoal;
        ProbNoCoal = ProbNoCoal * Branch.P22[i];
        
        P5NoCoal = P5NoCoal * Branch.P5NC[i];
    }
    //cout << "Case 2 = " << intervalSum << "\n";
    
    return intervalSum;
}

double GetCase3A(BranchNode Branch, int i ){
    
    
    double intervalSum =0.0;
    double ProbNoCoal = 1.0;
    
    if ( i !=0){
        
        for (int j =0; j < i; j++) {
            
            
            double C3ai = 0.0;
            
            C3ai = 2.0/5.0 * (1-Branch.P5NC[j] - Branch.P5P2NC[j]);
            intervalSum = intervalSum + C3ai * ProbNoCoal;
            ProbNoCoal = ProbNoCoal * Branch.P5NC[j];
            
        }
        
    }else{ intervalSum = 0.0;}
    
    
    
    
    return intervalSum;
    
    
}

double GetCase3B2(BranchNode Branch, int j, int i){
    
    double intervalSum =0.0;
    double P2NoCoal = 1.0;
    
    for (int k = j+1; k < i; k++) {
        double C3b2i = 0.0;
        
        C3b2i = 1.0/2.0 *(1 - Branch.P2NC[k]);
        intervalSum = intervalSum + C3b2i * P2NoCoal;
        P2NoCoal = P2NoCoal * Branch.P2NC[k];
        
        
    }
    intervalSum = intervalSum * Branch.T2P[i] + 1.0/2.0 * P2NoCoal *(Branch.T2P[i] - Branch.T2PP2[i]);
    
    
    
    return intervalSum;
    
    
    
}

double GetCase3B( BranchNode Branch, int i ){
    
    
    double intervalSum =0.0;
    double ProbNoCoal = 1.0;
    
    if ( i !=0){
        
        for (int j =0; j < i; j++) {
            //double Nj = Branch.Ne[j];
            //double tjp1 = Branch.Times[j+1];
            
            double C3bi = 0.0;
            double C3B2 = 0.0;
            C3B2 = GetCase3B2( Branch, j, i);
            
            C3bi = 4.0/5.0 * Branch.P5P2NC[j] * C3B2;
            intervalSum = intervalSum + C3bi * ProbNoCoal;
            ProbNoCoal = ProbNoCoal * Branch.P5NC[j];
            
        }
        
    }else{ intervalSum = 0.0;}
    
    
    return intervalSum;
}


double GetCase3(BranchNode Branch){
    
    
    double intervalSum =0.0;
    double ProbNoCoal = 1.0;
    double P5NoCoal = 1.0;
    
    for (int i =0; i < Branch.Counter; i++) {
        // double Ni = Branch.Ne[i];
        // double ti = Branch.Times[i];
        double C3i = 0.0;
        double C3A = GetCase3A(Branch, i);
        double C3B = GetCase3B(Branch, i);
        
        C3i = P5NoCoal * ( 2.0/5.0 * Branch.T2P[i] - 2.0/3.0 * Branch.T2PP2[i] + (2.0/3.0 - 2.0/5.0) * Branch.T2PP5[i]) + Branch.T2P[i] * C3A + C3B;
        
        
        
        intervalSum = intervalSum + C3i * ProbNoCoal;
        ProbNoCoal = ProbNoCoal * Branch.P22[i];
        P5NoCoal = P5NoCoal * Branch.P5NC[i];
    }
    
    //cout << "Case 3 = " << intervalSum << "\n";
    
    return intervalSum;
    
}


double GetTerm1(BranchNode Branch){
    
    double intervalSum = 0.0;
    double ProbNoCoal = 1.0;
    double P5NoCoal = 1.0;
    double T1i = 0;
    
    for(int i = 0; i < Branch.Counter; i++){
        
        
        T1i = P5NoCoal* Branch.TPP5[i];
        
        intervalSum = intervalSum + T1i * ProbNoCoal;
        ProbNoCoal = ProbNoCoal * Branch.P22[i];
        P5NoCoal = P5NoCoal * Branch.P5NC[i];
        //Conveniently P5nocoal * Pnocoal works out the same as p6nocoal.
        
    }
    
    return intervalSum;
}

double GetTerm2_1(BranchNode Branch, int i){
    
    
    // This is the probability that the first coalescence event happens in an interval before the present. Assuming all lineages have survived to final branch.
    
    
    double intervalSum =0.0;
    
    double P5NoCoal  = 1.0 ;
    
    if(i!=0){
        for(int k =0; k< i; k++){
            
            intervalSum = intervalSum + 1.0/5.0 * P5NoCoal * (1.0 - Branch.P5NC[k]);
            P5NoCoal = P5NoCoal* Branch.P5NC[k];
            
            
        }
    }else{
        intervalSum = 0.0;
    }
    
    return intervalSum;
    
}

double GetTerm2(BranchNode Branch){
    
    
    double intervalSum = 0.0;
    double ProbNoCoal = 1.0;
    double P5NoCoal = 1.0;
    double T2_1 = 0;
    
    
    for(int i =0; i < Branch.Counter; i++){
        
        //Probability no other pair( 5 of 6) has coalesced from interval 0 to i on this branch.
        T2_1 = GetTerm2_1(Branch, i);
        
        intervalSum = intervalSum + ProbNoCoal * (Branch.TP[i] * T2_1 + 1.0/5.0 * P5NoCoal *(Branch.TP[i] - Branch.TPP5[i]));
        
        ProbNoCoal = ProbNoCoal * Branch.P22[i];
        P5NoCoal = P5NoCoal * Branch.P5NC[i];
        
        
    }
    
    
    return intervalSum;
}


double GetTerm3_1(BranchNode Branch, int i){
    
    
    
    
    double intervalSum =0.0;
    
    double P2NoCoal  = 1.0 ;
    
    if(i!=0){
        for(int k =0; k< i; k++){
            
            intervalSum = intervalSum + 1.0/2.0 * P2NoCoal * (1.0 - Branch.P2NC[k]);
            P2NoCoal = P2NoCoal* Branch.P2NC[k];
            
            
        }
    }else{
        intervalSum = 0.0;
    }
    
    return intervalSum;
    
}

double GetTerm3(BranchNode Branch){
    
    double intervalSum = 0.0;
    double ProbNoCoal = 1.0;
    double P2NoCoal = 1.0;
    double T3_1;
    
    for(int i =0; i < Branch.Counter; i++){
        T3_1= GetTerm3_1(Branch, i);
        
        intervalSum = intervalSum + ProbNoCoal * (Branch.TP[i] *T3_1 + P2NoCoal* 0.5*(Branch.TP[i] - Branch.TPP2[i]) );
        ProbNoCoal = ProbNoCoal * Branch.P22[i];
        P2NoCoal = P2NoCoal * Branch.P2NC[i];
        
    }
    
    return intervalSum;
}

double GetTerm4(BranchNode Branch){
    
    double intervalSum = 0.0;
    double ProbNoCoal = 1.0;
    double P2NoCoal = 1.0;
    
    for(int i = 0; i < Branch.Counter ; i++){
        
        intervalSum = intervalSum + ProbNoCoal * P2NoCoal * Branch.TPP2[i];
        ProbNoCoal = ProbNoCoal * Branch.P22[i];
        P2NoCoal = P2NoCoal * Branch.P2NC[i];
        
    }
    
    return intervalSum;
    
}

double GetTerm5i(BranchNode Branch, BranchNode BranchOther, int i){
    
    double intervalSum = 0.0;
    double ProbNoCoal = 1.0; // The BranchOther probnocoal.
    double ProbNoCoali = 1.0; // The cur branch no coal.
    double Ni =Branch.Ne[i];
    double Nj = 0;
    
    double TPP22N;
    //First, find the interval where t_i is in BranchOther. This could be an exact match (0)
    double ti = Branch.Times[i];
    int counter_ti = -1;
    while(ti >BranchOther.Times[counter_ti+1]){
        counter_ti = counter_ti+1;
    }
    //if counter_ti = -1, means ti exists before the branchother starts. ti < D_other.
    
    //Second, find the interval where t_{i+1} is in branch other. this could be an exact match. (DMRCA)
    double tip1 = Branch.Times[i+1];
    int counter_tip1 = -1;
    while(tip1 >BranchOther.Times[counter_tip1+1]){
        counter_tip1 = counter_tip1+1;
    }
    //if counter_tip1 = -1, means ti exists before the branchother starts. ti < D_other.
    
    
    if( counter_ti == -1){
        
        // ti occurs before the start of BranchOther
        
        if(counter_tip1 == -1){
            //tip1 also occurs before the start of BranchOther
            
            //Probability 23 dont coal first = 1. Easy.
            intervalSum = intervalSum + Branch.TP[i];
            
            
        }else{
            
            //ti occurs before the start of BranchOther, but tip1 is inside.
            
            //First get the interval time before the start of BranchOther, where if our pair coalesces, it is autmoatically the minimum.
            intervalSum = intervalSum -(2.0 * Branch.Ne[i] + BranchOther.Times[0])* exp(-(BranchOther.Times[0] - Branch.Times[i])/(2.0*Branch.Ne[i])) + 2.0 * Branch.Ne[i] + Branch.Times[i];
            
            ProbNoCoali = ProbNoCoali *exp(-(BranchOther.Times[0] - Branch.Times[i])/(2.0*Branch.Ne[i]));
            //Next, get all terms that are in BranchOther (less than counter_tip1 < tip1)
            
            for(int j = 0; j < counter_tip1; j++){
                
                Nj = BranchOther.Ne[j];
                
                TPP22N = Nj * ( -(BranchOther.Times[j+1]*(Nj+Ni) + 2.0 * Nj*Ni)* exp(-(BranchOther.Times[j+1]- BranchOther.Times[j])*(Nj+ Ni)/(2.0*Nj*Ni)) + BranchOther.Times[j]*(Nj+Ni) + 2.0*Nj*Ni  )/pow((Nj+Ni),2.0);
                
                
                intervalSum = intervalSum + ProbNoCoal *ProbNoCoali* TPP22N;
                ProbNoCoal = ProbNoCoal * BranchOther.P22[j];
                ProbNoCoali = ProbNoCoali *  exp(-(BranchOther.Times[j+1] - BranchOther.Times[j])/(2.0*Ni));
                
            }
            //Lastly, get the interval from counter_tip1 to tip1
            Nj = BranchOther.Ne[counter_tip1];
            TPP22N = Nj * ( -(tip1*(Nj+Ni) + 2.0 * Nj*Ni)* exp(-(tip1- BranchOther.Times[counter_tip1])*(Nj+ Ni)/(2.0*Nj*Ni)) + BranchOther.Times[counter_tip1]*(Nj+Ni) + 2.0*Nj*Ni  )/pow((Nj+Ni),2.0);
            intervalSum = intervalSum + ProbNoCoal *ProbNoCoali* TPP22N;
            //ProbNoCoal = ProbNoCoal * BranchOther.P22[counter_tip1];
            
            
        }
        
        
    }else{
        
        //ti occurs either at the start or within BranchOther.
        
        if( counter_ti == counter_tip1){
            
            //ti is fully contained in a single interval on (otherbranch) (NEEDS TO BE ADDED TO TERM 6. )
            
            if(ti > BranchOther.Times[0]){
                for(int j = 0; j < counter_ti; j++){
                    
                    ProbNoCoal = ProbNoCoal * BranchOther.P22[j];
                    
                }
                ProbNoCoal = ProbNoCoal * exp(-(ti - BranchOther.Times[counter_ti])/(2.0*BranchOther.Ne[counter_ti]));
            }
            
            Nj = BranchOther.Ne[counter_ti];
            
            TPP22N = Nj * ( -(tip1*(Nj+Ni) + 2.0 * Nj*Ni)* exp(-(tip1- ti)*(Nj+ Ni)/(2.0*Nj*Ni)) + ti*(Nj+Ni) + 2.0*Nj*Ni  )/pow((Nj+Ni),2.0);
            intervalSum = intervalSum + ProbNoCoal *ProbNoCoali* TPP22N;
            
            
            
        }else{
            //ti spans at least two j intervals, and is fully contained within otherbranch.
            if( counter_tip1 <= counter_ti){
                
            }
            //Get the probability that other branch pair coal has not occured before ti.
            if(ti > BranchOther.Times[0]){
                for(int j = 0; j < counter_ti; j++){
                    
                    ProbNoCoal = ProbNoCoal * BranchOther.P22[j];
                    
                }
                ProbNoCoal = ProbNoCoal * exp(-(ti - BranchOther.Times[counter_ti])/(2.0*BranchOther.Ne[counter_ti]));
            }
            
            // Now get the term from ti to counter_ti + 1, HEY THIS COULD BE WRONG, WHAT IF ti+1  < time counter_ti+1 ?? I think I corrected this with the above equality check.
            Nj = BranchOther.Ne[counter_ti];
            TPP22N = Nj * ( -(BranchOther.Times[counter_ti+1]*(Nj+Ni) + 2.0 * Nj*Ni)* exp(-(BranchOther.Times[counter_ti+1]- ti)*(Nj+ Ni)/(2.0*Nj*Ni)) + ti*(Nj+Ni) + 2.0*Nj*Ni  )/pow((Nj+Ni),2.0);
            intervalSum = intervalSum + ProbNoCoal *ProbNoCoali* TPP22N;
            ProbNoCoal = ProbNoCoal * exp(-(BranchOther.Times[counter_ti+1] - ti)/(2.0*BranchOther.Ne[counter_ti]));
            ProbNoCoali = ProbNoCoali * exp(-(BranchOther.Times[counter_ti+1] - ti)/(2.0*Ni));
            
            //Now get all of the terms from counter_ti+1 to counter_tip1
            for(int j = counter_ti+1; j < counter_tip1; j++){
                Nj = BranchOther.Ne[j];
                TPP22N = Nj * ( -(BranchOther.Times[j+1]*(Nj+Ni) + 2.0 * Nj*Ni)* exp(-(BranchOther.Times[j+1]- BranchOther.Times[j])*(Nj+ Ni)/(2.0*Nj*Ni)) + BranchOther.Times[j]*(Nj+Ni) + 2.0*Nj*Ni  )/pow((Nj+Ni),2.0);
                intervalSum = intervalSum + ProbNoCoal *ProbNoCoali* TPP22N;
                ProbNoCoal = ProbNoCoal * BranchOther.P22[j];
                ProbNoCoali = ProbNoCoali* exp(-(BranchOther.Times[j+1] - BranchOther.Times[j])/(2.0*Branch.Ne[i]));
                
                
            }
            
            //Lastly, get the interval from counter_tip1 to tip1 (could be length zero but thats okay)
            Nj = BranchOther.Ne[counter_tip1];
            TPP22N = Nj * ( -(tip1*(Nj+Ni) + 2.0 * Nj*Ni)* exp(-(tip1- BranchOther.Times[counter_tip1])*(Nj+ Ni)/(2.0*Nj*Ni)) + BranchOther.Times[counter_tip1]*(Nj+Ni) + 2.0*Nj*Ni  )/pow((Nj+Ni),2.0);
            intervalSum = intervalSum + ProbNoCoal *ProbNoCoali* TPP22N;
            //ProbNoCoal = ProbNoCoal * BranchOther.P22[counter_tip1];
            
        } }
    
    return intervalSum;
}



double GetTerm5(BranchNode Branch, BranchNode BranchOther){
    
    // E(T12|T12min, coal in branch 12) P(coal in branch12 * T12 min)
    
    double intervalSum = 0.0;
    double ProbNoCoal = 1.0;
    double OtherBranchTerm = 0;
    
    for(int i =0; i < Branch.Counter; i++ ){
        OtherBranchTerm = GetTerm5i(Branch, BranchOther, i);
        
        
        intervalSum = intervalSum + ProbNoCoal * OtherBranchTerm;
        ProbNoCoal = ProbNoCoal * Branch.P22[i];
        
    }
    
    return intervalSum;
}



double GetTerm6i(BranchNode Branch, BranchNode BranchOther, int i){
    
    double intervalSum = 0.0;
    double ProbNoCoal = 1.0; // The BranchOther probnocoal.
    double ProbNoCoalPre_ti =1.0; // The BranchOther prob no coal before time ti.
    double ProbNoCoalPost_ti =1.0; // The BranchOther prob of no coal between ti and tj.
    double ProbNoCoali = 1.0; // The cur branch no coal.
    
    double Ni =Branch.Ne[i];
    double Nj = 0;
    double TPP1N;
    double TPP22N; //for use in the TPP1N integral.
    
    
    
    //First, find the interval where t_i is in BranchOther. This could be an exact match (0)
    
    
    double ti = Branch.Times[i];
    int counter_ti = -1;
    
    while(ti > BranchOther.Times[counter_ti+1]){
        
        counter_ti = counter_ti+1;
    }
    
    //if counter_ti = -1, means ti exists before the branchother starts. ti < D_other.
    
    //Second, find the interval where t_{i+1} is in branch other. this could be an exact match. (DMRCA)
    double tip1 = Branch.Times[i+1];
    int counter_tip1 = -1;
    while(tip1 >BranchOther.Times[counter_tip1+1]){
        counter_tip1 = counter_tip1+1;
    }
    //    if(tip1 == BranchOther.Times[counter_tip1 +1]){
    //        counter_tip1 = counter_tip1 +1;
    //    }
    
    //if counter_tip1 = -1, means ti exists before the branchother starts. ti < D_other.
    
    
    if( counter_ti == -1){
        
        // ti occurs before the start of BranchOther
        
        if(counter_tip1 == -1){
            //tip1 also occurs before the start of BranchOther
            
            //Probability 23 dont coal first = 1. So Prob 23 coal first = 0.
            intervalSum = intervalSum + 0.0;
            
        }else{
            
            //ti occurs before the start of BranchOther, but tip1 is inside.
            
            //First get the interval time before the start of BranchOther, where Prob 23 coal first = 0.
            intervalSum = intervalSum +0.0;
            
            ProbNoCoali = ProbNoCoali *exp(-(BranchOther.Times[0] - Branch.Times[i])/(2.0*Ni));
            
            //Next, get all terms that are in BranchOther (less than counter_tip1 < tip1)
            
            for(int j = 0; j < counter_tip1; j++){
                Nj = BranchOther.Ne[j];
                
                //E(t|other no coal before) we dont DIRECTLY use this calc, only indirect.
                TPP22N = Nj * ( -(BranchOther.Times[j+1]*(Nj+Ni) + 2.0 * Nj*Ni)* exp(-(BranchOther.Times[j+1]- BranchOther.Times[j])*(Nj+ Ni)/(2.0*Nj*Ni)) + BranchOther.Times[j]*(Nj+Ni) + 2.0*Nj*Ni  )/pow((Nj+Ni),2.0);
                
                
                TPP1N = -(BranchOther.Times[j+1] + 2*Ni)* exp(-(BranchOther.Times[j+1]- BranchOther.Times[j])/(2*Ni)) + BranchOther.Times[j] + 2*Ni - TPP22N;
                
                intervalSum = intervalSum + ProbNoCoali*(ProbNoCoalPost_ti * (TPP1N) + (1.0 - ProbNoCoalPost_ti)* (TPP1N+TPP22N));
                
                
                
                
                ProbNoCoal = ProbNoCoal * BranchOther.P22[j];
                ProbNoCoalPost_ti = ProbNoCoalPost_ti * BranchOther.P22[j];
                ProbNoCoali = ProbNoCoali* exp(-(BranchOther.Times[j+1] - BranchOther.Times[j])/(2.0*Branch.Ne[i]));
                
                //Since we start at the very base of branch23, dont worry about Pnocoal in (start 23, ti)
                
            }
            
            //Lastly, get the interval from counter_tip1 to tip1
            Nj = BranchOther.Ne[counter_tip1];
            TPP22N = Nj * ( -(tip1*(Nj+Ni) + 2.0 * Nj*Ni)* exp(-(tip1- BranchOther.Times[counter_tip1])*(Nj+ Ni)/(2.0*Nj*Ni)) + BranchOther.Times[counter_tip1]*(Nj+Ni) + 2.0*Nj*Ni  )/pow((Nj+Ni),2.0);
            TPP1N = -(tip1 + 2*Ni)* exp(-(tip1- BranchOther.Times[counter_tip1])/(2*Ni)) + BranchOther.Times[counter_tip1] + 2*Ni - TPP22N;
            
            
            intervalSum = intervalSum + ProbNoCoali*(ProbNoCoal * (TPP1N) + (1.0 - ProbNoCoal)* (TPP1N+TPP22N));
            
            //ProbNoCoal = ProbNoCoal * BranchOther.P22[counter_tip1];
            
            
        }
        
        
    }else{
        
        //ti occurs either at the start or within BranchOther.
        
        if( counter_ti == counter_tip1){
            
            //ti is fully contained in a single interval on (otherbranch)
            
            if(ti > BranchOther.Times[0]){
                for(int j = 0; j < counter_ti; j++){
                    
                    ProbNoCoal = ProbNoCoal * BranchOther.P22[j];
                    
                }
                
                
                ProbNoCoal = ProbNoCoal * exp(-(ti - BranchOther.Times[counter_ti])/(2.0*BranchOther.Ne[counter_ti]));
            }
            
            ProbNoCoalPre_ti = ProbNoCoal;
            
            // Probability the minumim happens before we even get to ti. * E(T|ti,ti+1)
            intervalSum = intervalSum + (1.0 - ProbNoCoalPre_ti) * Branch.TP[i];
            
            
            Nj = BranchOther.Ne[counter_ti];
            
            TPP22N = Nj * ( -(tip1*(Nj+Ni) + 2.0 * Nj*Ni)* exp(-(tip1- ti)*(Nj+ Ni)/(2.0*Nj*Ni)) + ti*(Nj+Ni) + 2.0*Nj*Ni  )/pow((Nj+Ni),2.0);
            TPP1N = -(tip1 + 2.0*Ni)* exp(-(tip1- ti)/(2.0*Ni)) + ti + 2.0*Ni - TPP22N;
            
            intervalSum = intervalSum + ProbNoCoalPre_ti * ProbNoCoalPost_ti * TPP1N + ProbNoCoalPre_ti * (1.0- ProbNoCoalPost_ti) * (TPP1N + TPP22N);
            
            
        }else{
            
            //ti spans at least two j intervals, and is fully contained within otherbranch.
            
            
            //Get the probability that 23 has not occured before ti.
            if(ti > BranchOther.Times[0]){
                for(int j = 0; j < counter_ti; j++){
                    
                    ProbNoCoal = ProbNoCoal * BranchOther.P22[j];
                    
                }
                ProbNoCoal = ProbNoCoal * exp(-(ti - BranchOther.Times[counter_ti])/(2.0*BranchOther.Ne[counter_ti]));
            }
            
            ProbNoCoalPre_ti = ProbNoCoal;
            
            // Term for 23 coal before t_i = P(23 coal before ti) * E(t|(ti,i+1))
            intervalSum = intervalSum + (1.0 - ProbNoCoalPre_ti) * Branch.TP[i];
            
            // Now get the term from ti to counter_ti + 1
            
            Nj = BranchOther.Ne[counter_ti];
            TPP22N = Nj * ( -(BranchOther.Times[counter_ti+1]*(Nj+Ni) + 2.0 * Nj*Ni)* exp(-(BranchOther.Times[counter_ti+1]- ti)*(Nj+ Ni)/(2.0*Nj*Ni)) + ti*(Nj+Ni) + 2.0*Nj*Ni  )/pow((Nj+Ni),2.0);
            TPP1N = -(BranchOther.Times[counter_ti+1] + 2*Ni)* exp(-(BranchOther.Times[counter_ti+1]- ti)/(2*Ni)) + ti + 2*Ni - TPP22N;
            
            intervalSum = intervalSum + ProbNoCoali*(ProbNoCoalPre_ti*ProbNoCoalPost_ti * TPP1N + ProbNoCoalPre_ti * (1.0- ProbNoCoalPost_ti) * (TPP1N + TPP22N));
            
            
            ProbNoCoal = ProbNoCoal * exp(-(BranchOther.Times[counter_ti+1] - ti)/(2.0*BranchOther.Ne[counter_ti]));
            ProbNoCoalPost_ti = ProbNoCoalPost_ti * exp(-(BranchOther.Times[counter_ti+1] - ti)/(2.0*BranchOther.Ne[counter_ti]));
            ProbNoCoali = ProbNoCoali * exp(-(BranchOther.Times[counter_ti+1] - ti)/(2.0*Ni));
            
            
            //Now get all of the terms from counter_ti+1 to counter_tip1
            for(int j = counter_ti+1; j < counter_tip1; j++){
                Nj = BranchOther.Ne[j];
                TPP22N = Nj * ( -(BranchOther.Times[j+1]*(Nj+Ni) + 2.0 * Nj*Ni)* exp(-(BranchOther.Times[j+1]- BranchOther.Times[j])*(Nj+ Ni)/(2.0*Nj*Ni)) + BranchOther.Times[j]*(Nj+Ni) + 2.0*Nj*Ni  )/pow((Nj+Ni),2.0);
                TPP1N = -(BranchOther.Times[j+1] + 2*Ni)* exp(-(BranchOther.Times[j+1]- BranchOther.Times[j])/(2*Ni)) + BranchOther.Times[j] + 2*Ni - TPP22N;
                
                intervalSum = intervalSum + ProbNoCoali*(ProbNoCoalPre_ti*ProbNoCoalPost_ti * TPP1N + ProbNoCoalPre_ti * (1.0- ProbNoCoalPost_ti) * (TPP1N + TPP22N));
                
                
                ProbNoCoal = ProbNoCoal * BranchOther.P22[j];
                ProbNoCoalPost_ti = ProbNoCoalPost_ti * BranchOther.P22[j];
                ProbNoCoali = ProbNoCoali* exp(-(BranchOther.Times[j+1] - BranchOther.Times[j])/(2.0*Branch.Ne[i]));
                
                
            }
            
            //Lastly, get the interval from counter_tip1 to tip1 (could be length zero but thats okay)
            Nj = BranchOther.Ne[counter_tip1];
            TPP22N = Nj * ( -(tip1*(Nj+Ni) + 2.0 * Nj*Ni)* exp(-(tip1- BranchOther.Times[counter_tip1])*(Nj+ Ni)/(2.0*Nj*Ni)) + BranchOther.Times[counter_tip1]*(Nj+Ni) + 2.0*Nj*Ni  )/pow((Nj+Ni),2.0);
            TPP1N = -(tip1 + 2*Ni)* exp(-(tip1- BranchOther.Times[counter_tip1])/(2*Ni)) + BranchOther.Times[counter_tip1] + 2*Ni - TPP22N;
            
            intervalSum = intervalSum + ProbNoCoali*(ProbNoCoalPre_ti*ProbNoCoalPost_ti * TPP1N + ProbNoCoalPre_ti * (1.0- ProbNoCoalPost_ti) * (TPP1N + TPP22N));
            
        }}
    
    return intervalSum;
}



double GetTerm6(BranchNode Branch, BranchNode BranchOther){
    
    // E(T12|T34 min, coal in branch 12) P(coal in branch12 * T34 min)
    
    double intervalSum = 0.0;
    double ProbNoCoal = 1.0;
    double OtherBranchTerm = 0;
    
    for(int i =0; i < Branch.Counter; i++ ){
        OtherBranchTerm = GetTerm6i(Branch, BranchOther, i);
        
        intervalSum = intervalSum + ProbNoCoal * OtherBranchTerm;
        ProbNoCoal = ProbNoCoal * Branch.P22[i];
        
    }
    
    return intervalSum;
}

double GetTerm7(BranchNode Branch){
    
    //Asym Case:  E(T02|T02 min, 01 no coal in4) P(w.e.)
    
    double intervalSum = 0.0;
    double ProbNoCoal = 1.0;
    double P2pNoCoal = 1.0;
    
    for(int i =0; i < Branch.Counter; i++ ){
        
        
        intervalSum = intervalSum + ProbNoCoal *P2pNoCoal *Branch.TPP2[i];
        ProbNoCoal = ProbNoCoal * Branch.P22[i];
        P2pNoCoal = P2pNoCoal * Branch.P2NC[i];
    }
    
    return intervalSum;
}

double GetTerm8(BranchNode Branch){
    
    //Asym Case:  E(T02|T02 in branch)
    
    double intervalSum = 0.0;
    double ProbNoCoal = 1.0;
    
    for(int i =0; i < Branch.Counter; i++ ){
        
        
        intervalSum = intervalSum + ProbNoCoal * Branch.TP[i];
        ProbNoCoal = ProbNoCoal * Branch.P22[i];
    }
    
    return intervalSum;
}


void TreeExtender(SpeciesTreeType SpeciesTree, SpeciesTreeType &ExtendedTree, int n, int n_indivs, Matrix<int, Dynamic, 2> SpeciesMap){
    
    // For use when Theory Function (with covariance) is needed
    //The extended tree will have the number of branches as if there were n_indivs species. This is 2*n_indivs -1.0 as opposed to 2*n-1.
    
    //If a specific species has two or more individuals, we will utilize the zero branch length tactic.
    
    //The first n_indivs branches will all be of zero length, population size shouldnt matter.
    
    //int NextBranch = n_indivs;
    
    //We will use the variable NextBranch, to assign labels to the internal 0 length branches if there are 3 or more indivs in a species.
    
    
    VectorXi mappedIndivs;
    int Branchdif = 2.0 *(n_indivs - n);
    
    //Assign the top of the tree
    ExtendedTree.Branches[2*n-2 +Branchdif].Ne = SpeciesTree.Branches[2*n-2].Ne;
    ExtendedTree.Branches[2*n-2 +Branchdif].Times[1] = SpeciesTree.Branches[2*n-2].Times[1];
    ExtendedTree.Branches[2*n-2 +Branchdif].Times[0] = SpeciesTree.Branches[2*n-2].Times[0];
    ExtendedTree.Branches[2*n-2 +Branchdif].Parent = 0;
    ExtendedTree.Branches[2*n-2 +Branchdif].Offspring[0] = SpeciesTree.Branches[2*n-2].Offspring[0]+Branchdif;
    ExtendedTree.Branches[2*n-2 +Branchdif].Offspring[1] = SpeciesTree.Branches[2*n-2].Offspring[1]+Branchdif;
    
    // cout<<"Top Branch = "  << 2*n-2 +Branchdif <<"\n";
    // Quick loop through branches n -> 2n-3 and assigning their parameters.
    for(int b = n; b < 2*n-2; b++){
        ExtendedTree.Branches[b +Branchdif].Ne = SpeciesTree.Branches[b].Ne;
        ExtendedTree.Branches[b +Branchdif].Times[1] = SpeciesTree.Branches[b].Times[1];
        ExtendedTree.Branches[b +Branchdif].Times[0] = SpeciesTree.Branches[b].Times[0];
        ExtendedTree.Branches[b +Branchdif].Parent = SpeciesTree.Branches[b].Parent+Branchdif;
        ExtendedTree.Branches[b +Branchdif].Offspring[0] = SpeciesTree.Branches[b].Offspring[0]+Branchdif;
        ExtendedTree.Branches[b +Branchdif].Offspring[1] = SpeciesTree.Branches[b].Offspring[1]+Branchdif;
        // cout<<" Internal Branch = "<< b <<" = "  << b +Branchdif <<"\n";
        
    }
    
    
    // Let's loop through the leaf branches, and subtend the appropriate branches to that.
    //cout <<"Branchdif = " <<Branchdif <<"\n";
    int NewBranch;
    int CurIndiv =0;
    int Cur_NewBranch = Branchdif;
    for(int b = 0; b < n; b++){
        VectorXi mappedIndivs;
        
        for(int x=0; x < SpeciesMap.rows(); x++){
            if(SpeciesMap(x,1)== b){
                
                mappedIndivs.conservativeResize(mappedIndivs.size()+1);
                mappedIndivs(mappedIndivs.size()-1) = SpeciesMap(x,0);
                //cout <<C1indivs(C1indivs.size()-1) <<"\n";
            }}
        
        int n_NewBranches = mappedIndivs.size()-1.0;
        
        NewBranch = Branchdif+ b;
        
        //b represents the current species.
        //NewBranch = b + Branchdif; // This represents the original leaf label.
        //This is NOT the right labelling.
        
        //Get the appropiate label by equation 2*n_indivs-1 = 2*n-1
        ExtendedTree.Branches[NewBranch].Ne = SpeciesTree.Branches[b].Ne;
        ExtendedTree.Branches[NewBranch].Times[0] = SpeciesTree.Branches[b].Times[0];
        ExtendedTree.Branches[NewBranch].Times[1] = SpeciesTree.Branches[b].Times[1];
        ExtendedTree.Branches[NewBranch].Parent = SpeciesTree.Branches[b].Parent + Branchdif;
        //cout<<" Leaf Branch = "<< b <<" = "  << NewBranch <<"\n";
        
        
        
        //cout <<"Species " << b << " indivs = " << mappedIndivs.transpose() <<"\n";
        //Query how many branches will be needed. (interior besides the actual num indivs)
        //int Cur_NewBranch = NewBranch; // this DEFINITION of this variable is not used.
        int Cur_Parent = NewBranch;
        //int My_Cur_NewBranch = Cur_NewBranch;
        //cout <<"Cur New Branch = " << Cur_NewBranch <<"\n";
        int new_n_NewBranches;
        if(mappedIndivs.size() > 1){
            
            while(mappedIndivs.size() >1){//New
                //If there is more than 1 individual
                new_n_NewBranches = mappedIndivs.size() -1.0;
                CurIndiv = mappedIndivs(mappedIndivs.size()-1);
                //cout<<" Cur indiv = " << CurIndiv <<"\n";
                //Assign the last indiv in the list to a branch.
                ExtendedTree.Branches[CurIndiv].Parent = Cur_Parent;
                ExtendedTree.Branches[Cur_Parent].Offspring[0] = CurIndiv;
                ExtendedTree.Branches[CurIndiv].Ne = 1.0;
                ExtendedTree.Branches[CurIndiv].Times[0]= 0.0;
                ExtendedTree.Branches[CurIndiv].Times[1] = 0.0;
                
                //Remove the last individual.
                mappedIndivs.conservativeResize(mappedIndivs.size()-1);
                if(mappedIndivs.size() ==1){
                    //if there is only one species left, assign it directly as the second child of the Cur_Parent.
                    CurIndiv = mappedIndivs(mappedIndivs.size()-1);
                    //cout<<" Cur indiv last = " << CurIndiv <<"\n";
                    ExtendedTree.Branches[CurIndiv].Parent = Cur_Parent;
                    ExtendedTree.Branches[Cur_Parent].Offspring[1] = CurIndiv;
                    ExtendedTree.Branches[CurIndiv].Ne = 1.0;
                    ExtendedTree.Branches[CurIndiv].Times[0]= 0.0;
                    ExtendedTree.Branches[CurIndiv].Times[1] = 0.0;
                    
                }else{
                    //if there is more than 1 individual remaining, add on a new branch that they will build off of.
                    Cur_NewBranch = Cur_NewBranch- 1.0;
                    //cout <<"Cur New Branch = "<< Cur_NewBranch <<"\n";
                    ExtendedTree.Branches[Cur_NewBranch].Parent = Cur_Parent;
                    ExtendedTree.Branches[Cur_Parent].Offspring[1] = Cur_NewBranch;
                    ExtendedTree.Branches[Cur_NewBranch].Ne = 1.0;
                    ExtendedTree.Branches[Cur_NewBranch].Times[0]= 0.0;
                    ExtendedTree.Branches[Cur_NewBranch].Times[1] = 0.0;
                    
                    Cur_Parent = Cur_NewBranch; //New
                }
                
                
            }//New
            
            //NewBranch = NewBranch + n_NewBranches;
            
            
            
        }//else do nothing
        
        
    }
    
    
}


void TheoryFunction(SpeciesTreeType SpeciesTree, Matrix<double , Dynamic, 1 > &TheoryMean,  Matrix<double, Dynamic, Dynamic> &TheoryCov, int n, Matrix<double, Dynamic, Dynamic> &TheorySharedBranchLength){
    
    // 012919 This does not work for multiple individuals in a species. Need to restructure this whole thing to do that.
    
    //changed theorymean to dynamic 4/17
    int nc2 = n*(n-1)/2.0;
    
    //int nc2i = n_indivs*(n_indivs-1)/2.0;
    
    //Temporary until I add it into the function call
    
    
    //For each individual, retrieve their history chain (parental ancestors) Requires indivs labeled 0 - (n-1)
    for(int i = 0; i < n; i++){
        int k = i;
        int counter = 0;
        
        //cout << "Ancestral Lineage of indiv " << i << " = \n" ; // CAN DELETE
        
        while(SpeciesTree.Branches[k].Parent !=0) {
            //SpeciesTree.Branches[i].Ancestors.resize(counter+1);
            
            // NOTE! THIS SETTING REQUIRES THE LAST BRANCH PARENT TO EQUAL 0.
            //SpeciesTree.Branches[i].Ancestors(counter) = k = SpeciesTree.Branches[k].Parent;
            k = SpeciesTree.Branches[k].Parent;
            //cout << k <<", \n";
            SpeciesTree.Branches[i].Ancestors.push_back(k);
            //cout << SpeciesTree.Branches[i].Ancestors(counter) << ", ";  // CAN DELETE
            counter = counter +1;
            //k = SpeciesTree.Branches[k].Parent;
            
            
        }
        
        // Note the Ancestral Lineage is now of each species, not individual.
        
        //cout <<"\n"; // CAN DELETE
        
    }
    
    
    int num_quartet = 0;
    //LOOP HERE over all set of 4 indivs.
    for (int i1 =0; i1< (n-1); i1 ++) {
        for(int i2 = i1+1; i2 <(n-1); i2++){
            for(int i3 = i2+1; i3 <(n-1);i3++){
                for(int i4 = i3+1; i4 < n; i4 ++){
                    
                    //cout << "Indivs: " << i1 << ", " << i2 << ", " << i3 <<", "<<i4 <<"\n";
                    
                    
                    // Given 4 individuals, get quartet tree.
                    int indivs [4] = {i1,i2,i3,i4};
                    
                    Quartet FourTree [4];
                    for (int m =0; m < 4; m++) {
                        FourTree[m].indiv=indivs[m];
                        //cout << "indiv = "<< FourTree[m].indiv << "\n";
                        double lenAnc = SpeciesTree.Branches[indivs[m]].Ancestors.size();
                        FourTree[m].Ancestry.clear();
                        for (int jh = 0 ; jh < lenAnc; jh++) {
                            
                            FourTree[m].Ancestry.insert(SpeciesTree.Branches[indivs[m]].Ancestors[jh]);
                            //FourTree[m].Ancestry[jh] = SpeciesTree.Branches[indivs[m]].Ancestors(jh);
                        }
                        //cout << "lenAnc = " << lenAnc << ", n = "<< n << "\n";
                        //FourTree[m].Ancestry = {SpeciesTree.Branches[indivs[m]].Ancestors(0),SpeciesTree.Branches[indivs[m]].Ancestors(0) + lenAnc} ; // OLD, doesn't work for small amounts of n_genes, surprisingly....
                    }
                    
                    
                    BranchNode Branches [7] ;
                    
                    //cout << "indiv 1 ancestry : \n";
                    ////print indiv 1 ancestry
                    //set<int>:: iterator itt;
                    //for (itt=  FourTree[1].Ancestry.begin(); itt != FourTree[1].Ancestry.end(); ++itt) {
                    //int ans = *itt;
                    //cout << ans<< "\n";
                    //}
                    
                    // Determine the branches common in all 4.
                    
                    // Turn everything into a set? Maybe not be the most memory efficient.
                    
                    
                    std::vector<int> s01(n);
                    std::vector<int> s02(n);
                    std::vector<int> s03(n);
                    
                    
                    std::vector<int> c01(n);
                    std::vector<int> c02(n);
                    std::vector<int> c03(n);
                    std::vector<int> c12(n);
                    std::vector<int> c13(n);
                    std::vector<int> c23(n);
                    
                    
                    std::vector<int>::iterator vit1;
                    std::vector<int>::iterator vit2;
                    std::vector<int>::iterator vit3;
                    
                    // std::set<int>::iterator it1;
                    
                    std::vector<int>::iterator vit4;
                    std::vector<int>::iterator vit5;
                    std::vector<int>::iterator vit6;
                    std::vector<int>::iterator vit7;
                    std::vector<int>::iterator vit8;
                    std::vector<int>::iterator vit9;
                    
                    
                    
                    
                    
                    vit1 = std::set_intersection (FourTree[0].Ancestry.begin(), FourTree[0].Ancestry.end(), FourTree[1].Ancestry.begin(), FourTree[1].Ancestry.end(), s01.begin());
                    s01.resize(vit1-s01.begin());
                    
                    vit2 = std::set_intersection (s01.begin(), s01.end(), FourTree[2].Ancestry.begin(), FourTree[2].Ancestry.end(), s02.begin());
                    s02.resize(vit2-s02.begin());
                    
                    vit3 = std::set_intersection (s02.begin(), s02.end(), FourTree[3].Ancestry.begin(), FourTree[3].Ancestry.end(), s03.begin());
                    s03.resize(vit3-s03.begin());
                    
                    
                    
                    // So s03 contains the information needed to build branch 6. Just ignore the first element, which should always be zero.
                    
                    
                    
                    Branches[6].Counter = s03.size() - 1.0 + 1.0;
                    Branches[6].Times[s03.size() - 1 + 1] = 400.0;
                    
                    for( int i = 0; i < s03.size(); i++){
                        Branches[6].Ne[i] = SpeciesTree.Branches[s03[i]].Ne;
                        Branches[6].Times[i] = SpeciesTree.Branches[s03[i]].Times[0];
                    }
                    //                    // Now to get branches 4 and 5. What do we need?
                    //                    cout << "S01 = \n";
                    //                    for(int i =0; i < s01.size(); i++){
                    //                        cout << s01[i] <<"\n";
                    //                    }
                    //                    cout << "S02 = \n";
                    //                    for(int i =0; i < s02.size(); i++){
                    //                        cout << s02[i] <<"\n";
                    //                    }
                    //
                    //                    cout << "S03 = \n";
                    //                    cout <<" s03 size = "<< s03.size() <<"\n";
                    
                    
                    // copy the 0th set ancestry
                    std::set<int> copiedset = FourTree[0].Ancestry;
                    
                    // Remove all elements which map to branch 6.
                    for(int i = 0; i < s03.size();i++){
                        //cout << s03[i] <<"\n";
                        FourTree[0].Ancestry.erase (s03[i]);
                        FourTree[1].Ancestry.erase (s03[i]);
                        FourTree[2].Ancestry.erase (s03[i]);
                        FourTree[3].Ancestry.erase(s03[i]);
                    }
                    // Sometimes s03 has or doesnt have 0, remove zero.
                    FourTree[0].Ancestry.erase(0);
                    FourTree[1].Ancestry.erase(0);
                    FourTree[2].Ancestry.erase(0);
                    FourTree[3].Ancestry.erase(0);
                    
                    
                    
                    //                    //print s03 which is the things that will go in branch 6
                    //                    cout <<" branch 6 components: \n";
                    //                    for (int i = 0; i <s03.size(); i++){
                    //                        cout<< s03[i] <<"\n";
                    //                    }
                    //
                    //                   // check
                    //                    cout <<"Here's whats left of the 2nd Ancestry list \n";
                    //                    for (it1 = FourTree[2].Ancestry.begin(); it1 !=FourTree[2].Ancestry.end(); ++it1) {
                    //                        cout <<*it1 <<"\n";
                    //                    }
                    
                    
                    // Now get the set similarities between all 4 sets (so 6 PAIRWISE comparisons).
                    vit4 = std::set_intersection(FourTree[0].Ancestry.begin(), FourTree[0].Ancestry.end(), FourTree[1].Ancestry.begin(), FourTree[1].Ancestry.end(), c01.begin());
                    c01.resize(vit4-c01.begin());
                    
                    vit5 = std::set_intersection(FourTree[0].Ancestry.begin(), FourTree[0].Ancestry.end(), FourTree[2].Ancestry.begin(), FourTree[2].Ancestry.end(), c02.begin());
                    c02.resize(vit5-c02.begin());
                    
                    vit6 = std::set_intersection(FourTree[0].Ancestry.begin(), FourTree[0].Ancestry.end(), FourTree[3].Ancestry.begin(), FourTree[3].Ancestry.end(), c03.begin());
                    c03.resize(vit6-c03.begin());
                    
                    vit7 = std::set_intersection(FourTree[1].Ancestry.begin(), FourTree[1].Ancestry.end(), FourTree[2].Ancestry.begin(), FourTree[2].Ancestry.end(), c12.begin());
                    c12.resize(vit7-c12.begin());
                    
                    vit8 = std::set_intersection(FourTree[1].Ancestry.begin(), FourTree[1].Ancestry.end(), FourTree[3].Ancestry.begin(), FourTree[3].Ancestry.end(), c13.begin());
                    c13.resize(vit8-c13.begin());
                    
                    vit9 = std::set_intersection(FourTree[2].Ancestry.begin(), FourTree[2].Ancestry.end(), FourTree[3].Ancestry.begin(), FourTree[3].Ancestry.end(), c23.begin());
                    c23.resize(vit9-c23.begin());
                    
                    //                    cout << "c01 size = "<< c01.size() <<"\n";
                    //                    cout << "c01 = " << c01[0] <<c01[1] << "\n";
                    //
                    //                    cout << "c02 size = "<< c02.size()<<"\n";
                    //                    cout << "c02 = " << c02[0] << "\n";
                    //
                    //                    cout << "c03 size = "<< c03.size()<<"\n";
                    //                    cout << "c03 = " << c03[0] <<"\n";
                    //                    cout << "c12 size = "<< c12.size()<<"\n";
                    //                    cout << "c12 = " << c12[0] <<c12[1] << "\n";
                    //
                    //                    cout << "c13 size = "<< c13.size()<<"\n";
                    //                    cout << "c13 = " << c13[0] <<c13[1] << "\n";
                    //                    cout << "c23 size = "<< c23.size()<<"\n";
                    //                    cout << "c23 = " << c23[0] <<c23[1] << "\n";
                    //
                    //
                    // Count the number of empty comparisons. (Should expect 3 for asym tree, 4 for sym).
                    int NumEmpty = c01.empty() + c02.empty()+c03.empty()+c12.empty()+c13.empty()+c23.empty();
                    // cout << "Number of empty pairwise comparisons = "<< NumEmpty <<"\n";
                    
                    int QuartMap [4];
                    
                    // TO DO: Now that we have num empty, we know the type of tree we need to build, so build it in each of the two cases.
                    
                    if(NumEmpty == 3){
                        
                        // We are in the asymmetric Case
                        Branches[6].Offspring[0] = 3;
                        Branches[6].Offspring[1] = 5;
                        Branches[5].Offspring[0] = 2;
                        Branches[5].Offspring[1] = 4;
                        Branches[4].Offspring[0] = 0;
                        Branches[4].Offspring[1] = 1;
                        Branches[4].Parent = 5;
                        Branches[5].Parent = 6;
                        
                        
                        
                        // First determine the individual
                        
                        
                        
                        // so 3 checks to do, if (1,2,3) empty, assign 0 to branch 5. if (1,4,5) then 1, if (2,4,6) then 2, and if (3,4,6) then 3.
                        
                        //If 0 is the last joining indiv
                        if( c01.empty() + c02.empty() + c03.empty() == 3){
                            
                            // Then set the ZEROTH person to be mapped to the THIRD in the tree.
                            
                            QuartMap[3] = indivs[0];
                            //Then determine which intersection of 12, 13, 23, contains the most.
                            
                            
                            if(c12.size() > c13.size() && c12.size() > c23.size()){
                                // if c12 is the longest, then they make up branch 4, and third indiv will be mapped to second.
                                QuartMap[2] = indivs[3];
                                QuartMap[0] = indivs[1];
                                QuartMap[1] = indivs[2];
                                
                                
                                // 12 - 13 is the set we want for branch 4, and  13 is the branch 5.
                                Branches[5].Counter = int(c13.size());
                                Branches[5].Times[c13.size()] = Branches[6].Times[0];
                                for( int i = 0; i < c13.size(); i++){
                                    Branches[5].Ne[i] = SpeciesTree.Branches[c13[i]].Ne;
                                    Branches[5].Times[i] = SpeciesTree.Branches[c13[i]].Times[0];
                                    // cout << Branches[5].Ne[i] <<"\n";
                                    //cout <<Branches[5].Times[i] <<"\n";
                                }
                                
                                // Now fill branch 4.
                                Branches[4].Counter = int(c12.size()) - int(c13.size());
                                Branches[4].Times[Branches[4].Counter] = Branches[5].Times[0];
                                for(int i =0; i < Branches[4].Counter; i++){
                                    Branches[4].Ne[i] = SpeciesTree.Branches[c12[i]].Ne;
                                    Branches[4].Times[i] = SpeciesTree.Branches[c12[i]].Times[0];
                                }
                                
                                
                                
                            }
                            
                            //            if c13.size() is the largest
                            if(c13.size() > c12.size() && c13.size() > c23.size()){
                                // if c12 is the longest, then they make up branch 4, and third indiv will be mapped to second.
                                QuartMap[2] = indivs[2];
                                QuartMap[0] = indivs[1];
                                QuartMap[1] = indivs[3];
                                
                                
                                // 13-12 is the set we want for branch 4, and  12 is the branch 5.
                                Branches[5].Counter = int(c12.size());
                                Branches[5].Times[c12.size()] = Branches[6].Times[0];
                                for( int i = 0; i < c12.size(); i++){
                                    Branches[5].Ne[i] = SpeciesTree.Branches[c12[i]].Ne;
                                    Branches[5].Times[i] = SpeciesTree.Branches[c12[i]].Times[0];
                                    // cout << Branches[5].Ne[i] <<"\n";
                                    // cout <<Branches[5].Times[i] <<"\n";
                                }
                                
                                // Now fill branch 4.
                                Branches[4].Counter = int(c13.size()) - int(c12.size());
                                Branches[4].Times[Branches[4].Counter] = Branches[5].Times[0];
                                for(int i =0; i < Branches[4].Counter; i++){
                                    Branches[4].Ne[i] = SpeciesTree.Branches[c13[i]].Ne;
                                    Branches[4].Times[i] = SpeciesTree.Branches[c13[i]].Times[0];
                                }
                                
                                
                                
                            }
                            
                            
                            
                            //            if(c23.size() == maxsize){
                            if(c23.size() > c13.size() && c23.size() > c12.size()){
                                // if c12 is the longest, then they make up branch 4, and third indiv will be mapped to second.
                                QuartMap[2] = indivs[1];
                                QuartMap[0] = indivs[2];
                                QuartMap[1] = indivs[3];
                                
                                
                                // 23-13 is the set we want for branch 4, and  13 is the branch 5.
                                Branches[5].Counter = int(c13.size());
                                Branches[5].Times[c13.size()] = Branches[6].Times[0];
                                for( int i = 0; i < c13.size(); i++){
                                    Branches[5].Ne[i] = SpeciesTree.Branches[c13[i]].Ne;
                                    Branches[5].Times[i] = SpeciesTree.Branches[c13[i]].Times[0];
                                    // cout << Branches[5].Ne[i] <<"\n";
                                    // cout <<Branches[5].Times[i] <<"\n";
                                }
                                
                                // Now fill branch 4.
                                Branches[4].Counter = int(c23.size()) - int(c13.size());
                                Branches[4].Times[Branches[4].Counter] = Branches[5].Times[0];
                                for(int i =0; i < Branches[4].Counter; i++){
                                    Branches[4].Ne[i] = SpeciesTree.Branches[c23[i]].Ne;
                                    Branches[4].Times[i] = SpeciesTree.Branches[c23[i]].Times[0];
                                }
                                
                                
                                
                            }
                            
                            
                        }
                        // if 1 is the last joining indiv
                        if(c01.empty() + c12.empty() + c13.empty() == 3){
                            
                            //Then set the FIRST element person to be mapped to the THIRD in the tree.
                            QuartMap[3] = indivs[1];
                            
                            // Then assign the pairs 02 03 23.
                            
                            // if 02 is the longest.
                            if(c02.size() > c03.size() && c02.size() > c23.size()){
                                QuartMap[2] = indivs[3];
                                QuartMap[0] = indivs[0];
                                QuartMap[1] = indivs[2];
                                
                                // 02 - 03 is the set we want for branch 4, and  03 is the branch 5.
                                Branches[5].Counter = int(c03.size());
                                Branches[5].Times[c03.size()] = Branches[6].Times[0];
                                for( int i = 0; i < c03.size(); i++){
                                    Branches[5].Ne[i] = SpeciesTree.Branches[c03[i]].Ne;
                                    Branches[5].Times[i] = SpeciesTree.Branches[c03[i]].Times[0];
                                    //  cout << Branches[5].Ne[i] <<"\n";
                                    //  cout <<Branches[5].Times[i] <<"\n";
                                }
                                
                                // Now fill branch 4.
                                Branches[4].Counter = int(c02.size()) - int(c03.size());
                                Branches[4].Times[Branches[4].Counter] = Branches[5].Times[0];
                                for(int i =0; i < Branches[4].Counter; i++){
                                    Branches[4].Ne[i] = SpeciesTree.Branches[c02[i]].Ne;
                                    Branches[4].Times[i] = SpeciesTree.Branches[c02[i]].Times[0];
                                }
                            } //close if 02 is longest
                            
                            // if 03 is the longest.
                            if(c03.size() > c02.size() && c03.size() > c23.size()){
                                QuartMap[2] = indivs[2];
                                QuartMap[0] = indivs[0];
                                QuartMap[1] = indivs[3];
                                
                                // 03 - 02 is the set we want for branch 4, and  02 is the branch 5.
                                Branches[5].Counter = int(c02.size());
                                Branches[5].Times[c02.size()] = Branches[6].Times[0];
                                for( int i = 0; i < c02.size(); i++){
                                    Branches[5].Ne[i] = SpeciesTree.Branches[c02[i]].Ne;
                                    Branches[5].Times[i] = SpeciesTree.Branches[c02[i]].Times[0];
                                    // cout << Branches[5].Ne[i] <<"\n";
                                    // cout <<Branches[5].Times[i] <<"\n";
                                }
                                
                                // Now fill branch 4.
                                Branches[4].Counter = int(c03.size()) - int(c02.size());
                                Branches[4].Times[Branches[4].Counter] = Branches[5].Times[0];
                                for(int i =0; i < Branches[4].Counter; i++){
                                    Branches[4].Ne[i] = SpeciesTree.Branches[c03[i]].Ne;
                                    Branches[4].Times[i] = SpeciesTree.Branches[c03[i]].Times[0];
                                }
                            } //close if 03 is longest
                            
                            // if 23 is the longest.
                            if(c23.size() > c02.size() && c23.size() > c03.size()){
                                QuartMap[2] = indivs[0];
                                QuartMap[0] = indivs[2];
                                QuartMap[1] = indivs[3];
                                
                                // 23 - 02 is the set we want for branch 4, and  02 is the branch 5.
                                Branches[5].Counter = int(c02.size());
                                Branches[5].Times[c02.size()] = Branches[6].Times[0];
                                for( int i = 0; i < c02.size(); i++){
                                    Branches[5].Ne[i] = SpeciesTree.Branches[c02[i]].Ne;
                                    Branches[5].Times[i] = SpeciesTree.Branches[c02[i]].Times[0];
                                    // cout << Branches[5].Ne[i] <<"\n";
                                    // cout <<Branches[5].Times[i] <<"\n";
                                }
                                
                                // Now fill branch 4.
                                Branches[4].Counter = int(c23.size()) - int(c02.size());
                                Branches[4].Times[Branches[4].Counter] = Branches[5].Times[0];
                                for(int i =0; i < Branches[4].Counter; i++){
                                    Branches[4].Ne[i] = SpeciesTree.Branches[c23[i]].Ne;
                                    Branches[4].Times[i] = SpeciesTree.Branches[c23[i]].Times[0];
                                }
                            } //close if 23 is longest
                            
                            
                        } // close if 1 is last joining indiv
                        
                        //if 2 is the last joining indiv
                        if(c02.empty() + c12.empty() + c23.empty() == 3){
                            //Then set the SECOND element person to be mapped to the THIRD in the tree.
                            QuartMap[3] = indivs[2];
                            // Then determine which 01 03 13 is the longest.
                            
                            //if 01 is longest
                            if(c01.size() > c13.size() && c01.size() > c03.size()){
                                // if c01 is the longest, then they make up branch 4, and third indiv will be mapped to second.
                                QuartMap[2] = indivs[3];
                                QuartMap[0] = indivs[0];
                                QuartMap[1] = indivs[1];
                                
                                
                                // 01 - 13 is the set we want for branch 4, and  13 is the branch 5.
                                Branches[5].Counter = int(c13.size());
                                Branches[5].Times[c13.size()] = Branches[6].Times[0];
                                for( int i = 0; i < c13.size(); i++){
                                    Branches[5].Ne[i] = SpeciesTree.Branches[c13[i]].Ne;
                                    Branches[5].Times[i] = SpeciesTree.Branches[c13[i]].Times[0];
                                    // cout << Branches[5].Ne[i] <<"\n";
                                    // cout <<Branches[5].Times[i] <<"\n";
                                }
                                
                                // Now fill branch 4.
                                Branches[4].Counter = int(c01.size()) - int(c13.size());
                                Branches[4].Times[Branches[4].Counter] = Branches[5].Times[0];
                                for(int i =0; i < Branches[4].Counter; i++){
                                    Branches[4].Ne[i] = SpeciesTree.Branches[c01[i]].Ne;
                                    Branches[4].Times[i] = SpeciesTree.Branches[c01[i]].Times[0];
                                }
                                
                            }// close if 01 longest
                            
                            //if 03 is longest
                            if(c03.size() > c13.size() && c03.size() > c01.size()){
                                // if c01 is the longest, then they make up branch 4, and third indiv will be mapped to second.
                                QuartMap[2] = indivs[1];
                                QuartMap[0] = indivs[0];
                                QuartMap[1] = indivs[3];
                                
                                
                                // 03 - 13 is the set we want for branch 4, and  13 is the branch 5.
                                Branches[5].Counter = int(c13.size());
                                Branches[5].Times[c13.size()] = Branches[6].Times[0];
                                for( int i = 0; i < c13.size(); i++){
                                    Branches[5].Ne[i] = SpeciesTree.Branches[c13[i]].Ne;
                                    Branches[5].Times[i] = SpeciesTree.Branches[c13[i]].Times[0];
                                    //  cout << Branches[5].Ne[i] <<"\n";
                                    // cout <<Branches[5].Times[i] <<"\n";
                                }
                                
                                // Now fill branch 4.
                                Branches[4].Counter = int(c03.size()) - int(c13.size());
                                Branches[4].Times[Branches[4].Counter] = Branches[5].Times[0];
                                for(int i =0; i < Branches[4].Counter; i++){
                                    Branches[4].Ne[i] = SpeciesTree.Branches[c03[i]].Ne;
                                    Branches[4].Times[i] = SpeciesTree.Branches[c03[i]].Times[0];
                                }
                                
                            }// close if 03 longest
                            
                            //if 13 is longest
                            if(c13.size() > c03.size() && c13.size() > c01.size()){
                                // if c01 is the longest, then they make up branch 4, and third indiv will be mapped to second.
                                QuartMap[2] = indivs[0];
                                QuartMap[0] = indivs[1];
                                QuartMap[1] = indivs[3];
                                
                                
                                // 13 - 03 is the set we want for branch 4, and  03 is the branch 5.
                                Branches[5].Counter = int(c03.size());
                                Branches[5].Times[c03.size()] = Branches[6].Times[0];
                                for( int i = 0; i < c03.size(); i++){
                                    Branches[5].Ne[i] = SpeciesTree.Branches[c03[i]].Ne;
                                    Branches[5].Times[i] = SpeciesTree.Branches[c03[i]].Times[0];
                                    // cout << Branches[5].Ne[i] <<"\n";
                                    // cout <<Branches[5].Times[i] <<"\n";
                                }
                                
                                // Now fill branch 4.
                                Branches[4].Counter = int(c13.size()) - int(c03.size());
                                Branches[4].Times[Branches[4].Counter] = Branches[5].Times[0];
                                for(int i =0; i < Branches[4].Counter; i++){
                                    Branches[4].Ne[i] = SpeciesTree.Branches[c13[i]].Ne;
                                    Branches[4].Times[i] = SpeciesTree.Branches[c13[i]].Times[0];
                                }
                                
                            }// close if 03 longest
                            
                            
                        } //close if 2 is last joining indiv
                        
                        
                        //if 3 is the last joining indiv
                        if (c03.empty() + c13.empty() + c23.empty() == 3) {
                            
                            //Then set the THIRD element person to be mapped to the THIRD in the tree.
                            QuartMap[3] = indivs[3];
                            
                            // Then determine which 01 02 12 is the longest.
                            
                            //if 01 is longest
                            if(c01.size() > c12.size() && c01.size() > c02.size()){
                                // if c01 is the longest, then they make up branch 4, and third indiv will be mapped to second.
                                QuartMap[2] = indivs[2];
                                QuartMap[0] = indivs[0];
                                QuartMap[1] = indivs[1];
                                
                                
                                // 01 - 12 is the set we want for branch 4, and  12 is the branch 5.
                                Branches[5].Counter = int(c12.size());
                                Branches[5].Times[c12.size()] = Branches[6].Times[0];
                                for( int i = 0; i < c12.size(); i++){
                                    Branches[5].Ne[i] = SpeciesTree.Branches[c12[i]].Ne;
                                    Branches[5].Times[i] = SpeciesTree.Branches[c12[i]].Times[0];
                                    // cout << Branches[5].Ne[i] <<"\n";
                                    //  cout <<Branches[5].Times[i] <<"\n";
                                }
                                
                                // Now fill branch 4.
                                Branches[4].Counter = int(c01.size()) - int(c12.size());
                                Branches[4].Times[Branches[4].Counter] = Branches[5].Times[0];
                                for(int i =0; i < Branches[4].Counter; i++){
                                    Branches[4].Ne[i] = SpeciesTree.Branches[c01[i]].Ne;
                                    Branches[4].Times[i] = SpeciesTree.Branches[c01[i]].Times[0];
                                }
                                
                            }// close if 01 longest
                            
                            //if 02 is longest
                            if(c02.size() > c12.size() && c02.size() > c01.size()){
                                // if c01 is the longest, then they make up branch 4, and third indiv will be mapped to second.
                                QuartMap[2] = indivs[1];
                                QuartMap[0] = indivs[0];
                                QuartMap[1] = indivs[2];
                                
                                
                                // 02 - 12 is the set we want for branch 4, and  12 is the branch 5.
                                Branches[5].Counter = int(c12.size());
                                Branches[5].Times[c12.size()] = Branches[6].Times[0];
                                for( int i = 0; i < c12.size(); i++){
                                    Branches[5].Ne[i] = SpeciesTree.Branches[c12[i]].Ne;
                                    Branches[5].Times[i] = SpeciesTree.Branches[c12[i]].Times[0];
                                    //cout << Branches[5].Ne[i] <<"\n";
                                    //  cout <<Branches[5].Times[i] <<"\n";
                                }
                                
                                // Now fill branch 4.
                                Branches[4].Counter = int(c02.size()) - int(c12.size());
                                Branches[4].Times[Branches[4].Counter] = Branches[5].Times[0];
                                for(int i =0; i < Branches[4].Counter; i++){
                                    Branches[4].Ne[i] = SpeciesTree.Branches[c02[i]].Ne;
                                    Branches[4].Times[i] = SpeciesTree.Branches[c02[i]].Times[0];
                                }
                                
                            }// close if 02 longest
                            
                            //if 12 is longest
                            if(c12.size() > c02.size() && c12.size() > c01.size()){
                                // if c01 is the longest, then they make up branch 4, and third indiv will be mapped to second.
                                QuartMap[2] = indivs[0];
                                QuartMap[0] = indivs[1];
                                QuartMap[1] = indivs[2];
                                
                                
                                // 12 - 02 is the set we want for branch 4, and  02 is the branch 5.
                                Branches[5].Counter = int(c02.size());
                                Branches[5].Times[c02.size()] = Branches[6].Times[0];
                                for( int i = 0; i < c02.size(); i++){
                                    Branches[5].Ne[i] = SpeciesTree.Branches[c02[i]].Ne;
                                    Branches[5].Times[i] = SpeciesTree.Branches[c02[i]].Times[0];
                                    // cout << Branches[5].Ne[i] <<"\n";
                                    //  cout <<Branches[5].Times[i] <<"\n";
                                }
                                
                                // Now fill branch 4.
                                Branches[4].Counter = int(c12.size()) - int(c02.size());
                                Branches[4].Times[Branches[4].Counter] = Branches[5].Times[0];
                                for(int i =0; i < Branches[4].Counter; i++){
                                    Branches[4].Ne[i] = SpeciesTree.Branches[c12[i]].Ne;
                                    Branches[4].Times[i] = SpeciesTree.Branches[c12[i]].Times[0];
                                }
                                
                            }// close if 02 longest
                            
                            
                            
                        }// close if 3 is last joining indiv
                        
                        
                        //
                        //Now print the asym tree details so we know if we did it right.
                        //                        cout << "Quart Map \n";
                        //                        for(int i = 0; i < 4; i++){
                        //                            cout << QuartMap[i] <<"\n";
                        //                        }
                        //                        cout << "Asym Tree Check: \n";
                        //                        cout << "Branch 4 Counter = " << Branches[4].Counter<<"\n";
                        //                        for (int i =0; i < Branches[4].Counter + 1; i++) {
                        //                            cout << "B4 Time "<<i << " = "<< Branches[4].Times[i]<<"\n";
                        //                        }
                        //                        for (int i =0; i < Branches[4].Counter ; i++) {
                        //                            cout << "B4 Ne "<<i << " = "<< Branches[4].Ne[i]<<"\n";
                        //                        }
                        //                        cout << "Branch 5 Counter = " << Branches[5].Counter<<"\n";
                        //                        for (int i =0; i < Branches[5].Counter + 1; i++) {
                        //                            cout << "B5 Time "<<i << " = "<< Branches[5].Times[i]<<"\n";
                        //                        }
                        //                        for (int i =0; i < Branches[5].Counter ; i++) {
                        //                            cout << "B5 Ne "<<i << " = "<< Branches[5].Ne[i]<<"\n";
                        //                        }
                        //                        cout << "Branch 6 Counter = " << Branches[6].Counter<<"\n";
                        //                        for (int i =0; i < Branches[6].Counter + 1; i++) {
                        //                            cout << "B6 Time "<<i << " = "<< Branches[6].Times[i]<<"\n";
                        //                        }
                        //                        for (int i =0; i < Branches[6].Counter ; i++) {
                        //                            cout << "B6 Ne "<<i << " = "<< Branches[6].Ne[i]<<"\n";
                        //                        }
                        
                        
                    }else{
                        if(NumEmpty ==4){
                            
                            //In the symmetric Case
                            Branches[6].Offspring[0] = 4;
                            Branches[6].Offspring[1] = 5;
                            Branches[4].Parent = 6;
                            Branches[5].Parent = 6;
                            Branches[4].Offspring[0] = 0;
                            Branches[4].Offspring[1] = 1;
                            Branches[5].Offspring[0] = 2;
                            Branches[5].Offspring[1] = 3;
                            
                            
                            //check if 01 = 02 or 01 = 03 or 02 = 03
                            if ( c01 == c02){
                                // Then 03 are together, and 12 are together.
                                QuartMap[0] = indivs[0];
                                QuartMap[1] = indivs[3];
                                QuartMap[2] = indivs[1];
                                QuartMap[3] = indivs[2];
                                
                                Branches[4].Counter = int(c03.size());
                                Branches[5].Counter = int(c12.size());
                                
                                Branches[4].Times[Branches[4].Counter] = Branches[6].Times[0];
                                Branches[5].Times[Branches[5].Counter] = Branches[6].Times[0];
                                
                                for(int i = 0; i < Branches[4].Counter; i++){
                                    Branches[4].Ne[i] = SpeciesTree.Branches[c03[i]].Ne;
                                    Branches[4].Times[i] = SpeciesTree.Branches[c03[i]].Times[0];
                                    
                                }
                                
                                for(int i = 0; i < Branches[5].Counter; i++){
                                    Branches[5].Ne[i] = SpeciesTree.Branches[c12[i]].Ne;
                                    Branches[5].Times[i] = SpeciesTree.Branches[c12[i]].Times[0];
                                    
                                }
                                
                            } // close if 03 are together
                            
                            if ( c01 == c03){
                                // Then 02 are together, and 13 are together.
                                QuartMap[0] = indivs[0];
                                QuartMap[1] = indivs[2];
                                QuartMap[2] = indivs[1];
                                QuartMap[3] = indivs[3];
                                
                                Branches[4].Counter = int(c02.size());
                                Branches[5].Counter = int(c13.size());
                                
                                Branches[4].Times[Branches[4].Counter] = Branches[6].Times[0];
                                Branches[5].Times[Branches[5].Counter] = Branches[6].Times[0];
                                
                                for(int i = 0; i < Branches[4].Counter; i++){
                                    Branches[4].Ne[i] = SpeciesTree.Branches[c02[i]].Ne;
                                    Branches[4].Times[i] = SpeciesTree.Branches[c02[i]].Times[0];
                                    
                                }
                                
                                for(int i = 0; i < Branches[5].Counter; i++){
                                    Branches[5].Ne[i] = SpeciesTree.Branches[c13[i]].Ne;
                                    Branches[5].Times[i] = SpeciesTree.Branches[c13[i]].Times[0];
                                    
                                }
                                
                            } // close if 02 are together
                            
                            if ( c02 == c03){
                                // Then 01 are together, and 23 are together.
                                QuartMap[0] = indivs[0];
                                QuartMap[1] = indivs[1];
                                QuartMap[2] = indivs[2];
                                QuartMap[3] = indivs[3];
                                
                                Branches[4].Counter = int(c01.size());
                                Branches[5].Counter = int(c23.size());
                                
                                Branches[4].Times[Branches[4].Counter] = Branches[6].Times[0];
                                Branches[5].Times[Branches[5].Counter] = Branches[6].Times[0];
                                
                                for(int i = 0; i < Branches[4].Counter; i++){
                                    Branches[4].Ne[i] = SpeciesTree.Branches[c01[i]].Ne;
                                    Branches[4].Times[i] = SpeciesTree.Branches[c01[i]].Times[0];
                                    
                                }
                                
                                for(int i = 0; i < Branches[5].Counter; i++){
                                    Branches[5].Ne[i] = SpeciesTree.Branches[c23[i]].Ne;
                                    Branches[5].Times[i] = SpeciesTree.Branches[c23[i]].Times[0];
                                    
                                }
                                
                            } // close if 01 are together
                            
                            //Now print the asym tree details so we know if we did it right.
                            //                            cout << "Quart Map \n";
                            //                            for(int i = 0; i < 4; i++){
                            //                                cout << QuartMap[i] <<"\n";
                            //                            }
                            //                            cout << "Sym Tree Check: \n";
                            //                            cout << "Branch 4 Counter = " << Branches[4].Counter<<"\n";
                            //                            for (int i =0; i < Branches[4].Counter + 1; i++) {
                            //                                cout << "B4 Time "<<i << " = "<< Branches[4].Times[i]<<"\n";
                            //                            }
                            //                            for (int i =0; i < Branches[4].Counter ; i++) {
                            //                                cout << "B4 Ne "<<i << " = "<< Branches[4].Ne[i]<<"\n";
                            //                            }
                            //                            cout << "Branch 5 Counter = " << Branches[5].Counter<<"\n";
                            //                            for (int i =0; i < Branches[5].Counter + 1; i++) {
                            //                                cout << "B5 Time "<<i << " = "<< Branches[5].Times[i]<<"\n";
                            //                            }
                            //                            for (int i =0; i < Branches[5].Counter ; i++) {
                            //                                cout << "B5 Ne "<<i << " = "<< Branches[5].Ne[i]<<"\n";
                            //                            }
                            //                            cout << "Branch 6 Counter = " << Branches[6].Counter<<"\n";
                            //                            for (int i =0; i < Branches[6].Counter + 1; i++) {
                            //                                cout << "B6 Time "<<i << " = "<< Branches[6].Times[i]<<"\n";
                            //                            }
                            //                            for (int i =0; i < Branches[6].Counter ; i++) {
                            //                                cout << "B6 Ne "<<i << " = "<< Branches[6].Ne[i]<<"\n";
                            //                            }
                            //
                            //
                            
                            
                        }else{
                            if(NumEmpty < 3 || NumEmpty > 4){
                                cout << "Error, number empty was not equal to 3 or 4!! \n";
                                cout << "Indivs: " << i1 << ", " << i2 << ", " << i3 <<", "<<i4 <<"\n";                            cout << "    Num Empty = " << NumEmpty <<"\n";
                                cout << c01.empty() <<", "<< c02.empty() <<", " << c03.empty() << ", "<< c12.empty() <<", "<< c13.empty()  << ", " << c23.empty() <<"\n";
                                
                                cout << "c01 size = "<< c01.size() <<"\n";
                                cout << "c01 = " << c01[0] <<c01[1] << "\n";
                                
                                cout << "c02 size = "<< c02.size()<<"\n";
                                cout << "c02 = " << c02[0] << "\n";
                                
                                cout << "c03 size = "<< c03.size()<<"\n";
                                cout << "c03 = " << c03[0] <<"\n";
                                cout << "c12 size = "<< c12.size()<<"\n";
                                cout << "c12 = " << c12[0] <<c12[1] << "\n";
                                
                                cout << "c13 size = "<< c13.size()<<"\n";
                                cout << "c13 = " << c13[0] <<c13[1] << "\n";
                                cout << "c23 size = "<< c23.size()<<"\n";
                                cout << "c23 = " << c23[0] <<c23[1] << "\n";
                                
                                
                                
                                
                                // Now to get branches 4 and 5. What do we need?
                                cout << "S01 size =  "<< s01.size() <<"\n";
                                for(int i =0; i < s01.size(); i++){
                                    cout << s01[i] <<"\n";
                                }
                                cout << "S02 size = "<< s02.size() <<" \n";
                                for(int i =0; i < s02.size(); i++){
                                    cout << s02[i] <<"\n";
                                }
                                cout << "S03 = "<< s03.size() <<"\n";
                                //cout <<" s03 size = "<< s03.size() <<"\n";
                                
                                // Remove all elements which map to branch 6.
                                for(int i = 0; i < s03.size();i++){
                                    cout << s03[i] <<"\n";}
                                
                                
                                //print s03 which is the things that will go in branch 6
                                cout <<" branch 6 components: \n";
                                for (int i = 0; i <s03.size(); i++){
                                    cout<< s03[i] <<"\n";
                                    
                                }
                                cout <<"Here's ORIGINAL 0th Ancestry list \n";
                                cout << " more orig size = " <<  SpeciesTree.Branches[indivs[0]].Ancestors.size() <<"\n";
                                
                                
                                cout << " size = " << copiedset.size() <<"\n";
                                for (std::set<int>::iterator it1 = copiedset.begin(); it1 !=copiedset.end(); ++it1) {
                                    cout <<*it1 <<"\n";
                                }
                                cout <<"Here's whats left of the 0th Ancestry list \n";
                                cout << " size = " << FourTree[0].Ancestry.size() <<"\n";
                                for (std::set<int>::iterator it1 = FourTree[0].Ancestry.begin(); it1 !=FourTree[0].Ancestry.end(); ++it1) {
                                    cout <<*it1 <<"\n";
                                }
                                cout <<"Here's whats left of the 1st Ancestry list \n";
                                for (std::set<int>::iterator it1 = FourTree[1].Ancestry.begin(); it1 !=FourTree[1].Ancestry.end(); ++it1) {
                                    cout <<*it1 <<"\n";
                                }
                                
                                
                                cout <<"Here's whats left of the 2nd Ancestry list \n";
                                for (std::set<int>::iterator it1 = FourTree[2].Ancestry.begin(); it1 !=FourTree[2].Ancestry.end(); ++it1) {
                                    cout <<*it1 <<"\n";
                                }
                                cout <<"Here's whats left of the 3rd Ancestry list \n";
                                for (std::set<int>::iterator it1 = FourTree[3].Ancestry.begin(); it1 !=FourTree[3].Ancestry.end(); ++it1) {
                                    cout <<*it1 <<"\n";
                                    
                                    
                                }
                                
                                
                                //return -99999;
                            }
                        }
                        
                    }
                    
                    
                    
                    
                    
                    //  cout<< "The intersection has " << s03.size() << " elements:\n";
                    // for( int i = 0; i < s03.size(); i++){
                    //    cout << s03[i] << "\n";
                    // }
                    
                    // Drop all branches that occur in only 1 ancestors list. They don't contribute to that specific quartet calculation.
                    
                    // Determine the branches that are common with all 4, so ones that appear in all 4 Ancestors lists. This will be at least one branch. Assign this to branch 6.
                    
                    //Second, determine, of the REMAINING branches, is there a set contained as a subset of the other? (i.e. 6,7,8 and 7,8 ) Which will tell you it is a asymmetric tree, or not.
                    
                    
                    
                    
                    
                    
                    // Determine the clustering of the branches into the THREE Quartet branches.
                    // How do I do this?
                    
                    
                    
                    
                    //    //Declare the Quartet Tree,
                    //   // BranchNode Branches [7];
                    //
                    //    int WhichTree = 0; // Set to 1 for symmetric species tree, 0 for asymmetric.
                    //
                    //    if( WhichTree == 1){
                    //    //  _____________________________________________________________________
                    //    // ______________________________________________________________________
                    //    // ______EXAMPLE Symmetric TREE TO WORK WITH_____________________________
                    //
                    //
                    //
                    //    Branches[0].Parent = Branches[1].Parent = 4;
                    //    Branches[2].Parent = Branches[3].Parent = 5;
                    //
                    //    Branches[4].Offspring[0] = 0; Branches[4].Offspring[1] = 1;
                    //    Branches[4].Parent = 6;
                    //    Branches[4].Counter = 5;
                    //    Branches[4].Times[0] = 1.0; Branches[4].Times[1] = 1.5; Branches[4].Times[2] = 2.0;Branches[4].Times[3] = 2.25; Branches[4].Times[4] = 2.5; Branches[4].Times[5] = 3.0;
                    //    Branches[4].Ne[0] = 1; Branches[4].Ne[1] = 1; Branches[4].Ne[2] = 1; Branches[4].Ne[3] = 1; Branches[4].Ne[4] = 1;
                    //
                    //
                    //    Branches[5].Offspring[0] = 2; Branches [5].Offspring[1] = 3;
                    //    Branches[5].Parent = 6;
                    //    Branches[5].Counter = 3;
                    //    Branches[5].Times[0] = 2.0;Branches[5].Times[1] = 2.25; Branches[5].Times[2] = 2.5; Branches[5].Times[3] = 3.0;
                    //    Branches[5].Ne[0] = 1; Branches[5].Ne[1] = 1; Branches[5].Ne[2] = 1;
                    //
                    //    Branches[6].Offspring[0] = 4; Branches[6].Offspring[1] = 5;
                    //    Branches[6].Counter = 3;
                    //    Branches[6].Times[0] = 3.0;
                    //    Branches[6].Times[1] = 3.5;
                    //    Branches[6].Times[2] = 4.0;
                    //    Branches[6].Times[3] = 400.0; //Never used. Approx infinity
                    //    Branches[6].Ne[0] = 1;
                    //    Branches[6].Ne[1] = 1;
                    //    Branches[6].Ne[2] = 1;
                    //    } else{
                    //    //  _____________________________________________________________________
                    //    // ______________________________________________________________________
                    //    // ______EXAMPLE Asymmetric TREE TO WORK WITH_____________________________
                    //
                    //
                    //
                    //    Branches[0].Parent = Branches[1].Parent = 4;
                    //    Branches[2].Parent = 5;
                    //    Branches[3].Parent = 6;
                    //
                    //    Branches[4].Offspring[0] = 0; Branches [4].Offspring[1] = 1;
                    //    Branches[4].Parent = 5;
                    //    Branches[4].Counter = 3;
                    //    Branches[4].Times[0] = 1.0; Branches[4].Times[1] = 1.25; Branches[4].Times[2] = 1.5;Branches[4].Times[3] = 2.0;
                    //    Branches[4].Ne[0] = 1; Branches[4].Ne[1] = 1; Branches[4].Ne[2] = 1;
                    //
                    //
                    //    Branches[5].Offspring[0] = 2; Branches [5].Offspring[1] = 4;
                    //    Branches[5].Parent = 6;
                    //    Branches[5].Counter = 2;
                    //    Branches[5].Times[0] = 2.0;Branches[5].Times[1] = 2.5; Branches[5].Times[2] = 3.0;
                    //    Branches[5].Ne[0] = 1; Branches[5].Ne[1] = 1;
                    //
                    //    Branches[6].Offspring[0] = 3; Branches[6].Offspring[1] = 5;
                    //    Branches[6].Counter = 2;
                    //    Branches[6].Times[0] = 3.0;
                    //    Branches[6].Times[1] = 4.0;
                    //    Branches[6].Times[2] = 400.0; //Never used. Approx infinity
                    //    Branches[6].Ne[0] = 1;
                    //    Branches[6].Ne[1] = 1;
                    //
                    //    }
                    //
                    
                    
                    //================| Begin Calculation of Probabilities |====================
                    
                    // Calculating the useful probabilities for each interval on each branch.
                    for (int i=4; i < 7 ; i++) {
                        for (int j =0 ; j < Branches[i].Counter  ; j++) {
                            //Calculate T = (t_{i+1}-t_i)/(2N_i)
                            double T = (Branches[i].Times[j+1]-Branches[i].Times[j])/(2.0*Branches[i].Ne[j]);
                            Branches[i].P22[j] = exp(-T);
                            Branches[i].P22A[j] = exp(-T);
                            Branches[i].P22B[j] = exp( - 3.0 * T); // rate 2/3 N
                            Branches[i].P22C[j] = exp(-6.0 * T); // Rate 1/3
                            Branches[i].P22H[j] = exp( - 2.0 * T);
                            
                            if( (i ==6) && (j == Branches[6].Counter - 1)){
                                Branches[i].P22[j] = 0.0;
                                Branches[i].P22A[j] = 0.0;
                                Branches[i].P22B[j] = 0.0;
                                Branches[i].P22C[j] = 0.0;
                                Branches[i].P22H[j] = 0.0;
                                
                            }
                            
                            
                            Branches[i].P33[j] = exp(-3*T);
                            Branches[i].P32[j] = (3.0/2.0) * (exp(-T) - exp(-3*T));
                            //Branches[i].P44[j] = exp(-6*T);
                            //Branches[i].P43[j] = 2*(exp(-3*T)-exp(-6*T));
                            //Branches[i].P42[j] = 9/5*exp(-T)-3*exp(-3*T)+6/5*exp(-6*T);
                            Branches[i].TP[j] = -(2.0 * Branches[i].Ne[j] + Branches[i].Times[j+1])* Branches[i].P22[j] + 2.0* Branches[i].Ne[j] + Branches[i].Times[j];
                            
                            Branches[i].TPA[j] = -(2.0 * Branches[i].Ne[j] + Branches[i].Times[j+1])* Branches[i].P22A[j] + 2.0* Branches[i].Ne[j] + Branches[i].Times[j];
                            
                            Branches[i].TPH[j] = -(Branches[i].Ne[j] + Branches[i].Times[j+1])* Branches[i].P22H[j] + Branches[i].Ne[j] + Branches[i].Times[j];
                            
                            Branches[i].TPB[j] = -(2.0/3.0 * Branches[i].Ne[j] + Branches[i].Times[j+1])* Branches[i].P22B[j] + 2.0/3.0* Branches[i].Ne[j] + Branches[i].Times[j];
                            
                            Branches[i].TPC[j] = -(1.0/3.0 * Branches[i].Ne[j] + Branches[i].Times[j+1])* Branches[i].P22C[j] + 1.0/3.0* Branches[i].Ne[j] + Branches[i].Times[j];
                            
                            Branches[i].TPP22[j] = (pow(Branches[i].Times[j+1],2) - pow(Branches[i].Times[j],2)) / ( 4.0 * Branches[i].Ne[j]) * Branches[i].P22[j];
                            Branches[i].T2P[j] = - (8 * pow(Branches[i].Ne[j],2) + 4 * Branches[i].Ne[j] * Branches[i].Times[j+1] + pow(Branches[i].Times[j+1],2))* Branches[i].P22[j] + 8 * pow(Branches[i].Ne[j],2) + 4 * Branches[i].Ne[j] * Branches[i].Times[j] + pow(Branches[i].Times[j],2);
                            Branches[i].PP22[j] = T * Branches[i].P22[j];
                            Branches[i].T2PP22[j] = (pow(Branches[i].Times[j+1],3) - pow(Branches[i].Times[j],3)) / ( 6.0 * Branches[i].Ne[j]) * Branches[i].P22[j];
                            Branches[i].PP22BA[j] = 3.0 / 2.0 *(Branches[i].P22A[j] - Branches[i].P22B[j]);
                            Branches[i].PP22CB[j] = 2.0 * ( Branches[i].P22B[j] - Branches[i].P22C[j]);
                            Branches[i].PP22CA[j] = 6.0/5.0 * (Branches[i].P22A[j] - Branches[i].P22C[j]);
                            Branches[i].PP22CBA[j] = 3.0/2.0 *(Branches[i].PP22CA[j] - Branches[i].PP22CB[j]);
                            
                        }
                    }
                    // Need to correct the calculations for the last interval.
                    // I THINK THIS IS REDUNDANT NOW, i think i have properly fixed this.
                    Branches[6].P22[Branches[6].Counter - 1] = 0;
                    Branches[6].TP[Branches[6].Counter - 1] =  2 * Branches[6].Ne[Branches[6].Counter - 1] + Branches[6].Times[Branches[6].Counter - 1];
                    
                    
                    // NEW CALCULATIONS for the last branch, 6.
                    
                    //IMPORTANT, need to go through and correct for last interval. NEED These for 5 too. if asymmetric.
                    
                    for (int i=0; i < Branches[6].Counter; i++) {
                        double Ni = Branches[6].Ne[i];
                        double T = (Branches[6].Times[i+1]-Branches[6].Times[i])/(2.0*Ni);
                        
                        
                        Branches[6].P5NC[i] = exp(-5*T);
                        Branches[6].P2NC[i] = exp(-2*T);
                        Branches[6].P3NC[i] = exp(-3*T);
                        Branches[6].P6NC[i] = exp(-6*T);
                        
                        if ( i == Branches[6].Counter -1){
                            
                            Branches[6].P5NC[i] = 0.0;
                            Branches[6].P2NC[i] = 0.0;
                            Branches[6].P3NC[i] = 0.0;
                            Branches[6].P6NC[i] = 0.0;
                            
                            
                        }
                        
                        
                        Branches[6].TP5[i] =  -(2.0*Ni/5.0 + Branches[6].Times[i+1]) * Branches[6].P5NC[i] + 2.0* Ni/5.0 + Branches[6].Times[i];
                        Branches[6].TPP5[i] = 1.0/6.0 *(- (2.0*Ni/ 6.0 + Branches[6].Times[i+1]) * Branches[6].P6NC[i] + 2.0*Ni/6.0 + Branches[6].Times[i]);
                        //cout << "TPP5 = " << Branches[6].TPP5[i] <<"\n";
                        Branches[6].T2PP5[i] = 1.0/108.0 * (-(4 * pow(Ni,2) + 12 * Ni * Branches[6].Times[i+1] + 18 * pow(Branches[6].Times[i+1],2)) * Branches[6].P6NC[i] + 4 * pow(Ni,2) + 12 * Ni * Branches[6].Times[i] + 18 * pow(Branches[6].Times[i],2));
                        Branches[6].TPP2[i] = 1.0/3.0 * ( -(2*Ni /3.0 + Branches[6].Times[i+1])* Branches[6].P3NC[i] + 2.0 * Ni / 3.0 + Branches[6].Times[i]);
                        Branches[6].T2PP2[i] = 1.0/27.0 * ( -(8 * pow(Ni,2) + 12 * Ni * Branches[6].Times[i+1] + 9 * pow(Branches[6].Times[i+1],2))*Branches[6].P3NC[i] + 8 * pow(Ni,2) + 12 * Ni * Branches[6].Times[i] + 9 * pow(Branches[6].Times[i],2));
                        Branches[6].P5P2NC[i] = 5.0/3.0 * ( Branches[6].P2NC[i] - Branches[6].P5NC[i]);
                        Branches[6].TPP2P22[i] = (Ni + Branches[6].Times[i])/2.0 * Branches[6].P2NC[i] - (Ni + Branches[6].Times[i+1])/2.0 * Branches[6].P3NC[i];
                        Branches[6].TPP5P22[i] = 1.0/5.0 * ( ((2.0*Ni)/5.0 + Branches[6].Times[i])* Branches[6].P22[i] - ((2.0 * Ni)/5.0 + Branches[6].Times[i+1]) * Branches[6].P6NC[i]);
                        Branches[6].TPP5P3NC[i] = 5.0/2.0*((Ni + Branches[6].Times[i]) * Branches[6].P3NC[i] - (Ni + Branches[6].Times[i+1]) * Branches[6].P5NC[i]);
                        
                        Branches[6].TPP5NCP3NC[i] = 1.0/9.0*(3* Branches[6].Times[i] + 2*Ni)*Branches[6].P3NC[i] - 1.0/9.0*(3* Branches[6].Times[i+1] + 2*Ni) * Branches[6].P6NC[i];
                        
                        Branches[6].P3P22[i] = 3.0/2.0 * ( Branches[6].P22[i] - Branches[6].P3NC[i]);
                        
                        Branches[6].TP2[i] = -(Ni + Branches[6].Times[i+1])* Branches[6].P2NC[i] + Ni + Branches[6].Times[i];
                        Branches[6].TP3[i] =-(2.0/3.0 *Ni + Branches[6].Times[i+1])* Branches[6].P3NC[i] + 2.0/3.0* Ni + Branches[6].Times[i];
                    }
                    
                    // Is the species tree symmetric? If so, need to pass information back and forth.
                    
                    //  determine if in symmetric or assymetric quartet case
                    bool isSim = (4 != Branches[5].Offspring[0]) && (4 != Branches[5].Offspring[1]);
                    // cout << " This equals one if symmetric: "<< isSim << "\n";
                    
                    // Declare superbranch outside of the if statement just so it exists.
                    BranchNode SuperBranch;
                    
                    if(isSim != 1){
                        
                        // Get the new calcs for branch 5 as well.
                        for (int i=0; i < Branches[5].Counter; i++) {
                            double Ni = Branches[5].Ne[i];
                            double T = (Branches[5].Times[i+1]-Branches[5].Times[i])/(2.0*Ni);
                            
                            
                            Branches[5].P5NC[i] = exp(-5*T);
                            Branches[5].P2NC[i] = exp(-2*T);
                            Branches[5].P3NC[i] = exp(-3*T);
                            Branches[5].P6NC[i] = exp(-6*T);
                            
                            
                            
                            
                            Branches[5].TP5[i] =  -(2.0*Ni/5.0 + Branches[5].Times[i+1]) * Branches[5].P5NC[i] + 2.0* Ni/5.0 + Branches[5].Times[i];
                            Branches[5].TPP5[i] = 1.0/6.0 *(- (2.0*Ni/ 6.0 + Branches[5].Times[i+1]) * Branches[5].P6NC[i] + 2.0*Ni/6.0 + Branches[5].Times[i]);
                            // cout << "TPP5 = " << Branches[5].TPP5[i] <<"\n";
                            Branches[5].T2PP5[i] = 1.0/108.0 * (-(4 * pow(Ni,2) + 12 * Ni * Branches[5].Times[i+1] + 18 * pow(Branches[5].Times[i+1],2)) * Branches[5].P6NC[i] + 4 * pow(Ni,2) + 12 * Ni * Branches[5].Times[i] + 18 * pow(Branches[5].Times[i],2));
                            Branches[5].TPP2[i] = 1.0/3.0 * ( -(2*Ni /3.0 + Branches[5].Times[i+1])* Branches[5].P3NC[i] + 2.0 * Ni / 3.0 + Branches[5].Times[i]);
                            Branches[5].T2PP2[i] = 1.0/27.0 * ( -(8 * pow(Ni,2) + 12 * Ni * Branches[5].Times[i+1] + 9 * pow(Branches[5].Times[i+1],2))*Branches[5].P3NC[i] + 8 * pow(Ni,2) + 12 * Ni * Branches[5].Times[i] + 9 * pow(Branches[5].Times[i],2));
                            Branches[5].P5P2NC[i] = 5.0/3.0 * ( Branches[5].P2NC[i] - Branches[5].P5NC[i]);
                            Branches[5].TPP2P22[i] = (Ni + Branches[5].Times[i])/2.0 * Branches[5].P2NC[i] - (Ni + Branches[5].Times[i+1])/2.0 * Branches[5].P3NC[i];
                            Branches[5].TPP5P22[i] = 1.0/5.0 * ( ((2.0*Ni)/5.0 + Branches[5].Times[i])* Branches[5].P22[i] - ((2.0 * Ni)/5.0 + Branches[5].Times[i+1]) * Branches[5].P6NC[i]);
                            Branches[5].TPP5P3NC[i] = 5.0/2.0*((Ni + Branches[5].Times[i]) * Branches[5].P3NC[i] - (Ni + Branches[5].Times[i+1]) * Branches[5].P5NC[i]);
                            
                            Branches[5].TPP5NCP3NC[i] = 1.0/9.0*(3* Branches[5].Times[i] + 2*Ni)*Branches[5].P3NC[i] - 1.0/9.0*(3* Branches[5].Times[i+1] + 2*Ni) * Branches[5].P6NC[i];
                            Branches[5].P3P22[i] = 3.0/2.0 * ( Branches[5].P22[i] - Branches[5].P3NC[i]);
                            
                            Branches[5].TP2[i] = -(Ni + Branches[5].Times[i+1])* Branches[5].P2NC[i] + Ni + Branches[5].Times[i];
                            Branches[5].TP3[i] =-(2.0/3.0 *Ni + Branches[5].Times[i+1])* Branches[5].P3NC[i] + 2.0/3.0* Ni + Branches[5].Times[i];
                        }
                        
                        
                        
                        //Fill SuperBranch
                        
                        int TotCount = Branches[5].Counter + Branches[6].Counter;
                        SuperBranch.Counter = TotCount;
                        
                        for (int i=0; i < Branches[5].Counter ; i++) {
                            
                            SuperBranch.Times[i] = Branches[5].Times[i];
                            SuperBranch.Ne[i] = Branches[5].Ne[i];
                            SuperBranch.P22[i] = Branches[5].P22[i];
                            SuperBranch.TP[i] = Branches[5].TP[i];
                            SuperBranch.TPP22[i] = Branches[5].TPP22[i];
                            SuperBranch.T2P[i] = Branches[5].T2P[i];
                            SuperBranch.T2PP22[i] = Branches[5].T2PP22[i];
                            SuperBranch.PP22[i] = Branches[5].PP22[i];
                            SuperBranch.EgtTi[i] = Branches[5].EgtTi[i];
                            SuperBranch.TP5[i] = Branches[5].TP5[i];
                            SuperBranch.P5NC[i] = Branches[5].P5NC[i];
                            SuperBranch.TPP5[i] = Branches[5].TPP5[i];
                            SuperBranch.T2PP5[i] = Branches[5].T2PP5[i];
                            SuperBranch.TPP2[i] = Branches[5].TPP2[i];
                            SuperBranch.T2PP2[i] = Branches[5].T2PP2[i];
                            SuperBranch.P5P2NC[i] = Branches[5].P5P2NC[i];
                            SuperBranch.TPP2P22[i] = Branches[5].TPP2P22[i];
                            SuperBranch.TPP5P22[i] = Branches[5].TPP5P22[i];
                            SuperBranch.TPP5P3NC[i] = Branches[5].TPP5P3NC[i];
                            SuperBranch.P3P22[i] = Branches[5].P3P22[i];
                            SuperBranch.TP2[i] = Branches[5].TP2[i];
                            SuperBranch.P2NC[i] = Branches[5].P2NC[i];
                            SuperBranch.TP3[i] = Branches[5].TP3[i];
                            SuperBranch.P3NC[i] = Branches[5].P3NC[i];
                            SuperBranch.P6NC[i] = Branches[5].P6NC[i];
                        }
                        
                        for (int i=Branches[5].Counter; i < TotCount; i++ ){
                            int j = i - Branches[5].Counter;
                            
                            SuperBranch.Times[i] = Branches[6].Times[j];
                            SuperBranch.Ne[i] = Branches[6].Ne[j];
                            SuperBranch.P22[i] = Branches[6].P22[j];
                            SuperBranch.TP[i] = Branches[6].TP[j];
                            SuperBranch.TPP22[i] = Branches[6].TPP22[j];
                            SuperBranch.T2P[i] = Branches[6].T2P[j];
                            SuperBranch.T2PP22[i] = Branches[6].T2PP22[j];
                            SuperBranch.PP22[i] = Branches[6].PP22[j];
                            SuperBranch.EgtTi[i] = Branches[6].EgtTi[j];
                            SuperBranch.TP5[i] = Branches[6].TP5[j];
                            SuperBranch.P5NC[i] = Branches[6].P5NC[j];
                            SuperBranch.TPP5[i] = Branches[6].TPP5[j];
                            SuperBranch.T2PP5[i] = Branches[6].T2PP5[j];
                            SuperBranch.TPP2[i] = Branches[6].TPP2[j];
                            SuperBranch.T2PP2[i] = Branches[6].T2PP2[j];
                            SuperBranch.P5P2NC[i] = Branches[6].P5P2NC[j];
                            SuperBranch.TPP2P22[i] = Branches[6].TPP2P22[j];
                            SuperBranch.TPP5P22[i] = Branches[6].TPP5P22[j];
                            SuperBranch.TPP5P3NC[i] = Branches[6].TPP5P3NC[j];
                            SuperBranch.P3P22[i] = Branches[6].P3P22[j];
                            SuperBranch.TP2[i] = Branches[6].TP2[j];
                            SuperBranch.P2NC[i] = Branches[6].P2NC[j];
                            SuperBranch.TP3[i] = Branches[6].TP3[j];
                            SuperBranch.P3NC[i] = Branches[6].P3NC[j];
                            SuperBranch.P6NC[i] = Branches[6].P6NC[j];
                            
                        }
                        
                    }
                    
                    
                    
                    
                    
                    // Calculating FullP22, FullP33, and FullP32 for branches 4 and 5.
                    Branches[4].FullP22 = CalcP22(Branches[4]);
                    Branches[5].FullP22 = CalcP22(Branches[5]);
                    Branches[5].FullP33 = CalcP33(Branches[5]);
                    Branches[5].FullP32 = CalcP32(Branches[5]);
                    //cout << "P22 on branch 5 = " << Branches[5].FullP22 << "\n";
                    
                    
                    
                    
                    
                    
                    
                    //START: Some Message Passing to spread Expectation(T | T > t_i)________
                    //and Expectation(T^2 | T > t_i)________
                    
                    // This is only useful for rates 2N for everything, not useful for variable rates.
                    
                    
                    // Initialize: Noting the top branch top interval expectation is just TP.
                    
                    Branches[6].EgtTi[Branches[6].Counter-1]= Branches[6].TP[Branches[6].Counter-1];
                    Branches[6].EgtT2i[Branches[6].Counter-1]= Branches[6].T2P[Branches[6].Counter-1];
                    
                    
                    // Pass through this branch
                    for( int j = (Branches[6].Counter- 2); j >= 0; j--){
                        Branches[6].EgtTi[j] = Branches[6].P22[j] * ( Branches[6].EgtTi[j+1] - (2* Branches[6].Ne[j] + Branches[6].Times[j+1]))+ 2* Branches[6].Ne[j] + Branches[6].Times[j];
                        Branches[6].EgtT2i[j] = Branches[6].P22[j] * Branches[6].EgtT2i[j+1] + Branches[6].T2P[j];
                        
                        //cout << " EgtTi = "<< Branches[6].EgtTi[j] << "\n";
                    }
                    
                    //Starting at the nodes after the top branch (so 5,4)
                    for (int i = 5; i > 3; i--) {
                        for( int j = Branches[i].Counter - 1; j >= 0; j--){
                            
                            
                            // if the top interval of the branch, pull value down from first interval of parent branch
                            if (j == (Branches[i].Counter - 1)){
                                
                                Branches[i].EgtTi[j]= Branches[i].P22[j]
                                * ( Branches[Branches[i].Parent].EgtTi[0]
                                   - (2* Branches[i].Ne[j] + Branches[i].Times[j+1]))
                                + 2 * Branches[i].Ne[j] + Branches[i].Times[j];
                                
                                Branches[i].EgtT2i[j] = Branches[i].P22[j] * Branches[Branches[i].Parent].EgtT2i[0] + Branches[i].T2P[j];
                                
                                
                                //Otherwise, just use the information from the interval directly above on the branch.
                            } else{
                                
                                Branches[i].EgtTi[j] = Branches[i].P22[j] *(Branches[i].EgtTi[j+1] -(2 * Branches[i].Ne[j] + Branches[i].Times[j+1])) + 2* Branches[i].Ne[j] + Branches[i].Times[j];
                                
                                Branches[i].EgtT2i[j] = Branches[i].P22[j] * Branches[i].EgtT2i[j+1] + Branches[i].T2P[j];
                                
                            }
                            //cout << " EgtTi = "<< Branches[i].EgtTi[j] << "\n";
                        }
                    }
                    //END: Some Message Passing to spread Expectation(T | T > t_i)________
                    
                    
                    
                    //___________________________________________________________________________________________________________________
                    // Calculating joint expectation in a single population.
                    
                    double  E0101 ;
                    double  E0102 ;
                    double  E0103 ;
                    double  E0112 ;
                    double  E0113 ;
                    double  E0123 ;
                    double  E0202 ;
                    double  E0203 ;
                    double  E0212 ;
                    double  E0213 ;
                    double  E0223 ;
                    double  E0303 ;
                    double  E0312 ;
                    double  E0313 ;
                    double  E0323 ;
                    double  E1212 ;
                    double  E1213 ;
                    double  E1223 ;
                    double  E1313 ;
                    double  E1323 ;
                    double  E2323 ;
                    double E01;
                    double E02;
                    double E03;
                    double E12;
                    double E13;
                    double E23;
                    
                    //Conditionals (Notation: T01gT23m is T01 given T23 is first coal, these include the probablity that T23 is the first coal, or whichever one.)
                    
                    double E01g01m;
                    double E01g23m;
                    double E02g02m;
                    double E02g13m;
                    double E03g03m;
                    double E03g12m;
                    double E12g12m;
                    double E12g03m;
                    double E13g13m;
                    double E13g02m;
                    double E23g23m;
                    double E23g01m;
                    
                    //Conditional Trios eg. E01g01m2 is E(T01 given 01 is min of trio 012)
                    //Trio 012
                    double E01g01m2;
                    double E01g02m;
                    double E01g12m;
                    double E02g02m1;
                    double E02g01m;
                    double E02g12m;
                    double E12g12m0;
                    double E12g01m;
                    double E12g02m;
                    //Trio 013
                    double E01g01m3;
                    double E01g03m;
                    double E01g13m;
                    double E03g03m1;
                    double E03g01m;
                    double E03g13m;
                    double E13g13m0;
                    double E13g01m;
                    double E13g03m;
                    //Trio 023
                    double E02g02m3;
                    double E02g03m;
                    double E02g23m;
                    double E03g03m2;
                    double E03g02m;
                    double E03g23m;
                    double E23g23m0;
                    double E23g02m;
                    double E23g03m;
                    //Trio 023
                    double E12g12m3;
                    double E12g13m;
                    double E12g23m;
                    double E13g13m2;
                    double E13g12m;
                    double E13g23m;
                    double E23g23m1;
                    double E23g12m;
                    double E23g13m;
                    
                    
                    if(isSim == 1){
                        
                        
                        double Case1 = GetCase1(Branches[6]);
                        double Case2 = GetCase2(Branches[6]);
                        double Case3 = GetCase3(Branches[6]);
                        double Case1alt = GetCase1alt(Branches[6]);
                        //double Case4 = GetCase4(Branches[6]);
                        //double Case5 = GetCase5(Branches[6]);
                        //double Case6 = GetCase6(Branches[6]);
                        
                        
                        double Expectation4ind = 2*Case1 + 2*Case2 +Case3 ;
                        double Expectation3ind =  2*Case1 + +Case1alt+ 1.5 * Case2 + 0.75 * Case3;
                        //cout<< "E(Tcd*Tab|Same Pop) = " << Expectation4ind<< "\n";
                        //cout << "E(Tab *Tac|Same Pop) = " << Expectation3ind << "\n";
                        
                        //conditional stuff.
                        double Term1 = GetTerm1(Branches[6]); // 4 indiv, given self min
                        double Term2 = GetTerm2(Branches[6]); // 4 indiv, given 2 others are min
                        double Term3 = GetTerm3(Branches[6]); // 3 indiv, given not min
                        double Term4 = GetTerm4(Branches[6]); // 3 indiv, given self min
                        double Term8_4 = GetTerm8(Branches[4]); //3 indiv given self min
                        double Term8_5 = GetTerm8(Branches[5]); //3 indiv given self min
                        
                        
                        
                        
                        // If Symmetric, easy to get the pairwise expectations as above.
                        
                        //Test
                        // double BigT01 = Branches[4].FullP22 * Expectation3ind;
                        //double SmallT01 = (Branches[4].EgtTi[0] - Branches[4].FullP22*Branches[6].EgtTi[0])*Branches[6].EgtTi[0];
                        //01 01
                        E0101 = Branches[4].EgtT2i[0];
                        // 01 02
                        
                        E0102 = Branches[4].FullP22 * Expectation3ind + Branches[6].EgtTi[0] * ( Branches[4].EgtTi[0] - Branches[4].FullP22 * Branches[6].EgtTi[0]);
                        //01 03
                        E0103 = E0102;
                        //01 12
                        E0112 = E0102;
                        //01 13
                        E0113 = E0102;
                        //01 23
                        E0123 = Branches[4].FullP22*Branches[5].FullP22 * Expectation4ind + Branches[6].EgtTi[0] * (Branches[4].EgtTi[0] - Branches[4].FullP22 * Branches[6].EgtTi[0])* Branches[5].FullP22 + Branches[6].EgtTi[0] *( Branches[5].EgtTi[0] - Branches[5].FullP22 * Branches[6].EgtTi[0])* Branches[4].FullP22 + (Branches[4].EgtTi[0] - Branches[4].FullP22 * Branches[6].EgtTi[0]) * (Branches[5].EgtTi[0] - Branches[5].FullP22 * Branches[6].EgtTi[0]);
                        //02 02
                        E0202 = Branches[6].EgtT2i[0];
                        //02 03
                        E0203 = (1- Branches[5].FullP22) * Branches[6].EgtT2i[0] + Branches[5].FullP22 * Expectation3ind;
                        //02 12
                        E0212 = (1-Branches[4].FullP22) * Branches[6].EgtT2i[0] + Branches[4].FullP22 * Expectation3ind;
                        //02 13
                        E0213 = Branches[4].FullP22 * Branches[5].FullP22 * Expectation4ind + (1- Branches[4].FullP22) * Branches[5].FullP22 * Expectation3ind + Branches[4].FullP22 * (1- Branches[5].FullP22) * Expectation3ind + (1-Branches[4].FullP22)* (1-Branches[5].FullP22) * Branches[6].EgtT2i[0];
                        //02 23
                        E0223 = Branches[5].FullP22 * Expectation3ind + Branches[6].EgtTi[0]*(Branches[5].EgtTi[0] - Branches[5].FullP22 * Branches[6].EgtTi[0]);
                        //03 03
                        E0303 = E0202;
                        //03 12
                        E0312 = E0213;
                        //03 13
                        E0313 = E0212;
                        //03 23
                        E0323 = E0223;
                        //12 12
                        E1212 = E0202;
                        //12 13
                        E1213 = E0203;
                        //12 23
                        E1223 = E0223;
                        //13 13
                        E1313 = E0202;
                        //13 23
                        E1323 = E0223;
                        //23 23
                        E2323 = Branches[5].EgtT2i[0];
                        
                        // 01
                        E01 = Branches[4].EgtTi[0];
                        // 02
                        E02 = Branches[6].EgtTi[0];
                        // 03
                        E03 = E02;
                        // 12
                        E12 = E02;
                        // 13
                        E13 = E02;
                        // 23
                        E23 = Branches[5].EgtTi[0];
                        
                        
                        //Expected Values given a certain pair is the minumum. For symmetric case.
                        E01g01m = Branches[5].FullP22* Branches[4].FullP22*Term1 + GetTerm5(Branches[4], Branches[5]);  //Verified simple case.
                        
                        E23g23m = Branches[4].FullP22*Branches[5].FullP22*Term1 + GetTerm5(Branches[5], Branches[4]);  //Verified simple case
                        
                        // ^ These two potentially (but unlikely) could be missing the term where the other doesnt coalesce in their branch . Term 5 probably handles it though.
                        
                        
                        E01g23m = Branches[5].FullP22 * Branches[4].FullP22 * Term2 +  GetTerm6(Branches[4],Branches[5]) + Branches[6].EgtTi[0] * (1.0 - Branches[5].FullP22)* Branches[4].FullP22; //Verified simple case
                        
                        
                        E23g01m = Branches[4].FullP22 *Branches[5].FullP22* Term2 +  GetTerm6(Branches[5],Branches[4]) + Branches[6].EgtTi[0] * (1.0 - Branches[4].FullP22)* Branches[5].FullP22; //Verified simple case
                        
                        E02g02m = Branches[4].FullP22*Branches[5].FullP22*Term1; //Verified simple case
                        E03g03m = E02g02m;
                        E12g12m = E02g02m;
                        E13g13m = E02g02m;
                        
                        E02g13m = Branches[4].FullP22*Branches[5].FullP22*Term2; //Verified simple case
                        E03g12m = E02g13m;
                        E12g03m = E02g13m;
                        E13g02m = E02g13m;
                        
                        //Conditional Trios (symmetric)
                        //Q0
                        E01g01m2 = Term8_4 +Branches[4].FullP22 * Term4;  //Verified simple case
                        E01g01m3 =E01g01m2;
                        //Q1
                        E01g02m = Branches[4].FullP22 * Term3; //Verified simple case
                        E01g12m = E01g02m;
                        E01g03m = E01g02m;
                        E01g13m = E01g02m;
                        //Q2
                        E02g02m1 = Branches[4].FullP22 * Term4; //Verified simple case
                        E12g12m0 = E02g02m1;
                        E03g03m1 = E02g02m1;
                        E13g13m0 = E02g02m1;
                        //Q3
                        E02g01m = Branches[4].FullP22 * Term3 + Branches[6].EgtTi[0] * (1.0 - Branches[4].FullP22); //Verified simple case
                        E12g01m = E02g01m;
                        E03g01m = E02g01m;
                        E13g01m = E02g01m;
                        //Q4
                        E02g12m = Branches[4].FullP22 * Term3; //Verified simple case
                        E12g02m = E02g12m;
                        E03g13m = E02g12m;
                        E13g03m = E02g12m;
                        //Q5
                        E02g02m3 = Branches[5].FullP22 * Term4; //Verified simple case
                        E03g03m2 = E02g02m3;
                        E12g12m3 = E02g02m3;
                        E13g13m2 = E02g02m3;
                        //Q6
                        E02g03m = Branches[5].FullP22 * Term3; //Verified simple case
                        E03g02m = E02g03m;
                        E12g13m = E02g03m;
                        E13g12m = E02g03m;
                        //Q7
                        E02g23m = Branches[5].FullP22 * Term3 + Branches[6].EgtTi[0] * (1.0 - Branches[5].FullP22); //Verified simple case
                        E03g23m = E02g23m;
                        E12g23m = E02g23m;
                        E13g23m = E02g23m;
                        //Q8
                        E23g23m0 = Term8_5 + Branches[5].FullP22 * Term4;//Verified simple case
                        E23g23m1 = E23g23m0;
                        //Q9
                        E23g02m = Branches[5].FullP22*Term3; //Verified simple case
                        E23g03m = E23g02m;
                        E23g12m = E23g02m;
                        E23g13m = E23g02m;
                        
                        
                    }else{
                        //If Asymmetric, harder to get pairwise expectations, requires more work.
                        
                        double Case1 = GetCase1(Branches[6]);
                        double Case2 = GetCase2(Branches[6]);
                        double Case3 = GetCase3(Branches[6]);
                        double Case1alt = GetCase1alt(Branches[6]);
                        //double Case4 = GetCase4(Branches[6]);
                        //double Case5 = GetCase5(Branches[6]);
                        //double Case6 = GetCase6(Branches[6]);
                        
                        
                        double Expectation4ind = 2*Case1 + 2*Case2 +Case3 ;
                        double Expectation3ind =  2*Case1 + +Case1alt+ 1.5 * Case2 + 0.75 * Case3;
                        
                        double Term1 = GetTerm1(Branches[6]); // 4 indiv, given self min
                        double Term2 = GetTerm2(Branches[6]); // 4 indiv, given 2 others are min
                        double Term3 = GetTerm3(Branches[6]); // 3 indiv, given not min
                        double Term4 = GetTerm4(Branches[6]); // 3 indiv, given self min
                        double Term7 = GetTerm7(Branches[5]); //4 indiv given self in branch 5.
                        double Term8 = GetTerm8(Branches[4]); // 4 indiv, given self in branch 4.
                        double Term4_5 = GetTerm4(Branches[5]); // 3 indiv, given self min
                        double Term3_5 = GetTerm3(Branches[5]); // 3 indiv, given not min
                        double Term8_5 = GetTerm8(Branches[5]); // 4 indiv, given self in branch 4.
                        
                        // cout<< "E(Tcd*Tab|Same Pop) = " << Expectation4ind<< "\n";
                        //cout << "E(Tab *Tac|Same Pop) = " << Expectation3ind << "\n";
                        //Need to create a superbranch of Branches 5 and 6.
                        // BranchNode SuperBranch = ;
                        // How do I get this superbranch??
                        // First need to have the approp quantities for branch 5. DONE
                        //How do I do this.... I Think done.
                        double Case1SB = GetCase1(SuperBranch);
                        double Case1altSB = GetCase1alt(SuperBranch);
                        double Case2SB = GetCase2(SuperBranch);
                        double Case3SB = GetCase3(SuperBranch);
                        // double Expectation4indSB = 2*Case1SB + 2*Case2SB +Case3SB ; // Dont think I'll ever need this one.
                        double Expectation3indSB =  2*Case1SB + Case1altSB + 1.5 * Case2SB + 0.75 * Case3SB; // This is a key quantity!!
                        //cout<< "E(Tcd*Tab|SB) = " << Expectation4indSB<< "\n";
                        //cout << "E(Tab *Tac|SB) = " << Expectation3indSB << "\n";
                        //01 01
                        E0101 = Branches[4].EgtT2i[0];
                        //01 02
                        E0102 = Branches[4].FullP22 * Expectation3indSB + Branches[5].EgtTi[0] * ( Branches[4].EgtTi[0] - Branches[4].FullP22 * Branches[5].EgtTi[0]);
                        //01 03
                        E0103 = Branches[4].FullP22 * Branches[5].FullP22 * Expectation3ind + Branches[6].EgtTi[0] *(Branches[4].EgtTi[0] - Branches[4].FullP22*Branches[5].FullP22* Branches[6].EgtTi[0]);
                        //01 12
                        E0112 = E0102;
                        //01 13
                        E0113 = E0103;
                        //01 23 **
                        E0123 = ((Branches[4].EgtTi[0] - Branches[4].FullP22 * Branches[5].EgtTi[0]) + Branches[4].FullP22 * (Branches[5].EgtTi[0] - Branches[5].FullP22 * Branches[6].EgtTi[0])) * Branches[6].EgtTi[0] + Branches[4].FullP22 *( 2.0/3.0 * Branches[5].FullP32 * Expectation3ind +Branches[5].FullP33 * Expectation4ind);
                        //02 02
                        E0202 = Branches[5].EgtT2i[0];
                        //02 03
                        E0203 = Branches[5].FullP22 * Expectation3ind + Branches[6].EgtTi[0] * ( Branches[5].EgtTi[0] - Branches[5].FullP22 * Branches[6].EgtTi[0]);
                        //02 12
                        E0212 = (1-Branches[4].FullP22) * Branches[5].EgtT2i[0] + Branches[4].FullP22 * Expectation3indSB;
                        //02 13
                        E0213 = (Branches[5].EgtTi[0] - Branches[5].FullP22* Branches[6].EgtTi[0]) * Branches[6].EgtTi[0] + ( (1- Branches[4].FullP22)* Branches[5].FullP22 + Branches[4].FullP22 * 2.0/3.0 * Branches[5].FullP32) * Expectation3ind + Branches[4].FullP22 * Branches[5].FullP33 * Expectation4ind;
                        //02 23
                        E0223 = E0203;
                        //03 03
                        E0303 = Branches[6].EgtT2i[0];
                        //03 12
                        E0312 = E0213;
                        //03 13
                        E0313 = ((1-Branches[4].FullP22) + Branches[4].FullP22 *(1-Branches[5].FullP22))* Branches[6].EgtT2i[0] + Branches[4].FullP22 * Branches[5].FullP22 * Expectation3ind;
                        //03 23
                        E0323 = (1-Branches[5].FullP22)* Branches[6].EgtT2i[0] + Branches[5].FullP22 * Expectation3ind;
                        //12 12
                        E1212 = E0202;
                        //12 13
                        E1213 = E0203;
                        //12 23
                        E1223 = E0203;
                        //13 13
                        E1313 = E0303;
                        //13 23
                        E1323 = E0323;
                        //23 23
                        E2323 = E0303;
                        //01
                        E01 = Branches[4].EgtTi[0];
                        //02
                        E02 = Branches[5].EgtTi[0];
                        //03
                        E03 = Branches[6].EgtTi[0];
                        //12
                        E12 = E02;
                        //13
                        E13 = E03;
                        //23
                        E23 = E03;
                        
                        //Expected Values given a certain pair is the minumum. For asymmetric case.
                        E03g03m = Branches[4].FullP22* Branches[5].FullP33 * Term1; //Verified simple case
                        E13g13m = E03g03m;
                        E23g23m = E03g03m;
                        
                        
                        E03g12m = Branches[4].FullP22 * (1.0/3.0 * (1-Branches[5].FullP33))* Branches[6].EgtTi[0] + Branches[4].FullP22* Branches[5].FullP33 * Term2; //Verified simple case
                        E13g02m = E03g12m;
                        
                        
                        E23g01m = (1.0 -Branches[4].FullP22) * Branches[6].EgtTi[0] + Branches[4].FullP22 *  (1.0/3.0 * (1.0 - Branches[5].FullP33)) * Branches[6].EgtTi[0] + Branches[4].FullP22 * Branches[5].FullP33 * Term2; //Verified simple case
                        
                        E02g13m = Branches[4].FullP22 * Branches[5].FullP33* Term2; //Verified simple case
                        E12g03m = E02g13m;
                        // Yes Fullp22 for branch 5 shoudl be there twice, dont want 01 or 12 to coal in branch 5, represents P2pNCoal on branch.
                        
                        E02g02m = Branches[4].FullP22 * Term7 + Branches[4].FullP22 * Branches[5].FullP33* Term1; //Verified simple case
                        E12g12m = E02g02m;
                        
                        E01g01m = Term8 + Branches[4].FullP22* Branches[5].FullP33* Term1 + Branches[4].FullP22* Term7; //Verified simple case
                        
                        E01g23m = Term2 *Branches[4].FullP22 * Branches[5].FullP33; //Verified simple case
                        
                        //Conditional Trios (symmetric)
                        //Q10
                        E01g01m2 = Term8 + Branches[4].FullP22*Term4_5 + Branches[5].FullP33* Branches[4].FullP22 * Term4; //Verified simple case
                        //Q11
                        E01g02m = Branches[4].FullP22 * Term3_5 + 0.5*(1.0 - Branches[5].FullP22*Branches[5].FullP22)*Branches[4].FullP22 * Branches[5].FullP22* Branches[6].EgtTi[0] + Branches[5].FullP33 * Branches[4].FullP22* Term3; //Verified simple case
                        E01g12m = E01g02m;
                        //Q12
                        E02g02m1 = Branches[4].FullP22 * Term4_5 + Branches[4].FullP22 * Branches[5].FullP33 * Term4; //Verfied simple case
                        E12g12m0 = E02g02m1;
                        //Q13
                        E02g01m = (1.0 - Branches[4].FullP22)* Branches[5].EgtTi[0] + Branches[4].FullP22 * Term3_5 + Branches[4].FullP22 * 0.5*(1.0 - Branches[5].FullP22*Branches[5].FullP22)* Branches[5].FullP22* Branches[6].EgtTi[0] + Branches[4].FullP22 * Branches[5].FullP33 * Term3; //Verified simple case
                        E12g01m = E02g01m;
                        //Q14
                        E02g12m = Branches[4].FullP22 * Term3_5 + Branches[4].FullP22 * 0.5*(1.0 - Branches[5].FullP22*Branches[5].FullP22)* Branches[5].FullP22* Branches[6].EgtTi[0] + Branches[4].FullP22 * Branches[5].FullP33 * Term3; //Verified simple case
                        E12g02m = E02g12m;
                        //Q15
                        E01g01m3 = Term8 + Branches[4].FullP22 * Term8_5 + Branches[4].FullP22 * Branches[5].FullP22 * Term4; //Verified simple case
                        //Q16
                        E01g03m = Branches[4].FullP22 * Branches[5].FullP22 * Term3; //Verified simple case
                        E01g13m = E01g03m;
                        //Q17
                        E03g03m1 = Branches[4].FullP22*Branches[5].FullP22 * Term4; //Verified simple case
                        E13g13m0 = E03g03m1;
                        //Q18
                        E03g01m = (1.0 - Branches[4].FullP22)* Branches[6].EgtTi[0] + Branches[4].FullP22 * (1.0 - Branches[5].FullP22) * Branches[6].EgtTi[0] + Branches[4].FullP22 * Branches[5].FullP22 * Term3; //Verified simple case
                        E13g01m = E03g01m;
                        //Q19
                        E03g13m = Branches[4].FullP22 * Branches[5].FullP22 * Term3; //Verified simple case
                        E13g03m = E01g13m;
                        //Q20
                        E02g02m3 = Term8_5 + Branches[5].FullP22 * Term4; //Verified simple case
                        E12g12m3 = E02g02m3;
                        //Q21
                        E02g03m = Branches[5].FullP22 * Term3; //Verified simple case
                        E02g23m = E02g03m;
                        E12g13m = E02g03m;
                        E12g23m = E02g03m;
                        //Q22
                        E13g13m2 = Branches[5].FullP22 * Term4; //Verified simple case
                        E23g23m1 = E13g13m2;
                        E03g03m2 = E13g13m2;
                        E23g23m0 = E13g13m2;
                        //Q23
                        E03g02m = (1.0 - Branches[5].FullP22) * Branches[6].EgtTi[0] + Branches[5].FullP22 * Term3; //Verified simple case
                        E23g02m = E03g02m;
                        E13g12m = E03g02m;
                        E23g12m = E03g02m;
                        //Q24
                        E03g23m =Branches[5].FullP22 * Term3; //Verified simple case
                        E23g03m = E03g23m;
                        E13g23m = E03g23m;
                        E23g13m = E03g23m;
                        
                    }//close asym
                    // cout <<"\n Pair " << QuartMap[0] <<", " <<QuartMap[2] << " | " << QuartMap[1] << ", " << QuartMap[2] << "\n";
                    // cout << E03g23m <<"\n";
                    
                    
                    
                    
                    
                    //
                    //                    cout << "ET01T01 = " << E0101 <<"\n";
                    //                    cout << "ET01T02 = " << E0102 <<"\n";
                    //                    cout << "ET01T03 = " << E0103 <<"\n";
                    //                    cout << "ET01T12 = " << E0112 <<"\n";
                    //                    cout << "ET01T13 = " << E0113 <<"\n";
                    //                    cout << "ET01T23 = " << E0123 <<"\n";
                    //                    cout << "ET02T02 = " << E0202 <<"\n";
                    //                    cout << "ET02T03 = " << E0203 <<"\n";
                    //                    cout << "ET02T12 = " << E0212 <<"\n";
                    //                    cout << "ET02T13 = " << E0213 <<"\n";
                    //                    cout << "ET02T23 = " << E0223 <<"\n";
                    //                    cout << "ET03T03 = " << E0303 <<"\n";
                    //                    cout << "ET03T12 = " << E0312 <<"\n";
                    //                    cout << "ET03T13 = " << E0313 <<"\n";
                    //                    cout << "ET03T23 = " << E0323 <<"\n";
                    //                    cout << "ET12T12 = " << E1212 <<"\n";
                    //                    cout << "ET12T13 = " << E1213 <<"\n";
                    //                    cout << "ET12T23 = " << E1223 <<"\n";
                    //                    cout << "ET13T13 = " << E1313 <<"\n";
                    //                    cout << "ET13T23 = " << E1323 <<"\n";
                    //                    cout << "ET23T23 = " << E2323 <<"\n";
                    //                    cout << "ET01 =" << E01 << "\n";
                    //                    cout << "ET02 =" << E02 << "\n";
                    //                    cout << "ET03 =" << E03 << "\n";
                    //                    cout << "ET12 =" << E12 << "\n";
                    //                    cout << "ET13 =" << E13 << "\n";
                    //                    cout << "ET23 =" << E23 << "\n";
                    
                    
                    //                    // Data type to store the covariance
                    //                     double CovQuartMat [6][6];
                    //                    // How to enumerate over pairs?
                    //                    // 21 unique entries (aka not the symmetric entry in the matrix)
                    //                    //Row/Column order = 01, 02, 03, 12, 13, 23
                    //
                    //                    // Now actually calculate the covariances, and put them in the correct matrix.
                    //                    //recall entries are 0:01 1:02 2:03 3:12 4:13 5:23 ( have a function now to solve this 1 to 1 entry)
                    //                    CovQuartMat[0][0] = E0101 - E01*E01;
                    //                    CovQuartMat[1][1] = E0202 - E02*E02;
                    //                    CovQuartMat[2][2] = E0303 - E03*E03;
                    //                    CovQuartMat[3][3] = E1212 - E12*E12;
                    //                    CovQuartMat[4][4] = E1313 - E13*E13;
                    //                    CovQuartMat[5][5] = E2323 - E23*E23;
                    //                    CovQuartMat[0][1] = CovQuartMat[1][0] = E0102 - E01 * E02;
                    //                    CovQuartMat[0][2] = CovQuartMat[2][0] = E0103 - E01 * E03;
                    //                    CovQuartMat[0][3] = CovQuartMat[3][0] = E0112 - E01 * E12;
                    //                    CovQuartMat[0][4] = CovQuartMat[4][0] = E0113 - E01 * E13;
                    //                    CovQuartMat[0][5] = CovQuartMat[5][0] = E0123 - E01 * E23;
                    //                    CovQuartMat[1][2] = CovQuartMat[2][1] = E0203 - E02 * E03;
                    //                    CovQuartMat[1][3] = CovQuartMat[3][1] = E0212 - E02 * E12;
                    //                    CovQuartMat[1][4] = CovQuartMat[4][1] = E0213 - E02 * E13;
                    //                    CovQuartMat[1][5] = CovQuartMat[5][1] = E0223 - E02 * E23;
                    //                    CovQuartMat[2][3] = CovQuartMat[3][2] = E0312 - E03 * E12;
                    //                    CovQuartMat[2][4] = CovQuartMat[4][2] = E0313 - E03 * E13;
                    //                    CovQuartMat[2][5] = CovQuartMat[5][2] = E0323 - E03 * E23;
                    //                    CovQuartMat[3][4] = CovQuartMat[4][3] = E1213 - E12 * E13;
                    //                    CovQuartMat[3][5] = CovQuartMat[5][3] = E1223 - E12 * E23;
                    //                    CovQuartMat[4][5] = CovQuartMat[5][4] = E1323 - E13 * E23;
                    //
                    //                    cout << "Pairwise Variances = \n";
                    //                    for (int i = 0; i < 6; i++) {
                    //                        cout << CovQuartMat[i][0] << "\n";
                    //                    }
                    //                    cout << "01 v 02 = " << CovQuartMat[0][1] << "\n";
                    //                    cout << "01 v 03 = " << CovQuartMat[0][2] << "\n";
                    //
                    
                    // Now, put them into the CORRECT places for the full covariance matrix
                    //Use the MatLabel function
                    int Zero = MatLabel(QuartMap[0], QuartMap[1], n);
                    int One = MatLabel(QuartMap[0], QuartMap[2], n);
                    int Two = MatLabel(QuartMap[0], QuartMap[3], n);
                    int Three = MatLabel(QuartMap[1], QuartMap[2], n);
                    int Four = MatLabel(QuartMap[1], QuartMap[3], n);
                    int Five = MatLabel(QuartMap[2], QuartMap[3], n);
                    
                    //                                            cout << QuartMap[0] <<", " << QuartMap[1] << " map to " << Zero <<"\n";
                    //                                            cout << QuartMap[0] <<", " << QuartMap[2] << " map to " << One <<"\n";
                    //                                            cout << QuartMap[0] <<", " << QuartMap[3] << " map to " << Two <<"\n";
                    //                                            cout << QuartMap[1] <<", " << QuartMap[2] << " map to " << Three <<"\n";
                    //                                            cout << QuartMap[1] <<", " << QuartMap[3] << " map to " << Four <<"\n";
                    //                                            cout << QuartMap[2] <<", " << QuartMap[3] << " map to " << Five <<"\n";
                    
                    
                    TheoryCov(Zero,Zero) = E0101 - E01*E01;
                    TheoryCov(One,One) = E0202 - E02*E02;
                    TheoryCov(Two,Two) = E0303 - E03*E03;
                    TheoryCov(Three,Three) = E1212 - E12*E12;
                    TheoryCov(Four,Four) = E1313 - E13*E13;
                    TheoryCov(Five,Five) = E2323 - E23*E23;
                    TheoryCov(Zero,One) = TheoryCov(One,Zero) = E0102 - E01 * E02;
                    TheoryCov(Zero,Two) = TheoryCov(Two,Zero) = E0103 - E01 * E03;
                    TheoryCov(Zero,Three) = TheoryCov(Three,Zero) = E0112 - E01 * E12;
                    TheoryCov(Zero,Four) = TheoryCov(Four,Zero) = E0113 - E01 * E13;
                    TheoryCov(Zero,Five) = TheoryCov(Five,Zero) = E0123 - E01 * E23;
                    TheoryCov(One,Two) = TheoryCov(Two,One) = E0203 - E02 * E03;
                    TheoryCov(One,Three) = TheoryCov(Three,One) = E0212 - E02 * E12;
                    TheoryCov(One,Four) = TheoryCov(Four,One) = E0213 - E02 * E13;
                    TheoryCov(One,Five) = TheoryCov(Five,One) = E0223 - E02 * E23;
                    TheoryCov(Two,Three) = TheoryCov(Three,Two) = E0312 - E03 * E12;
                    TheoryCov(Two,Four) = TheoryCov(Four,Two) = E0313 - E03 * E13;
                    TheoryCov(Two,Five) = TheoryCov(Five,Two) = E0323 - E03 * E23;
                    TheoryCov(Three,Four) = TheoryCov(Four,Three) = E1213 - E12 * E13;
                    TheoryCov(Three,Five) = TheoryCov(Five,Three) = E1223 - E12 * E23;
                    TheoryCov(Four,Five) = TheoryCov(Five,Four) = E1323 - E13 * E23;
                    //cout << " Cov(0,1) =  " << TheoryCov(0,1) <<"\n";
                    
                    TheoryMean(Zero,0) = E01;
                    TheoryMean(One,0) = E02;
                    TheoryMean(Two,0) = E03;
                    TheoryMean(Three,0) = E12;
                    TheoryMean(Four,0) = E13;
                    TheoryMean(Five,0) = E23;
                    
                    //cout << "Theory Mean  = " << TheoryMean.transpose()<<"\n";
                    //cout <<"Theory Cov = \n"<< TheoryCov <<"\n";
                    //Shared Branch length calculations
                    //Self
                    TheorySharedBranchLength(Zero,Zero) = 2.0*E01;
                    TheorySharedBranchLength(One,One) = 2.0* E02;
                    TheorySharedBranchLength(Two,Two) = 2.0 * E03;
                    TheorySharedBranchLength(Three, Three) = 2.0 * E12;
                    TheorySharedBranchLength(Four,Four) = 2.0 * E13;
                    TheorySharedBranchLength(Five,Five) = 2.0 * E23;
                    
                    //Triplets ( Check that this formula makes sense.)
                    TheorySharedBranchLength(Zero,One) =   2.0 * E01g12m - E12g12m0 + E01g01m2 + E02g02m1;
                    TheorySharedBranchLength(Zero,Two) =   2.0 * E01g13m - E13g13m0 + E01g01m3 + E03g03m1;
                    TheorySharedBranchLength(Zero,Three) = 2.0 * E01g02m - E02g02m1 + E01g01m2 + E12g12m0;
                    TheorySharedBranchLength(Zero,Four) =  2.0 * E01g03m - E03g03m1 + E01g01m3 + E13g13m0;
                    TheorySharedBranchLength(One,Two) =    2.0 * E02g23m - E23g23m0 + E02g02m3 + E03g03m2;
                    TheorySharedBranchLength(One,Three) =  2.0 * E02g01m - E01g01m2 + E02g02m1 + E12g12m0;
                    TheorySharedBranchLength(One,Five) =   2.0 * E02g03m - E03g03m2 + E02g02m3 + E23g23m0;
                    TheorySharedBranchLength(Two,Four) =   2.0 * E03g01m - E01g01m3 + E03g03m1 + E13g13m0;
                    TheorySharedBranchLength(Two,Five) =   2.0 * E03g02m - E02g02m3 + E03g03m2 + E23g23m0;
                    TheorySharedBranchLength(Three,Four) = 2.0 * E12g23m - E23g23m1 + E12g12m3 + E13g13m2;
                    TheorySharedBranchLength(Three,Five) = 2.0 * E12g13m - E13g13m2 + E12g12m3 + E23g23m1;
                    TheorySharedBranchLength(Four,Five) =  2.0 * E13g12m - E12g12m3 + E13g13m2 + E23g23m1;
                    
                    TheorySharedBranchLength(One,Zero) = TheorySharedBranchLength(Zero,One);
                    TheorySharedBranchLength(Two,Zero) = TheorySharedBranchLength(Zero,Two);
                    TheorySharedBranchLength(Three,Zero) = TheorySharedBranchLength(Zero,Three);
                    TheorySharedBranchLength(Four,Zero) = TheorySharedBranchLength(Zero,Four);
                    TheorySharedBranchLength(Two,One) = TheorySharedBranchLength(One,Two);
                    TheorySharedBranchLength(Three,One) = TheorySharedBranchLength(One,Three);
                    TheorySharedBranchLength(Five,One) = TheorySharedBranchLength(One,Five);
                    TheorySharedBranchLength(Four,Two) = TheorySharedBranchLength(Two,Four);
                    TheorySharedBranchLength(Five,Two) = TheorySharedBranchLength(Two,Five);
                    TheorySharedBranchLength(Four,Three) = TheorySharedBranchLength(Three,Four);
                    TheorySharedBranchLength(Five,Three) = TheorySharedBranchLength(Three,Five);
                    TheorySharedBranchLength(Five,Four) = TheorySharedBranchLength(Four,Five);
                    
                    
                    //Two pairs of unique indiviudals.
                    
                    
                    double EMinCoal = E01g01m + E23g01m + E01g23m + E23g23m  + E02g02m + E13g02m + E02g13m + E13g13m + E03g03m + E12g03m + E03g12m + E12g12m;
                    //cout <<"EMin Coal = " << EMinCoal <<"\n";
                    
                    TheorySharedBranchLength(Zero,Five) = E01 + E23 - EMinCoal;
                    TheorySharedBranchLength(One,Four) = E02 + E13 - EMinCoal;
                    TheorySharedBranchLength(Two,Three) = E03 + E12 - EMinCoal;
                    
                    TheorySharedBranchLength(Five,Zero) = TheorySharedBranchLength(Zero,Five);
                    TheorySharedBranchLength(Four,One) = TheorySharedBranchLength(One,Four);
                    TheorySharedBranchLength(Three,Two) = TheorySharedBranchLength(Two,Three);
                    
                    //cout <<" Shared Branch Length Matrix \n" << TheorySharedBranchLength <<"\n ";
                    // Stuff for the LogNormal bias correction:
                    

                    
                    
                    


                    
                    
                }}}} // close the for loop through all quartets.
    
    
    
} // END THEORY FUNCTION.


void Process_Input( string i_file, int & n_species, SpeciesTreeType & Species_Tree){
    string line_;
    int counter = 0;
    int branch = 0;
    double time = 0;
    int joiner = 0;
    int joinee = 0;
    Eigen::VectorXi locations;
    
    //if input file exists, open it, and parse the contents
    ifstream file_ (i_file);
    
    if(file_.is_open()){
        std::string line_;
        while(getline(file_, line_)){
            
            if(line_[0] == '#' || line_.empty())
                continue;
            
            if(line_.find("=") != string::npos){
               // cout << line_ <<"\n";

                //remove spaces
            line_.erase(std::remove_if(line_.begin(), line_.end(), ::isspace), line_.end());
            auto delimiterPos = line_.find("=");
            auto value = line_.substr(delimiterPos + 1);
            n_species= std::stoi(value);
            //cout  << value << " " << n_species <<"\n";
                // declare the right amount of branches
            Species_Tree.Branches.resize(2*n_species -1);
            locations.resize(n_species);
                continue;
                
            }
            if (line_.find("-n")!= string::npos){
                counter = counter + 1;
                char *cstr = &line_[0u];
                int b = 0;
                while(isblank(cstr[b])==true){
                    b++;
                }
                while(isblank(cstr[b])==false){
                    b++;
                }
                while(isblank(cstr[b])==true){
                    b++;
                }
                branch = atoi(&cstr[b]) -1;
                //cout <<"branch "<< branch <<"\n";
                locations[branch] = branch;
                b=b+1;

                while(isblank(cstr[b])==true){
                    b++;
                }

                Species_Tree.Branches[branch].Ne = atof(&cstr[b]);
                Species_Tree.Branches[branch].Times[0] = 0;

                //cout << "pop size = " << atof(&cstr[b]) <<"\n";
                continue;
                
            }
            if (line_.find("-t") != string::npos){

                
                //cout <<counter <<"\n";
                
                char *cstr = &line_[0u];
                int b = 0;
                while(isblank(cstr[b])==true){
                    b++;
                }
                while(isblank(cstr[b])==false){
                    b++;
                }
                while(isblank(cstr[b])==true){
                    b++;
                }
                Species_Tree.Branches[counter].Times[0] =  atof(&cstr[b]);
                time =atof(&cstr[b]);
                
                b=b+1;
                while(isblank(cstr[b])==false){
                    b++;
                }
                while(isblank(cstr[b])==true){
                    b++;
                }
                joiner = atoi(&cstr[b])-1;
                
                b=b+1;
                
                while(isblank(cstr[b])==true){
                    b++;
                }
                joinee = atoi(&cstr[b])-1;
                
                b=b+1;
                
                while(isblank(cstr[b])==true){
                    b++;
                }
                Species_Tree.Branches[counter].Ne =  atof(&cstr[b]);

            }
           // cout <<"joiner " << joiner  <<" " << locations[joiner] <<"\n";
            //cout <<"joinee " << joinee  <<" " << locations[joinee] <<"\n";

            Species_Tree.Branches[locations[joiner]].Parent= counter;
            Species_Tree.Branches[locations[joinee]].Parent= counter;
            Species_Tree.Branches[locations[joiner]].Times[1]= time;
            Species_Tree.Branches[locations[joinee]].Times[1]= time;
            Species_Tree.Branches[counter].Offspring[0] = locations[joiner];
            Species_Tree.Branches[counter].Offspring[1] = locations[joinee];

            locations[joiner] = counter;
            locations[joinee] = counter;
            if(counter+1 == 2*n_species -1){
                Species_Tree.Branches[counter].Parent = 0;

            }
            counter = counter +1;

        }
    }
    else{
        std::cerr<<"Couldn't open file for reading.\n";
        
    }
}




Matrix<int, Dynamic, 2> MapMaker(int n_species){
    
    Matrix<int, Dynamic, 2 > TheMap;
    TheMap.resize(4*n_species, 2);

            for(int l = 0; l < n_species; l++){
                for(int i =0; i < 4; i++){
                    TheMap(l*4+i,1)  = l;

                }
            }
    
        for(int l =0; l < (n_species * 4); l++){
            TheMap(l, 0) = l;
            
        }

    return TheMap;
}

//Final Tree Printer, a way to view the tree.
void TreePrinter(SpeciesTreeType MyTree, int n){
    
    cout <<"================================================\n";
    for(int branch = 2*n-2; branch >=0; branch--){
        cout <<"----------------------------------------\n";
        
        cout<< "Branch " << branch <<":\n";
        
        if(branch != 2*n-2){
        cout << "    Parent : " << MyTree.Branches[branch].Parent <<"\n";
        cout<< "    End Time : " << MyTree.Branches[branch].Times[1] <<"\n";
        }
        cout<< "    Start Time : " << MyTree.Branches[branch].Times[0] <<"\n";
        
        cout << "    Pop Size : "<< MyTree.Branches[branch].Ne<<"\n";
        if(branch >= n){
            cout << "    Offspring : " << MyTree.Branches[branch].Offspring[0] << ", " << MyTree.Branches[branch].Offspring[1]<< "\n";
        }
        
        
    }
    cout <<"================================================\n";
    
}


void get_EigentoData(Matrix<double, Dynamic, Dynamic> src, char* pathAndName)
{
    ofstream fichier(pathAndName, ios::out | ios::trunc);
    if(fichier)
    {
        // instructions
        fichier << src << "\n";
        fichier.close();
    }
    else  // sinon
    {
        cerr << "Unable to write matrix to file." << endl;
    }
}

// Store Vectors into files.
void get_EigenVectoData(Matrix<double, Dynamic, 1> src, char* pathAndName)
{
    //changed src to dynamic from n? 4/17
    ofstream fichier(pathAndName, ios::out | ios::trunc);
    if(fichier)
    {
        // instructions
        fichier << src << "\n";
        fichier.close();
    }
    else  // sinon
    {
        cerr << "Unable to write vector to file." << endl;
    }
}

string IndivLabel(int indiv, int n_species){
    
    string lab;
    
    double indD = indiv/4.0;
    int species = 0;
    while(indD >= species+1 ){
        species = species+1;
    }
    if(indD-species ==0){
        lab = to_string(species) + "_a";
    }else if(indD-species == 0.25){
        lab = to_string(species) + "_b";
    }else if(indD-species == 0.5){
        lab = to_string(species) + "_c";
    }else if(indD-species == 0.75){
        lab = to_string(species) + "_d";
    }
    
    
    return lab;
}

void VectorLabels(int n_species){
    
    vector<string> pair_labels;
    string pair_lab;
    
    int n_indivs = 4* n_species;
    for(int i1 =0; i1 < n_indivs-1; i1++){
        for(int i2 = i1+1 ; i2 < n_indivs; i2++){
            pair_lab =to_string(MatLabel(i1, i2, n_indivs)+1)+"\t"+ IndivLabel(i1, n_species) + "," + IndivLabel(i2, n_species);
            pair_labels.push_back(pair_lab);
            
        }
    }
    char LabelNames[15] = "Labels.txt";
    std::ofstream OutFile(LabelNames);
    for(const auto &e : pair_labels) OutFile << e<< "\n";
    
    
}



int main(int argc, char **argv){
    
    // We will take in one command from the user, the name of an executable file. This will basically just have a tree command that we need to parse. Parsing it should not be too hard.
    if (argc !=2){
        cout << "Please provide the location and name of the input control file. Only one file can be read at a time. \n";
        
    }
    
    string InputFile = argv[1];
    
    char TheoryMeanName[15] = "Means.txt";
    char TheoryCovName[15] = "Covariance.txt";
    char TheorySBLName[25] = "SharedBranchLength.txt";

    
    int n_species;
    
    SpeciesTreeType Tree;
    
    Process_Input(InputFile, n_species, Tree);
    
    cout << "Inputted Tree: \n";
    TreePrinter(Tree, n_species);
    
    Matrix<int, Dynamic, 2> MapFile = MapMaker(n_species);
        
    int n_indivs = 4* n_species;
    int nc2 = n_indivs * (n_indivs -1)/2.0;
    
    SpeciesTreeType ExtendedTree;
    ExtendedTree.Branches.resize(2*n_indivs-1);
    TreeExtender(Tree, ExtendedTree, n_species, n_indivs, MapFile);
    Matrix<double, Dynamic,1> TheoryMean;
    Matrix<double, Dynamic,Dynamic> TheoryCov;
    Matrix<double, Dynamic, Dynamic> TheorySharedBranchLength;
    TheorySharedBranchLength.resize(nc2,nc2);
    TheoryMean.resize(nc2,1);
    TheoryCov.resize(nc2,nc2);
    TheoryFunction(ExtendedTree, TheoryMean, TheoryCov, n_indivs, TheorySharedBranchLength);
    
    get_EigenVectoData(TheoryMean, TheoryMeanName);
    get_EigentoData(TheoryCov, TheoryCovName);
    get_EigentoData(TheorySharedBranchLength, TheorySBLName);
    
    VectorLabels(n_species);
    
    
    return 0;
    
}



