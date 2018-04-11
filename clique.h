/*#ifndef 3CLIQUE_H
#define 3CLIQUE_H*/

#include <map>
#include <cassert>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <fstream>
#include <set>
#include <ctime>
#include <list>
#include <deque>

#include <memory.h>
#include <math.h>

using namespace std;


#define AVE 1
int const NC = 39142;
int  const NN = 39142;
int  const N=95506;
int const k=3;
int const pow_cut = 65;
double alpha =1.8;

double DIS[21]={
   0,
0.0484,
0.04803,
0.0615,
0.06797,
0.0754,
0.06635,
0.06437,
0.06153,
0.05795,
0.06392,
0.05077,
0.04821,
0.04365,
0.04152,
0.0407,
0.03632,
0.0333,
0.03152,
0.02909,
0.02951
};

double Pro[21];
int dis_nodes[N];


double S1[NN];
double M[NN];
int LABEL_[N];
vector <int > ptr;
vector <int> cliquenodes;
typedef  std::vector<std::vector <int> > Edge;

#include "mt19937ar.h"
extern mt19937 rand_num;
typedef  std::vector<std::set <int> > Net; // needs laboratory's library

int findroot(int i);
int kCLIQUEFIND(int v1,int v2);
void percolation_32(int *labels_1);
int CLUSTER(int newnum,int v1,int v2);
void ClusterCombine(int r1,int r2);

int percolate_edge_author(FILE *fp);


//#endif // 3CLIQUE_H
