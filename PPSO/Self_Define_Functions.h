#pragma once
#pragma once
#ifndef SELF_DEFINE_FUNCTIONS_H_INCLUDED
#define SELF_DEFINE_FUNCTIONS_H_INCLUDED


#include <time.h>
#include <cstdio>
#include "unistd.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <math.h>
#include <string>
#include <string.h>

// The following is the library of random number generators in boost
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/cauchy_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>



//the settings of the program parameters
const int dim = 30;//the dimension size of the problem to be optimized
const int MAX_FV = 10000*dim;//the maximum number of fitness evaluations

//the settings of parameters in PSO
const int Population_Size = 256;
const double weight = 1.0;


void cec14_test_func(double* x, double* f, int nx, int mx, int func_num);
double cal_glbest(int index, double** gbest, double* c, int N, double** population, int** subpopulation, int h,int t,int ***gbest_index);
void replaceinfo(double* results, int** Buffer, int* subpopulation, int hierarchical, double* subpopulation_fitness, double* buffer_field_fitness, int subpopulationsize, int length, double* temp_buffer_field_result, double* fitness2);
int search_index(double x, double* fitness, int length);
int search_max_num(double* goal, int length);
#endif // SELF_DEFINE_FUNCTIONS_H_INCLUDED
