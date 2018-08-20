#ifndef INCLUDES_H
#define INCLUDES_H

using namespace std;

// External includes
#include <algorithm>
#include <iostream>
#include <map>
#include <math.h>
#include <stdio.h>
#include <string>
#include <thread>
#include <utility>
#include <thread>
#include <vector>
#include <fstream>
#include <functional>
#include <omp.h>

// Internal includes
#include "functions.h"
#include "numberparameters.h"


#define MAXIMIZATION 1
#define MINIMIZATION 0


// Typedefs
typedef vector<bool> bit_vector;
typedef function <double(vector<double>)> opt_function;
#endif // INCLUDES_H
