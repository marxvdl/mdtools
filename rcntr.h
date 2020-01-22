/*
 * hmmpred_parse.h
 *
 *  Created on: 01/03/2011
 *      Author: marx
 */

#ifndef RCNTR_H_
#define RCNTR_H_

#include "misc.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

vector< vector<real> > parseHmmPredOutput(string filename);
vector<real> pickMatrixLine(vector< vector<real> > & matrix, uint col);

#endif /* RCNTR_H_ */
