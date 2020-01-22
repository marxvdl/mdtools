/*
 * Copyright 2010-2020, Marx Gomes van der Linden
 *                      marx.linden@ifb.edu.br
 * 
 * This file is part of HmmPred.
 * 
 * HmmPred is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * HmmPred is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with HmmPred.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "Protein.h"
#include "rcntr.h"

#include <cassert>
#include <algorithm>

/**
 * Generates the "rcntr" file.
 */
void Protein::generateRcntr(int layers, bool quad, ostream& output, string referenceAtom){
	int i = 1;
		
	if(layers == 0){
		output << atoms.size() << " distances from center " << centerX << " " << centerY << " " << centerZ << endl;
		output << "1.000 1.000 2000" << endl;

		for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it)
			if((referenceAtom == "") or (it->type == referenceAtom)) {
				if(quad)
					output << it->index << " " << it->distanceToCenter << " 100.00" << endl;
				else
					output << it->index << " " << it->distanceToCenter << " 0.000 1.000" << endl;
			}
			
	}
	else if(layers == 1){
		output << atoms.size() << " distances from center " << centerX << " " << centerY << " " << centerZ << endl;
		output << "1.000 1.000 2000" << endl;
		
		real mult = //2.7 * cuberoot( (real) atoms[atoms.size()-1].residueIndex );
		            2.7 * cuberoot(        atoms[atoms.size()-1].residueIndex);

		for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it)
			if((referenceAtom == "") or (it->type == referenceAtom)) {
				if(quad)
					output << it->index << " " << it->distanceToCenter/mult << " 100.00" << endl;
				else
					output << it->index << " " << it->distanceToCenter/mult << " 0.000 1.000" << endl;
			}
			
	}
	else{
		calculateLayerInfo(layers);

		stringstream ss;
		for(uint j=1; j<layerBoundaries.size(); j++)
			ss << " " << layerBoundaries[j];

		output << atoms.size() << " distances from center " << centerX << " " << centerY << " " << centerZ << "   |   boundaries:" << ss.str() << "  |  size: " << toLetter().size() << endl;
		output << "1.000 1.000 2000" << endl;

		for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it){
			int layer = getLayerNumberFromBurial(it->distanceToCenter);
			
			if((referenceAtom == "") or (it->type == referenceAtom)) {				
				if(quad){
					if(layer == 0)
						output << it->index << " " << 0.0 << " " << quadLimit(layerDeltas[layer] * 2) << " 1.000" << endl;
					else
						output << it->index << " " << layerDistances[layer] << " " << quadLimit(layerDeltas[layer]) << " 1.000" << endl;
				}
				else {
					output << it->index << " " << layerDistances[layer] << " " << layerDeltas[layer] << " 1.000" << endl;
				}
			}
			
		}
	}
}

/**
 * Generates the "rcntr" file according to predicted probabilities.
 */
void Protein::generateRcntr(int layers, bool quad, vector< vector< double > >& probs, ostream& output, string referenceAtom){
	real maxDist =
		2.7 * pow( (real)probs.size(), 1.0/3.0);
		//2.7 * cuberoot(atoms[atoms.size()-1].residueIndex);

	output << probs.size() << " distances from center " << centerX << " " << centerY << " " << centerZ << "    # radius: " << maxDist << endl;
	output << "1.000 1.000 2000" << endl;

	vector<real> prob1 = pickMatrixLine(probs, 1);
	vector<real> prob2 = pickMatrixLine(probs, 2);
	vector<real> prob3 = pickMatrixLine(probs, 3);
	uint i = 1;
	for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it){
// 		real regression = 0.488 + 0.332*prob1[i-1] + 0.586*prob2[i-1] + 0.996*prob3[i-1];  //CA
		real regression = 0.499 + 0.335*prob1[i-1] + 0.566*prob2[i-1] + 0.949*prob3[i-1]; //CB
		
		if(it->type == referenceAtom || (referenceAtom == "CB" && it->residueType == "GLY" && it->type == "CA")){
			real dist = regression * maxDist;
			if(dist < 0)
				dist = 0;
			
			if(quad)
				output << it->index << " " << dist << " 100.00" << endl;
			else
				output << it->index << " " << dist << " 0.000 1.000" << endl;
			
			i++;
		}
	}

	assert(i-1 == probs.size());


}

/**
 * Reads the first entry in a HmmPred results file with probabilities
 * and returns it as a matrix of probabilities.
 */
vector< vector<real> > parseHmmPredOutput(string filename){
	ifstream input(filename.c_str());

	vector< vector<real> > probs;

	if(input.fail()){
		cout << "Load from HmmPred result fail!" << endl;
		return probs;
	}

	string line;
	while( getline (input, line) ){
		if(line[0] == '>')
			break;
	}

	string priSeq;
	getline (input, priSeq);
	getline (input, line);

	while( getline (input, line) ){
		if(line.size() == 0)
			break;

		istringstream iss(line);

		vector<real> probReals;
		real sub;
		while (iss >> sub)
			probReals.push_back(sub);

		probs.push_back(probReals);
	}

	assert(priSeq.length() == probs.size());

	return probs;
}


vector<real> pickMatrixLine(vector< vector<real> > & matrix, uint line){
	vector<real> vec(matrix.size(), 0.0);
	for(uint i=0; i<matrix.size(); i++)
		vec[i] = matrix[i][line];
	return vec;
}


/**
 * Generates the "gar" file.
 */
void Protein::generateGar(int layers, ostream& output){
	if(layers == 0){
		for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it)
			output << "burial\t" << it->residueIndex << " " << it->type << "\t" << it->distanceToCenter << " 0 0\n";
	}
	else{
		calculateLayerInfo(layers);

		stringstream ss;
		for(uint j=1; j<layerBoundaries.size(); j++)
			ss << " " << layerBoundaries[j];

		for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it){
			int layer = getLayerNumberFromBurial(it->distanceToCenter);
			output << "burial\t" << it->residueIndex << " " << it->type << "\t" << layerDistances[layer] << " " << layerDeltas[layer] << " " << layerDeltas[layer] << "\n";
		}
	}
}
