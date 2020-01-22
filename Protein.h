/*
 * PDB2Levels.h
 *
 *  Created on: Mar 19, 2010
 *      Author: marx
 */

#ifndef PDB2LEVELS_H_
#define PDB2LEVELS_H_

#include "Atom.h"

#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <iostream>
#include <map>
#include <boost/concept_check.hpp>

using namespace std;

class Protein {
private:
	//varying layer width
	real midLayerWidth;
	real bottomLayerWidth;
	//end
	
	bool keepCenterAtOrigin;
	real centerX, centerY, centerZ;

	vector<real> layerBoundaries;
	vector<real> layerDistances;
	vector<real> layerDeltas;

	void calculateLayerInfo(int layers);
	char getLayerCharFromBurial(real burial);
	int getLayerNumberFromBurial(real burial);

	int unkBuffer;
	bool firstIteratorRound;
	int lastResidue;
	unsigned int residueIteratorAtom;
	void initResidueIterator();
	string nextResidue();

	string sequenceMapping(map<string, char>& mapping);
	real getRadiusOfGyration();
	bool loadLayerInfoFromBData(int nlayers, string referenceFile);
	
	string addRamaLine(ostream& output, vector< map<string,uint> >& atomMap, uint x, string a1, uint y, string a2, uint z, string a3, uint k, string a4, const string STR );
	real quadLimit(real delta){
		return 2.0 / pow2(delta);
	}

public:
	vector<Atom> atoms;
	real gyrRadius;

	void loadFromMDRead(istream& input);
	void loadFromRg5DB(istream& input);
	void loadFromPDB(istream& input);

	void generateRcntr(int layers, bool quad, ostream& output, string referenceAtom);
	void generateRcntr(int layers, bool quad, vector< vector<real> >& probs, ostream& output, string referenceAtom);
	void generateDat(int layers, ASM atomSizes, real scale, ostream& output);
	void generateGar(int layers, ostream& output);
	void generatePDB(int layers, ostream& output);
	void generateHbonds(ostream& output, bool all);
	void generateChiral(ostream& output);
	void generateRama(ostream& output);
	void generateHistogram(ostream& output, const string& referenceAtom);
	void generateBData(ostream& output, string referenceFile="");
	void generateClashReport(ostream& output);

	string toBurialSequence(int layers, string referenceAtom, bool normalize, bool rg5);
	string toBurialSequenceCAdir(int layers=2, bool invert=false);
	string toBurialSequenceAllAtoms(int layers, bool normalize);
	string toRelativeBurialSequence(int layers, string referenceAtom);

	string to3Letters();
	string toHP();
	string toHPN();
	string toHPNAC();
	string toLetter();
	string toHPNhp();
	string toHPNhpn();
	
	Protein(bool k, real mlwidth, real bwwidth) : 
		keepCenterAtOrigin(k), 
		midLayerWidth(mlwidth),
		bottomLayerWidth(bwwidth)
		{}
	
};

#endif /* PDB2LEVELS_H_ */

