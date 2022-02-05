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
#include "DatFile.h"
#include "ClashFinder.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <limits>
#include <cassert>
#include <set>

/**
 * Loads a protein structure from and MD file (generated by the READ program).
 */
void Protein::loadFromMDRead(istream& input){
	if(input.fail()){
		cout << "Load from MD fail!" << endl;
		return;
	}

	string line;

	real totalX = 0, totalY = 0, totalZ = 0;

	while( getline (input, line) ) {
		istringstream linein(line);

		if(line[12] == 'N' || line[12] == 'C' || line[12] == 'O' || line[12] == 'S'){
			Atom atom;
			linein >> atom.index >> atom.residueIndex >> atom.type >> atom.residueType
			>> atom.x >> atom.y >> atom.z;
			totalX += atom.x;
			totalY += atom.y;
			totalZ += atom.z;
			atoms.push_back(atom);
		}
	}

	if(keepCenterAtOrigin){
		centerX = 0.0;
		centerY = 0.0;
		centerZ = 0.0;
	}
	else{
		centerX = totalX / atoms.size();
		centerY = totalY / atoms.size();
		centerZ = totalZ / atoms.size();
	}

	for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it){
		it->calculateRealBurial(centerX, centerY, centerZ);
	}
}

/**
 * Loads a protein structure from a ".rg5" database file.
 */
void Protein::loadFromRg5DB(istream& input){
	if(input.fail()){
		cout << "Load from RG5 fail!" << endl;
		return;
	}

	string line;
	//    1 N   ARG   578 14.6316 1.19765 1.26261
	while( getline (input, line) ) {
		istringstream linein(line);

		Atom atom;
		atom.index = 0;
		linein >> atom.index;
		if(atom.index != 0){
			linein >> atom.type >> atom.residueType >> atom.residueIndex
			>> atom.distanceToCenter >> atom.burial;
			atom.x = atom.y = atom.z = 0;
			atoms.push_back(atom);
		}
	}
}



/**
 * Loads a protein structure from a PDB file.
 */
void Protein::loadFromPDB(istream& input){
	if(input.fail()){
		cout << "Load from PDB fail!" << endl;
		return;
	}

	string line;

	real totalX = 0, totalY = 0, totalZ = 0;

	while( getline (input, line) ) {

		//		COLUMNS        DATA  TYPE    FIELD        DEFINITION
		//		-------------------------------------------------------------------------------------
		//		 1 -  6        Record name   "ATOM  "
		//		 7 - 11        Integer       serial       Atom  serial number.
		//		13 - 16        Atom          name         Atom name.
		//		17             Character     altLoc       Alternate location indicator.
		//		18 - 20        Residue name  resName      Residue name.
		//		22             Character     chainID      Chain identifier.
		//		23 - 26        Integer       resSeq       Residue sequence number.
		//		27             AChar         iCode        Code for insertion of residues.
		//		31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
		//		39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
		//		47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
		//		55 - 60        Real(6.2)     occupancy    Occupancy.
		//		61 - 66        Real(6.2)     tempFactor   Temperature  factor.
		//		77 - 78        LString(2)    element      Element symbol, right-justified.
		//		79 - 80        LString(2)    charge       Charge  on the atom.

		if(line[0] == 'A' && line[1] == 'T' && line[2] == 'O' && line[3] == 'M'){
			Atom atom;
			atom.index = string2int( line.substr(6,5) );
			atom.type = trim( line.substr(12,4) );
			atom.residueType = line.substr(17,3);
			atom.residueIndex = string2int( line.substr(22,4) );
			atom.x = string2real( line.substr(30,8) );
			atom.y = string2real( line.substr(38,8) );
			atom.z = string2real( line.substr(46,8) );
			totalX += atom.x;
			totalY += atom.y;
			totalZ += atom.z;
			atoms.push_back(atom);
		}
	}

	if(keepCenterAtOrigin){
		centerX = 0.0;
		centerY = 0.0;
		centerZ = 0.0;
	}
	else{
		centerX = totalX / atoms.size();
		centerY = totalY / atoms.size();
		centerZ = totalZ / atoms.size();
	}

	for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it){
		it->calculateRealBurial(centerX, centerY, centerZ);
	}
}

real Protein::getRadiusOfGyration(){
	real n = atoms.size();
	real avgX = 0;
	real avgY = 0;
	real avgZ = 0;		
	for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it){
		avgX += it->x;
		avgY += it->y;
		avgZ += it->z;			
	}
	avgX /= n;
	avgY /= n;
	avgZ /= n;		
	
	real sum = 0;
	for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it)
		sum += pow2(it->x - avgX) + pow2(it->y - avgY) + pow2(it->z - avgZ);
	return sqrt(sum / n);
}

/**
 * Returns a sequence of burial layers for all atoms.
 */
string Protein::toBurialSequenceAllAtoms(int layers, bool normalize){
	stringstream seq;
	
	if(layers == 1){
		//real radgyr = getRadiusOfGyration();
		real radgyr = 2.7 * cuberoot(atoms[atoms.size()-1].residueIndex);
		
		for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it)
			it->distanceToCenter /= radgyr;
	}	
	else if(layers == 0 && normalize){
		real max = -1;
		for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it)
			if(it->distanceToCenter > max)
				max = it->distanceToCenter;
			assert(max != -1);
		for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it)
			it->distanceToCenter /= max;
	}
	
	if(layers == 0 || layers == 1){
		for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it)
			seq << it->distanceToCenter << "\n";
	}
	else{	
		calculateLayerInfo(layers);
		
		for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it)
			seq << getLayerCharFromBurial(it->distanceToCenter);
	}
	
	return seq.str();
}


/**
 * Returns a sequence of burial layers of atoms of type 'referenceAtom'.
 */
string Protein::toBurialSequence(int layers, string referenceAtom, bool normalize, bool rg5){
	stringstream seq;
	
	if(referenceAtom == "all"){
		if(rg5)
			return "Option not supported for RG5\n";
		else
			return toBurialSequenceAllAtoms(layers, normalize);
	}
	
	if(!rg5 && layers == 1){
		//real radgyr = getRadiusOfGyration();
		real radgyr = 2.7 * cuberoot(atoms[atoms.size()-1].residueIndex);
		
		for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it)
			it->distanceToCenter /= radgyr;
	}	

	if(layers == 0 && normalize){
		real max = -1;
		for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it)
			if(it->type == referenceAtom && it->distanceToCenter > max)
				max = it->distanceToCenter;
		assert(max != -1);
		for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it)
			it->distanceToCenter /= max;
	}
	
	if(layers == 0 || layers == 1){
		bool foundRefAtom = false;

		int lastResidue = atoms[0].residueIndex;
		for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it){
			if(it->residueIndex != lastResidue){
				if(!foundRefAtom)
					seq << "-\n";
				foundRefAtom = false;

				int gap = it->residueIndex - lastResidue - 1;
				for(int i=0; i<gap; i++)
					seq << "-\n";

				lastResidue = it->residueIndex;

			}

			if(!foundRefAtom && 
				( it->type == referenceAtom 
				|| 
					//GLY defaults to CA
					(it->residueType == "GLY" && it->type == "CA")
				)
			){
				if(rg5 && layers == 1)
					seq << it->burial << endl;
				else
					seq << it->distanceToCenter << endl;
				foundRefAtom = true;
			}
		}

		if(!foundRefAtom)
			seq << "-\n";
	}
	else{

		//seq << layerBoundaries[0] << " " << layerBoundaries[1] << " " << layerBoundaries[2] << endl;
		calculateLayerInfo(layers);

		bool foundRefAtom = false;

		int lastResidue = atoms[0].residueIndex;
		for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it){
			if(it->residueIndex != lastResidue){
				if(!foundRefAtom)
					seq << 'X';
				foundRefAtom = false;

				int gap = it->residueIndex - lastResidue - 1;
				for(int i=0; i<gap; i++)
					seq << 'X';

				lastResidue = it->residueIndex;

			}

		if(!foundRefAtom && 
			( it->type == referenceAtom 
			|| 
			//GLY defaults to CA
			(it->residueType == "GLY" && it->type == "CA")
			)
			){
				seq << getLayerCharFromBurial(it->distanceToCenter);
				foundRefAtom = true;
			}
		}

		if(!foundRefAtom)
			seq << 'X';
	}

	return seq.str();
}

/**
 * Returns a sequence of burial layers (2 layers) of atoms of type CA, plus the direction (given by CB).
 */
string Protein::toBurialSequenceCAdir(int layers, bool invert){
	if(layers == 0)
		layers = 2;
	
	vector<real> CAburial;
	vector<real> CBburial;	
	
	calculateLayerInfo(layers);
	
	///////////
	///  CA  //
	///////////
	bool foundCA = false;

	int lastResidue = atoms[0].residueIndex;
	for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it){
		if(it->residueIndex != lastResidue){
				if(!foundCA)
				CAburial.push_back(-1);
			foundCA = false;

			int gap = it->residueIndex - lastResidue - 1;
			for(int i=0; i<gap; i++)
				CAburial.push_back(-1);

			lastResidue = it->residueIndex;

		}

		if(!foundCA && it->type == "CA"){
				CAburial.push_back(it->distanceToCenter);
			foundCA = true;
		}
	}
	if(!foundCA)
		CAburial.push_back(-1);
	
	///////////
	///  CB  //
	///////////
	bool foundCB = false;
	
	lastResidue = atoms[0].residueIndex;
	for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it){
		if(it->residueIndex != lastResidue){
			if(!foundCB)
				CBburial.push_back(-1);
			foundCB = false;
			
			int gap = it->residueIndex - lastResidue - 1;
			for(int i=0; i<gap; i++)
				CBburial.push_back(-1);
			
			lastResidue = it->residueIndex;
			
		}
		
		if(!foundCB && it->type == "CB"){
			CBburial.push_back(it->distanceToCenter);
			foundCB = true;
		}
	}
	if(!foundCB)
		CBburial.push_back(-1);
	
	///////////////////////
	uint size = CBburial.size();
	assert(size = CAburial.size());
	if(CBburial[size-1] == -1 && CAburial[size-1] == -1){
		CAburial.pop_back();
		CBburial.pop_back();		
	}

	///////////////////////
		
	if(invert){
		//Assemble CB + direction of CA
		stringstream seq;
		for(uint i=0; i<CAburial.size(); i++){
			if(CAburial[i] == -1){
				seq << 'X';
			}
			else{
				if(CBburial[i] == -1)
					if(layers == 2)
						seq << getLayerCharFromBurial(CAburial[i]);
					else
						switch(getLayerCharFromBurial(CAburial[i])){
							case '0': seq << 'a'; break;
							case '1': seq << 'c'; break;
							case '2': seq << 'e'; break;
							case '3': seq << 'g'; break;
							case '4': seq << 'i'; break;
							case '5': seq << 'k'; break;
						}
						else{
							bool outwards = CAburial[i] >= CBburial[i];
							char cblayer = getLayerCharFromBurial(CBburial[i]);
							
							if(layers == 2){
								if(cblayer == '0')
									seq << (outwards? '^' : '0');
								else
									seq << (outwards? '1' : 'v');
							}
							else{
								switch(cblayer){
									case '0': seq << (outwards? 'b' : 'a'); break;
									case '1': seq << (outwards? 'd' : 'c'); break;
									case '2': seq << (outwards? 'f' : 'e'); break;
									case '3': seq << (outwards? 'h' : 'g'); break;
									case '4': seq << (outwards? 'j' : 'i'); break;
									case '5': seq << (outwards? 'l' : 'k'); break;
								}
							}
							
						}			
			}
			
		}	
		if(!foundCA && !foundCB)
			seq << 'X';
		return seq.str();
	}
	else{
		//Assemble CA + direction of CB
		stringstream seq;
		for(uint i=0; i<CAburial.size(); i++){
			if(CAburial[i] == -1){
				seq << 'X';
			}
			else{
				if(CBburial[i] == -1)
					if(layers == 2)
						seq << getLayerCharFromBurial(CAburial[i]);
					else
						switch(getLayerCharFromBurial(CAburial[i])){
							case '0': seq << 'a'; break;
							case '1': seq << 'c'; break;
							case '2': seq << 'e'; break;
							case '3': seq << 'g'; break;
							case '4': seq << 'i'; break;
							case '5': seq << 'k'; break;
						}
						else{
							bool outwards = CBburial[i] >= CAburial[i];
							char calayer = getLayerCharFromBurial(CAburial[i]);
							
							if(layers == 2){
								if(calayer == '0')
									seq << (outwards? '^' : '0');
								else
									seq << (outwards? '1' : 'v');
							}
							else{
								switch(calayer){
									case '0': seq << (outwards? 'b' : 'a'); break;
									case '1': seq << (outwards? 'd' : 'c'); break;
									case '2': seq << (outwards? 'f' : 'e'); break;
									case '3': seq << (outwards? 'h' : 'g'); break;
									case '4': seq << (outwards? 'j' : 'i'); break;
									case '5': seq << (outwards? 'l' : 'k'); break;
								}
							}
							
						}			
			}
			
		}	
		if(!foundCA && !foundCB)
			seq << 'X';
		return seq.str();
	}	
	
}


/**
 * Returns a relative sequence of burial layers of CA atoms.
 * (deprecated -- do not use)
 */
string Protein::toRelativeBurialSequence(int layers, string referenceAtom){
	stringstream seq;

	int lastLayer = -1;
	if(layers == 0){
		cerr << "Error: Number of layers should be >0" << endl;
		return "Error: Number of layers should be >0";
	}
	else{
		calculateLayerInfo(layers);

		for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it){
			if(it->type == referenceAtom){
				int layer = getLayerCharFromBurial(it->distanceToCenter);
				if(lastLayer == -1){
					//	seq << layer;
				}
				else {
					char c = (layer == lastLayer)?
							'.'
							:
							( (layer > lastLayer)? '+' : '-' );
					seq << c;
				}
				lastLayer = layer;
			}
		}
	}

	return seq.str();
}

inline char Protein::getLayerCharFromBurial(real burial){
	for(unsigned int i=1; i<layerBoundaries.size(); i++){
		if(burial <= layerBoundaries[i])
			return i-1 + '0';
	}
	if(burial < 0)
		return 'X';
	else
		return layerBoundaries.size() - 2 + '0';
}

int Protein::getLayerNumberFromBurial(real burial){
	for(unsigned int i=1; i<layerBoundaries.size(); i++){
		if(burial <= layerBoundaries[i])
			return i-1;
	}
	return layerBoundaries.size()-2;
}


bool atomComparator (Atom a, Atom b) {
	return a.distanceToCenter < b.distanceToCenter;
}

void Protein::calculateLayerInfo(int layers){	
	layerBoundaries.resize(layers+1);

	//
	// 1. Calculate the boundaries between different layers
	//

	// 1.1. Splitting by residue
	// (each layer has the same number of CA atoms)
	if(splitLayersByRes){

		//filter CA atoms
		vector<Atom> caAtoms;		
		copy_if(
			atoms.begin(), 
			atoms.end(), 
			back_inserter(caAtoms), 
			[](Atom a){
				return a.type == "CA";
			} 
		);

		sort(caAtoms.begin(), caAtoms.end(), atomComparator);
		layerBoundaries[0] = 0;
		real layerSize = (real)caAtoms.size() / (real)layers;
		for(int i=1; i<layers; i++){
			layerBoundaries[i] = caAtoms[(int)(i*layerSize)].distanceToCenter;
		}
		layerBoundaries[layers] = caAtoms[caAtoms.size() - 1].distanceToCenter;
	}

	// 1.2. Splitting by other criteria
	else{	
		vector<Atom> sortedAtoms = atoms;
		sort(sortedAtoms.begin(), sortedAtoms.end(), atomComparator);

		//default: each layer has the same number of atoms
		if(midLayerWidth == -1){
			layerBoundaries[0] = 0;
			real layerSize = (real)sortedAtoms.size() / (real)layers;
			for(int i=1; i<layers; i++){
				layerBoundaries[i] = sortedAtoms[(int)(i*layerSize)].distanceToCenter;
			}
			layerBoundaries[layers] = sortedAtoms[sortedAtoms.size() - 1].distanceToCenter;
		}
		
		//special case: mid layer width specified by user
		else{
			if(layers != 3)
				throw "mid-layer-size only makes sense with 3 layers!";
						
			real layer1Size;
			real layer0Size;
			real layer2Size;		
			
			//varying only middle layer (by proportion to original mid-layer width)
			if(bottomLayerWidth == -1){
				real layerSize = (real)sortedAtoms.size() / (real)layers;
				
				layer1Size = layerSize * midLayerWidth;
				layer0Size = ((real)sortedAtoms.size() - layer1Size) / 2.0;
				layer2Size = layer0Size;			
			}
			
			//varying middle and bottom layer (by proportion with total width)
			else{
				real totalSize = (real)sortedAtoms.size();
				
				layer0Size = bottomLayerWidth * totalSize;
				layer1Size = midLayerWidth * totalSize;
				layer2Size = totalSize - layer0Size - layer1Size;		
			}
			
			layerBoundaries[1] = sortedAtoms[(int)  layer0Size               ].distanceToCenter;
			layerBoundaries[2] = sortedAtoms[(int) (layer0Size + layer1Size) ].distanceToCenter;
			layerBoundaries[3] = sortedAtoms[       sortedAtoms.size() - 1   ].distanceToCenter;
		}
	}

	//
	// 2. Use the boundaries to calculate distances and deltas
	//
	layerDistances.resize(layers);
	layerDeltas.resize(layers);
	for(int i=0; i<layers; i++){
		layerDeltas[i] = (layerBoundaries[i+1] - layerBoundaries[i]) / 2.0;
		layerDistances[i] = layerDeltas[i] + layerBoundaries[i];
	}
}

//These are used in generateHbonds to represent hydrogen bond donors and acceptors
struct donor {
	int O, C;
	donor(){ reset(); }
	void reset(){ O=-1; C=-1; }
	bool isValid(){ return O!=-1 && C!=-1; }
	string str(){
		stringstream ss;
		ss << "O=" << O << ", C=" << C;
		return ss.str();
	}
};
struct acceptor {
	int N, C1, C2;
	acceptor(){ reset(); }
	void reset(){ N=-1; C1=-1; C2 = -1; }
	bool isValid(){ return N!=-1 && C1!=-1 && C2!=-1; }
	string str(){
		stringstream ss;
		ss << "N=" << N << ", C1=" << C1 << ", C2=" << C2;
		return ss.str();
	}
};

/**
 * Generates the "hbonds" file.
 */
void Protein::generateHbonds(ostream& output, bool all){
	int size = atoms[atoms.size()-1].residueIndex + 1;
	vector<int> N(size,-1), CA(size,-1), C(size,-1), O(size,-1);

	for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it){
		if(it->type == "N")
			N[it->residueIndex] = it->index;
		else if(it->type == "CA")
			CA[it->residueIndex] = it->index;
		else if(it->type == "C")
			C[it->residueIndex] = it->index;
		else if(it->type == "O")
			O[it->residueIndex] = it->index;
	}

	list<string> lines;

	//Pairs each N-CA-C group with the O's and C's from every residue
	//except from their immediate neighbours
	for(int i=1; i<size - 1; i++){
		if(
		   N[ i ] == -1 || CA[ i ] == -1 || C[i] == -1 ||
		   N[i+1] == -1 || CA[i+1] == -1			
		) 
			continue;
		if(
				atoms[N[i+1]].residueType == "PRO" ||
				atoms[CA[i+1]].residueType == "PRO" ||
				atoms[C[i]].residueType == "PRO"
		)
			continue;

		stringstream N_CA_Css;
		N_CA_Css << N[i+1] << " " << CA[i+1] << " " << C[i];
		string N_CA_C = N_CA_Css.str();

		for(int j=1; j<size; j++){
			//if(O[j] == -1 || C[j] == -1) continue;

			if(j >= i && j < i+3) continue; //skip neighbours
			stringstream line;
			line << N_CA_C << " " << O[j] << " " << C[j] << endl;
			lines.push_back(line.str());
		}
	}

	//Generate bonds from side chains
	if(all){
		donor don1,don2;
		acceptor acc1,acc2,acc3;

		vector<donor> donors;
		vector<acceptor> acceptors;

		int lastRes = -1;
		for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it){
			if(it->residueIndex != lastRes){
				if(don1.isValid()) donors.push_back(don1);
				if(don2.isValid()) donors.push_back(don2);

				if(acc1.isValid()) acceptors.push_back(acc1);
				if(acc2.isValid()) acceptors.push_back(acc2);
				if(acc3.isValid()) acceptors.push_back(acc3);

				don1.reset(); don2.reset();
				acc1.reset(); acc2.reset(); acc3.reset();
			}

			//			Resíduo   Doadores   Aceptores
			//			 ARG                 NH1,NE,NH2
			//			 ASN        OD1         ND2
			//			 ASP      OD1,OD2     OD1,OD2
			//			 GLN        OE1         NE2
			//			 GLU      OE1,OE2     OE1,OE2
			//			 HIS                  ND1,ND2
			//			 LYS                    NZ
			//			 SER        OG          OG
			//			 THR        OG1         OG1
			//			 TYR        OH          OH

			//			if(it->type == "") don1.O = it->index;
			//			if(it->type == "") don1.C = it->index;
			//
			//			if(it->type == "") acc1.C1 = it->index;
			//			if(it->type == "") acc1.N  = it->index;
			//			if(it->type == "") acc1.C2 = it->index;

			if(it->residueType == "ARG"){
				//NH1, NH2?
				if(it->type == "CD") acc1.C1 = it->index;
				if(it->type == "NE") acc1.N  = it->index;
				if(it->type == "CZ") acc1.C2 = it->index;
			}
			else if(it->residueType == "ASN"){
				//ND2?
				if(it->type == "OD1") don1.O = it->index;
				if(it->type == "CG")  don1.C = it->index;
			}
			else if(it->residueType == "ASP"){
				//OD1, OD2?
				if(it->type == "OD1") don1.O = it->index;
				if(it->type == "CG")  don1.C = it->index;

				if(it->type == "OD2") don2.O = it->index;
				if(it->type == "CG")  don2.C = it->index;
			}
			else if(it->residueType == "GLN"){
				//NE2?
				if(it->type == "OE1") don1.O = it->index;
				if(it->type == "CD")  don1.C = it->index;
			}
			else if(it->residueType == "GLU"){
				//OE1, OE2?
				if(it->type == "OE1") don1.O = it->index;
				if(it->type == "CD")  don1.C = it->index;

				if(it->type == "OE2") don2.O = it->index;
				if(it->type == "CD")  don2.C = it->index;
			}
			else if(it->residueType == "HIS"){
				if(it->type == "CG")  acc1.C1 = it->index;
				if(it->type == "ND1") acc1.N  = it->index;
				if(it->type == "CE1") acc1.C2 = it->index;

				if(it->type == "CD2") acc2.C1 = it->index;
				if(it->type == "NE2") acc2.N  = it->index;
				if(it->type == "CE1") acc2.C2 = it->index;
			}
			else if(it->residueType == "LYS"){
				//NZ?
			}
			else if(it->residueType == "SER"){
				//OG?
				if(it->type == "OG") don1.O = it->index;
				if(it->type == "CB") don1.C = it->index;
			}
			else if(it->residueType == "THR"){
				//OG1?
				if(it->type == "OG1") don1.O = it->index;
				if(it->type == "CB")  don1.C = it->index;
			}
			else if(it->residueType == "TYR"){
				//OH?
				if(it->type == "OH") don1.O = it->index;
				if(it->type == "CZ") don1.C = it->index;
			}
			lastRes = it->residueIndex;
		}
		if(don1.isValid()) donors.push_back(don1);
		if(don2.isValid()) donors.push_back(don2);

		if(acc1.isValid()) acceptors.push_back(acc1);
		if(acc2.isValid()) acceptors.push_back(acc2);
		if(acc3.isValid()) acceptors.push_back(acc3);

		//		cout << "Donors: " << endl;
		//		for(vector<donor>::iterator don = donors.begin(); don != donors.end(); ++don){
		//			cout << don->str() << endl;
		//		}
		//		cout << endl << "Acceptors: " << endl;
		//		for(vector<acceptor>::iterator acc = acceptors.begin(); acc != acceptors.end(); ++acc){
		//			cout << acc->str() << endl;
		//		}

		for(vector<donor>::iterator don = donors.begin(); don != donors.end(); ++don){
			int donRes = atoms[don->O - 1].residueIndex;

			for(int i=1; i<size - 1; i++){
				if(
					N[ i ] == -1 || CA[ i ] == -1 || C[i] == -1 ||
					N[i+1] == -1 || CA[i+1] == -1			
				) 
					continue;
				if(
						atoms[N[i+1]].residueType == "PRO" ||
						atoms[CA[i+1]].residueType == "PRO" ||
						atoms[C[i]].residueType == "PRO"
				)
					continue;

				int accRes = atoms[C[i] - 1].residueIndex;

				if((donRes >= accRes-1) && (donRes <= accRes+2)) continue; //skip neighbours

				stringstream line;
				line << N[i+1] << " " << CA[i+1] << " " << C[i] << " " << don->O << " " << don->C << endl;
				lines.push_back(line.str());
			}
		}

	}
	//
	real nres = (real) atoms[atoms.size()-1].residueIndex;	
	
	output << lines.size() << " 40\n0.000 1.000 2000" << endl;	
// 	output << "3.0 0.5 0.7 " << cuberoot(nres) * 3.5287429016 << " 100.0 100.0 100.0 10.0 5.0" << endl;
	output << "3.0 0.5 0.7 " << (cuberoot(nres) * 4.23856314) - 2.0  << " 100.0 100.0 100.0 10.0 5.0" << endl;
	for(list<string>::iterator it = lines.begin(); it != lines.end(); ++it){
		output << *it;
	}
}

/**
 * Generate the "chiral" file.
 */
void Protein::generateChiral(ostream& output){
	int size = atoms[atoms.size()-1].residueIndex + 1;
	vector<int> N(size,-1), CA(size,-1), C(size,-1), CB(size,-1);

	for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it){
		if(it->type == "N")
			N[it->residueIndex] = it->index;
		else if(it->type == "CA")
			CA[it->residueIndex] = it->index;
		else if(it->type == "C")
			C[it->residueIndex] = it->index;
		else if(it->type == "CB")
			CB[it->residueIndex] = it->index;
	}

	list<string> lines;

	for(int i=1; i<size; i++){
		if(N[i] == -1 || CA[i] == -1 || C[i] == -1 || CB[i] == -1) continue;

		stringstream line;
		line.width(4);
		line << N[i];
		line.width(5);
		line << CA[i];
		line.width(5);
		line << C[i];
		line.width(5);
		line << CB[i] << "   2.500  20.000 " << endl;
		lines.push_back(line.str());
	}

	output << lines.size() << " chiral centers" << endl;
	for(list<string>::iterator it = lines.begin(); it != lines.end(); ++it){
		output << *it;
	}
}

/**
 * Should be called before any nextResidue() calls.
 */
inline void Protein::initResidueIterator(){
	lastResidue = std::numeric_limits<int>::min();
	residueIteratorAtom = 0;
	firstIteratorRound = true;
	unkBuffer = 0;
}

/**
 * Each call to this function returns the following residue as a 3-letter code.
 */
inline string Protein::nextResidue(){
	if(unkBuffer){
		unkBuffer--;
		return "UNK";
	}

	while(atoms[residueIteratorAtom].residueIndex == lastResidue){
		if(++residueIteratorAtom >= atoms.size())
			return "end";
	}

	if(firstIteratorRound || atoms[residueIteratorAtom].residueIndex == lastResidue + 1){
		lastResidue = atoms[residueIteratorAtom].residueIndex;
		firstIteratorRound = false;
		return atoms[residueIteratorAtom].residueType;
	}
	else{
		unkBuffer = atoms[residueIteratorAtom].residueIndex - lastResidue - 2;
		lastResidue = atoms[residueIteratorAtom].residueIndex - 1;
		return "UNK";
	}
}

void Protein::generateHistogram(ostream& output, const std::string& referenceAtom){
	real radgyr = getRadiusOfGyration();
	real slotSize = 0.1;
	uint numberOfSlots = ceil(2.0 / slotSize);
	
	for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it)
		it->distanceToCenter /= radgyr;
	
	vector<uint> hgram(numberOfSlots, 0);
	for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it){
		
		if(referenceAtom == "" or it->type == referenceAtom){		
			real burial = it->distanceToCenter;
			uint slot = burial / slotSize;
			hgram[slot]++;
		}
	}
	
	for(uint i=0; i<numberOfSlots; i++)
		output << i*slotSize << " " << hgram[i] << "\n";	
}

/**
 * Returns the amino acid sequence as a series of 3-letter codes.
 */
string Protein::to3Letters(){
	initResidueIterator();
	string aa = nextResidue();
	stringstream seq;
	bool first = true;
	while(aa != "end"){
		if(first)
			seq << aa;
		else
			seq << " " << aa;
		aa = nextResidue();
		first = false;
	}
	return seq.str();
}

/**
 * Returns a sequence of letters that correspond to the amino acid sequence,
 * according to a given mapping.
 */
string Protein::sequenceMapping(map<string, char>& mapping){
	initResidueIterator();
	stringstream seq;
	string aa = nextResidue();
	while(aa != "end"){
		if(mapping.count(aa))
			seq << mapping[aa];
		else
			seq << 'X';
		aa = nextResidue();
	}
	return seq.str();
}

/**
 * Returns the amino acid sequence as a string of HP (Hidrophobic/Polar) codes.
 */
string Protein::toHP(){
	map<string, char> t;
	t["ALA"] = t["CYS"] = t["PHE"] = t["GLY"] = t["ILE"] = t["LEU"] = t["MET"] = t["VAL"] = t["TRP"] = t["TYR"] = 'H';
	t["ASP"] = t["GLU"] = t["HIS"] = t["LYS"] = t["ASN"] = t["PRO"] = t["GLN"] = t["ARG"] = t["SER"] = t["THR"] = 'P';
	t["UNK"] = 'X';
	return sequenceMapping(t);
}

/**
 * Returns the amino acid sequence as a string of HPN (Hidrophobic/Neutral/Polar) codes.
 */
string Protein::toHPN(){
	map<string, char> t;
	t["CYS"] = t["PHE"] = t["ILE"] = t["LEU"] = t["MET"] = t["VAL"] = t["TRP"] = t["TYR"] = 'H';
	t["ASP"] = t["GLU"] = t["LYS"] = t["ASN"] = t["PRO"] = t["GLN"] = t["ARG"] = 'P';
	t["ALA"] = t["GLY"] = t["HIS"] = t["SER"] = t["THR"] = 'N';
	t["UNK"] = 'X';
	return sequenceMapping(t);
}

/**
 * Returns the amino acid sequence as a string of HPNAC (Hidrophobic/Neutral/Polar/Aromatic/Charged) codes.
 */
string Protein::toHPNAC(){
	map<string, char> t;
	t["CYS"]=t["ILE"]=t["LEU"]=t["MET"]=t["VAL"]='H';
	t["ASN"]=t["PRO"]=t["GLN"]='P';
	t["ALA"]=t["GLY"]=t["HIS"]=t["SER"]=t["THR"]='N';
	t["PHE"]=t["TRP"]=t["TYR"]='A';
	t["ASP"]=t["GLU"]=t["LYS"]=t["ARG"]='C';
	t["UNK"] = 'X';
	return sequenceMapping(t);
}

/**
 * Returns the amino acid sequence as a string of 1-letter codes.
 */
string Protein::toLetter(){
	map<string, char> t;
	t["CYS"] = 'C'; t["ILE"] = 'I'; t["LEU"] = 'L'; t["MET"] = 'M';
	t["VAL"] = 'V'; t["ASN"] = 'N'; t["PRO"] = 'P'; t["GLN"] = 'Q';
	t["ALA"] = 'A'; t["GLY"] = 'G'; t["HIS"] = 'H'; t["SER"] = 'S';
	t["THR"] = 'T'; t["PHE"] = 'F'; t["TRP"] = 'W'; t["TYR"] = 'Y';
	t["ASP"] = 'D'; t["GLU"] = 'E'; t["LYS"] = 'K'; t["ARG"] = 'R';
	t["UNK"] = 'X';
	return sequenceMapping(t);
}

/**
 * Returns the amino acid sequence as a string of HPNhp codes.
 */
string Protein::toHPNhp(){
	map<string, char> t;
	t["PHE"] = t["ILE"] = t["LEU"] = t["VAL"] = t["TRP"] = 'H';
	t["ASP"] = t["GLU"] = t["LYS"] = 'P';
	t["ALA"] = t["GLY"] = t["HIS"] = t["SER"] = t["THR"] = 'N';
	t["CYS"] = t["MET"] = t["TYR"] = 'h';
	t["ASN"] = t["PRO"] = t["GLN"] = t["ARG"] = 'p';
	t["UNK"] = 'X';
	return sequenceMapping(t);
}

/**
 * Returns the amino acid sequence as a string of HPNhpn codes.
 */
string Protein::toHPNhpn(){
	map<string, char> t;
	t["PHE"] = t["ILE"] = t["LEU"] = t["VAL"] = t["TRP"] = 'H';
	t["ASP"] = t["GLU"] = t["LYS"] = 'P';
	t["ALA"] = t["HIS"] = 'N';
	t["GLY"] = t["SER"] = t["THR"] = 'n';
	t["CYS"] = t["MET"] = t["TYR"] = 'h';
	t["ASN"] = t["PRO"] = t["GLN"] = t["ARG"] = 'p';
	t["UNK"] = 'X';
	return sequenceMapping(t);
}

/**
 * Generates the ".dat" file with the simulation topology and parameters.
 */
void Protein::generateDat(int layers, ASM atomSizes, real scale, ostream& output){
	DatFile dat(atoms, atomSizes, scale);
	dat.write(layers, output);	
}

/**
 * Generates a file reporting the a list of pairs of atoms that are closest to each other in space.
 */ 
void Protein::generateClashReport(ostream& output){
	ClashFinder cf(atoms, 0, 0);
	cf.findClashes(output);
}

/**
 * Generates a PDB file.
 */
void Protein::generatePDB(int layers, ostream& output){

	for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it){
		char line[68];
		sprintf(line,
			"ATOM  %5d  %3s %3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
			it->index, it->type.c_str(), it->residueType.c_str(), it->residueIndex,
			it->x, it->y, it->z, it->distanceToCenter, it->distanceToCenter
		);
		output << line;
	}	
}

bool Protein::loadLayerInfoFromBData(int nlayers, string referenceFile){
	ifstream input(referenceFile.c_str());
	if(input.fail()){
		cout << "Error loading " << referenceFile << endl;
		return false;
	}
	
	string line;
	while( getline (input, line) ) {
		istringstream linein(line);
		
		string token;
		int n;

		linein >> token;
		if(token != "#lb")
			continue;
				
		linein >> n;
		if(n == nlayers){
			layerBoundaries.resize(n+1);
			for(int i=0; i<=n; i++)
				linein >> layerBoundaries[i];
			return true;
		}
	}
	
	cout << "Could not parse file " << referenceFile << endl;
	return false;
}



void Protein::generateBData(ostream& output, string referenceFile){
	
	real mult = 2.7 * cuberoot(atoms[atoms.size()-1].residueIndex);
	
	vector< vector<char> > layer(25, vector<char>());
	
	stringstream ss;
	
	//
	// Calculate layer boundaries
	//
 	for(int nl=2; nl<=24; nl++){
		
		if(referenceFile == "")
			calculateLayerInfo(nl);
		else{
			if( !loadLayerInfoFromBData(nl, referenceFile) )
				return;
		}
		
						
		layer[nl].reserve(atoms.size());
		
		for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it)
			layer[nl].push_back( getLayerCharFromBurial(it->distanceToCenter) );
		
		
		ss << "#lb " << nl << "\t";
		for(int i=0; i<layerBoundaries.size(); i++)
			ss << layerBoundaries[i] << " ";
		ss << "\n";
	}
	
	//
	// Print out bdata info
	//	
	//(atmN, *, type, restype, resN, burialAngstrom, burialRadgyr, l2, l3, l4, l5, l6)
	int i = 0;
	output << "#   1  2  3   4     5                 6                   7  8  9 10 11 12 13 14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29\n"
	       << "#                   burial_angs       burial_rg          l2 l3 l4 l5 l6 l7 l8 l9 l10 l11 l12 l13 l14 l15 l16 l17 l18 l19 l20 l21 l22 l23 l24\n";
	for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it){
		
		char line[500];
		sprintf(line,
			"%5d %3s %3s %-4d  %16.13f %16.13f     %c  %c  %c  %c  %c  %c  %c  %c  %2d  %2d  %2d  %2d  %2d  %2d  %2d  %2d  %2d  %2d  %2d  %2d  %2d  %2d  %2d\n",
			it->index, it->type.c_str(), it->residueType.c_str(), it->residueIndex,
			it->distanceToCenter, it->distanceToCenter / mult,
			layer[2][i], layer[3][i], layer[4][i], layer[5][i], layer[6][i],
			layer[7][i], layer[8][i], layer[9][i], layer[10][i]-48, layer[11][i]-48,
			layer[12][i]-48, layer[13][i]-48, layer[14][i]-48, layer[15][i]-48, layer[16][i]-48,
			layer[17][i]-48, layer[18][i]-48, layer[19][i]-48, layer[20][i]-48, layer[21][i]-48,
			layer[22][i]-48, layer[23][i]-48, layer[24][i]-48
		);
		
		output << line;
		i++;		
	}
	
	output << ss.str();	
}

inline string Protein::addRamaLine(ostream& output, vector< map<string,uint> >& atomMap, uint x, string a1, uint y, string a2, uint z, string a3, uint k, string a4, const string STR ){
	if(!atomMap[x].count(a1)) throw a1 + " not found in residue " + int2string(x+1);
	if(!atomMap[y].count(a2)) throw a2 + " not found in residue " + int2string(y+1);
	if(!atomMap[z].count(a3)) throw a3 + " not found in residue " + int2string(z+1);
	if(!atomMap[k].count(a4)) throw a4 + " not found in residue " + int2string(k+1);
	
	char* tmp = new char[16];
	sprintf(tmp, "%3d %3d %3d %3d", 
		atoms[ atomMap[x][a1] ].index , 
		atoms[ atomMap[y][a2] ].index , 
		atoms[ atomMap[z][a3] ].index , 
		atoms[ atomMap[k][a4] ].index 
	); 
	string s = tmp + STR;
	delete tmp;
	return s;
}

/**
 * Generates the "rama" file with phi/psi angle parameters.
 */
void Protein::generateRama(ostream& output){	
	//this copied from DatFile.cpp:
	vector< map<string,uint> > atomMap;
	int lastResidueIndex = -1;
	map<string,uint> tmp;
	for(uint i=0; i<atoms.size(); i++){
		if(atoms[i].residueIndex != lastResidueIndex){
			if(lastResidueIndex != -1)
				atomMap.push_back(tmp);
			tmp.clear();
			lastResidueIndex = atoms[i].residueIndex;
		}
		tmp[atoms[i].type] = i;
	}
	atomMap.push_back(tmp);
	//
	
// 	const string PHISTR  = "   0.00 0.00000000000   0.27 0.00000000000   0.42 0.00000000000   1 # phi\n";
// 	const string PSISTR  = "   0.45 3.14159265359   1.58 3.14159265359   0.55 3.14159265359   2 # psi\n";
// 	const string PHILSTR = "   2.00 0.00000000000   2.00 0.00000000000   0.40 0.00000000000   3 # phi'\n";
// 	const string PSILSTR = "   0.20 0.00000000000   0.20 0.00000000000   0.40 0.00000000000   4 # psi'\n";
	
	const string PHISTR  = "   1.50  1.04719755120  1.00 -1.04719755120  0.50 0.00000000000  0.50 3.14159265359 # phi\n";
 	const string PSISTR  = "   1.50 -2.09439510239  1.50  1.04719755120  1.00 0.00000000000  0.50 3.14159265359 # psi\n";

	list<string> lines;
	
	for(uint i=0; i<atomMap.size()-1; i++){
			if(atoms[atomMap[i]["CA"]].residueType == "GLY"){
				cout << "GLY skipping " << atoms[atomMap[i]["CA"]].index << "\n";
				continue;
			}
			
			try{
				if(atoms[atomMap[i+1]["CA"]].residueType == "PRO")
					cout << "PRO skipping " << atoms[atomMap[i+1]["CA"]].index << "\n";
				else
					lines.push_back( addRamaLine(output, atomMap, i,"C",  i+1,"N",  i+1,"CA", i+1,"C", PHISTR) ); //phi
			}
			catch(string e){ cerr << e << " (phi)\n"; }
			
			try{
				lines.push_back( addRamaLine(output, atomMap, i,"N",    i,"CA",   i,"C",  i+1,"N", PSISTR) ); //psi
			}
			catch(string e){ cerr << e << " (psi)\n"; }
			
// 			try{
// 				lines.push_back( addRamaLine(output, atomMap, i,"C",  i+1,"N",  i+1,"CA", i+1,"CB", PHILSTR) ); //phi'
// 			}
// 			catch(string e){ cerr << e << " (phi')\n"; }
// 			
// 			try{
// 				lines.push_back( addRamaLine(output, atomMap, i,"CB",   i,"CA",   i,"C",  i+1,"N",  PSILSTR) ); //psi'
// 			}
// 			catch(string e){ cerr << e << " (psi')\n"; }
			
	}
	
	output << lines.size() << " phi/psi angles\n";
	for(list<string>::iterator it = lines.begin(); it != lines.end(); it++)
		output << *it;
}


