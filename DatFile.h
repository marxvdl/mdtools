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

#ifndef DATFILE_H
#define DATFILE_H

#include "Protein.h"
#include "Atom.h"
#include <vector>


class DatFile
{
protected:
	vector<Atom> atoms;
	vector< map<string,uint> > atomMap;
	vector< vector<uint> > connected2;
	ASM atomSizes;
	real scale;
	
	void buildConnectedMap();
	
	void populateBondedPairs(list< string >& lines);	
	void populateAngles(list< string >& lines);
	void populateDihedrals(list< string >& lines);
	void populateRepulsivePais(list< string >& lines);
	void populateStructure(list< string >& lines);
	
	void addBond(list< string >& lines, uint i, string a1, uint j, string a2);
	void addAngle(list< string >& lines, uint x, string a1, uint y, string a2, uint z, string a3 );
	void addDihedral(list< string >& lines, uint x, string a1, uint y, string a2, uint z, string a3, uint k, string a4, real param );
	
	string writeBondString(uint a1, uint a2);
	string writeAngleString(uint a1, uint a2, uint a3);
	string writeDihedralString(uint a1, uint a2, uint a3, uint a4, real param);
	string addRepulsivePair(uint a1, uint a2);
	bool unboundedBy2(uint a1, uint a2);
	bool notNeighboringCA(uint a1, uint a2);
	bool notNeighboringOCA(uint a1, uint a2);
	bool dontFormDihedral(uint a1, uint a2);
	
	bool containedInDihedral(string a, string b, uint n, const string dihedrals[][4]);
	inline bool compare(uint i, uint j, string a, string b);
	
	void printSection(ostream& output, list<string>& section){
		for(list<string>::iterator it = section.begin(); it != section.end(); ++it)
			output << *it;
	}

public:
	DatFile(vector< Atom >& proteinAtoms, ASM asM, real scale);
	void write(int layers, ostream& output);
	
};

#endif // DATFILE_H
