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

#include "ClashFinder.h"

struct AtomPair{
	real dist;
	int i;
	int j;
	
	AtomPair(int _i, int _j, real _dist) : i(_i), j(_j), dist(_dist) {}
};

bool atomPairSort(const AtomPair & a, const AtomPair & b) { return a.dist < b.dist; }

void ClashFinder::findClashes(ostream& output){
// 	const uint OUTPUT_MAX_SIZE = 10;	
	
	buildConnectedMap();
	
	list<AtomPair> pairs;

	for(uint i=0; i<atoms.size(); i++)
		for(uint j=i+1; j<atoms.size(); j++)
			if( unboundedBy2(i,j) && notNeighboringCA(i,j) && notNeighboringOCA(i,j) && dontFormDihedral(i,j) ) {
				Atom& atm1 = atoms[i];
				Atom& atm2 = atoms[j];
				
				if(atm1.residueIndex == atm2.residueIndex)
					continue;
		
				real dist = atoms[i].distanceTo(atoms[j]);
				AtomPair ap(i,j,dist);
				pairs.push_back( ap );
			}
			
	pairs.sort(atomPairSort);
	
// 	uint n = 0;
	for(list<AtomPair>::iterator it = pairs.begin(); it != pairs.end(); it++){
// 		if(n++ == OUTPUT_MAX_SIZE)
// 			break;
		output << " " << it->dist << "\t" << it->i+1 << " " << it->j+1 << "\n";
	}
}

