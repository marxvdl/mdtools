#include "DatFile.h"

DatFile::DatFile(vector< Atom >& proteinAtoms, ASM asM, real sca){
	atomSizes = asM;
	scale = sca;
	
	//
	// 1. Remove hidrogens from the list of atoms
	//
	uint myIndex = 1;
	for(vector<Atom>::iterator it = proteinAtoms.begin(); it != proteinAtoms.end(); ++it)
// 		if(
// 			it->type == "NH1" || it->type == "NH2" || it->type == "OH" || it->type == "CH2" ||
// 			it->type.find_first_of('H') == string::npos
// 		) 
		{ 
			//does not containt 'H'
			Atom tmp = *it;
			tmp.index = myIndex++;
			atoms.push_back(tmp);
		}
// 		else{
// 			cout << "skipping atom " << it->index << " " << it->type << "\n";
// 		}
	
	//
	// 2. Make an atom map for internal reference
	//
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
}


void DatFile::write(int layers, ostream& output){		
	list<string> bondedPairs;
	populateBondedPairs(bondedPairs);	
	printSection(output, bondedPairs);
	
	list<string> angles;
	populateAngles(angles);
	printSection(output, angles);
	
	list<string> dihedrals;
	populateDihedrals(dihedrals);
	printSection(output, dihedrals);
	
	output << "           0 Contacts.\n";
	
	list<string> repulsivePairs;
	populateRepulsivePais(repulsivePairs);
	printSection(output, repulsivePairs);
		
	list<string> structure;
	populateStructure(structure);
	printSection(output, structure);
}


inline void DatFile::addBond(list< string >& lines, uint i, string a1, uint j, string a2){
	if(!atomMap[i].count(a1)) return;
	if(!atomMap[j].count(a2)) return;
	
	lines.push_back( writeBondString(atomMap[i][a1], atomMap[j][a2]) );
}

void DatFile::populateBondedPairs(list< string >& lines){
	for(uint i=0; i<atomMap.size(); i++){
		//Backbone
		addBond(lines,  i,"N", i,"CA");
		addBond(lines,  i,"CA", i,"C");
		addBond(lines,  i,"C", i,"O");
		if(i != atomMap.size()-1)
			addBond(lines,  i,"C", i+1,"N");
		
		//Sidechain
			string resType = atoms[atomMap[i]["CA"]].residueType;
			if(resType == "GLY")
				continue;
			
			addBond(lines,  i,"CB", i,"CA");
			
			if(resType == "ARG"){
				addBond(lines,  i,"CB", i,"CG");
				addBond(lines,  i,"CG", i,"CD");
				addBond(lines,  i,"CD", i,"NE");
				addBond(lines,  i,"NE", i,"CZ");
				addBond(lines,  i,"CZ", i,"NH1");
				addBond(lines,  i,"CZ", i,"NH2");
			}
			else if(resType == "ASN"){
				addBond(lines,  i,"CB", i,"CG");
				addBond(lines,  i,"CG", i,"OD1");
				addBond(lines,  i,"CG", i,"ND2");
			}
			else if(resType == "ASP"){
				addBond(lines,  i,"CB", i,"CG");
				addBond(lines,  i,"CG", i,"OD1");
				addBond(lines,  i,"CG", i,"OD2");
			}
			else if(resType == "CYS"){
				addBond(lines,  i,"CB", i,"SG");
			}
			else if(resType == "GLN"){
				addBond(lines,  i,"CB", i,"CG");
				addBond(lines,  i,"CG", i,"CD");
				addBond(lines,  i,"CD", i,"OE1");
				addBond(lines,  i,"CD", i,"NE2");
			}
			else if(resType == "GLU"){
				addBond(lines,  i,"CB", i,"CG");
				addBond(lines,  i,"CG", i,"CD");
				addBond(lines,  i,"CD", i,"OE1");
				addBond(lines,  i,"CD", i,"OE2");
			}
			else if(resType == "HIS"){
				addBond(lines,  i,"CB", i,"CG");
				addBond(lines,  i,"CG", i,"ND1");
				addBond(lines,  i,"CG", i,"CD2");
				addBond(lines,  i,"ND1", i,"CE1");
				addBond(lines,  i,"CD2", i,"NE2");
				addBond(lines,  i,"NE2", i,"CE1");
			}
			else if(resType == "ILE"){
				addBond(lines,  i,"CB", i,"CG1");
				addBond(lines,  i,"CG1", i,"CD1");
				addBond(lines,  i,"CB", i,"CG2");
			}
			else if(resType == "LEU"){
				addBond(lines,  i,"CB", i,"CG");
				addBond(lines,  i,"CG", i,"CD1");
				addBond(lines,  i,"CG", i,"CD2");
			}
			else if(resType == "LYS"){
				addBond(lines,  i,"CB", i,"CG");
				addBond(lines,  i,"CG", i,"CD");
				addBond(lines,  i,"CD", i,"CE");
				addBond(lines,  i,"CE", i,"NZ");
			}
			else if(resType == "MET"){
				addBond(lines,  i,"CB", i,"CG");
				addBond(lines,  i,"CG", i,"SD");
				addBond(lines,  i,"SD", i,"CE");
			}
			else if(resType == "PHE"){
				addBond(lines,  i,"CB", i,"CG");
				addBond(lines,  i,"CG", i,"CD1");
				addBond(lines,  i,"CG", i,"CD2");
				addBond(lines,  i,"CD1", i,"CE1");
				addBond(lines,  i,"CD2", i,"CE2");
				addBond(lines,  i,"CE1", i,"CZ");
				addBond(lines,  i,"CE2", i,"CZ");
			}
			else if(resType == "PRO"){
				addBond(lines,  i,"CB", i,"CG");
				addBond(lines,  i,"CG", i,"CD");
				addBond(lines,  i,"CD", i,"N");
			}
			else if(resType == "SER"){
				addBond(lines,  i,"CB", i,"OG");
			}
			else if(resType == "THR"){
				addBond(lines,  i,"CB", i,"OG1");
				addBond(lines,  i,"CB", i,"CG2");
			}
			else if(resType == "TRP"){
				addBond(lines,  i,"CB", i,"CG");
				addBond(lines,  i,"CG", i,"CD1");
				addBond(lines,  i,"CG", i,"CD2");
				addBond(lines,  i,"CD1", i,"NE1");
				addBond(lines,  i,"CD2", i,"CE2");
				addBond(lines,  i,"NE1", i,"CE2");
				addBond(lines,  i,"CE2", i,"CZ2");
				addBond(lines,  i,"CZ2", i,"CH2");
				addBond(lines,  i,"CZ3", i,"CH2");
				addBond(lines,  i,"CE3", i,"CZ2");
				addBond(lines,  i,"CD2", i,"CE3");
			}
			else if(resType == "TYR"){
				addBond(lines,  i,"CB", i,"CG");
				addBond(lines,  i,"CG", i,"CD1");
				addBond(lines,  i,"CG", i,"CD2");
				addBond(lines,  i,"CD1", i,"CE1");
				addBond(lines,  i,"CD2", i,"CE2");
				addBond(lines,  i,"CE1", i,"CZ");
				addBond(lines,  i,"CE2", i,"CZ");
				addBond(lines,  i,"CZ", i,"OH");
			}
			else if(resType == "VAL"){
				addBond(lines,  i,"CB", i,"CG1");
				addBond(lines,  i,"CB", i,"CG2");
			}
	}
	char header[36];
	sprintf(header, "%12d Hookean bonded pairs.\n", (int)lines.size());
	lines.push_front(header);
}

inline string DatFile::writeBondString(uint a1, uint a2){
	char line[31];
	sprintf(line, "%5d %5d %8.3f %8.3f\n", 
		atoms[a1].index,
	 atoms[a2].index,
	 atoms[a1].distanceTo(atoms[a2]),
		100.0		
	);
	return line;       
}

inline void DatFile::addAngle(list< string >& lines, uint x, string a1, uint y, string a2, uint z, string a3 ){
	if(!atomMap[x].count(a1)) return;
	if(!atomMap[y].count(a2)) return;
	if(!atomMap[z].count(a3)) return;
	
	lines.push_back( writeAngleString(atomMap[x][a1], atomMap[y][a2], atomMap[z][a3]) );
}

void DatFile::populateAngles(list< string >& lines){
	for(uint i=0; i<atomMap.size(); i++){	
		//Backbone
		addAngle(lines,  i,"N", i,"CA", i,"C");
		addAngle(lines,  i,"CA", i,"C", i,"O");
		if(i != atomMap.size()-1){
			addAngle(lines,  i,"CA", i,"C", i+1,"N");
			addAngle(lines,  i,"O", i,"C", i+1,"N");
			addAngle(lines,  i,"C", i+1,"N", i+1,"CA");
		}
		
		//Sidechain
		string resType = atoms[atomMap[i]["CA"]].residueType;
		if(resType == "GLY")
			continue;
		
		addAngle(lines,  i,"N", i,"CA", i,"CB");
		addAngle(lines,  i,"C", i,"CA", i,"CB");
		
		if(resType == "ARG"){                                                                                                               
			addAngle(lines,  i,"CA", i,"CB", i,"CG");
			addAngle(lines,  i,"CB", i,"CG", i,"CD");
			addAngle(lines,  i,"CG", i,"CD", i,"NE");                                                                                                          
			addAngle(lines,  i,"CD", i,"NE", i,"CZ");
			addAngle(lines,  i,"NE", i,"CZ", i,"NH1");                                                                                                         
			addAngle(lines,  i,"NE", i,"CZ", i,"NH2");
			addAngle(lines,  i,"NH1", i,"CZ", i,"NH2");                                                                                                        
		}
		else if(resType == "ASN"){                                                                                                                          
			addAngle(lines,  i,"CA", i,"CB", i,"CG");
			addAngle(lines,  i,"CB", i,"CG", i,"OD1");                                                                                                         
			addAngle(lines,  i,"CB", i,"CG", i,"ND2");
			addAngle(lines,  i,"OD1", i,"CG", i,"ND2");
		}
		else if(resType == "ASP"){
			addAngle(lines,  i,"CA", i,"CB", i,"CG");
			addAngle(lines,  i,"CB", i,"CG", i,"OD1");
			addAngle(lines,  i,"CB", i,"CG", i,"OD2");
			addAngle(lines,  i,"OD1", i,"CG", i,"OD2");
		}
		else if(resType == "CYS"){
			addAngle(lines,  i,"CA", i,"CB", i,"SG");
		}
		else if(resType == "GLN"){
			addAngle(lines,  i,"CA", i,"CB", i,"CG");
			addAngle(lines,  i,"CB", i,"CG", i,"CD");
			addAngle(lines,  i,"CG", i,"CD", i,"OE1");
			addAngle(lines,  i,"CG", i,"CD", i,"NE2");
			addAngle(lines,  i,"OE1", i,"CD", i,"NE2");
		}
		else if(resType == "GLU"){
			addAngle(lines,  i,"CA", i,"CB", i,"CG");
			addAngle(lines,  i,"CB", i,"CG", i,"CD");
			addAngle(lines,  i,"CG", i,"CD", i,"OE1");
			addAngle(lines,  i,"CG", i,"CD", i,"OE2");
			addAngle(lines,  i,"OE1", i,"CD", i,"OE2");
		}
		else if(resType == "HIS"){
			addAngle(lines,  i,"CA", i,"CB", i,"CG");
			addAngle(lines,  i,"CB", i,"CG", i,"ND1");
			addAngle(lines,  i,"CB", i,"CG", i,"CD2");
			addAngle(lines,  i,"ND1", i,"CG", i,"CD2");
			addAngle(lines,  i,"CG", i,"ND1", i,"CE1");
			addAngle(lines,  i,"CG", i,"CD2", i,"NE2");
			addAngle(lines,  i,"CD2", i,"NE2", i,"CE1");
			addAngle(lines,  i,"ND1", i,"CE1", i,"NE2");
		}
		else if(resType == "ILE"){
			addAngle(lines,  i,"CA", i,"CB", i,"CG1");
			addAngle(lines,  i,"CB", i,"CG1", i,"CD1");
			addAngle(lines,  i,"CA", i,"CB", i,"CG2");
			addAngle(lines,  i,"CG1", i,"CB", i,"CG2");
		}
		else if(resType == "LEU"){
			addAngle(lines,  i,"CA", i,"CB", i,"CG");
			addAngle(lines,  i,"CB", i,"CG", i,"CD1");
			addAngle(lines,  i,"CB", i,"CG", i,"CD2");
			addAngle(lines,  i,"CD1", i,"CG", i,"CD2");
		}
		else if(resType == "LYS"){
			addAngle(lines,  i,"CA", i,"CB", i,"CG");
			addAngle(lines,  i,"CB", i,"CG", i,"CD");
			addAngle(lines,  i,"CG", i,"CD", i,"CE");
			addAngle(lines,  i,"CD", i,"CE", i,"NZ");
		}
		else if(resType == "MET"){
			addAngle(lines,  i,"CA", i,"CB", i,"CG");
			addAngle(lines,  i,"CB", i,"CG", i,"SD");
			addAngle(lines,  i,"CG", i,"SD", i,"CE");
		}
		else if(resType == "PHE"){
			addAngle(lines,  i,"CA", i,"CB", i,"CG");
			addAngle(lines,  i,"CB", i,"CG", i,"CD1");
			addAngle(lines,  i,"CB", i,"CG", i,"CD2");
			addAngle(lines,  i,"CD1", i,"CG", i,"CD2");
			addAngle(lines,  i,"CG", i,"CD1", i,"CE1");
			addAngle(lines,  i,"CG", i,"CD2", i,"CE2");
			addAngle(lines,  i,"CD1", i,"CE1", i,"CZ");
			addAngle(lines,  i,"CD2", i,"CE2", i,"CZ");
			addAngle(lines,  i,"CE1", i,"CZ", i,"CE2");
		}
		else if(resType == "PRO"){
			addAngle(lines,  i,"CA", i,"CB", i,"CG");
			addAngle(lines,  i,"CB", i,"CG", i,"CD");
			addAngle(lines,  i,"CG", i,"CD", i,"N");
			addAngle(lines,  i,"CD", i,"N", i,"CA");
			if(i != 0)
				addAngle(lines,  i,"CD", i,"N", i-1,"C");
		}
		else if(resType == "SER"){
			addAngle(lines,  i,"CA", i,"CB", i,"OG");
		}
		else if(resType == "THR"){
			addAngle(lines,  i,"CA", i,"CB", i,"OG1");
			addAngle(lines,  i,"CA", i,"CB", i,"CG2");
			addAngle(lines,  i,"OG1", i,"CB", i,"CG2");
		}
		else if(resType == "TRP"){
			addAngle(lines,  i,"CA", i,"CB", i,"CG");
			addAngle(lines,  i,"CB", i,"CG", i,"CD1");
			addAngle(lines,  i,"CB", i,"CG", i,"CD2");
			addAngle(lines,  i,"CD1", i,"CG", i,"CD2");
			addAngle(lines,  i,"CG", i,"CD2", i,"CE2");
			addAngle(lines,  i,"CD2", i,"CE2", i,"NE1");
			addAngle(lines,  i,"CE2", i,"NE1", i,"CD1");
			addAngle(lines,  i,"NE1", i,"CD1", i,"CG");
			addAngle(lines,  i,"CD2", i,"CE3", i,"CZ3");
			addAngle(lines,  i,"CE3", i,"CZ3", i,"CH2");
			addAngle(lines,  i,"CZ3", i,"CH2", i,"CZ2");
			addAngle(lines,  i,"CH2", i,"CZ2", i,"CE2");
			addAngle(lines,  i,"CZ2", i,"CE2", i,"CD2");
			addAngle(lines,  i,"CE2", i,"CD2", i,"CE3");
			addAngle(lines,  i,"CG", i,"CD2", i,"CE3");
			addAngle(lines,  i,"NE1", i,"CE2", i,"CZ2");
		}
		else if(resType == "TYR"){
			addAngle(lines,  i,"CA", i,"CB", i,"CG");
			addAngle(lines,  i,"CB", i,"CG", i,"CD1");
			addAngle(lines,  i,"CB", i,"CG", i,"CD2");
			addAngle(lines,  i,"CD1", i,"CG", i,"CD2");
			addAngle(lines,  i,"CG", i,"CD1", i,"CE1");
			addAngle(lines,  i,"CG", i,"CD2", i,"CE2");
			addAngle(lines,  i,"CD1", i,"CE1", i,"CZ");
			addAngle(lines,  i,"CD2", i,"CE2", i,"CZ");
			addAngle(lines,  i,"CE1", i,"CZ", i,"CE2");
			addAngle(lines,  i,"CE1", i,"CZ", i,"OH");
			addAngle(lines,  i,"CE2", i,"CZ", i,"OH");
		}
		else if(resType == "VAL"){
			addAngle(lines,  i,"CA", i,"CB", i,"CG1");
			addAngle(lines,  i,"CA", i,"CB", i,"CG2");
			addAngle(lines,  i,"CG1", i,"CB", i,"CG2");
		}		
	}
	char header[32];
	sprintf(header, "%12d Bond angle pairs.\n", (int)lines.size());
	lines.push_front(header);
	
	
}

inline string DatFile::writeAngleString(uint a1, uint a2, uint a3){
	char line[37];
	sprintf(line, "%5d %5d %5d %8.3f %8.3f\n", 
		atoms[a1].index,
	 atoms[a2].index,
	 atoms[a3].index,
	 atoms[a1].angleTo(atoms[a2], atoms[a3]),
		20.0		
	);
	return line; 
}

inline void DatFile::addDihedral(list< string >& lines, uint x, string a1, uint y, string a2, uint z, string a3, uint k, string a4, real param ){
	if(!atomMap[x].count(a1)) return;
	if(!atomMap[y].count(a2)) return;
	if(!atomMap[z].count(a3)) return;
	if(!atomMap[k].count(a4)) return;
	
	lines.push_back( writeDihedralString(atomMap[x][a1], atomMap[y][a2], atomMap[z][a3], atomMap[k][a4], param) );
}

void DatFile::populateDihedrals(list< string >& lines){	
// 	//CA dihedrals
// 	for(uint i=0; i<atomMap.size()-3; i++)
// 		addDihedral(lines,  i,"CA",  i+1,"CA",  i+2,"CA",  i+3,"CA", 2.0);
	
	//Peptide bond
	for(uint i=0; i<atomMap.size()-1; i++){
		addDihedral(lines,  i,"O",  i,"C",  i+1,"N",  i+1,"CA", 5.0);
		addDihedral(lines,  i,"CA",  i,"C",  i+1,"N",  i+1,"CA", 5.0);

		//Side chain cyclic groups
		string resType = atoms[atomMap[i]["CA"]].residueType;
		if(resType == "TYR"){
			addDihedral(lines, i,"CB", i,"CG", i,"CD1", i,"CE1", 10.0);
			addDihedral(lines, i,"CB", i,"CG", i,"CD2", i,"CE2", 10.0);
			addDihedral(lines, i,"CG", i,"CD1", i,"CE1", i,"CZ", 10.0);
			addDihedral(lines, i,"CD1", i,"CE1", i,"CZ", i,"CE2", 10.0);
			addDihedral(lines, i,"CE1", i,"CZ", i,"CE2", i,"CD2", 10.0);
			addDihedral(lines, i,"CZ", i,"CE2", i,"CD2", i,"CG", 10.0);
			addDihedral(lines, i,"CE2", i,"CD2", i,"CG", i,"CD1", 10.0);
			addDihedral(lines, i,"CD2", i,"CG", i,"CD1", i,"CE1", 10.0);
			addDihedral(lines, i,"CD1", i,"CE1", i,"CZ", i,"OH", 10.0);
			addDihedral(lines, i,"CD2", i,"CE2", i,"CZ", i,"OH", 10.0);			
		}
		else if(resType == "PHE"){
			addDihedral(lines, i,"CB", i,"CG", i,"CD1", i,"CE1", 10.0);
			addDihedral(lines, i,"CB", i,"CG", i,"CD2", i,"CE2", 10.0);
			addDihedral(lines, i,"CG", i,"CD1", i,"CE1", i,"CZ", 10.0);
			addDihedral(lines, i,"CD1", i,"CE1", i,"CZ", i,"CE2", 10.0);
			addDihedral(lines, i,"CE1", i,"CZ", i,"CE2", i,"CD2", 10.0);
			addDihedral(lines, i,"CZ", i,"CE2", i,"CD2", i,"CG", 10.0);
			addDihedral(lines, i,"CE2", i,"CD2", i,"CG", i,"CD1", 10.0);
			addDihedral(lines, i,"CD2", i,"CG", i,"CD1", i,"CE1", 10.0);
		}
		else if(resType == "HIS"){
			addDihedral(lines, i,"CB", i,"CG", i,"ND1", i,"CE1", 10.0);
			addDihedral(lines, i,"CB", i,"CG", i,"CD2", i,"NE2", 10.0);
			addDihedral(lines, i,"CG", i,"ND1", i,"CE1", i,"NE2", 10.0);
			addDihedral(lines, i,"ND1", i,"CE1", i,"NE2", i,"CD2", 10.0);
			addDihedral(lines, i,"CE1", i,"NE2", i,"CD2", i,"CG", 10.0);
			addDihedral(lines, i,"NE2", i,"CD2", i,"CG", i,"ND1", 10.0);
			addDihedral(lines, i,"CD2", i,"CG", i,"ND1", i,"CE1", 10.0);
		}
		else if(resType == "PRO"){
			addDihedral(lines, i,"CA", i,"CB", i,"CG", i,"CD", 10.0);
			addDihedral(lines, i,"CB", i,"CG", i,"CD", i,"N", 10.0);
			addDihedral(lines, i,"CG", i,"CD", i,"N", i,"CA", 10.0);
			addDihedral(lines, i,"CD", i,"N", i,"CA", i,"CB", 10.0);
		}
	}
	
	char header[30];
	sprintf(header, "%12d Dihedral pairs.\n", (int)lines.size());
	lines.push_front(header);
}

inline string DatFile::writeDihedralString(uint a1, uint a2, uint a3, uint a4, real param){
	char line[65];
	sprintf(line, "%5d %5d %5d %5d %9.5f %9.5f %9.5f %9.5f\n", 
		atoms[a1].index,
		atoms[a2].index,
		atoms[a3].index,
		atoms[a4].index,
		atoms[a1].dihedralTo(atoms[a2], atoms[a3], atoms[a4]),
		param, 0.0, 0.0
	);
	return line; 	
}

inline bool isBackBoneAtom(const string & atom){
	if(atom == "N" || atom == "CA" || atom == "C" || atom == "O")
		return true;
	else
		return false;
}


void DatFile::populateRepulsivePais(list< string >& lines){
	buildConnectedMap();
	
	if(atomSizes == 400){
		for(uint i=0; i<atoms.size(); i++)
			for(uint j=i+1; j<atoms.size(); j++)
				if( isBackBoneAtom(atoms[i].type) && isBackBoneAtom(atoms[j].type) && unboundedBy2(i,j) && notNeighboringCA(i,j) && notNeighboringOCA(i,j) && dontFormDihedral(i,j) )
					lines.push_back( addRepulsivePair(i,j) );
	}
	else{
		for(uint i=0; i<atoms.size(); i++)
			for(uint j=i+1; j<atoms.size(); j++)
				if( unboundedBy2(i,j) && notNeighboringCA(i,j) && notNeighboringOCA(i,j) && dontFormDihedral(i,j) )
					lines.push_back( addRepulsivePair(i,j) );
	}
	
			
	char header[30];
	sprintf(header, "%12d Repulsive pairs\n", (int)lines.size());
	lines.push_front(header);	
}



inline bool DatFile::compare(uint i, uint j, string a, string b){
	return 	(atoms[i].type == a  && atoms[j].type == b) || ( atoms[i].type == b && atoms[j].type == a);
}

inline real pow2scale(real scale, real x){
	return pow2( 2.5 + scale * (x - 2.5) );
}


string DatFile::addRepulsivePair(uint a1, uint a2){
	
	if( (abs(atoms[a1].residueIndex - atoms[a2].residueIndex) > 2) && (atoms[a1].type[0] == 'C' && atoms[a2].type[0] == 'C') )  {
		
		Atom& atm1 = atoms[a1];
		Atom& atm2 = atoms[a2];
		
		real dist = atoms[a1].distanceTo(atoms[a2]);
// 		if(dist < 4.0)
// 			cerr << "Warning: atoms " << atm1.index << " (" << atm1.type << " in residue " << atm1.residueIndex << ") " 
// 			     << "and " << atm2.index << " (" << atm2.type << " in residue " << atm2.residueIndex << ") are too close! (" << dist << " Å)."
// 			     << endl;		
	}
	
	real repSigma = pow2(2.5);
	
	switch(atomSizes){
		case ASM_standard:
			if( compare(a1,a2, "CB","O") )
// 				repSigma = pow2scale(scale, 3.00);
				repSigma = pow2(3.00);
			break;
			
		case ASM_bigNCB:
			if( compare(a1,a2, "CB","O") || compare(a1,a2, "CB","N") )
				repSigma = pow2scale(scale, 3.00);
			
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2scale(scale, 3.0);
			break;
			
		case ASM_ramaNormal:			
			if     ( compare(a1,a2, "C","C") || compare(a1,a2, "C","CA") || compare(a1,a2, "C","CB") || compare(a1,a2, "CA","CA") || compare(a1,a2, "CA","CB") || compare(a1,a2, "CB","CB"))
				repSigma = pow2scale(scale, 3.20);
			
			else if( compare(a1,a2, "C","O") || compare(a1,a2, "CA","O") || compare(a1,a2, "CB","O") )
				repSigma = pow2scale(scale, 2.80);
			
			else if( compare(a1,a2, "C","N") || compare(a1,a2, "CA","N") || compare(a1,a2, "CB","N") )
				repSigma = pow2scale(scale, 2.90);
			
			else if( compare(a1,a2, "O","O") )
				repSigma = pow2scale(scale, 2.80);
			
			else if( compare(a1,a2, "O","N") )
				repSigma = pow2scale(scale, 2.70);
			
			else if( compare(a1,a2, "N","N") )
				repSigma = pow2scale(scale, 2.70);
			
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2scale(scale, 2.80);
			break;
			
		case ASM_ramaOuter:			
			if     ( compare(a1,a2, "C","C") || compare(a1,a2, "C","CA") || compare(a1,a2, "C","CB") || compare(a1,a2, "CA","CA") || compare(a1,a2, "CA","CB") || compare(a1,a2, "CB","CB"))
				repSigma = pow2scale(scale, 3.00);
			
			else if( compare(a1,a2, "C","O") || compare(a1,a2, "CA","O") || compare(a1,a2, "CB","O") )
				repSigma = pow2scale(scale, 2.70);
			
			else if( compare(a1,a2, "C","N") || compare(a1,a2, "CA","N") || compare(a1,a2, "CB","N") )
				repSigma = pow2scale(scale, 2.80);
			
			else if( compare(a1,a2, "O","O") )
				repSigma = pow2scale(scale, 2.70);
			
			else if( compare(a1,a2, "O","N") )
				repSigma = pow2scale(scale, 2.60);
			
			else if( compare(a1,a2, "N","N") )
				repSigma = pow2scale(scale, 2.60);
			
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2scale(scale, 2.70);
			break;
			
		case ASM_ramaCBCNO:
			if     ( compare(a1,a2, "C","CB") )
				repSigma = pow2scale(scale, 3.20);
			
			else if( compare(a1,a2, "CB","O") )
				repSigma = pow2scale(scale, 2.80);
			
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2scale(scale, 2.90);

			if( compare(a1,a2, "CB","O") )
				repSigma = pow2scale(scale, 2.80);
			break;
			
		case ASM_ramaCBCNO2:
			if     ( compare(a1,a2, "C","CB") )
				repSigma = pow2scale(scale, 3.20);
			
			else if( compare(a1,a2, "CB","O") )
				repSigma = pow2scale(scale, 2.80);
			
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2scale(scale, 2.90);
			
			else if( compare(a1,a2, "N","O") )
				repSigma = pow2scale(scale, 2.70);

			if( compare(a1,a2, "CB","O") )
				repSigma = pow2scale(scale, 2.80);
			break;
			
		
		
		#include "repsizes.cpp" // Generated by generate-combinatory-repsizes.pl
			
		
		case 260:
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2scale(scale, 2.60);
			break;
		case 270:
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2scale(scale, 2.70);
			break;
		case 280:
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2scale(scale, 2.80);
			break;
		case 290:
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2scale(scale, 2.90);
			break;

		
		case 300: //ramaNormal
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2(2.80);
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.90);
			break;
			
		case 301: //ramaOuter
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2(2.70);
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.80);
			break;
			
		case 310: //ramaNormal
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2(2.80);
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.90);
			else if( compare(a1,a2, "CB","C") )
				repSigma = pow2(3.20);
			else if( compare(a1,a2, "O","O") )
				repSigma = pow2(2.8);
			break;
			
		case 311: //ramaOuter
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2(2.70);
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.80);
			else if( compare(a1,a2, "CB","C") )
				repSigma = pow2(3.00);
			else if( compare(a1,a2, "O","O") )
				repSigma = pow2(2.7);
			break;
			
		case 327:
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2(2.70);			
			break;
			
		case 328:
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2(2.80);			
			break;
			
		case 330:
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2(3.00);
			break;
			
		case 350:
			//all 2.5
			break;
			
		
			
		//////////////////////
		// CB-O  2.8
		//////////////////////	
		case 365:
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2(2.80);
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.50);
			break;
		
		case 366:
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2(2.80);
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.60);
			break;
			
		case 367:
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2(2.80);
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.70);
			break;
			
		case 368:
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2(2.80);
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.80);
			break;
			
		case 369:
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2(2.80);
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.90);
			break;
			
		//////////////////////
		// CB-O  3.0
		//////////////////////
		case 375:
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2(3.00);
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.50);
			break;
		
		case 376:
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2(3.00);
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.60);
			break;
			
		case 377:
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2(3.00);
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.70);
			break;
			
		case 378:
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2(3.00);
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.80);
			break;
			
		case 379:
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2(3.00);
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.90);
			break;
			
		
		//////////////////////
		// CB-O  2.8/3.0
		//////////////////////
		case 385:
			if( compare(a1,a2, "CB","O") ){
				if(atoms[a1].residueIndex == atoms[a2].residueIndex)
					repSigma = pow2(2.80);
				else
					repSigma = pow2(3.00);
			}			
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.50);
			break;
		
		case 386:
			if( compare(a1,a2, "CB","O") ){
				if(atoms[a1].residueIndex == atoms[a2].residueIndex)
					repSigma = pow2(2.80);
				else
					repSigma = pow2(3.00);
			}
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.60);
			break;
			
		case 387:
			if( compare(a1,a2, "CB","O") ){
				if(atoms[a1].residueIndex == atoms[a2].residueIndex)
					repSigma = pow2(2.80);
				else
					repSigma = pow2(3.00);
			}
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.70);
			break;
			
		case 388:
			if( compare(a1,a2, "CB","O") ){
				if(atoms[a1].residueIndex == atoms[a2].residueIndex)
					repSigma = pow2(2.80);
				else
					repSigma = pow2(3.00);
			}
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.80);
			break;
			
		case 389:
			if( compare(a1,a2, "CB","O") ){
				if(atoms[a1].residueIndex == atoms[a2].residueIndex)
					repSigma = pow2(2.80);
				else
					repSigma = pow2(3.00);
			}
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.90);
			break; 
	
			
		//////////////////////
		// CB-O  2.5/3.0
		//////////////////////
		case 395:
			if( compare(a1,a2, "CB","O") ){
				if(atoms[a1].residueIndex == atoms[a2].residueIndex)
					repSigma = pow2(2.50);
				else
					repSigma = pow2(3.00);
			}			
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.50);
			break;
		
		case 396:
			if( compare(a1,a2, "CB","O") ){
				if(atoms[a1].residueIndex == atoms[a2].residueIndex)
					repSigma = pow2(2.50);
				else
					repSigma = pow2(3.00);
			}
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.60);
			break;
			
		case 397:
			if( compare(a1,a2, "CB","O") ){
				if(atoms[a1].residueIndex == atoms[a2].residueIndex)
					repSigma = pow2(2.50);
				else
					repSigma = pow2(3.00);
			}
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.70);
			break;
			
		case 398:
			if( compare(a1,a2, "CB","O") ){
				if(atoms[a1].residueIndex == atoms[a2].residueIndex)
					repSigma = pow2(2.50);
				else
					repSigma = pow2(3.00);
			}
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.80);
			break;
			
		case 399:
			if( compare(a1,a2, "CB","O") ){
				if(atoms[a1].residueIndex == atoms[a2].residueIndex)
					repSigma = pow2(2.50);
				else
					repSigma = pow2(3.00);
			}
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.90);
			break; 
			
		case 400:
			break;
			
		///////////////////////////////////////
		// Ho BK, Thomas A, Brasseur R (2003)
		///////////////////////////////////////	
		case 500: //100%
			if     ( compare(a1,a2, "CB","O") )
				repSigma = pow2(3.15);
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(3.30);
			else if( compare(a1,a2, "O","O") )
				repSigma = pow2(2.80);
			else if( compare(a1,a2, "O","N") )
				repSigma = pow2(2.95);
			break;
		
		case 595: //95%
			if     ( compare(a1,a2, "CB","O") )
				repSigma = pow2(3.15 * 0.95);
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(3.30 * 0.95);
			else if( compare(a1,a2, "O","O") )
				repSigma = pow2(2.80 * 0.95);
			else if( compare(a1,a2, "O","N") )
				repSigma = pow2(2.95 * 0.95);
			break;
			
		case 590: //90%
			if     ( compare(a1,a2, "CB","O") )
				repSigma = pow2(3.15 * 0.90);
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(3.30 * 0.90);
			else if( compare(a1,a2, "O","O") )
				repSigma = pow2(2.80 * 0.90);
			else if( compare(a1,a2, "O","N") )
				repSigma = pow2(2.95 * 0.90);
			break;
			
		case 585: //85%
			if     ( compare(a1,a2, "CB","O") )
				repSigma = pow2(3.15 * 0.85);
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(3.30 * 0.85);
			else if( compare(a1,a2, "O","O") )
				repSigma = pow2(2.80 * 0.85);
			else if( compare(a1,a2, "O","N") )
				repSigma = pow2(2.95 * 0.85);			
			break;
			
		///////////////////////////////////////
		// CB-O = 3.0
		///////////////////////////////////////	
		case 501:
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2(3.00);
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.50);
			break;
			
		case 502:
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2(3.00);
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.80);
			break;
			
		case 503:
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2(3.00);
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(2.90);
			break;
			
		case 504:
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2(3.00);
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(3.00);
			break;
			
		case 505:
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2(3.00);
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(3.10);
			break;
			
		case 506:
			if( compare(a1,a2, "CB","O") )
				repSigma = pow2(3.00);
			else if( compare(a1,a2, "CB","N") )
				repSigma = pow2(3.20);
			break;
			
		default:
			cerr << "Error: unknown atomSizes number " << atomSizes << "\n";			
			return "";
		
	}

	real distanceSquare = atoms[a1].distanceSquareTo(atoms[a2]);
	
	real sigma = (repSigma < distanceSquare)? repSigma : distanceSquare;
		
	char line[35];
	sprintf(line, "%5d %5d  %9.5f  %9.5f\n", 
		atoms[a1].index,
		atoms[a2].index,
		sigma,
		1.0
	);
	return line; 	
}

void DatFile::buildConnectedMap(){
	vector< vector<uint> > connected( atoms.size(), vector<uint>(0) );
	//              n-1           n           n+1
	//           ::::::::::   ::::::::::   ::::::::::
	// ... - C - N - CA - C - N - CA - C - N - CA - C - ...
	//       |       |    |       |    |       |    |
	//       O      ...   O      ...   O      ...   O	
	
	#define CONNECT(a,b)  if( atomMap[a].count(b) && atomMap[a][b]>i ){ connected[i].push_back(atomMap[a][b]); }
	
	for(uint i=0; i<atoms.size(); i++){
		
		Atom& atm = atoms[i];
		uint n = atoms[i].residueIndex - 1;
		uint lastN = atoms[atoms.size()-1].residueIndex - 1;
		
		//Backbone
		if(atm.type == "N"){	
			if(n != 0) 
				CONNECT(n-1,"C")
			CONNECT(n,"CA")
		}
		else if(atm.type == "CA"){
			CONNECT(n,"N")
			CONNECT(n,"C")	
		}
		else if(atm.type == "C"){
			CONNECT(n,"CA")
			CONNECT(n,"O")
			if(n != lastN) CONNECT(n+1,"N")
		}
		else if(atm.type == "O"){
			CONNECT(n,"C")
		}
		
		//Side chain
		if(atm.residueType != "GLY"){
			if(atm.type == "CA")
				CONNECT(n,"CB")
				else if(atm.type == "CB")
					CONNECT(n,"CA")
		}
		
		if(atm.residueType == "ARG"){
			if(atm.type == "CB"){
				CONNECT(n,"CG")
			}
			else if(atm.type == "CD"){
				CONNECT(n,"CG")
				CONNECT(n,"NE")
			}
			else if(atm.type == "CG"){
				CONNECT(n,"CB")
				CONNECT(n,"CD")
			}
			else if(atm.type == "CZ"){
				CONNECT(n,"NE")
				CONNECT(n,"NH1")
				CONNECT(n,"NH2")
			}
			else if(atm.type == "NE"){
				CONNECT(n,"CD")
				CONNECT(n,"CZ")
			}
			else if(atm.type == "NH1"){
				CONNECT(n,"CZ")
			}
			else if(atm.type == "NH2"){
				CONNECT(n,"CZ")
			}
		}
		else if(atm.residueType == "ASN"){
			if(atm.type == "CB"){
				CONNECT(n,"CG")
			}
			else if(atm.type == "CG"){
				CONNECT(n,"CB")
				CONNECT(n,"OD1")
				CONNECT(n,"ND2")
			}
			else if(atm.type == "ND2"){
				CONNECT(n,"CG")
			}
			else if(atm.type == "OD1"){
				CONNECT(n,"CG")
			}
		}
		else if(atm.residueType == "ASP"){
			if(atm.type == "CB"){
				CONNECT(n,"CG")
			}
			else if(atm.type == "CG"){
				CONNECT(n,"CB")
				CONNECT(n,"OD1")
				CONNECT(n,"OD2")
			}
			else if(atm.type == "OD1"){
				CONNECT(n,"CG")
			}
			else if(atm.type == "OD2"){
				CONNECT(n,"CG")
			}
		}
		else if(atm.residueType == "CYS"){
			if(atm.type == "CB"){
				CONNECT(n,"SG")
			}
			else if(atm.type == "SG"){
				CONNECT(n,"CB")
			}
		}
		else if(atm.residueType == "GLN"){
			if(atm.type == "CB"){
				CONNECT(n,"CG")
			}
			else if(atm.type == "CD"){
				CONNECT(n,"CG")
				CONNECT(n,"OE1")
				CONNECT(n,"NE2")
			}
			else if(atm.type == "CG"){
				CONNECT(n,"CB")
				CONNECT(n,"CD")
			}
			else if(atm.type == "NE2"){
				CONNECT(n,"CD")
			}
			else if(atm.type == "OE1"){
				CONNECT(n,"CD")
			}
		}
		else if(atm.residueType == "GLU"){
			if(atm.type == "CB"){
				CONNECT(n,"CG")
			}
			else if(atm.type == "CD"){
				CONNECT(n,"CG")
				CONNECT(n,"OE1")
				CONNECT(n,"OE2")
			}
			else if(atm.type == "CG"){
				CONNECT(n,"CB")
				CONNECT(n,"CD")
			}
			else if(atm.type == "OE1"){
				CONNECT(n,"CD")
			}
			else if(atm.type == "OE2"){
				CONNECT(n,"CD")
			}
		}
		else if(atm.residueType == "HIS"){
			if(atm.type == "CB"){
				CONNECT(n,"CG")
			}
			else if(atm.type == "CD2"){
				CONNECT(n,"CG")
			}
			else if(atm.type == "CG"){
				CONNECT(n,"CB")
				CONNECT(n,"ND1")
				CONNECT(n,"CD2")
			}
			else if(atm.type == "ND1"){
				CONNECT(n,"CG")
			}
		}
		else if(atm.residueType == "ILE"){
			if(atm.type == "CB"){
				CONNECT(n,"CG1")
				CONNECT(n,"CG2")
			}
			else if(atm.type == "CD1"){
				CONNECT(n,"CG1")
			}
			else if(atm.type == "CG1"){
				CONNECT(n,"CB")
				CONNECT(n,"CD1")
			}
			else if(atm.type == "CG2"){
				CONNECT(n,"CB")
			}
		}
		else if(atm.residueType == "LEU"){
			if(atm.type == "CB"){
				CONNECT(n,"CG")
			}
			else if(atm.type == "CD1"){
				CONNECT(n,"CG")
			}
			else if(atm.type == "CD2"){
				CONNECT(n,"CG")
			}
			else if(atm.type == "CG"){
				CONNECT(n,"CB")
				CONNECT(n,"CD1")
				CONNECT(n,"CD2")
			}
		}
		else if(atm.residueType == "LYS"){
			if(atm.type == "CB"){
				CONNECT(n,"CG")
			}
			else if(atm.type == "CD"){
				CONNECT(n,"CG")
				CONNECT(n,"CE")
			}
			else if(atm.type == "CE"){
				CONNECT(n,"CD")
				CONNECT(n,"NZ")
			}
			else if(atm.type == "CG"){
				CONNECT(n,"CB")
				CONNECT(n,"CD")
			}
			else if(atm.type == "NZ"){
				CONNECT(n,"CE")
			}
		}
		else if(atm.residueType == "MET"){
			if(atm.type == "CB"){
				CONNECT(n,"CG")
			}
			else if(atm.type == "CE"){
				CONNECT(n,"SD")
			}
			else if(atm.type == "CG"){
				CONNECT(n,"CB")
				CONNECT(n,"SD")
			}
			else if(atm.type == "SD"){
				CONNECT(n,"CG")
				CONNECT(n,"CE")
			}
		}
		else if(atm.residueType == "PHE"){
			if(atm.type == "CB"){
				CONNECT(n,"CG")
			}
			else if(atm.type == "CD1"){
				CONNECT(n,"CG")
			}
			else if(atm.type == "CD2"){
				CONNECT(n,"CG")
			}
			else if(atm.type == "CG"){
				CONNECT(n,"CB")
				CONNECT(n,"CD1")
				CONNECT(n,"CD2")
			}
		}
		else if(atm.residueType == "PRO"){
			if(atm.type == "CB"){
				CONNECT(n,"CG")
			}
			else if(atm.type == "CD"){
				CONNECT(n,"CG")
				CONNECT(n,"N")
			}
			else if(atm.type == "CG"){
				CONNECT(n,"CB")
				CONNECT(n,"CD")
			}
			else if(atm.type == "N"){
				CONNECT(n,"CA")
				CONNECT(n,"CD")
			}
		}
		else if(atm.residueType == "SER"){
			if(atm.type == "CB"){
				CONNECT(n,"OG")
			}
			else if(atm.type == "OG"){
				CONNECT(n,"CB")
			}
		}
		else if(atm.residueType == "THR"){
			if(atm.type == "CB"){
				CONNECT(n,"OG1")
				CONNECT(n,"CG2")
			}
			else if(atm.type == "CG2"){
				CONNECT(n,"CB")
			}
			else if(atm.type == "OG1"){
				CONNECT(n,"CB")
			}
		}
		else if(atm.residueType == "TRP"){
			if(atm.type == "CB"){
				CONNECT(n,"CG")
			}
			else if(atm.type == "CD1"){
				CONNECT(n,"CG")
			}
			else if(atm.type == "CD2"){
				CONNECT(n,"CG")
			}
			else if(atm.type == "CG"){
				CONNECT(n,"CB")
				CONNECT(n,"CD1")
				CONNECT(n,"CD2")
			}
		}
		else if(atm.residueType == "TYR"){
			if(atm.type == "CB"){
				CONNECT(n,"CG")
			}
			else if(atm.type == "CD1"){
				CONNECT(n,"CG")
			}
			else if(atm.type == "CD2"){
				CONNECT(n,"CG")
			}
			else if(atm.type == "CG"){
				CONNECT(n,"CB")
				CONNECT(n,"CD1")
				CONNECT(n,"CD2")
			}
		}
		else if(atm.residueType == "VAL"){
			if(atm.type == "CB"){
				CONNECT(n,"CG1")
				CONNECT(n,"CG2")
			}
			else if(atm.type == "CG1"){
				CONNECT(n,"CB")
			}
			else if(atm.type == "CG2"){
				CONNECT(n,"CB")
			}
		}
		
	}
	
	//map atoms connected up to 2 bonds
	connected2.resize( atoms.size(), vector<uint>(0) );
	for(uint i=0; i<atoms.size(); i++){
		for(uint j=0; j<connected[i].size(); j++){
			uint aj = connected[i][j];
			//add atom connected by 1 bond
			connected2[i].push_back(aj);
			
			//add atoms connected by 2 bonds (upstream)
			for(uint k=0; k<connected[aj].size(); k++)
				connected2[i].push_back(connected[aj][k]);
			
			//add atoms connected by 2 bonds (downstream)
			for(uint k=0; k<connected[i].size(); k++){
				if(j==k) continue;
	 
				uint ak = connected[i][k];
				if(ak > aj)
					connected2[aj].push_back(ak);
			}
		}
	}
}

/**
 * Determines whether atoms a1 and a2 are not connected 
 * by at most 2 covalent bonds.
 */
bool DatFile::unboundedBy2(uint a1, uint a2){
	for(uint i=0; i<connected2[a1].size(); i++)
		if(connected2[a1][i] == a2)
			return false;
		return true;
}

/**
 * Determines whether atoms a1 and a2 are not CA atoms connected
 * to the same peptide bond.
 */
bool DatFile::notNeighboringCA(uint a1, uint a2){
	Atom& atm1 = atoms[a1];
	Atom& atm2 = atoms[a2];
	if(atm1.type != "CA" || atm2.type != "CA")
		return true;
	return abs(atm1.residueIndex - atm2.residueIndex) > 1;
}

/**
 * Determines whether atoms a1 and a2 are not CA and O atoms connected
 * in the same peptide bond.
 */
bool DatFile::notNeighboringOCA(uint a1, uint a2){
	Atom& atm1 = atoms[a1];
	Atom& atm2 = atoms[a2];
	
	if(atm1.type == "O" && atm2.type == "CA")
		return atm1.residueIndex != atm2.residueIndex - 1;
	
	if(atm2.type == "O" && atm1.type == "CA")
		return atm2.residueIndex != atm1.residueIndex - 1;
	
	return true;
}

bool DatFile::containedInDihedral(string a, string b, uint n, const string dihedrals[][4]){
	for(uint i=0; i<n; i++){
		bool a_in = false;
		bool b_in = false;
		
		for(uint j=0; j<4; j++)
			if(dihedrals[i][j] == a){
				a_in = true;
				break;
			}
			
		if(a_in == false)
			continue;
		
		for(uint j=0; j<4; j++)
			if(dihedrals[i][j] == b){
				b_in = true;
				break;
			}
		
		if(b_in == true)
			return true;
	}
	return false;
}

bool DatFile::dontFormDihedral(uint a1, uint a2){
	static const string TYR_DIHEDRALS[][4] = {
		{"CB", "CG", "CD1", "CE1"},
		{"CB", "CG", "CD2", "CE2"},
		{"CG", "CD1", "CE1", "CZ"},
		{"CD1", "CE1", "CZ", "CE2"},
		{"CE1", "CZ", "CE2", "CD2"},
		{"CZ", "CE2", "CD2", "CG"},
		{"CE2", "CD2", "CG", "CD1"},
		{"CD2", "CG", "CD1", "CE1"},
		{"CD1", "CE1", "CZ", "OH"},
		{"CD2", "CE2", "CZ", "OH"}
	};
	static const string PHE_DIHEDRALS[][4] = {
		{"CB", "CG", "CD1", "CE1"},
		{"CB", "CG", "CD2", "CE2"},
		{"CG", "CD1", "CE1", "CZ"},
		{"CD1", "CE1", "CZ", "CE2"},
		{"CE1", "CZ", "CE2", "CD2"},
		{"CZ", "CE2", "CD2", "CG"},
		{"CE2", "CD2", "CG", "CD1"},
		{"CD2", "CG", "CD1", "CE1"}
	};
	static const string HIS_DIHEDRALS[][4] = {
		{"CB", "CG", "ND1", "CE1"},
		{"CB", "CG", "CD2", "NE2"},
		{"CG", "ND1", "CE1", "NE2"},
		{"ND1", "CE1", "NE2", "CD2"},
		{"CE1", "NE2", "CD2", "CG"},
		{"NE2", "CD2", "CG", "ND1"},
		{"CD2", "CG", "ND1", "CE1"}
	};
	static const string PRO_DIHEDRALS[][4] = {
		{"CA", "CB", "CG", "CD"},
		{"CB", "CG", "CD", "N"},
		{"CG", "CD", "N", "CA"},
		{"CD", "N", "CA", "CB"}
	};
	
	Atom& atm1 = atoms[a1];
	Atom& atm2 = atoms[a2];
	
	if(atm1.residueIndex != atm2.residueIndex)
		return true;
	
	if(atm1.residueType == "TYR")
		return !containedInDihedral(atm1.type, atm2.type, 10, TYR_DIHEDRALS);
	
	if(atm1.residueType == "PHE")
		return !containedInDihedral(atm1.type, atm2.type, 8, PHE_DIHEDRALS);
	
	if(atm1.residueType == "HIS")
		return !containedInDihedral(atm1.type, atm2.type, 7, HIS_DIHEDRALS);
	
	if(atm1.residueType == "PRO")
		return !containedInDihedral(atm1.type, atm2.type, 4, PRO_DIHEDRALS);
	
	return true;
}

void DatFile::populateStructure(list< string >& lines){
	char header1[30];
	char header2[49];
	sprintf(header1, "%12d Atom Positions.\n", (int)atoms.size());
	sprintf(header2, "%12d is the length of chain           1\n", (int)atoms.size());	
	
	lines.push_back(header1);	
	lines.push_back("           1 chains\n");
	lines.push_back(header2);	
	
	uint index = 1;
	for(vector<Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it){
		char line[58];
		sprintf(line, "%5d %4d  %-3s %s %8.3f %8.3f %8.3f %8.3f\n",
			index++, it->residueIndex, 
			it->type.c_str(), it->residueType.c_str(),
			it->x, it->y, it->z,
			1.0
		);
		lines.push_back(line);
	}
	
}
