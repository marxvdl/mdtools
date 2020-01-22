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

#include <iostream>
#include <fstream>
#include <string>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include "Protein.h"
#include "rcntr.h"

using namespace std;

namespace po = boost::program_options;

enum Operation { RCNTR, HBONDS, CHIRAL, RAMA, SEQ, GAR, DAT, PDB, HGRAM, BDATA, CLASH };

void usage(po::options_description& desc);

/**
 * Parses the command line arguments and invokes the appropriate tool.
 */
int main(int argc, char **argv) {
	
	//
	// 1. Parse the arguments
	//

	int layers;
	string seqType;
	string inputFilename;
	string outputFilename;
	string inputFormat;
	string referenceAtom;
	string hmmpredFilename;
	string atomSizes;
	string bdataRefFile;
	real scale;
	real midLayerWidth;
	real bottomLayerWidth;

	po::positional_options_description p;
	po::options_description desc("");
	desc.add_options()
	    ("help", "Displays this help and exit ")
	    ("operation",         po::value<string>())
	    ("all,a", "")
	    ("input,i",           po::value<string>(&inputFilename)    -> default_value("")        )
	    ("output,o",          po::value<string>(&outputFilename)   -> default_value("")        )
	    ("layers,l",          po::value<int>   (&layers)           -> default_value(0)         )
	    ("seqtype,t",         po::value<string>(&seqType)          -> default_value("Letter")  )
	    ("input-format,f",    po::value<string>(&inputFormat)      -> default_value("")        )
	    ("reference-atom,r",  po::value<string>(&referenceAtom)    -> default_value("")        )
	    ("hmmpred,h",         po::value<string>(&hmmpredFilename)  -> default_value("")        )
	    ("atom-sizes,s",      po::value<string>(&atomSizes)        -> default_value("standard"))
	    ("scale,c",           po::value<real>  (&scale)            -> default_value(1)         )
	    ("bdata-ref-file,e",  po::value<string>(&bdataRefFile)     -> default_value("")        )
	    ("normalize,n", "")
	    ("keep-center,k", "")
	    ("mid-layer-width,w", po::value<real>  (&midLayerWidth)    -> default_value(-1)        )
	    ("bot-layer-width,y", po::value<real>  (&bottomLayerWidth) -> default_value(-1)        )
	    ("quadratic,q", "")
	;

	p.add("operation", 1);
	p.add("input", 1);
	p.add("output", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).
			options(desc).positional(p).run(), vm);
	po::notify(vm);
	
	if ( vm.count("help")  || !vm.count("operation") ) {
		usage(desc);
	    return 1;
	}

	bool normalize = vm.count("normalize");
	bool quadratic = vm.count("quadratic");

	Operation operation;
	string opstring = vm["operation"].as<string>();
	if(opstring == "rcntr")
			operation = RCNTR;
	else if(opstring == "gar")
		operation = GAR;
	else if(opstring == "dat")
		operation = DAT;
	else if(opstring == "pdb")
		operation = PDB;
	else if(opstring == "hbonds")
		operation = HBONDS;
	else if(opstring == "chiral")
		operation = CHIRAL;
	else if(opstring == "rama")
		operation = RAMA;
	else if(opstring == "seq")
		operation = SEQ;
	else if(opstring == "hgram")
		operation = HGRAM;
	else if(opstring == "bdata")
		operation = BDATA;
	else if(opstring == "clash")
		operation = CLASH;
	else {
		cout << "Error: operation must be one of \"rcntr gar dat pdb hbonds chiral rama seq hgram bdata\"." << endl;
		cout << endl;
		usage(desc);
		return 1;
	}

	if(vm.count("input"))
		inputFilename = vm["input"].as<string>();

	if(vm.count("output"))
		outputFilename = vm["output"].as<string>();

	//
	// 2. Invoke appropriate tool
	//
	Protein protein(vm.count("keep-center"), midLayerWidth, bottomLayerWidth);
	if(inputFormat == ""){
		if(operation == SEQ)
			inputFormat = "rg5";
		else if (operation == DAT)
			inputFormat = "pdb";
		else
			inputFormat = "md";
	}
	else {
		if(inputFormat != "rg5" && inputFormat != "md" && inputFormat != "pdb") {
			cout << "Error: invalid input format \"" << inputFormat << "\"." << endl;
			return 2;
		}
	}


	if(inputFilename == "") {
		if(inputFormat == "rg5")
			protein.loadFromRg5DB(cin);
		else if(inputFormat == "md")
			protein.loadFromMDRead(cin);
		else if(inputFormat == "pdb")
			protein.loadFromPDB(cin);
	}
	else{
		ifstream inputStream(inputFilename.c_str());
		if(inputFormat == "rg5")
			protein.loadFromRg5DB(inputStream);
		else if(inputFormat == "md")
			protein.loadFromMDRead(inputStream);
		else if(inputFormat == "pdb")
			protein.loadFromPDB(inputStream);
	}

	try{
		switch(operation){
			case RCNTR:
				if(hmmpredFilename == ""){
					if(outputFilename == "")
							protein.generateRcntr(layers, quadratic, cout, referenceAtom);	
					else {
						ofstream outputStream(outputFilename.c_str());
						protein.generateRcntr(layers, quadratic, outputStream, referenceAtom);
					}
				}
				else{
					vector< vector<real> > probs = parseHmmPredOutput(hmmpredFilename);
					if(outputFilename == "")
						protein.generateRcntr(layers, quadratic, probs, cout, referenceAtom);

					else {
						ofstream outputStream(outputFilename.c_str());
						protein.generateRcntr(layers, quadratic, probs, outputStream, referenceAtom);
					}
				}
				break;

			case GAR:
				if(outputFilename == "")
						protein.generateGar(layers, cout);
					else {
						ofstream outputStream(outputFilename.c_str());
						protein.generateGar(layers, outputStream);
					}
				break;
				
			case DAT:
				
				ASM asM;
				if(atomSizes == "standard")
					asM = ASM_standard;
				else if(atomSizes == "bigNCB")
					asM = ASM_bigNCB;
				else if(atomSizes == "ramaNormal")
					asM = ASM_ramaNormal;
				else if(atomSizes == "ramaOuter")
					asM = ASM_ramaOuter;
				else if(atomSizes == "ramaCBCNO")
					asM = ASM_ramaCBCNO;
				else if(atomSizes == "ramaCBCNO2")
					asM = ASM_ramaCBCNO2;
				else if(atomSizes[0] == 'n')
					asM = string2int(atomSizes.substr(1));
				else{
					cerr << "Error: Atom size should be one of: standard bigNCB ramaNormal ramaOuter ramaCBCNO ramaCBCNO2\n";
					return 2;
				}
				
				if(scale < 0 or scale > 1){
					cerr << "Error: Scale should be a number between 0 and 1.\n";
					return 2;
				}
				
				try{
					if(outputFilename == "")
						protein.generateDat(layers, asM, scale, cout);
					else {
						ofstream outputStream(outputFilename.c_str());
						protein.generateDat(layers, asM, scale, outputStream);
					}
				}
				catch(char const* msg){
					cerr << "Error: " << msg << "\n";
					return 2;
				}
				break;
				
			case PDB:
				try{
					if(outputFilename == "")
						protein.generatePDB(layers, cout);
					else {
						ofstream outputStream(outputFilename.c_str());
						protein.generatePDB(layers, outputStream);
					}
				}
				catch(char const* msg){
					cerr << "Error: " << msg << "\n";
					return 2;
				}
				break;	

			case HBONDS:
				if(outputFilename == "")
					protein.generateHbonds(cout, vm.count("all"));
				else {
					ofstream outputStream(outputFilename.c_str());
					protein.generateHbonds(outputStream, vm.count("all"));
				}
				break;

			case CHIRAL:
				if(outputFilename == "")
					protein.generateChiral(cout);
				else {
					ofstream outputStream(outputFilename.c_str());
					protein.generateChiral(outputStream);
				}
				break;
				
			case RAMA:
				if(outputFilename == "")
					protein.generateRama(cout);
				else {
					ofstream outputStream(outputFilename.c_str());
					protein.generateRama(outputStream);
				}
				break;
				
			case HGRAM:
				if(outputFilename == "")
					protein.generateHistogram(cout, referenceAtom);
				else {
					ofstream outputStream(outputFilename.c_str());
					protein.generateHistogram(outputStream, referenceAtom);
				}
				break;
				
			case BDATA:
			if(outputFilename == "")
					protein.generateBData(cout, bdataRefFile);
				else {
					ofstream outputStream(outputFilename.c_str());
					protein.generateBData(outputStream, bdataRefFile);
				}
				break;
			
			case CLASH:
			if(outputFilename == "")
					protein.generateClashReport(cout);
				else {
					ofstream outputStream(outputFilename.c_str());
					protein.generateClashReport(outputStream);
				}
				break;

			case SEQ:
				if(referenceAtom == "")
					referenceAtom = "CA";
				
				string seq;
				if(seqType == "HP")
					seq = protein.toHP();
				else if(seqType == "HPN")
					seq = protein.toHPN();
				else if(seqType == "HPNAC")
					seq = protein.toHPNAC();
				else if(seqType == "Letter")
					seq = protein.toLetter();
				else if(seqType == "HPNhp")
					seq = protein.toHPNhp();
				else if(seqType == "HPNhpn")
					seq = protein.toHPNhpn();
				else if(seqType == "burial"){
					if(referenceAtom == "CAdir")
						seq = protein.toBurialSequenceCAdir(layers);
					else if(referenceAtom == "CBdir")
						seq = protein.toBurialSequenceCAdir(layers, true);
					else
						seq = protein.toBurialSequence(layers, referenceAtom, normalize, inputFormat=="rg5");
				}
				else if(seqType == "rburial")
					seq = protein.toRelativeBurialSequence(layers, referenceAtom);
				else {
					cerr << "Invalid mapping: " << seqType << endl;
					return 1;
				}

				if(outputFilename == "")
					cout << seq << endl;
				else {
					ofstream outputStream(outputFilename.c_str());
					outputStream << seq << endl;
				}

				break;
				

		}
	}
	catch(char const* s){
		cerr << "Error: " << s << endl;
		return 1;
	}
	return 0;

}

void usage(po::options_description& desc){
	cout << "Usage: MDTools OPERATION [options...] [INPUT] [OUTPUT]\n"
	     << "       MDTools OPERATION [options...] [-i INPUT] [-o OUTPUT]\n"
	     << "\n"
	     << "Where the valid options depend on the OPERATION performed.\n"
	     << "OPERATION should be one of:\n"
	     << "\n"
	     << "  rcntr  -> Generates the \"rcntr\" file with the residue burial info.\n"
	     << "            Input file should be an MD file generated by Read.exe\n"
	     << "            Options:\n"
	     << "               -l, --layers=N        Splits the protein into N layers\n"
	     << "                                     N=0: use exact distances (default).\n"
	     << "                                     N=1: burials in radii of gyration.\n"
	     << "\n"
	     << "               -h, --hmmpred=FILE    Uses informations from an HmmPred results\n"
	     << "                                     file with probabilities.\n"
	     << "                                     (Will grab the first result).\n"
	     << "\n"
	     << "               -r, --reference-atom  Atom predicted in HmmPred results.\n"
	     << "\n"
	     << "               -q, --quadratic       Produces output in the quadratic potential format.\n"
	     << "\n"
	     << "  gar    -> Generates the \".gar\" file with the residue burial info formatted\n"
	     << "            as burial restraints for GAPF.\n"
	     << "            Options: Same as above.\n"
	     << "\n"
	     << "  dat    -> Generates the \".dat\" file with the topologies for MD.\n"
	     << "\n"
	     << "               -s, --atom-sizes      Atom radii for non-contacts. Should be one of:\n"
	     << "\n"
	     << "                                     standard   :  CB-O = 3.0; all others = 2.5 (default)\n"
	     << "                                     bigNCB     :  Same as above, but 3.0 for CB-N\n"
	     << "                                     ramaNormal :  Ramachandran's \"Normally allowed\"\n"
	     << "                                     ramaOuter  :  Ramachandran's \"Outer limit\"\n"
	     << "                                     ramaCBCNO  :  Like ramaNormal, but only for CB-N, CB-C, CB-O\n"
	     << "\n"
	     << "               -c, --scale           Scale atom sizes in relation to standard (0 to 1).\n"
	     << "\n"
	     << "  pdb    -> Generates a PDB file.\n"
	     << "\n"
	     << "  chiral -> Generates the \"chiral\" file with chirality data.\n"
	     << "            By default, input file should be an MD file generated by Read.exe\n"
	     << "\n"
	     << "  hbonds -> Generates the \"hbonds\" file with hydrogen bonding data from the.\n"
	     << "            protein backbone\n"
	     << "            By default, input file should be an MD file generated by Read.exe\n"
	     << "            Options:\n"
	     << "               -a, --all             Also generates hydrogen bond pairs from\n"
	     << "                                     sidechains.\n"
	     << "  rama   -> Generates the \"rama\" file with phi/psi angle parameters according\n"
	     << "            to the CHARMM22* potential (Piana, Lindorff-Larsen, Shaw (2011)).\n"
	     << "\n"
	     << "  seq    -> Converts a structure file into a sequence file.\n"
	     << "            By default, input file should be a database \".rg5\" file\n"
	     << "            Options:\n"
	     << "               -t, --seqtype=TYPE    Choose mapping to be used on output.\n"
	     << "                                     TYPE should be one of:\n"
	     << "                                     * HP\n"
	     << "                                     * HPN\n"
	     << "                                     * HPNAC\n"
	     << "                                     * Letter (default)\n"
	     << "                                     * HPNhp\n"
	     << "                                     * HPNhpn\n"
	     << "                                     * burial    \\ These options also require the\n"
	     << "                                     * rburial   / -l switch (see \"rcntr\" above)\n"
	     << "\n"
	     << "               -n, --normalize       If \"burial\" or \"rburial\" are chosen as\n"
	     << "                                     the sequence type and exact distances are used\n"
	     << "                                     (--layers=0), normalize distances from 0 to 1.\n"
	     << "\n"
	     << "               -r, --reference-atom  If \"burial\" or \"rburial\" are chosen as\n"
	     << "                                     the sequence type, specifies the atom whose\n"
	     << "                                     burial is to be calculated.\n"
	     << "                                     (default: CA)\n"
	     << "                                     (special values: all, CAdir)\n"
	     << "\n"
	     << "  hgram  -> Plots a burial histogram.\n"
	     << "            Options:\n"
	     << "               -r, --reference-atom  If present, plots only the burial of given\n"
	     << "                                     atom.\n"
	     << "\n"
	     << "  bdata  -> Generates a file that stores relevant protein burial data.\n"
	     << "               -e, --bdata-ref-file  If present, uses burial layer limits from\n"
	     << "                                     given bdata file.\n"
	     << "\n"
	     << "  clash  -> Generates a atom-clash report for the protein.\n"
	     << "\n"
	     << "In all cases:\n"
	     << "   -f, --input-format   Changes the input file type.\n"
	     << "                        Should be one of: md, rg5, pdb\n"
	     << "\n"
     	     << "   -k, --keep-center    Calculate all burials in relation to (0,0,0).\n"	     
	     << "\n"
	     << "if INPUT or OUTPUT are omitted, stdin and stdout are used, respectively.\n"
	     << "\n"	
	     << "MDTools was compiled on " << __DATE__ << ", " << __TIME__
	     #ifdef __INTEL_COMPILER
	     << " using the Intel C++ Compiler version " << 	__INTEL_COMPILER << "." << __INTEL_COMPILER_BUILD_DATE
	     #else
	     #ifdef __GNUC__
	     << " using the GNU C++ compiler version " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__
	     #endif
	     #endif
	     #ifdef __LP64__
	     << " (64-bit)"
	     #else
	     << " (32-bit)"
	     #endif
	     ".\n"; 	     
}

