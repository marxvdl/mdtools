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

#ifndef MISC_H_
#define MISC_H_

#include <sstream>
#include <string>
#include <cmath>

//typedef double real;
#define real double

typedef unsigned int uint;
typedef int ASM;
const real PI = 3.1415926535897932384626433832795;

using namespace std;

enum AtomSizesModel { ASM_standard, ASM_bigNCB, ASM_ramaNormal, ASM_ramaOuter, ASM_ramaCBCNO, ASM_ramaCBCNO2 };

inline string trim(const string &str){
    size_t s = str.find_first_not_of(" \n\r\t");
    size_t e = str.find_last_not_of (" \n\r\t");

    if(( string::npos == s) || ( string::npos == e))
        return "";
    else
        return str.substr(s, e-s+1);
}


inline int string2int(string s){
	int result;
	stringstream(s) >> result;
	return result;
}

inline string int2string(int n){
	stringstream ss;
	ss << n;
	return ss.str();
}

inline real string2real(string s){
	real result;
	stringstream(s) >> result;
	return result;
}

inline bool stringStartsWith(string big, string small){
	if(big.length() < small.length())
		return false;

	for(uint i=0; i<small.length(); i++)
		if(big[i] != small[i])
			return false;

	return true;
}

inline string real2string(real n){
	stringstream ss;
	ss << n;
	return ss.str();
}

inline real pow2(real n){
	return n*n;
}

inline real cuberoot(real n){
	return pow(n, 1.0/3.0);
}


inline bool stringIsIn(string& s, const string list[], uint n){
	for(uint i=0; i<n; i++)
		if(list[i] == s)
			return true;
		return false;
}

#endif /* MISC_H_ */

