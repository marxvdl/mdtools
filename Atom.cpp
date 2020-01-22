/*
 * Atom.cpp
 *
 *  Created on: Mar 22, 2010
 *      Author: marx
 */

#include "Atom.h"

void Atom::calculateBurial(int levels, real cx, real cy, real cz){
	real dx = x - cx;
	real dy = y - cy;
	real dz = z - cz;
	distanceToCenter = sqrt(dx*dx + dy*dy + dz*dz);
}

void Atom::calculateRealBurial(real cx, real cy, real cz){
	real dx = x - cx;
	real dy = y - cy;
	real dz = z - cz;
	distanceToCenter = sqrt(dx*dx + dy*dy + dz*dz);
}
