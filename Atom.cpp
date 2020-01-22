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
