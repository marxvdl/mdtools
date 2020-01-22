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

// C 1.3
// O 1.3
// N 1.3
case 100:
	if     ( compare(a1,a2, "C","C") || compare(a1,a2, "C","CA") || compare(a1,a2, "C","CB") || compare(a1,a2, "CA","CA") || compare(a1,a2, "CA","CB") || compare(a1,a2, "CB","CB"))
		repSigma = pow2(2.6);
	
	else if( compare(a1,a2, "C","O") || compare(a1,a2, "CA","O") || compare(a1,a2, "CB","O") )
		repSigma = pow2(2.6);
	
	else if( compare(a1,a2, "C","N") || compare(a1,a2, "CA","N") || compare(a1,a2, "CB","N") )
		repSigma = pow2(2.6);
	
	else if( compare(a1,a2, "O","O") )
		repSigma = pow2(2.6);
	
	else if( compare(a1,a2, "O","N") )
		repSigma = pow2(2.6);
	
	else if( compare(a1,a2, "N","N") )
		repSigma = pow2(2.6);
	
	if( compare(a1,a2, "CB","O") )
		repSigma = pow2(2.6);
break;		

// C 1.3
// O 1.35
// N 1.3
case 101:
	if     ( compare(a1,a2, "C","C") || compare(a1,a2, "C","CA") || compare(a1,a2, "C","CB") || compare(a1,a2, "CA","CA") || compare(a1,a2, "CA","CB") || compare(a1,a2, "CB","CB"))
		repSigma = pow2(2.6);
	
	else if( compare(a1,a2, "C","O") || compare(a1,a2, "CA","O") || compare(a1,a2, "CB","O") )
		repSigma = pow2(2.65);
	
	else if( compare(a1,a2, "C","N") || compare(a1,a2, "CA","N") || compare(a1,a2, "CB","N") )
		repSigma = pow2(2.6);
	
	else if( compare(a1,a2, "O","O") )
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "O","N") )
		repSigma = pow2(2.65);
	
	else if( compare(a1,a2, "N","N") )
		repSigma = pow2(2.6);
	
	if( compare(a1,a2, "CB","O") )
		repSigma = pow2(2.65);
break;		

// C 1.3
// O 1.3
// N 1.35
case 102:
	if     ( compare(a1,a2, "C","C") || compare(a1,a2, "C","CA") || compare(a1,a2, "C","CB") || compare(a1,a2, "CA","CA") || compare(a1,a2, "CA","CB") || compare(a1,a2, "CB","CB"))
		repSigma = pow2(2.6);
	
	else if( compare(a1,a2, "C","O") || compare(a1,a2, "CA","O") || compare(a1,a2, "CB","O") )
		repSigma = pow2(2.6);
	
	else if( compare(a1,a2, "C","N") || compare(a1,a2, "CA","N") || compare(a1,a2, "CB","N") )
		repSigma = pow2(2.65);
	
	else if( compare(a1,a2, "O","O") )
		repSigma = pow2(2.6);
	
	else if( compare(a1,a2, "O","N") )
		repSigma = pow2(2.65);
	
	else if( compare(a1,a2, "N","N") )
		repSigma = pow2(2.7);
	
	if( compare(a1,a2, "CB","O") )
		repSigma = pow2(2.6);
break;		

// C 1.3
// O 1.35
// N 1.35
case 103:
	if     ( compare(a1,a2, "C","C") || compare(a1,a2, "C","CA") || compare(a1,a2, "C","CB") || compare(a1,a2, "CA","CA") || compare(a1,a2, "CA","CB") || compare(a1,a2, "CB","CB"))
		repSigma = pow2(2.6);
	
	else if( compare(a1,a2, "C","O") || compare(a1,a2, "CA","O") || compare(a1,a2, "CB","O") )
		repSigma = pow2(2.65);
	
	else if( compare(a1,a2, "C","N") || compare(a1,a2, "CA","N") || compare(a1,a2, "CB","N") )
		repSigma = pow2(2.65);
	
	else if( compare(a1,a2, "O","O") )
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "O","N") )
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "N","N") )
		repSigma = pow2(2.7);
	
	if( compare(a1,a2, "CB","O") )
		repSigma = pow2(2.65);
break;		

// C 1.35
// O 1.3
// N 1.3
case 104:
	if     ( compare(a1,a2, "C","C") || compare(a1,a2, "C","CA") || compare(a1,a2, "C","CB") || compare(a1,a2, "CA","CA") || compare(a1,a2, "CA","CB") || compare(a1,a2, "CB","CB"))
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "C","O") || compare(a1,a2, "CA","O") || compare(a1,a2, "CB","O") )
		repSigma = pow2(2.65);
	
	else if( compare(a1,a2, "C","N") || compare(a1,a2, "CA","N") || compare(a1,a2, "CB","N") )
		repSigma = pow2(2.65);
	
	else if( compare(a1,a2, "O","O") )
		repSigma = pow2(2.6);
	
	else if( compare(a1,a2, "O","N") )
		repSigma = pow2(2.6);
	
	else if( compare(a1,a2, "N","N") )
		repSigma = pow2(2.6);
	
	if( compare(a1,a2, "CB","O") )
		repSigma = pow2(2.65);
break;		

// C 1.35
// O 1.35
// N 1.3
case 105:
	if     ( compare(a1,a2, "C","C") || compare(a1,a2, "C","CA") || compare(a1,a2, "C","CB") || compare(a1,a2, "CA","CA") || compare(a1,a2, "CA","CB") || compare(a1,a2, "CB","CB"))
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "C","O") || compare(a1,a2, "CA","O") || compare(a1,a2, "CB","O") )
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "C","N") || compare(a1,a2, "CA","N") || compare(a1,a2, "CB","N") )
		repSigma = pow2(2.65);
	
	else if( compare(a1,a2, "O","O") )
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "O","N") )
		repSigma = pow2(2.65);
	
	else if( compare(a1,a2, "N","N") )
		repSigma = pow2(2.6);
	
	if( compare(a1,a2, "CB","O") )
		repSigma = pow2(2.7);
break;		

// C 1.35
// O 1.3
// N 1.35
case 106:
	if     ( compare(a1,a2, "C","C") || compare(a1,a2, "C","CA") || compare(a1,a2, "C","CB") || compare(a1,a2, "CA","CA") || compare(a1,a2, "CA","CB") || compare(a1,a2, "CB","CB"))
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "C","O") || compare(a1,a2, "CA","O") || compare(a1,a2, "CB","O") )
		repSigma = pow2(2.65);
	
	else if( compare(a1,a2, "C","N") || compare(a1,a2, "CA","N") || compare(a1,a2, "CB","N") )
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "O","O") )
		repSigma = pow2(2.6);
	
	else if( compare(a1,a2, "O","N") )
		repSigma = pow2(2.65);
	
	else if( compare(a1,a2, "N","N") )
		repSigma = pow2(2.7);
	
	if( compare(a1,a2, "CB","O") )
		repSigma = pow2(2.65);
break;		

// C 1.35
// O 1.35
// N 1.35
case 107:
	if     ( compare(a1,a2, "C","C") || compare(a1,a2, "C","CA") || compare(a1,a2, "C","CB") || compare(a1,a2, "CA","CA") || compare(a1,a2, "CA","CB") || compare(a1,a2, "CB","CB"))
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "C","O") || compare(a1,a2, "CA","O") || compare(a1,a2, "CB","O") )
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "C","N") || compare(a1,a2, "CA","N") || compare(a1,a2, "CB","N") )
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "O","O") )
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "O","N") )
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "N","N") )
		repSigma = pow2(2.7);
	
	if( compare(a1,a2, "CB","O") )
		repSigma = pow2(2.7);
break;		

// C 1.4
// O 1.3
// N 1.3
case 108:
	if     ( compare(a1,a2, "C","C") || compare(a1,a2, "C","CA") || compare(a1,a2, "C","CB") || compare(a1,a2, "CA","CA") || compare(a1,a2, "CA","CB") || compare(a1,a2, "CB","CB"))
		repSigma = pow2(2.8);
	
	else if( compare(a1,a2, "C","O") || compare(a1,a2, "CA","O") || compare(a1,a2, "CB","O") )
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "C","N") || compare(a1,a2, "CA","N") || compare(a1,a2, "CB","N") )
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "O","O") )
		repSigma = pow2(2.6);
	
	else if( compare(a1,a2, "O","N") )
		repSigma = pow2(2.6);
	
	else if( compare(a1,a2, "N","N") )
		repSigma = pow2(2.6);
	
	if( compare(a1,a2, "CB","O") )
		repSigma = pow2(2.7);
break;		

// C 1.4
// O 1.35
// N 1.3
case 109:
	if     ( compare(a1,a2, "C","C") || compare(a1,a2, "C","CA") || compare(a1,a2, "C","CB") || compare(a1,a2, "CA","CA") || compare(a1,a2, "CA","CB") || compare(a1,a2, "CB","CB"))
		repSigma = pow2(2.8);
	
	else if( compare(a1,a2, "C","O") || compare(a1,a2, "CA","O") || compare(a1,a2, "CB","O") )
		repSigma = pow2(2.75);
	
	else if( compare(a1,a2, "C","N") || compare(a1,a2, "CA","N") || compare(a1,a2, "CB","N") )
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "O","O") )
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "O","N") )
		repSigma = pow2(2.65);
	
	else if( compare(a1,a2, "N","N") )
		repSigma = pow2(2.6);
	
	if( compare(a1,a2, "CB","O") )
		repSigma = pow2(2.75);
break;		

// C 1.4
// O 1.3
// N 1.35
case 110:
	if     ( compare(a1,a2, "C","C") || compare(a1,a2, "C","CA") || compare(a1,a2, "C","CB") || compare(a1,a2, "CA","CA") || compare(a1,a2, "CA","CB") || compare(a1,a2, "CB","CB"))
		repSigma = pow2(2.8);
	
	else if( compare(a1,a2, "C","O") || compare(a1,a2, "CA","O") || compare(a1,a2, "CB","O") )
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "C","N") || compare(a1,a2, "CA","N") || compare(a1,a2, "CB","N") )
		repSigma = pow2(2.75);
	
	else if( compare(a1,a2, "O","O") )
		repSigma = pow2(2.6);
	
	else if( compare(a1,a2, "O","N") )
		repSigma = pow2(2.65);
	
	else if( compare(a1,a2, "N","N") )
		repSigma = pow2(2.7);
	
	if( compare(a1,a2, "CB","O") )
		repSigma = pow2(2.7);
break;		

// C 1.4
// O 1.35
// N 1.35
case 111:
	if     ( compare(a1,a2, "C","C") || compare(a1,a2, "C","CA") || compare(a1,a2, "C","CB") || compare(a1,a2, "CA","CA") || compare(a1,a2, "CA","CB") || compare(a1,a2, "CB","CB"))
		repSigma = pow2(2.8);
	
	else if( compare(a1,a2, "C","O") || compare(a1,a2, "CA","O") || compare(a1,a2, "CB","O") )
		repSigma = pow2(2.75);
	
	else if( compare(a1,a2, "C","N") || compare(a1,a2, "CA","N") || compare(a1,a2, "CB","N") )
		repSigma = pow2(2.75);
	
	else if( compare(a1,a2, "O","O") )
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "O","N") )
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "N","N") )
		repSigma = pow2(2.7);
	
	if( compare(a1,a2, "CB","O") )
		repSigma = pow2(2.75);
break;		

// C 1.45
// O 1.3
// N 1.3
case 112:
	if     ( compare(a1,a2, "C","C") || compare(a1,a2, "C","CA") || compare(a1,a2, "C","CB") || compare(a1,a2, "CA","CA") || compare(a1,a2, "CA","CB") || compare(a1,a2, "CB","CB"))
		repSigma = pow2(2.9);
	
	else if( compare(a1,a2, "C","O") || compare(a1,a2, "CA","O") || compare(a1,a2, "CB","O") )
		repSigma = pow2(2.75);
	
	else if( compare(a1,a2, "C","N") || compare(a1,a2, "CA","N") || compare(a1,a2, "CB","N") )
		repSigma = pow2(2.75);
	
	else if( compare(a1,a2, "O","O") )
		repSigma = pow2(2.6);
	
	else if( compare(a1,a2, "O","N") )
		repSigma = pow2(2.6);
	
	else if( compare(a1,a2, "N","N") )
		repSigma = pow2(2.6);
	
	if( compare(a1,a2, "CB","O") )
		repSigma = pow2(2.75);
break;		

// C 1.45
// O 1.35
// N 1.3
case 113:
	if     ( compare(a1,a2, "C","C") || compare(a1,a2, "C","CA") || compare(a1,a2, "C","CB") || compare(a1,a2, "CA","CA") || compare(a1,a2, "CA","CB") || compare(a1,a2, "CB","CB"))
		repSigma = pow2(2.9);
	
	else if( compare(a1,a2, "C","O") || compare(a1,a2, "CA","O") || compare(a1,a2, "CB","O") )
		repSigma = pow2(2.8);
	
	else if( compare(a1,a2, "C","N") || compare(a1,a2, "CA","N") || compare(a1,a2, "CB","N") )
		repSigma = pow2(2.75);
	
	else if( compare(a1,a2, "O","O") )
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "O","N") )
		repSigma = pow2(2.65);
	
	else if( compare(a1,a2, "N","N") )
		repSigma = pow2(2.6);
	
	if( compare(a1,a2, "CB","O") )
		repSigma = pow2(2.8);
break;		

// C 1.45
// O 1.3
// N 1.35
case 114:
	if     ( compare(a1,a2, "C","C") || compare(a1,a2, "C","CA") || compare(a1,a2, "C","CB") || compare(a1,a2, "CA","CA") || compare(a1,a2, "CA","CB") || compare(a1,a2, "CB","CB"))
		repSigma = pow2(2.9);
	
	else if( compare(a1,a2, "C","O") || compare(a1,a2, "CA","O") || compare(a1,a2, "CB","O") )
		repSigma = pow2(2.75);
	
	else if( compare(a1,a2, "C","N") || compare(a1,a2, "CA","N") || compare(a1,a2, "CB","N") )
		repSigma = pow2(2.8);
	
	else if( compare(a1,a2, "O","O") )
		repSigma = pow2(2.6);
	
	else if( compare(a1,a2, "O","N") )
		repSigma = pow2(2.65);
	
	else if( compare(a1,a2, "N","N") )
		repSigma = pow2(2.7);
	
	if( compare(a1,a2, "CB","O") )
		repSigma = pow2(2.75);
break;		

// C 1.45
// O 1.35
// N 1.35
case 115:
	if     ( compare(a1,a2, "C","C") || compare(a1,a2, "C","CA") || compare(a1,a2, "C","CB") || compare(a1,a2, "CA","CA") || compare(a1,a2, "CA","CB") || compare(a1,a2, "CB","CB"))
		repSigma = pow2(2.9);
	
	else if( compare(a1,a2, "C","O") || compare(a1,a2, "CA","O") || compare(a1,a2, "CB","O") )
		repSigma = pow2(2.8);
	
	else if( compare(a1,a2, "C","N") || compare(a1,a2, "CA","N") || compare(a1,a2, "CB","N") )
		repSigma = pow2(2.8);
	
	else if( compare(a1,a2, "O","O") )
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "O","N") )
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "N","N") )
		repSigma = pow2(2.7);
	
	if( compare(a1,a2, "CB","O") )
		repSigma = pow2(2.8);
break;		

// C 1.5
// O 1.3
// N 1.3
case 116:
	if     ( compare(a1,a2, "C","C") || compare(a1,a2, "C","CA") || compare(a1,a2, "C","CB") || compare(a1,a2, "CA","CA") || compare(a1,a2, "CA","CB") || compare(a1,a2, "CB","CB"))
		repSigma = pow2(3);
	
	else if( compare(a1,a2, "C","O") || compare(a1,a2, "CA","O") || compare(a1,a2, "CB","O") )
		repSigma = pow2(2.8);
	
	else if( compare(a1,a2, "C","N") || compare(a1,a2, "CA","N") || compare(a1,a2, "CB","N") )
		repSigma = pow2(2.8);
	
	else if( compare(a1,a2, "O","O") )
		repSigma = pow2(2.6);
	
	else if( compare(a1,a2, "O","N") )
		repSigma = pow2(2.6);
	
	else if( compare(a1,a2, "N","N") )
		repSigma = pow2(2.6);
	
	if( compare(a1,a2, "CB","O") )
		repSigma = pow2(2.8);
break;		

// C 1.5
// O 1.35
// N 1.3
case 117:
	if     ( compare(a1,a2, "C","C") || compare(a1,a2, "C","CA") || compare(a1,a2, "C","CB") || compare(a1,a2, "CA","CA") || compare(a1,a2, "CA","CB") || compare(a1,a2, "CB","CB"))
		repSigma = pow2(3);
	
	else if( compare(a1,a2, "C","O") || compare(a1,a2, "CA","O") || compare(a1,a2, "CB","O") )
		repSigma = pow2(2.85);
	
	else if( compare(a1,a2, "C","N") || compare(a1,a2, "CA","N") || compare(a1,a2, "CB","N") )
		repSigma = pow2(2.8);
	
	else if( compare(a1,a2, "O","O") )
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "O","N") )
		repSigma = pow2(2.65);
	
	else if( compare(a1,a2, "N","N") )
		repSigma = pow2(2.6);
	
	if( compare(a1,a2, "CB","O") )
		repSigma = pow2(2.85);
break;		

// C 1.5
// O 1.3
// N 1.35
case 118:
	if     ( compare(a1,a2, "C","C") || compare(a1,a2, "C","CA") || compare(a1,a2, "C","CB") || compare(a1,a2, "CA","CA") || compare(a1,a2, "CA","CB") || compare(a1,a2, "CB","CB"))
		repSigma = pow2(3);
	
	else if( compare(a1,a2, "C","O") || compare(a1,a2, "CA","O") || compare(a1,a2, "CB","O") )
		repSigma = pow2(2.8);
	
	else if( compare(a1,a2, "C","N") || compare(a1,a2, "CA","N") || compare(a1,a2, "CB","N") )
		repSigma = pow2(2.85);
	
	else if( compare(a1,a2, "O","O") )
		repSigma = pow2(2.6);
	
	else if( compare(a1,a2, "O","N") )
		repSigma = pow2(2.65);
	
	else if( compare(a1,a2, "N","N") )
		repSigma = pow2(2.7);
	
	if( compare(a1,a2, "CB","O") )
		repSigma = pow2(2.8);
break;		

// C 1.5
// O 1.35
// N 1.35
case 119:
	if     ( compare(a1,a2, "C","C") || compare(a1,a2, "C","CA") || compare(a1,a2, "C","CB") || compare(a1,a2, "CA","CA") || compare(a1,a2, "CA","CB") || compare(a1,a2, "CB","CB"))
		repSigma = pow2(3);
	
	else if( compare(a1,a2, "C","O") || compare(a1,a2, "CA","O") || compare(a1,a2, "CB","O") )
		repSigma = pow2(2.85);
	
	else if( compare(a1,a2, "C","N") || compare(a1,a2, "CA","N") || compare(a1,a2, "CB","N") )
		repSigma = pow2(2.85);
	
	else if( compare(a1,a2, "O","O") )
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "O","N") )
		repSigma = pow2(2.7);
	
	else if( compare(a1,a2, "N","N") )
		repSigma = pow2(2.7);
	
	if( compare(a1,a2, "CB","O") )
		repSigma = pow2(2.85);
break;		

