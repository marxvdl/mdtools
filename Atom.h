/*
 * Atom.h
 *
 *  Created on: Mar 19, 2010
 *      Author: marx
 */

#ifndef ATOM_H_
#define ATOM_H_

#include <sstream>
#include <string>
#include <cmath>

#include "misc.h"

using namespace std;

class Atom {
private:
	static const real RADTODEG = 57.295779579;
	real distance(real ax,real ay,real az, real bx,real by,real bz){
		real dx = ax - bx;
		real dy = ay - by;
		real dz = az - bz;
		
		return sqrt(pow2(dx) + pow2(dy) + pow2(dz));
	}
	
	real distanceSquare(real ax,real ay,real az, real bx,real by,real bz){
		real dx = ax - bx;
		real dy = ay - by;
		real dz = az - bz;
		
		return pow2(dx) + pow2(dy) + pow2(dz);
	}
public:
	int index;
	string type;
	string residueType;
	int residueIndex;
	real x, y, z;

	real distanceToCenter;
	real burial;
	void calculateBurial(int levels, real cx, real cy, real cz);
	void calculateRealBurial(real cx, real cy, real cz);

	string toString(){
		ostringstream r;
	    r << index << ": " << type << ", " << residueType << " (" << residueIndex << "), pos: " << x << ", " << y << ", " << z << "; Distance = " << distanceToCenter << "; Burial = " << burial;
		return r.str();
	}
	
	
	/**
	 * Returns the distance between this atom and b.
	 */ 
	real distanceTo(Atom &b){
		return distance(x,y,z, b.x,b.y,b.z);
	}
	
	real distanceSquareTo(Atom &b){
		return distanceSquare(x,y,z, b.x,b.y,b.z);
	}

	/**
	 * Returns the angle formed between vectors b->[this atom] and b->c.
	 */
	real angleTo(const Atom& b, const Atom& c){
		real m_ba = distance(  x,  y,  z,   b.x,b.y,b.z);
		real m_bc = distance(b.x,b.y,b.z,   c.x,c.y,c.z);
		
		real v_ba[] = {b.x -  x,  b.y -  y, b.z -  z};
		real v_bc[] = {b.x - c.x, b.y -c.y, b.z - c.z};
			
		real dotp = v_ba[0]*v_bc[0] + v_ba[1]*v_bc[1] + v_ba[2]*v_bc[2];
			
		return acos ( dotp / (m_ba * m_bc) );	
	}
	
	
	/**
	 * Returns the dihedral angle formed by this atom and atoms b,c,d.
	 * Adapted from:
	 * http://www.bio.net/bionet/mm/xtal-log/1997-February/002899.html
	 */
	real dihedralTo(const Atom& b, const Atom& c, const Atom& d){
		real m_ab = distance(  x,  y,  z,   b.x,b.y,b.z);
		real m_ac = distance(  x,  y,  z,   c.x,c.y,c.z);
		real m_ad = distance(  x,  y,  z,   d.x,d.y,d.z);
		real m_bc = distance(b.x,b.y,b.z,   c.x,c.y,c.z);
		real m_bd = distance(b.x,b.y,b.z,   d.x,d.y,d.z);
		real m_cd = distance(c.x,c.y,c.z,   d.x,d.y,d.z);
		
		real p = pow2(m_ab) * (  pow2(m_bc) + pow2(m_cd) - pow2(m_bd) ) +
		         pow2(m_bc) * ( -pow2(m_bc) + pow2(m_cd) + pow2(m_bd) ) +
			 pow2(m_ac) * (  pow2(m_bc) - pow2(m_cd) + pow2(m_bd) ) -
			 2 * pow2(m_bc) * pow2(m_ad);
				
		real q = (m_ab + m_bc + m_ac) * ( m_ab + m_bc - m_ac) *
		         (m_ab - m_bc + m_ac) * (-m_ab + m_bc + m_ac ) *
		         (m_bc + m_cd + m_bd) * ( m_bc + m_cd - m_bd ) *
		         (m_bc - m_cd + m_bd) * (-m_bc + m_cd + m_bd );
			 
		real z = p / sqrt(q);
		
		if(abs(z-1.0) < 0.000001)
			return 0;
		else 
			return acos(z);
	}

};


#endif /* ATOM_H_ */
