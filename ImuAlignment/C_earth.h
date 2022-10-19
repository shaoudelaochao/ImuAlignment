#pragma once
#ifndef _EARTHPAR_H
#define _EARTHPAR_H
////////WGS_84坐标系的椭球参数
class C_earth
{
public:
	double Re, ff, Rp, e, wie;
	double GM, ge, gp;
	
	double beta1, beta2, beta3, beta4, beta5, beta6;
	double g0, deg, hour, degh, degh2, ug, msh2;
	C_earth();
	

};

#endif