#include"stdafx.h"

class C_SINSTools
{
public:
	C_SINSTools();
	~C_SINSTools();
public:

	void GravityUpdate(double *gn);////重力更新

	double Cnb2Att(double Cnb[], double *Att);//姿态矩阵转为欧拉角
	void Att2Cnb(  double *att, double *Cnb);

	void Rov(double *C,double *faiw);///////旋转矢量转化为变化矩阵

	int Diagonal(int nx, double *Pk, double f, ...);/////对角线元素设置

	void StaticUpdate(  double *Dw1,   double *Dw2,   double *Dw3,   double *Dw4,
		   double *Df1,    double *Df2,    double *Df3,   double *Df4,   double Dt,
		double *Cnb, double *Vn, double *P);//////////姿态，速度，位置更新算法
	 
	




};