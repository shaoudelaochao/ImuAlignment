#include"stdafx.h"

class C_SINSTools
{
public:
	C_SINSTools();
	~C_SINSTools();
public:

	void GravityUpdate(double *gn);////��������

	double Cnb2Att(double Cnb[], double *Att);//��̬����תΪŷ����
	void Att2Cnb(  double *att, double *Cnb);

	void Rov(double *C,double *faiw);///////��תʸ��ת��Ϊ�仯����

	int Diagonal(int nx, double *Pk, double f, ...);/////�Խ���Ԫ������

	void StaticUpdate(  double *Dw1,   double *Dw2,   double *Dw3,   double *Dw4,
		   double *Df1,    double *Df2,    double *Df3,   double *Df4,   double Dt,
		double *Cnb, double *Vn, double *P);//////////��̬���ٶȣ�λ�ø����㷨
	 
	




};