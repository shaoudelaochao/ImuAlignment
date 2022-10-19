#include"stdafx.h"
////初始对准
class C_AlignCoarse
{
public:
	C_AlignCoarse();
	~C_AlignCoarse();
public:

	void InitAlignC(void);////初始化
	void AlignC( double *euler ,double *Cnb, double *eb, double *db);////初始对准函数
	
	double gn[3];
    double var1[3], var2[3], var3[3];

	double fbib[3], wbib[3], fbib_1[3], wbib_1[3]; 
	double vbr1[3], vbr2[3], vbr3[3];
	double skew_vbr1[9], skew_vbr2[9];
	double Att[3];

};