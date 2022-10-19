#pragma once


#include "stdio.h"
class C_Matrix
{
public:
	C_Matrix();
	~C_Matrix();
public:

	/*向量赋值*/
	int CVect3(double p[], double a, double b, double c);

	/*矩阵求逆*/
	int MatrixInverse(double a[], int n);

	/*矩阵相乘 a是m行n列 b是n行k列  c为m行k列*/
	void MatrixMult(double a[], double b[], double c[], int m, int n, int k);

	/*矩阵转置a[m][n]*/
	void MatrixTranslate(double a[], double b[], int m, int n);

	/*矩阵相加*/
	void MatrixAdd(double a[], double b[], double c[], int n, int m);

	/*矩阵相减*/
	void MatrixSub(double a[], double b[], double c[], int n, int m);

	/*乔里斯基分解*/
	void Cholesky(double a[], double u[], int n);

	/*I*单位矩阵*/
	void uintMatrix(double E[], int n);

	/*Tr*/////对角线元素相加
	double Trace(double a[], int n);

	/*0 Mat*///赋值为0矩阵
	void zero(double a[], int n, int m);

	/*Cross production of two vectors*/
	void crossprod(double a[], double b[], double c[]);

	/*skew*///反对称矩阵
	int skew(double a[], double *skew_a);

	int MatrixMult_quat(double q1[], double q2[], double *q3);

	int MatrixMult_3_6_to_33(double v3[], double v6[], double *matrix33);

	//取整函数
	double round(double val);
	//矩阵扩大倍数X=ff*X_1;
	double  MatrixEnlarge(double ff, double *X_1, double *X, int nx);

	//矩阵对矩阵赋值b[i]=a[i];
	double MatrixAssign(double *a, double *b, int nx);

	///中位数计算
	double MatrixMedian(double pp[], int n, bool IsAbs);
	
};

