#pragma once


#include "stdio.h"
class C_Matrix
{
public:
	C_Matrix();
	~C_Matrix();
public:

	/*������ֵ*/
	int CVect3(double p[], double a, double b, double c);

	/*��������*/
	int MatrixInverse(double a[], int n);

	/*������� a��m��n�� b��n��k��  cΪm��k��*/
	void MatrixMult(double a[], double b[], double c[], int m, int n, int k);

	/*����ת��a[m][n]*/
	void MatrixTranslate(double a[], double b[], int m, int n);

	/*�������*/
	void MatrixAdd(double a[], double b[], double c[], int n, int m);

	/*�������*/
	void MatrixSub(double a[], double b[], double c[], int n, int m);

	/*����˹���ֽ�*/
	void Cholesky(double a[], double u[], int n);

	/*I*��λ����*/
	void uintMatrix(double E[], int n);

	/*Tr*/////�Խ���Ԫ�����
	double Trace(double a[], int n);

	/*0 Mat*///��ֵΪ0����
	void zero(double a[], int n, int m);

	/*Cross production of two vectors*/
	void crossprod(double a[], double b[], double c[]);

	/*skew*///���Գƾ���
	int skew(double a[], double *skew_a);

	int MatrixMult_quat(double q1[], double q2[], double *q3);

	int MatrixMult_3_6_to_33(double v3[], double v6[], double *matrix33);

	//ȡ������
	double round(double val);
	//����������X=ff*X_1;
	double  MatrixEnlarge(double ff, double *X_1, double *X, int nx);

	//����Ծ���ֵb[i]=a[i];
	double MatrixAssign(double *a, double *b, int nx);

	///��λ������
	double MatrixMedian(double pp[], int n, bool IsAbs);
	
};

