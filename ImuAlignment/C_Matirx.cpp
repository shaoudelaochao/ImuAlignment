/*************************************************
Copyright (C),2017.12.20
File name:		C_matrix.cpp
Author:			CGF	   Version:           Date: 2017.12.20
Description:		Matrix tools
Others:			// 其它内容的说明
Function List:	// 主要函数列表，每条记录应包括函数名及功能简要说明

1 CVect3
2 MatrixInverse
3 MatrixMult
4 MatrixTranslate
5 MatrixAdd
6 MatrixSub
7 Cholesky(double a[], double u[], int n);
8 uintMatrix
9 Trace
10 zero
11 crossprod
12 skew
13 MatrixMult_quat
14 int MatrixMult_3_6_to_33
15 doubel rond
History:      // 修改历史记录列表，每条修改记录应包括修改日期、修改
// 者及修改内容简述
1. Date:
Author:
Modification:
2. ...
*************************************************/

#include "stdafx.h"



C_Matrix::C_Matrix()
{
}


C_Matrix::~C_Matrix()
{
}

int C_Matrix::CVect3(double *p, double a, double b, double c)
{
	memset(p, 0, sizeof(double)* 3 * 1);
	p[0] = a;
	p[1] = b;
	p[2] = c;
	return 1;
}

int C_Matrix::MatrixInverse(double a[], int n)
{
	int		i, j, k, l, u, v;
	int		is[24], js[24];
	double	d, p;

	for (k = 0; k <= n - 1; k++)
	{
		d = 0.0;
		for (i = k; i <= n - 1; i++)
		{
			for (j = k; j <= n - 1; j++)
			{
				l = i*n + j;
				p = fabs(a[l]);
				if (p>d)
				{
					d = p;
					is[k] = i;
					js[k] = j;
				}
			}
		}
		if (d<1.0e-9 && d>-1.0e-9)
		{
			return(-1);
		}
		if (is[k] != k)
		{
			for (j = 0; j <= n - 1; j++)
			{
				u = k*n + j;
				v = is[k] * n + j;
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
		}
		if (js[k] != k)
		{
			for (i = 0; i <= n - 1; i++)
			{
				u = i*n + k;
				v = i*n + js[k];
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
		}
		l = k*n + k;
		a[l] = 1.0 / a[l];
		for (j = 0; j <= n - 1; j++)
		{
			if (j != k)
			{
				u = k*n + j;
				a[u] = a[u] * a[l];
			}
		}
		for (i = 0; i <= n - 1; i++)
		{
			if (i != k)
			{
				for (j = 0; j <= n - 1; j++)
				{
					if (j != k)
					{
						u = i*n + j;
						a[u] = a[u] - a[i*n + k] * a[k*n + j];
					}
				}
			}
		}
		for (i = 0; i <= n - 1; i++)
		{
			if (i != k)
			{
				u = i*n + k;
				a[u] = -a[u] * a[l];
			}
		}
	}
	for (k = n - 1; k >= 0; k--)
	{
		if (js[k] != k)
		{
			for (j = 0; j <= n - 1; j++)
			{
				u = k*n + j;
				v = js[k] * n + j;
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
		}
		if (is[k] != k)
		{
			for (i = 0; i <= n - 1; i++)
			{
				u = i*n + k;
				v = i*n + is[k];
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
		}
	}
	return(0);
}

void C_Matrix::MatrixMult(double a[], double b[], double c[], int m, int n, int k)
{
	int i, j, l;
	for (i = 0; i<m; i++)
	for (j = 0; j<k; j++)
	{
		c[i*k + j] = 0.0;
		for (l = 0; l<n; l++)
		{
			c[i*k + j] = c[i*k + j] + a[i*n + l] * b[l*k + j];
		}
	}
	return;
}

void C_Matrix::MatrixTranslate(double a[], double b[], int m, int n)
{
	int i, j;
	for (i = 0; i<m; i++)
	for (j = 0; j<n; j++)
		b[j*m + i] = a[i*n + j];
	return;
}

void C_Matrix::MatrixAdd(double a[], double b[], double c[], int n, int m)
{
	int i = 0, j = 0;
	for (i = 0; i<n; i++)
	{
		for (j = 0; j<m; j++)
		{
			c[i*m + j] = a[i*m + j] + b[i*m + j];
		}
	}
	return;
}

void C_Matrix::MatrixSub(double a[], double b[], double c[], int n, int m)
{
	int i = 0, j = 0;
	for (i = 0; i<n; i++)
	{
		for (j = 0; j<m; j++)
		{
			c[i*m + j] = a[i*m + j] - b[i*m + j];
		}
	}
	return;
}
///得到开弓号的上三角矩阵
void C_Matrix::Cholesky(double a[], double u[], int n)
{
	int i = 0, j = 0, k = 0, k2 = 0, f = 0;
	double Uki = 0, Uij = 0;
	u[0] = sqrt(a[0]);
	for (f = 0; f<n*n; f++)
	{
		u[f] = 0;
	}
	for (i = 0; i<n - 1; i++)
	{
		{
			for (j = i + 1; j<n; j++)
			{
				Uki = 0;
				Uij = 0;
				for (k = 0; k<i; k++)
				{
					Uki = Uki + u[k*n + i] * u[k*n + i];
					Uij = Uij + u[k*n + i] * u[k*n + j];
				}
				u[i*n + i] = sqrt(a[i*n + i] - Uki);
				u[i*n + j] = (a[i*n + j] - Uij) / u[i*n + i];
			}
		}
	}
	Uki = 0;
	for (k2 = 0; k2<n; k2++)
	{
		Uki = Uki + u[k2*n + i] * u[k2*n + i];
	}
	u[i*n + i] = sqrt(a[i*n + i] - Uki);
	return;
}



void C_Matrix::uintMatrix(double E[], int n)
{
	int i = 0, j = 0;
	for (i = 0; i<n; i++)
	for (j = 0; j<n; j++)
	{
		E[i*n + j] = 0.0;
		if (i == j)
			E[i*n + j] = 1.0;
	}
}
//对角线元素相加
double C_Matrix::Trace(double a[], int n)
{
	int i, j;
	double t = 0;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (i == j)
				t += *(a + i*n + j);
		}
	}
	return t;
}

void C_Matrix::zero(double a[], int m, int n)
{
	int i, j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
			a[i*n + j] = 0;
	}
}

void C_Matrix::crossprod(double a[], double b[], double c[])
{
	zero(c, 3, 1);
	c[0] = a[1] * b[2] - b[1] * a[2];
	c[1] = b[0] * a[2] - a[0] * b[2];
	c[2] = a[0] * b[1] - b[0] * a[1];
}

int C_Matrix::skew(double a[], double *skew_a)
{
	skew_a[0] = 0;    skew_a[1] = -a[2]; skew_a[2] = a[1];
	skew_a[3] = a[2]; skew_a[4] = 0;     skew_a[5] = -a[0];
	skew_a[6] = -a[1]; skew_a[7] = a[0];  skew_a[8] = 0;
	return 1;
}

int C_Matrix::MatrixMult_quat(double q1[], double q2[], double *q3)
{
	C_Matrix m_matrix;

	int i;
	//double s1, s2;
	double v1[3];
	double v2[3];
	double v1_t[3];
	double v1_tv2[1];
	double v1v2[3];
	double vect3[3];
	double vect3_f[3];
	memset(v1, 0, sizeof(double)* 3 * 1);
	memset(v2, 0, sizeof(double)* 3 * 1);
	memset(v1_t, 0, sizeof(double)* 1 * 3);
	memset(v1v2, 0, sizeof(double)* 3 * 1);
	memset(v1_tv2, 0, sizeof(double)* 1 * 1);
	memset(vect3, 0, sizeof(double)* 3 * 1);
	memset(vect3_f, 0, sizeof(double)* 3 * 1);

	for (i = 0; i < 3; i++)
	{
		v1[i] = q1[1 + i];
		v2[i] = q2[1 + i];
	}
	m_matrix.MatrixTranslate(v1, v1_t, 3, 1);
	m_matrix.MatrixMult(v1_t, v2, v1_tv2, 1, 3, 1);
	q3[0] = q1[0] * q2[0] - v1_tv2[0];

	m_matrix.crossprod(v1, v2, v1v2);
	for (i = 0; i < 3; i++)
	{
		v1[i] *= q2[0];
		v2[i] *= q1[0];
	}

	m_matrix.MatrixAdd(v1, v2, vect3, 3, 1);
	m_matrix.MatrixAdd(vect3, v1v2, vect3_f, 3, 1);
	for (i = 0; i < 3; i++)
		q3[1 + i] = vect3_f[i];
	if (q3[0] < 0)
	{
		for (i = 0; i < 4; i++)
			q3[i] = -q3[i];
	}

	//delete v1, v2, v1_t, v1_tv2, v1v2, vect3, vect3_f;

	return 1;
}

int C_Matrix::MatrixMult_3_6_to_33(double v3[], double v6[], double *matrix33)
{
	matrix33[0] = v3[0], matrix33[1] = v6[0], matrix33[2] = v6[1];
	matrix33[3] = v6[2], matrix33[4] = v3[1], matrix33[5] = v6[3];
	matrix33[6] = v6[4], matrix33[7] = v6[5], matrix33[8] = v3[2];
	/*fun(a, b) = (a1 b1 b2)
	(b3 a2 b4)
	(b5 b6 a3)*/
	return 1;
}


double C_Matrix::round(double val)
{
	return (val> 0.0) ? floor(val + 0.5) : ceil(val - 0.5);
}



//矩阵扩大倍数nx行ny列
double  C_Matrix::MatrixEnlarge(double ff, double *X_1, double *X, int nx)
{


	for (int j = 0; j < nx; j++)

	{
		X[j] = ff*X_1[j];

	}
	return 1;

}
////矩阵赋值 把a矩阵赋值给B矩阵
double C_Matrix::MatrixAssign(double *a, double *b, int nx)
{

	for (int j = 0; j < nx; j++)

	{
		b[j] = a[j];

	}
	return 1;

}


double C_Matrix::MatrixMedian(double pp[], int n, bool IsAbs)
{
	double *p = new double[n];
	if (IsAbs)
	{
		for (int i = 0; i < n; i++) p[i] = fabs(pp[i]);
	}
	else
	{
		for (int i = 0; i < n; i++) p[i] = pp[i];
	}


	int k = int(n / 2.0);
	while (k>0)
	{
		for (int j = k; j <= n - 1; j++)
		{
			double t = p[j];
			int i = j - k;

			while ((i >= 0) && (p[i] > t))
			{
				p[i + k] = p[i];
				i = i - k;

			}
			p[i + k] = t;

		}
		k = k / 2;
	}
	double mean = (n % 2 == 1) ? p[n / 2] : (p[n / 2] + p[n / 2 - 1]) / 2.0;

	delete[]p;
	return mean;

}



