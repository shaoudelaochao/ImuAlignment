#include "stdafx.h"


C_KalmanFilter::C_KalmanFilter()
{
}


C_KalmanFilter::~C_KalmanFilter()
{
}


C_KalmanFilter::C_KalmanFilter(int nx0, int nz0)
{
	nx = nx0;
	nz = nz0;

}


//
//int C_KalmanFilter::SetPk(double *Pk, double f, ...)
//{
//	int i;
//	va_list vl;
//	va_start(vl, f);
//	memset(Pk, 0, sizeof(double)* nx * nx);
//	for (i = 0; i < nx; i++)
//	{
//		Pk[i * nx + i] = f*f;
//		f = va_arg(vl, double);
//	}//PK主对角线元素赋值
//	va_end(vl);
//	return 1;
//}
//
//int C_KalmanFilter::SetRk(double *Rk, double f, ...)
//{
//	int i;
//	va_list vl;
//	va_start(vl, f);
//	memset(Rk, 0, sizeof(double)* nz * nz);
//	for (i = 0; i < nz; i++)
//	{
//		Rk[i * nz + i] = f * f;
//		f = va_arg(vl, double);
//	}//PK主对角线元素赋值
//	va_end(vl);
//	return 1;
//}

void C_KalmanFilter::SetHk(double *Hk)
{
	Hk[7*0+3] = Hk[7*1+4] = 1.0;
}



void C_KalmanFilter::SetFt(  double *Pos0,   double g, double *Ft)
{
	C_earth cglv;
	double wN = cglv.wie*cos(Pos0[0]);
	double wU = cglv.wie*sin(Pos0[0]);

	///F状态转移矩阵
	Ft[7*0+1] = wU;      Ft[7*0+2] = -wN;
	Ft[7*1+0] = -wU;     Ft[7*1+5] = -1;
	Ft[7*2+0] = wN;     Ft[7*2+6] = -1;
	Ft[7*3+1] = -g;     Ft[7*4+0] = g;


}

void C_KalmanFilter::SetPkRkQt(double *Pk, double *Rk, double *Qt)
{


	C_Matrix matrix;
	C_SINSTools tool;
	C_earth cglv;
	//////////////////////////////////////////////////////
	///////   设置Xk Pk Qt  Rk的值  ///////////////////////
	/////// Xk=[失准角（ENU） 速度误差（EN）  陀螺随机常值漂移（NU）]
	double eb = 0.1*cglv.degh;////陀螺常值漂移（deg/h）
	double db = 50 * cglv.ug;/////角速度常值漂移（ug）

	/////std
	double std_att[3] = {    0.1*cglv.deg, 0.1*cglv.deg,0.2*cglv.deg };
	double std_vn[3] = { 0.02,0.02,0.02};
	double std_eb[3] = { 0.1*cglv.degh, 0.1*cglv.degh,0.1*cglv.degh };


	/////角度随机游走(deg/sqrt(h))/////////////////
	double web2 = 0.0003*cglv.degh2;
	double web1 = 0.0005*cglv.degh2;
	double web3 = 0.0003*cglv.degh2;

	/////速度随机游走(m/s/sqrt(h))////////////////
	double wdb1 = 4.5e-4*cglv.msh2;
	double wdb2 = 7.2e-4*cglv.msh2;
	double wdb3 = 1.5e-3 *cglv.msh2;

	/////状态协方差矩阵（先验协方差）
	tool.Diagonal(7, Pk, std_att[0], std_att[1],  std_att[2],
		std_vn[0], std_vn[1], std_eb[0], std_eb[1]);////////状态向量对应的误差平方


	/////量测噪声设置
	tool.Diagonal(2, Rk, 0.01, 0.01);
	///////Qt为功率谱密度矩阵，一般根据加速度计和陀螺仪的白噪声，零偏驱动噪声，
	///////比例因子误差误差驱动白噪声水平给定，一般由传感器标定参数给出
	tool.Diagonal(6, Qt, web1, web2, web3, wdb1, wdb2, wdb3);
	/////////////////////////Xk PK Qt Rk设置完毕///////////////////

}



//////以e为底的矩阵函数///////////////
///近似处理eA=I+A+A_2/2!+A_3/3!+A_4/4!+A_5/5!+A_5/6!+.......
void C_KalmanFilter::expmFk(  double *Ft,   double Dt, double *Fk)
{
	C_Matrix matrix;
	double Ft_1[7 * 7] = { 0 }, Ft_2[7 * 7] = { 0 }, Ft_3[7 * 7] = { 0 };
	double	Ft_4[7 * 7] = { 0 }, Ft_5[7 * 7] = { 0 }, Ft_6[7 * 7] = { 0 };
	double	Ft_7[7 * 7] = { 0 }, Ft_8[7 * 7] = { 0 }, Ft_9[7 * 7] = { 0 }, I_7[7 * 7] = { 0 };

	for (int i = 0; i < 7 * 7; i++)  Ft_1[i] = Ft[i] * Dt;
	
	matrix.MatrixMult(Ft_1, Ft_1, Ft_2, 7, 7, 7);
	matrix.MatrixMult(Ft_2, Ft_1, Ft_3, 7, 7, 7);
	matrix.MatrixMult(Ft_3, Ft_1, Ft_4, 7, 7, 7);
	matrix.MatrixMult(Ft_4, Ft_1, Ft_5, 7, 7, 7);
	matrix.MatrixMult(Ft_6, Ft_1, Ft_6, 7, 7, 7);
	matrix.MatrixMult(Ft_6, Ft_1, Ft_7, 7, 7, 7);
	matrix.MatrixMult(Ft_7, Ft_1, Ft_8, 7, 7, 7);
	matrix.MatrixMult(Ft_8, Ft_1, Ft_9, 7, 7, 7);
	matrix.uintMatrix(I_7, 7);
	for (int i = 0; i < 7 * 7; i++)
	{
		Fk[i] = I_7[i] + Ft_1[i] + (Ft_2[i] * 1.0) / 2 + (Ft_3[i] * 1.0) / (2 * 3) + (Ft_4[i] * 1.0) / (2 * 3 * 4)
			+ (Ft_5[i] * 1.0) / (2 * 3 * 4 * 5) + (Ft_6[i] * 1.0) / (2 * 3 * 4 * 5 * 6) + (Ft_7[i] * 1.0) / (2 * 3 * 4 * 5 * 6 * 7)
			+ (Ft_8[i] * 1.0) / (2 * 3 * 4 * 5 * 6 * 7 * 8) + (Ft_9[i] * 1.0) / (2 * 3 * 4 * 5 * 6 * 7 * 8 * 9);
	}
}

void C_KalmanFilter::setQk(  double *Cnb,   double *Qt,   double *Fk,
	       double Dt, double *Gt,double *Qk )
{
	C_Matrix matrix;
	//////////////////////////////////////////////////////////////
	////////// Qk = 0.5*Dt*(Fk*Gt*Qt*Gt'*Fk' + Gt*Qt*Gt')/////////
	///G噪声分配矩阵
	Gt[0] = -Cnb[0];    Gt[1] = -Cnb[1];   Gt[2] = -Cnb[2];
	Gt[6] = -Cnb[3];    Gt[7] = -Cnb[4];   Gt[8] = -Cnb[5];
	Gt[12] = -Cnb[6];   Gt[13] = -Cnb[7];  Gt[14] = -Cnb[8];
	Gt[21] = Cnb[0];    Gt[22] = Cnb[1];   Gt[23] = Cnb[2];
	Gt[27] = Cnb[3];    Gt[28] = Cnb[4];   Gt[29] = Cnb[5];

	double FG[7 * 6] = { 0 }, FGQ[7 * 6] = { 0 }, Gt_T[6 * 7] = { 0 };
	double FGQGt_T[7 * 7] = {0.0}, Fk_T[7 * 7] = { 0.0 },FGQGt_TFk_T[7 * 7] = { 0.0 };
	matrix.MatrixMult(Fk, Gt, FG, 7, 7, 6);//Fk*Gt
	matrix.MatrixMult(FG, Qt, FGQ, 7, 6, 6);///Fk*Gt*Qt
	matrix.MatrixTranslate(Gt, Gt_T, 7, 6);///
	matrix.MatrixMult(FGQ, Gt_T, FGQGt_T, 7, 6, 7);//Fk*Gt*Qt*Gt'
	matrix.MatrixTranslate(Fk, Fk_T, 7, 7);//Fk'
	matrix.MatrixMult(FGQGt_T, Fk_T, FGQGt_TFk_T, 7, 7, 7);//Fk*Gt*Qt*Gt'*Fk'

	////////////////////////////////////////////////////////////////
	////////////////////Gt*Qt*Gt'//////////////////////////////////
	double GQ[7 * 6] = { 0 }, GQGt_T[7 * 7] = { 0.0 };
	matrix.MatrixMult(Gt, Qt, GQ, 7, 6, 6);
	matrix.MatrixMult(GQ, Gt_T, GQGt_T, 7, 6, 7);

	for (int i = 0; i < 7 * 7; i++)  Qk[i] = 0.5*Dt*(FGQGt_TFk_T[i] + GQGt_T[i]);

}





int C_KalmanFilter::Timeupdate(double *Ft, double *Xk, double *Pk, double *Qt,int nx)
{
	C_Matrix m_matrix;
	int i;
	/************************************************************************/
	//double *Xk_tu = new double[nx * 1];
	double *Xk_tu = (double *)malloc(nx * 1 * sizeof(double));
	memset(Xk_tu, 0, sizeof(double)*nx * 1);
	m_matrix.MatrixMult(Ft, Xk, Xk_tu, nx, nx, 1);
	for (i = 0; i < nx; i++)
		Xk[i] = Xk_tu[i];////////Xk=Ft*Xk;
	/***************************Xk预测状态向量更新完毕*************************/

	/************************************************************************/
	//double *FP = new double[nx * nx];
	double *FP = (double *)malloc(nx*nx*sizeof(double));
	memset(FP, 0, sizeof(double)*nx*nx);

	m_matrix.MatrixMult(Ft, Pk, FP, nx, nx, nx);//return Ft*Pk


	//double *Ft_T = new double[nx * nx];
	double *Ft_T = (double *)malloc(nx*nx*sizeof(double));//Ft Translate
	memset(Ft_T, 0, sizeof(double)*nx*nx);
	m_matrix.MatrixTranslate(Ft, Ft_T, nx, nx);


	//double *FPF = new double[nx * nx]; 
	double *FPF = (double *)malloc(nx*nx*sizeof(double));
	memset(FPF, 0, sizeof(double)*nx*nx);
	m_matrix.MatrixMult(FP, Ft_T, FPF, nx, nx, nx);//return Ft*Pk*(~Ft)


	//double *Pk_update = new double[nx * nx];
	double *Pk_update = (double *)malloc(nx*nx*sizeof(double));
	memset(Pk_update, 0, sizeof(double)*nx*nx);
	m_matrix.MatrixAdd(FPF, Qt, Pk_update, nx, nx);/////////Pk=Ft*Pk*(~Ft)+Qks

	for (i = 0; i < nx*nx; i++)
		Pk[i] = Pk_update[i];
	/***************************Pk预测状态协方差矩阵更新完毕*****************************/

	double sigma_p[3] = {  sqrt(Pk[0 * 7 + 0])*kRAD2DEG,sqrt(Pk[1 * 7 +1])*kRAD2DEG, sqrt(Pk[2 * 7 + 2])*kRAD2DEG};
	for (i = 0; i < 3; i++)
	{
		if (sigma_p[i] < 0.0001)   Pk[i * 7 + i] = pow2( 0.0001*kDEG2RAD);
	}
	double sigma_vn[3] = { sqrt(Pk[3 * 7 + 3]),sqrt(Pk[4 * 7 + 4]), 0 };
	for (i = 0; i < 2; i++)
	{
		if (sigma_vn[i] < 0.001)   Pk[(i+3)* 7 + i+3] = pow2(0.001);
	}



	free(Xk_tu), free(FP), free(Ft_T), free(FPF), free(Pk_update);
	return 1;

}

int C_KalmanFilter::MeasUpdate(double *Pk, double *Rk, double *Hk, double *Xk, double *Zk,int nx,int nz)
{
	C_Matrix m_matrix;
	int i;


	/**********************************local parameters of K***************************************/
	double *K = (double *)malloc(nx*nz*sizeof(double));//K
	memset(K, 0, sizeof(double)*nx*nz);

	double *K_T = (double *)malloc(nz*nx*sizeof(double));//K_T
	memset(K_T, 0, sizeof(double)*nz*nx);

	double *H_T = (double *)malloc(nx*nz*sizeof(double));//H_T
	memset(H_T, 0, sizeof(double)*nx*nz);

	double *PH_T = (double *)malloc(nx*nz*sizeof(double)); //PH_T
	memset(PH_T, 0, sizeof(double)*nx*nz);

	double *HP = (double *)malloc(nz*nx*sizeof(double)); //HP
	memset(HP, 0, sizeof(double)*nz*nx);

	double *HPH_T = (double *)malloc(nz*nz*sizeof(double));//HPH_T
	memset(HPH_T, 0, sizeof(double)*nz*nz);

	double *HPH_TpR = (double *)malloc(nz*nz*sizeof(double));//HPH_T+R
	memset(HPH_TpR, 0, sizeof(double)*nz*nz);

	/*-------------K-------------*/
	m_matrix.MatrixTranslate(Hk, H_T, nz, nx);        //return H_T
	m_matrix.MatrixMult(Pk, H_T, PH_T, nx, nx, nz);   //return P*(~H)
	m_matrix.MatrixMult(Hk, Pk, HP, nz, nx, nx);      //return H*P
	m_matrix.MatrixMult(HP, H_T, HPH_T, nz, nx, nz); //return H*P*(~H)
	m_matrix.MatrixAdd(HPH_T, Rk, HPH_TpR, nz, nz);   //return H*P*(~H)+R
	m_matrix.MatrixInverse(HPH_TpR, nz);              //return 1/(H*P*(~H)+R)
	m_matrix.MatrixMult(PH_T, HPH_TpR, K, nx, nz, nz);//return K=P*(~H)/(H*P*(~H)+R)
	m_matrix.MatrixTranslate(K, K_T, nx, nz);        //reurn (~K)

	/******************增益滤波矩阵更新完毕****************************************************/

	/*------------------local parameters of P------------------------------------------------*/
	double *P_update = (double *)malloc(nx*nx*sizeof(double));//P_update
	memset(P_update, 0, sizeof(double)*nx*nx);

	double *I = (double *)malloc(nx*nx*sizeof(double));//I
	memset(I, 0, sizeof(double)*nx*nx);
	m_matrix.uintMatrix(I, nx);

	double *KH = (double *)malloc(nx*nx*sizeof(double));//K*H
	memset(KH, 0, sizeof(double)*nx*nx);

	double *ImKH = (double *)malloc(nx*nx*sizeof(double));//I-K*H
	memset(ImKH, 0, sizeof(double)*nx*nx);

	double *ImKH_T = (double *)malloc(nx*nx*sizeof(double));//(I-K*H)_T
	memset(ImKH_T, 0, sizeof(double)*nx*nx);

	double *ImKHP = (double *)malloc(nx*nx*sizeof(double));//(I-K*H)*P
	memset(ImKHP, 0, sizeof(double)*nx*nx);

	double *ImKHPImKH_T = (double *)malloc(nx*nx*sizeof(double));//(I-K*H)*P*(I-K*H)_T
	memset(ImKHPImKH_T, 0, sizeof(double)*nx*nx);

	double *KR = (double *)malloc(nx*nz*sizeof(double));//K*R
	memset(KR, 0, sizeof(double)*nx*nz);

	double *KRK_T = (double *)malloc(nx*nx*sizeof(double));//K*R*K_T
	memset(KRK_T, 0, sizeof(double)*nx*nx);

	/*-------------P_update-------------*/
	m_matrix.MatrixMult(K, Hk, KH, nx, nz, nx);//return KH
	m_matrix.MatrixSub(I, KH, ImKH, nx, nx);//return I-KH
	m_matrix.MatrixTranslate(ImKH, ImKH_T, nx, nx);//return (I-K*H)_T
	m_matrix.MatrixMult(ImKH, Pk, ImKHP, nx, nx, nx);//return (I - K*H)*P
	m_matrix.MatrixMult(ImKHP, ImKH_T, ImKHPImKH_T, nx, nx, nx);//return (I-K*H)*P*(I-K*H)_T
	m_matrix.MatrixMult(K, Rk, KR, nx, nz, nz);//return K*R
	m_matrix.MatrixMult(KR, K_T, KRK_T, nx, nz, nx);//return K*R*K_T
	m_matrix.MatrixAdd(ImKHPImKH_T, KRK_T, P_update, nx, nx);


	memset(Pk, 0, sizeof(double)*nx*nx);
	for (i = 0; i < nx*nx; i++)
		Pk[i] = P_update[i];
	/********************************新的状态协方差矩阵***************************************/


	/*-----------------------------------------local parameters of X-------------------------*/
	double *X_update = (double *)malloc(nx * 1 * sizeof(double));//X_update
	memset(X_update, 0, sizeof(double)*nx * 1);

	double *HX = (double *)malloc(nz * 1 * sizeof(double));//H*X
	memset(HX, 0, sizeof(double)*nz * 1);

	double *ZmHX = (double *)malloc(nz * 1 * sizeof(double));//Z-HX
	memset(ZmHX, 0, sizeof(double)*nz * 1);

	double *KZmHX = (double *)malloc(nx * 1 * sizeof(double));//K*(Z-HX)
	memset(KZmHX, 0, sizeof(double)*nx * 1);

	/*-------------X-------------*/
	m_matrix.MatrixMult(Hk, Xk, HX, nz, nx, 1);//return H*X
	m_matrix.MatrixSub(Zk, HX, ZmHX, nz, 1);//return Z-H*X///为新息向量

	m_matrix.MatrixMult(K, ZmHX, KZmHX, nx, nz, 1);//return K*(Z-H*X)
	m_matrix.MatrixAdd(Xk, KZmHX, X_update, nx, 1);//return X+K*(Z-H*X)
	memset(Xk, 0, sizeof(double)*nx * 1);
	for (i = 0; i < nx; i++)
		Xk[i] = X_update[i];



	/*------------------------------------新的状态估计完成------------------------------------------*/
	free(K), free(K_T), free(HP); free(H_T), free(PH_T), free(HPH_T), free(HPH_TpR);
	free(P_update), free(I), free(KH), free(ImKH), free(ImKH_T), free(ImKHP), free(ImKHPImKH_T), free(KR), free(KRK_T);
	free(KZmHX), free(ZmHX), free(HX), free(X_update);

	return 1;
}


int C_KalmanFilter::FilterBack(double *Xk, double *Vn, double *Cnb,double *eb)
 {
	/////理想的导航坐标系N系到载体坐标系B系的姿态矩阵Cnb
	////与计算坐标系之间有偏差Cn'b。
	////他们之间的偏差为n'系和n系之间的偏差
	////Cn'b=Cn'n*Cnb
	/////捷联惯导算法与组合导航技术，严恭敏，78p
	C_Matrix matrix;
	int m= 1;
	double n =0.0;
	double  XK_[3 * 1] = { 0 }, skew_XK_[3 * 3] = {0.0};
	double  I3[3 * 3] = { 0 }, Cnn[3 * 3] = { 0.0 };
	double  Vn2[3 * 1] = { 0 }, Cnb_2[3 * 3] = { 0.0 };
	
	/////------失准角反馈----------------------------
	XK_[0] = Xk[0];XK_[1] = Xk[1];XK_[2] = Xk[2];
	matrix.skew(XK_, skew_XK_);
	matrix.uintMatrix(I3, 3);
	matrix.MatrixAdd(I3,skew_XK_,Cnn,3,3);//eye(3) + Askew(m*Xk(1:3,1));
	matrix.MatrixMult(Cnn,Cnb,Cnb_2,3,3,3);
	for (int i = 0; i < 9; i++)
	{
		Cnb[i] = Cnb_2[i];
	}


	//-------------速度误差反馈--------------------------------------
	Vn2[0] = Vn[0] - m*Xk[3];Vn2[1] = Vn[1] - m*Xk[4];
	Vn[0] = Vn2[0]; Vn[1] = Vn2[1];

	//-------------零偏反馈------------------------------------------
    eb[0] = eb[0]- Xk[5]; eb[1] = eb[1]- Xk[6];

	//-----------误差状态初始化-------------------------------------
	for (int i = 0; i < 5; i++)  Xk[i] = 0.0;

	return 1;
}