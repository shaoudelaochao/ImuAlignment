#include"stdafx.h"
#include "PSINS.h"
using namespace std;
C_SINSTools::C_SINSTools()
{

}

C_SINSTools::~C_SINSTools()
{

}

void C_SINSTools::Att2Cnb( double *att,double *Cnb)
{
	double	si = sin(att[0]), ci = cos(att[0]),
		sj = sin(att[1]), cj = cos(att[1]),
		sk = sin(att[2]), ck = cos(att[2]);
	
	Cnb[0] = cj*ck - si*sj*sk;	Cnb[1] = -ci*sk;	Cnb[2] = sj*ck + si*cj*sk;
	Cnb[3] = cj*sk + si*sj*ck;	Cnb[4] = ci*ck;	Cnb[5] = sj*sk - si*cj*ck;
	Cnb[6] = -ci*sj;				Cnb[7] = si;		Cnb[8] = ci*cj;
	
}


///参考：捷联惯导算法与组合导航讲义，严恭敏////
///姿态矩阵转化为欧拉角
////载体坐标系到当地水平坐标系的装换可以通过绕Z,X,Y三次旋转得到
////Rnb=RZ(y)RX(-p)RY(-r)
double C_SINSTools::Cnb2Att(double Cnb[], double *Att)
{   
	C_Matrix matrix;

	double C32 = Cnb[7];     double CC = sqrt(1 - C32*C32);
	double C31 = -Cnb[6]*1.0 / CC; double C33 = Cnb[8] *1.0 / CC;
	double C12 = Cnb[1]*1.0  / CC;   double C22 = Cnb[4]*1.0  / CC;

	
	matrix.CVect3(Att, atan2(C31,C33),asin(C32),atan2(C12,C22));
                      //r 横滚角（-180~180）//ｐ　俯仰角(-90~90)//ｙ　航向角(0~360)
	if (Att[2] < 0)
		Att[2] = Att[2] + 2 * PI;
	
	
	return 1;
}

////捷联惯导算法与组合导航讲义，严恭敏45~50
void  C_SINSTools::GravityUpdate(double *gn)////重力更新
{
	C_Matrix matrix;
	C_earth glv;
	
	double ge = glv.ge;
	double beta1 = glv.beta1;double beta2 = glv.beta2;
    double beta3 = glv.beta3;double beta4 = glv.beta4;
    double beta5 = glv.beta5;double beta6 = glv.beta6;
	double L = L0*glv.deg;

	 
	double gL = ge*(1 + beta1*sin(L)*sin(L) + beta2*sin(2 * L)*sin(2 * L));
	double gLh = gL - (beta3 + beta4*sin(L)*sin(L)) *H0*1e-6 + beta5*H0*H0*1e-12;

	double theta = beta6*H0*sin(2 * L);
	///考虑重力矢量投影到地理坐标系
	matrix.CVect3(gn,0,-theta,-gLh);
}



///////////////////////////////////////////////////////
///////旋转矢量转化为变化矩阵，捷联惯导算法与组合导航讲义，严恭敏 p19////////////////////////////
void C_SINSTools::Rov(double *C, double *faiw)
{ 
	C_Matrix matrix;
	double fai = sqrt(faiw[0] * faiw[0] + faiw[1] * faiw[1] + faiw[2] * faiw[2]);
	double fai_2 = fai*fai;
	int nx = 3;
    double a = 0;
	double b = 0;
	if (fai_2 < 1.e-8)
	{
		a = 1 - fai_2*(1 / 6.0 - fai_2 / 120.0);
		b = 0.5 - fai_2*(1 / 24.0 - fai_2 / 720.0);
	}


	else
	{
		a = sin(fai)*1.0 / fai;
	    b = (1 - cos(fai))*1.0 / (fai_2);
	}
		

	double *skew_faiw = (double *)malloc(nx * nx * sizeof(double));
	memset(skew_faiw , 0, sizeof(double)*nx * nx);
    double *skew_faiw_2 = (double *)malloc(nx * nx * sizeof(double));
	memset(skew_faiw_2 , 0, sizeof(double)*nx * nx);

    matrix.skew(faiw, skew_faiw);
	matrix.MatrixMult(skew_faiw,skew_faiw,skew_faiw_2,nx,nx,nx);

    double *I3 = (double *)malloc(nx * nx * sizeof(double));
	memset( I3, 0, sizeof(double)*nx * nx);

	matrix.uintMatrix(I3,3);

	for (int i = 0; i < nx*nx; i++)
	{
	
		C[i] = I3[i] + a*skew_faiw[i] + b*skew_faiw_2[i];
	
	}

	free(skew_faiw), free(skew_faiw_2), free(I3);


}





////对角线元素赋值
int C_SINSTools::Diagonal(int nx, double *Pk, double f, ...)
{
	int i;
	va_list vl;
	va_start(vl, f);
	memset(Pk, 0, sizeof(double)* nx * nx);
	for (i = 0; i < nx; i++)
	{
		Pk[i * nx + i] = f*f;
		f = va_arg(vl, double);
	}//PK主对角线元素赋值
	va_end(vl);
	return 1;
}

//////捷联惯导算法与组合导航讲义，严恭敏  p153-156  p69-p79
void C_SINSTools::StaticUpdate( double *Dw1,  double *Dw2,  double *Dw3,  double *Dw4,
	                                                 double *Df1,   double *Df2,   double *Df3,  double *Df4,  double Dt,
	                                             double *Cnb, double *Vn, double *P)
{




	C_earth cglv;
	C_Matrix matrix;
	C_SINSTools  tool;
	/////////////////////////////////////////////////////////
	///////////姿态更新//////////////////////////////////////
	///////Cb(m)n(m)=Cn(m-1)n(m)*Cb(m-1)n(m-1)*Cb(m)b(m-1)////
	
	///------------Cn(m-1)n(n)------------------------------
	/////Cn(m-1)n(n)可以由（Wnin=Wnie+Wnen，静止状态时Wnen=0）获得
	int nx = 3;
	double *wnie = (double *)malloc(nx * 1 * sizeof(double));
	memset(wnie, 0, sizeof(double)*nx * 1);
	matrix.CVect3(wnie,0,cglv.wie*cos(P[0]),cglv.wie*sin(P[0]));
	 
	///////////////////////////////////////////////////////////////////
	/////////////////////faiw=DT*wnie   Cnn=Rov(faiw)////////////////////////
   double *faiw= (double *)malloc(nx * 1 * sizeof(double));
	memset(faiw, 0, sizeof(double)*nx * 1);
	matrix.MatrixEnlarge(Dt,wnie,faiw,nx);

	double *Cnn =(double *)malloc(nx * nx * sizeof(double));
	memset(Cnn, 0, sizeof(double)*nx * nx);
    tool.Rov(Cnn,faiw);

    double *Cnn_t =(double *)malloc(nx * nx * sizeof(double));
	memset(Cnn_t, 0, sizeof(double)*nx * nx);
	matrix.MatrixTranslate(Cnn, Cnn_t, nx, nx);
	//------------------Cn(m-1)n(n)-求解完毕--------------------------------------
    
	////-----------------Cb(m)b(m-1)-------------------------------
	//////////////////Dw=Dw1+Dw2+Dw3+Dw4///////////////////////////
    double *Dw = (double *)malloc(nx * 1 * sizeof(double));
	memset(Dw, 0, sizeof(double)*nx * 1);
	for (int i = 0; i < nx; i++)
	{
		Dw[i] = Dw1[i] + Dw2[i] + Dw3[i] + Dw4[i];
	}
	
	///////////////////////////////////////////////////////////////
	////////Dwn=Dw+skew[(k1*Dw1+K2*Dw2+k3*Dw3)]*Dw4;///////////////
	////////////陀螺的四子样圆锥误差补偿//////////////////////////////////
	double *k_Dw123 = (double *)malloc(nx * 1* sizeof(double));
	memset(k_Dw123, 0, sizeof(double)*nx * 1);

    double *skew_k_Dw123 = (double *)malloc(nx * nx* sizeof(double));
	memset(skew_k_Dw123, 0, sizeof(double)*nx * nx);
    double k1 = 54.0 / 105, k2 = 92.0 / 105, k3 = 214.0 / 105;
	for (int i = 0; i < nx; i++)
	{
		k_Dw123[i] = k1*Dw1[i] + k2*Dw2[i] + k3*Dw3[i];
	
	}
    matrix.skew(k_Dw123,skew_k_Dw123);

    double *skew_Dw123_4 = (double *)malloc(nx * 1* sizeof(double));
	memset(skew_Dw123_4, 0, sizeof(double)*nx * 1);
	matrix.MatrixMult(skew_k_Dw123, Dw4, skew_Dw123_4, nx, nx, 1);


	double *Dwn = (double *)malloc(nx * 1 * sizeof(double));
	memset(Dwn, 0, sizeof(double)*nx * 1);
	matrix.MatrixAdd(Dw, skew_Dw123_4, Dwn, nx, 1);
	//////////////////////////////////////////////////////////////////
	/////////////////////////Cbb=Rov(Dwn)  ///////////////////////
	double *Cbb = (double *)malloc(nx * nx * sizeof(double));
	memset(Cbb, 0, sizeof(double)*nx * nx);
	tool.Rov(Cbb,Dwn);

	//--------------Cb(m)b(m-1)-----------------------------------

	////////////////////////Cnb2=cnn_1*Cnb1*Cbb////////////////////////
    double *Cnnb = (double *)malloc(nx * nx * sizeof(double));
	memset(Cnnb, 0, sizeof(double)*nx * nx);

    double *Cnb_2 = (double *)malloc(nx * nx * sizeof(double));
	memset(Cnb_2, 0, sizeof(double)*nx * nx);

	matrix.MatrixMult(Cnn_t,Cnb,Cnnb,nx,nx,nx);
	matrix.MatrixMult(Cnnb,Cbb,Cnb_2,nx,nx,nx);
	
	///-------------姿态更新完毕----------------------------------------



    ///----------------速度更新-----------------------------------
	////////gn=GravityUpdate(),Vncorg=gn*Dt;=////////////////////
	double *gn = (double *)malloc(nx * 1 * sizeof(double));
	memset(gn, 0, sizeof(double)*nx * 1);

	double *Vncorg = (double *)malloc(nx * 1 * sizeof(double));
	memset(Vncorg, 0, sizeof(double)*nx * 1);

	tool.GravityUpdate(gn);
	for (int i = 0; i < nx; i++)
	{
		Vncorg[i] = Dt*gn[i];
	}

	/////////Vncorg/////////////////////////////////////////////////////////
	//////////////////////Df=Df1+Df2+Df3+Df4///////////////////////
	double *Df = (double *)malloc(nx * 1 * sizeof(double));
	memset(Df, 0, sizeof(double)*nx * 1);
	for (int i = 0; i < nx; i++)
	{
		Df[i] = Df1[i] + Df2[i] + Df3[i] + Df4[i];
	}
	
	////////////////////////////////////////////////////////////////////////
	/////////////////Vbrot=0.5*Askew(Dw)*Df/////////////////////////////////
	//////////////////速度的旋转误差补偿/////////////////////////////////////
	double *skew_Dw = (double *)malloc(nx * nx* sizeof(double));
	memset(skew_Dw, 0, sizeof(double)*nx * nx);
    matrix.skew(Dw,skew_Dw);

 
    double *Vbrot = (double *)malloc(nx * 1* sizeof(double));
	memset(Vbrot, 0, sizeof(double)*nx * 1);

	matrix.MatrixMult(skew_Dw,Df,Vbrot,nx,nx,1);
	for (int i = 0; i < nx; i++)
	{
		Vbrot[i] = 0.5*Vbrot[i];
	
	}

	////////////////////////////////////////////////////////////////////////
	//////////////////Vbscul = Askew(k1*Dw1 + k2*Dw2 + k3*Dw3)*Df +...//////
	//////////////////Askew(k1*Df1 + k2*Df2 + k3*Df3)*Dw;///////////////////
	/////////////////////划桨误差补偿算法////////////////////////////////////

    double *skew_k_Dw_Df = (double *)malloc(nx * 1 * sizeof(double));
	memset(skew_k_Dw_Df, 0, sizeof(double)*nx * 1);
	matrix.MatrixMult(skew_k_Dw123, Df, skew_k_Dw_Df, nx, nx, 1);//Askew(k1*Dw1 + k2*Dw2 + k3*Dw3)*Df 


	double *k_Df123 = (double *)malloc(nx * 1 * sizeof(double));
	memset(k_Df123, 0, sizeof(double)*nx * 1);
    for (int i = 0; i < nx; i++)
	  {
		k_Df123[i] = k1*Df1[i] + k2*Df2[i] + k3*Df3[i];

	 }

	double *skew_k_Df123 = (double *)malloc(nx * nx* sizeof(double));
	memset(skew_k_Df123, 0, sizeof(double)*nx * nx);
	matrix.skew(k_Df123, skew_k_Df123);

	double *skew_k_Df_Dw = (double *)malloc(nx * 1 * sizeof(double));
	memset(skew_k_Df_Dw, 0, sizeof(double)*nx * 1);
	matrix.MatrixMult(skew_k_Df123, Dw, skew_k_Df_Dw, nx, nx, 1);//Askew(k1*Df1 + k2*Df2 + k3*Df3)*Dw


    double *Vbscul = (double *)malloc(nx * 1 * sizeof(double));
	memset(Vbscul, 0, sizeof(double)*nx * 1);

	matrix.MatrixAdd(skew_k_Dw_Df, skew_k_Df_Dw, Vbscul,nx,1);



	////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////



	////////////////////////////////////////////////////////////////////////
	////////////Vnn = Cnb1*(b1*Df1 +  b2*Df2 + b3*Df3 + b4*Df4)/////////////
	double b1 = 13.0 / 35;
    double b2 = 59.0 / 210;
    double b3 = -1.0 / 105;
    double b4 = 19.0 / 14;

	double *b_Df = (double *)malloc(nx * 1 * sizeof(double));
	memset(b_Df, 0, sizeof(double)*nx * 1);
	for (int i = 0; i < nx; i++)
	{
		b_Df[i] = b1*Df1[i] + b2*Df2[i] + b3*Df3[i] + b4*Df4[i];
	}

	double *Vnn = (double *)malloc(nx * 1 * sizeof(double));
	memset(Vnn, 0, sizeof(double)*nx * 1);
	matrix.MatrixMult(Cnb,b_Df,Vnn,nx,nx,1);

	////////////////////////////////////////////////////////////////////////
	///////////////////////Vnin = Askew(faiw)*Vnn///////////////////////////

	double *skew_faiw = (double *)malloc(nx * nx * sizeof(double));
	memset(skew_faiw, 0, sizeof(double)*nx * nx);
	matrix.skew(faiw,skew_faiw);

	double *Vnin = (double *)malloc(nx * 1 * sizeof(double));
	memset(Vnin, 0, sizeof(double)*nx * 1);

	matrix.MatrixMult(skew_faiw, Vnn, Vnin,nx,nx,1);


	////////////////////////////////////////////////////////////////////////
	///////////////////////Vnib = Cnb1*(Vbrot + Vbscul);///////////////////////////
	double *VV = (double *)malloc(nx * 1 * sizeof(double));
	memset(VV, 0, sizeof(double)*nx * 1);

    double *Vnib = (double *)malloc(nx * 1 * sizeof(double));
	memset(Vnib, 0, sizeof(double)*nx * 1);

	matrix.MatrixAdd(Vbrot,Vbscul,VV,nx,1);
	matrix.MatrixMult(Cnb,VV,Vnib,nx,nx,1);

	/////////////////////////////////////////////////////////////////////////////////
	///////////////Vnsf = Cnb1*Df - Vnin + Vnib;////////////////////////////////
	double *Cnb_Df = (double *)malloc(nx * 1 * sizeof(double));
	memset(Cnb_Df, 0, sizeof(double)*nx * 1);

    double *Vnsf = (double *)malloc(nx * 1 * sizeof(double));
	memset(Vnsf, 0, sizeof(double)*nx * 1);

	matrix.MatrixMult(Cnb,Df,Cnb_Df,nx,nx,1);

	for (int i = 0; i < nx; i++)
	{
		Vnsf[i] = Cnb_Df[i] + Vnib[i] - Vnin[i];
	}
	//////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////
	///////////////Vn=Vnsf+Vncorg;////////////////////////////////
	matrix.MatrixAdd(Vnsf, Vncorg,Vn,nx, 1);

	for (int i = 0; i < nx*nx; i++)
	{
		Cnb[i] = Cnb_2[i];

	}

	free(wnie), free(faiw), free(Cnn), free(Cnn_t), free(Dw), free(k_Dw123), free(skew_k_Dw123);
	free(skew_Dw123_4), free(Dwn), free(Cbb), free(Cnnb), free(Cnb_2), free(gn), free(Vncorg);
	free(Df), free(skew_Dw), free(Vbrot), free(skew_k_Dw_Df), free(k_Df123), free(skew_k_Df123), free(skew_k_Df_Dw), free(Vbscul);
	free(b_Df), free(Vnn), free(skew_faiw), free(Vnin), free(VV), free(Vnib), free(Cnb_Df), free(Vnsf);

}




