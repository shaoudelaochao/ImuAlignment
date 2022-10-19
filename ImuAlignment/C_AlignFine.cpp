#include"stdafx.h"

char enufilehead1[512] = "% (e/n/u-baseline=WGS84,Q=1:fix,2:float,3:sbas,4:dgps,5:single,6:ppp,ns=# of satellites)";
char enufilehead2[512] = "%  GPST                  e-baseline(m)  n-baseline(m)  u-baseline(m)   Q  ns   sde(m)   sdn(m)   sdu(m)  sden(m)  sdnu(m)  sdue(m) age(s)  ratio";



C_AlignFine::C_AlignFine()
{
}

C_AlignFine::~C_AlignFine()
{
}

void C_AlignFine::InitAlignFine(void)
{
	int i = 0;
	for (i = 0; i < 3; i++)
	{
        c_euler[i]=0;
        f_euler[i] = 0;
        Pos0[i] = 0;
        gn[i] = 0;
       Dw1_0[i] = Dw2_0[i] = Dw3_0[i] = Dw4_0[i] = 0;
       Df1_0[i] = Df2_0[ i] = Df3_0[i] =Df4_0[i] = 0;
       Dw1[i] = Dw2[i] = Dw3[i] = Dw4[i] = 0;
       Df1[i] = Df2[i] = Df3[i] = Df4[i] = 0;
	   Vn[i] = ACC[i] = AT[i] = 0;
	   eb[i] = db[i] = 0;
	}
	
	for (i = 0; i < 7*7; i++)
	{
      Ft[i] = 0; Fk[i] = 0; Qk[i] = 0; Pk[i] = 0;
	}

	for (i = 0; i < 3*3; i++)  Mne[i] = 0,Cnb[i]=0;
	for (i = 0; i < 7*6; i++) Gt[i] = 0;
	for (i = 0; i < 2*7; i++) Hk[i] = 0;
	for (i = 0; i < 7*1; i++) Xk[i] = 0;
	for (i = 0; i < 6*6; i++) Qt[i]=0;
	for (i = 0; i < 2*2; i++) Rk[i] = 0;
	for (i = 0; i < 2*1; i++) Zk[i] = 0;



}



static void OpenOutFile( char outfile[6][512] , FILE *openfile[6])
{
	int i = 0;
	for (i = 0; i < 6; i++)
	{
		strcpy(outfile[i], IMUFileLog);
	}

	strcat(outfile[0], "_att.pos");
	strcat(outfile[1], "_Patt.pos");
	strcat(outfile[2], "_vn.pos");
	strcat(outfile[3], "_Pvn.pos");
	strcat(outfile[4], "_eb.pos");
	strcat(outfile[5], "_Peb.pos");

	for (int i = 0; i < 6; i++)
	{
		if (!(openfile[i] = fopen(outfile[i], "wt")))
		{
			printf("loadopts: options file open error(%s)\n", openfile[i]);	
		}

		fprintf(openfile[i], "%s\n", enufilehead1);
		fprintf(openfile[i], "\n");
		fprintf(openfile[i], "%s\n", enufilehead2);
		fprintf(openfile[i], "\n");

	}
	
}




static void OutResult(FILE *openfile[6],double *att,double *vn,double *eb, double *Pk,  gtime_t  times)
{

	char time_s[63] = "";
	char sep[6] = "  ";
	time2str(times, time_s, 3);
	
	double Att[3] = { 0.0 }, Patt[3] = { 0.0 }, vel[3] = { 0.0 }, Pvel[3] = { 0.0 },eb_g[3] = { 0.0 }, Peb[3] = { 0.0 };


	Att[0] = att[0] *kRAD2DEG; Att[1] = att[1] *kRAD2DEG; Att[2] = att[2]*kRAD2DEG;
	vel[0] = vn[0]; vel[1] = vn[1]; vel[2] = vn[2];
	eb_g[0]=eb[0]* kRAD2DEG; eb_g[1] = eb[1] * kRAD2DEG; eb_g[2] = eb[2] * kRAD2DEG;


	Patt[0] = sqrt(Pk[0 * 7 + 0])*kRAD2DEG; 
	Patt[1] = sqrt(Pk[1 * 7 + 1]) * kRAD2DEG;
	Patt[2] = sqrt(Pk[2 * 7 + 2]) * kRAD2DEG;

	Pvel[0] = sqrt(Pk[3 * 7 + 3]);
	Pvel[1] = sqrt(Pk[4 * 7 + 4]);

	Peb[0] = sqrt(Pk[5 * 7 + 5])*kRAD2DEG;
	Peb[1] = sqrt(Pk[6 * 7 + 6])*kRAD2DEG;



	fprintf(openfile[0], "%s%s%13.8f%s%13.8f%s%13.8f%s%d\n",
				time_s, sep, Att[0] , sep, Att[ 1],  sep, Att[ 2] , sep, 4);

	fprintf(openfile[1], "%s%s%13.8f%s%13.8f%s%13.8f%s%d\n",
		time_s, sep, Patt[0] , sep, Patt[1] , sep, Patt[2] , sep, 4);

	fprintf(openfile[2], "%s%s%13.8f%s%13.8f%s%13.8f%s%d\n",
		time_s, sep, vel[0], sep, vel[1], sep, vel[2], sep, 4);

	fprintf(openfile[3], "%s%s%13.8f%s%13.8f%s%13.8f%s%d\n",
		time_s, sep, Pvel[0], sep, Pvel[1], sep, Pvel[2], sep, 4);

	fprintf(openfile[4], "%s%s%13.8f%s%13.8f%s%13.8f%s%d\n",
		time_s, sep, eb[0] * 3600, sep, eb[1] * 3600, sep, eb[2] * 3600, sep, 4);

	fprintf(openfile[5], "%s%s%13.8f%s%13.8f%s%13.8f%s%d\n",
		time_s, sep, Peb[0]*3600, sep, Peb[1] * 3600, sep, Peb[2] * 3600, sep, 4);
	
}












void C_AlignFine::AlignFine(  double *C_att,  double *Cnb, double *F_att, double *eb, double *db)
{
	C_Matrix matrix;
	C_SINSTools tool;
	C_earth cglv;
	C_AlignCoarse aligncoarse;
	C_KalmanFilter KF;
	gtime_t  times;

	/////输出文件---------------------------------
	char outfile[6][512] = { " " };
	FILE *openfile[6];

	OpenOutFile(outfile, openfile);
	
	//---------------------------------------------
	//////////////////////打开IMU文件/////////////////
	IMUObs_t    imu_data = {0},imuobs[4] = { 0 };
	CIMU();
	CSINS sin_;
	sin_.att.i = C_att[0]; sin_.att.j = C_att[1]; sin_.att.k = C_att[2];
	sin_.Cnb = a2mat(sin_.att); sin_.qnb = a2qua(sin_.att);
	sin_.pos.i = L0* cglv.deg; sin_.pos.j = B0* cglv.deg; sin_.pos.k = H0;///静基座下的机械编排位置不变
	sin_.eb.i = eb[0];	sin_.eb.j = eb[1];	sin_.eb.k = eb[2];

	//------------姿态和位置初始化
	Cnb[0] = sin_.Cnb.e00; Cnb[1] = sin_.Cnb.e01; Cnb[2] = sin_.Cnb.e02;
	Cnb[3] = sin_.Cnb.e10; Cnb[4] = sin_.Cnb.e11; Cnb[5] = sin_.Cnb.e12;
	Cnb[6] = sin_.Cnb.e20; Cnb[7] = sin_.Cnb.e21; Cnb[8] = sin_.Cnb.e22;
	matrix.CVect3(Pos0, L0*cglv.deg, B0* cglv.deg, H0);


	/////////   设置Ft  Gt Hk   ///////////////////////
	////////    严恭敏，156////////////////////////////
	tool.GravityUpdate(gn);
	double g = gn[2];

	KF.SetFt(Pos0, g, Ft);
	KF.SetHk(Hk);
	KF.SetPkRkQt(Pk, Rk, Qt);


	int lenNum = 0;
	int MeasUpdateNum = 0;

	double timespan = 1.0 / 200;//////200HZ的频率
	//////////////////////打开IMU文件/////////////////
	FILE *imu_file;
	if (!(imu_file = fopen(IMUFileLog, "rt")))
	{
		printf("loadopts: options file open error(%s)\n", imu_file);
	}

	while (!feof(imu_file))
	{

		fscanf(imu_file, "%lf   %lf  %lf   %lf   %lf   %lf   %lf   \n",
			&imu_data.gtime,
			&imu_data.gyro[0], &imu_data.gyro[1], &imu_data.gyro[2],
			&imu_data.acc[0], &imu_data.acc[1], &imu_data.acc[2]);
		    imu_data.time.time = int(imu_data.gtime);
			imu_data.time.sec = imu_data.gtime- int(imu_data.gtime);

			
			memcpy(&imuobs[lenNum++], &imu_data,  sizeof(IMUObs_t));
			

			if (lenNum == 4 )
			{


				lenNum = 0;
				double dt = 4 * timespan;
				double *Vn = (double *)malloc(2 * 1 * sizeof(double));
				memset(Vn, 0, sizeof(double) * 2 * 1);

				tool.StaticUpdate(imuobs[0].gyro, imuobs[1].gyro, imuobs[2].gyro, imuobs[3].gyro,
					imuobs[0].acc, imuobs[1].acc, imuobs[2].acc, imuobs[3].acc,
					dt, Cnb, Vn, Pos0);

#if 0
				double Dt = 0;
				sin_.vn.i = 0, sin_.vn.j = 0, sin_.vn.k =0;
				sin_.UpdateSins(imuobs,4, Dt);
				dt = Dt;
				//------------姿态和位置初始化
				Cnb[0] = sin_.Cnb.e00; Cnb[1] = sin_.Cnb.e01; Cnb[2] = sin_.Cnb.e02;
				Cnb[3] = sin_.Cnb.e10; Cnb[4] = sin_.Cnb.e11; Cnb[5] = sin_.Cnb.e12;
				Cnb[6] = sin_.Cnb.e20; Cnb[7] = sin_.Cnb.e21; Cnb[8] = sin_.Cnb.e22;
				Vn[0]= sin_.vn.i ,  Vn[1]=sin_.vn.j ,  Vn[2]= sin_.vn.k;
				if (dt < 0.0001) continue;
#endif



				KF.expmFk(Ft, dt, Fk);////Ft离散化
				KF.setQk(Cnb, Qt, Fk, dt, Gt, Qk);///Qk设置
				KF.Timeupdate(Fk, Xk, Pk, Qk, 7);///时间更新

				MeasUpdateNum++;

				if (MeasUpdateNum % 2 == 0)
				{
					Zk[0] = Vn[0]; Zk[1] = Vn[1];
					KF.MeasUpdate(Pk, Rk, Hk, Xk, Zk, 7, 2);
					KF.FilterBack(Xk, Vn, Cnb, eb);
				}



				sin_.vn.i = Vn[0], sin_.vn.j = Vn[1], sin_.vn.k = Vn[2];
				sin_.Cnb.e00 = Cnb[0]; sin_.Cnb.e01 = Cnb[1]; sin_.Cnb.e02 = Cnb[2];
				sin_.Cnb.e10 = Cnb[3]; sin_.Cnb.e11 = Cnb[4]; sin_.Cnb.e12 = Cnb[5];
				sin_.Cnb.e20 = Cnb[6]; sin_.Cnb.e21 = Cnb[7]; sin_.Cnb.e22 = Cnb[8];
				sin_.att = m2att(sin_.Cnb);
				sin_.qnb = m2qua(sin_.Cnb);


				double att[3] = { sin_.att.i,sin_.att.j,sin_.att.k };

				printf("%s  %lf %lf %lf \n", "att", att[0] * kRAD2DEG, att[1] * kRAD2DEG, att[2] * kRAD2DEG);
				times.time = int(imuobs[0].time.time);
				times.sec = imuobs[0].time.sec;
				OutResult(openfile, att, Vn, eb, Pk, times);	
				for (int i = 0; i < 4;i++)   memset(&imuobs[i], 0x0, sizeof(IMUObs_t));
				
			}

	}
	F_att[0] = sin_.att.i* kRAD2DEG, F_att[1] = sin_.att.j* kRAD2DEG, F_att[2] = sin_.att.k* kRAD2DEG;

}












void C_AlignFine::PrintData()
{
	
	C_AlignCoarse aligncoarse;
	C_AlignFine alignfine;

	///初始化
	aligncoarse.InitAlignC();
	alignfine.InitAlignFine();

	aligncoarse.AlignC(c_euler, Cnb, eb, db);///粗对准
	cout << "粗对准结果: 横滚  俯仰  航向 " << endl;
	cout << c_euler[0]*kRAD2DEG << "  " << c_euler[1] * kRAD2DEG << "  " << c_euler[2] * kRAD2DEG << endl;

	alignfine.AlignFine(c_euler,Cnb, f_euler, eb, db);///精对准
	cout << "精对准结果 横滚  俯仰  航向 " << endl;
    cout << f_euler[0]  << "  " << f_euler[1]  << "  " << f_euler[2]  << endl;


}


