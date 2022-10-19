#include"stdafx.h"

C_AlignCoarse::C_AlignCoarse()
{
}

C_AlignCoarse::~C_AlignCoarse()
{
}

void C_AlignCoarse::InitAlignC(void)////初始化
{
	int i = 0;
	for (i = 0; i < 3; i++)
	{
		gn[i] = var1[i] = var2[i] = var3[i] = fbib[i] = wbib[i] = fbib_1[i] = wbib_1[i] = 0;
		vbr1[i] = vbr2[i] = vbr3[i] =Att[i]= 0;
		skew_vbr1[9], skew_vbr2[9];
	}

	for (i = 0; i < 9; i++)
	{
		skew_vbr1[i] = skew_vbr2[i] = 0;
	}
	
}



/**************解析初对准的方法************************/
/////参考文献：捷联惯导系统粗对准方法比较    魏春岭 　 张洪钺
void C_AlignCoarse::AlignC( double *euler,double *Cnb, double *eb,double *db)
{

	C_Matrix matrix;
	C_SINSTools tool;
	C_earth cglv;
	IMUObs_t  imu_data = {0};
	int lenNum = 0;
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
		    lenNum++;
			for (int i = 0; i < 3; i++)
			{
				wbib[i] += imu_data.gyro[i];
				fbib[i] += imu_data.acc[i];
			}
	}
	///////////原始数据给出的是增量形式，200hz，
	//////////角增量转换为角速率形式
	for (int i = 0; i < 3; i++)
	{
		fbib[i] = 200 * fbib[i] / lenNum;
		wbib[i] = 200 * wbib[i] / lenNum;

		eb[i] = wbib[i] / 200.0;
		db[i] = fbib[i] / 200.0;
	}


	double latitude = L0*cglv.deg;
    double longitude = B0*cglv.deg;

	tool.GravityUpdate(gn);//重力更新

	double gy = gn[1];
	double gz = gn[2];

	double wN = cglv.wie*cos(latitude);
	double wU = cglv.wie*sin(latitude);

	matrix.CVect3(var1,0,-gy,-gz);
	matrix.CVect3(var2, (gz*wN - gy*wU), 0, 0);

    matrix.CVect3(var3,0,var2[0]*gz,-var2[0]*gy);

	double C1[9] = {var1[0],var1[1],var1[2],
                    var2[0],var2[1],var2[2],
                    var3[0],var3[1],var3[2]};

	
	///////////////////////////////////////////

	matrix.MatrixAssign(fbib,vbr1,3);

	matrix.skew(vbr1,skew_vbr1);//反对称矩阵

	matrix.MatrixMult(skew_vbr1,wbib,vbr2,3,3,1);

	matrix.skew(vbr2,skew_vbr2);//反对称矩阵

	matrix.MatrixMult(skew_vbr2, vbr1, vbr3, 3, 3, 1);

	double C2[9] = {vbr1[0],vbr1[1],vbr1[2],
	                vbr2[0],vbr2[1],vbr2[2],
	                vbr3[0],vbr3[1],vbr3[2]};
	                                         
	matrix.MatrixInverse(C1,3);

	matrix.MatrixMult(C1,C2, Cnb,3,3,3);
	
	CVect3 att(0.0);
	att = m2att(Cnb);
	euler[0] = att.i; euler[1] = att.j; euler[2] = att.k;

}
