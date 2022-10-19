#include"stdafx.h"

// time struct
// time struct
typedef struct
{
	int time;     // seconds from 1970-1-1 00:00:00 UTC
	double  sec;    // fraction of second under 1 s
} Time_gbt;

// imu data struct
typedef struct
{
	Time_gbt   time;
	double    gyro[3], acc[3];
	int imufreq;
	double temp;
	double gtime;
} IMUObs_t;


////精对准
class C_AlignFine
{
public:
	C_AlignFine();
	~C_AlignFine();
public:

	void InitAlignFine(void);////初始化
	void AlignFine( double *C_att, double *Cnb, double *F_att, double *eb, double *db);////精对准函数
	void PrintData();///输出函数

	double eb[3];
	double db[3];
	double c_euler[3];
    double f_euler[3];
	
	double Pos0[3];
	double Mne[9];////////NED到ENU的转换矩阵/////////////////////////////////////
	double Cnb[9];
	double gn[3];
	
	double Ft[7 * 7];
	double Fk[7 * 7];
	double Qk[7 * 7];
	double Gt[7 * 6];
	double Hk[2 * 7];
	double Xk[7 * 1];/////// Xk=[失准角（ENU） 速度误差（EN）  陀螺随机常值漂移（NU）]
	double Pk[7 * 7];
	double Qt[6 * 6];
	double Rk[2 * 2];
	double Zk[2 * 1];
	
	double Dw1_0[3], Dw2_0[3], Dw3_0[3], Dw4_0[3];
	double Df1_0[3], Df2_0[3], Df3_0[3], Df4_0[3];
	
	double Dw1[3], Dw2[3], Dw3[3], Dw4[3];
	double Df1[3], Df2[3], Df3[3], Df4[3];

	double Vn[3];
	double ACC[3];
	double AT[3];
	
};