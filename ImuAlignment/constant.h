//////常数参数设置
/*------------------------------------------------------*/
/*-------------------------------------------------------*/

#define IMUFileLog   "..//imu_data.txt"
#define PI   3.14159265358979
#define kPI  3.14159265358979



/*--------------Earth Information Initialization--------------*/

#ifndef PI
#define PI		3.14159265358979
#endif
#define PI_2	(PI/2.0)
#define PI_4	(PI/4.0)
#define _2PI	(2.0*PI)

#define sqrt2	1.414213562373095	// sqrt(2) ...
#define sqrt3	1.732050807568877
#define sqrt5	2.236067977499790
#define sqrt6	2.449489742783178
#define sqrt7	2.645751311064591
#define sqrt8	2.828427124746190

#define DEG		(PI/180.0)		// arcdeg
#define MIN		(DEG/60.0)		// arcmin
#define SEC		(MIN/60.0)		// arcsec
#define HUR		3600.0			// hur
#define SHUR	60.0			// sqrt(hur)
#define DPS		(DEG/1.0)		// deg/s
#define DPH		(DEG/HUR)		// deg/h
#define DPSH	(DEG/SHUR)		// deg/sqrt(h)
#define G0		9.7803267714
#define MG		(G0/1.0e3)
#define UG		(G0/1.0e6)		// ug
#define UGPSHZ	(UG/1)			// ug/sqrt(Hz)
#define RE		6378137.0
#define PPM		1.0e-6


#define kDEG2RAD         (kPI/180.0)          /* deg to rad */
#define kRAD2DEG         (180.0/kPI)          /* rad to deg */
#define kSC2RAD      3.1415926535898     /* semi-circle to radian (IS-GPS) */
#define kAU          149597870691.0      /* 1 AU (m) */
#define kAS2R        (kDEG2RAD/3600.0)        /* arc sec to radian */

		

//////初始对准时间段：纬度，经度，大地高，的平均值。

#define L0   42.0475181397189//	////latitude, longitude, high
#define B0   120.657165442596//
#define H0   145.577921603147//