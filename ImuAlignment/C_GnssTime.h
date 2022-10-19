#pragma once
#include"stdafx.h"
#ifdef WIN_DLL
#define EXPORT __declspec(dllexport) /* for Windows DLL */
#else
#define EXPORT
#endif
// time struct
typedef struct
{
	int time;     // seconds from 1970-1-1 00:00:00 UTC
	double sec;    // fraction of second under 1 s
} gtime_t;




	/* time and string functions -------------------------------------------------*/
	EXPORT double  str2num(const char *s, int i, int n);
	EXPORT int     str2time(const char *s, int i, int n, gtime_t *t);
	EXPORT void    time2str(gtime_t t, char *str, int n);
	EXPORT gtime_t epoch2time(const double *ep);
	EXPORT void    time2epoch(gtime_t t, double *ep);
	EXPORT void    time2epoch_n(gtime_t t, double *ep, int n);
	EXPORT gtime_t gpst2time(int week, double sec);
	EXPORT double  time2gpst(gtime_t t, int *week);
	EXPORT gtime_t gst2time(int week, double sec);
	EXPORT double  time2gst(gtime_t t, int *week);
	EXPORT gtime_t bdt2time(int week, double sec);
	EXPORT double  time2bdt(gtime_t t, int *week);
	EXPORT char    *time_str(gtime_t t, int n);

	EXPORT gtime_t timeadd(gtime_t t, double sec);
	EXPORT double  timediff(gtime_t t1, gtime_t t2);
	EXPORT gtime_t gpst2utc(gtime_t t);
	EXPORT gtime_t utc2gpst(gtime_t t);
	EXPORT gtime_t gpst2bdt(gtime_t t);
	EXPORT gtime_t bdt2gpst(gtime_t t);
	EXPORT gtime_t timeget(void);
	EXPORT void    timeset(gtime_t t);
	EXPORT void    timereset(void);
	EXPORT double  time2doy(gtime_t t);
	EXPORT double  utc2gmst(gtime_t t, double ut1_utc);
	EXPORT int read_leaps(const char *file);

	EXPORT int adjgpsweek(int week);
	EXPORT uint32_t tickget(void);
	EXPORT void sleepms(int ms);
	


