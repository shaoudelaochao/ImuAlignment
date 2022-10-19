
class C_KalmanFilter
{
public:
	C_KalmanFilter();
	~C_KalmanFilter();
public:
	int nx, nz;


	C_KalmanFilter(int nx0, int nz0);

	

	int Timeupdate(double *Ft, double *Xk, double *Pk, double *Qt,int nx);////时间更新
	int MeasUpdate(double *Pk, double *Rk, double *Hk, double *Xk, double *Zk,int nx,int nz);///量测更新
	int FilterBack(double *Xk1,double *Vn1,double *Cnb1,double *eb);////反馈



	void SetFt(double *Pos0, double g, double *Ft);
    void SetHk(double *Hk);

	void SetPkRkQt(double *Pk, double *Rk, double *Qt);

	void expmFk( double *Ft, double Dt, double *Fk );///Fk设置




	void setQk( double *Cnb,  double *Qt,  double *Fk,
	 double Dt, double *Gt, double *Qk);///Qk设置

	/*int SetPk(double *Pk, double f, ...);
	int SetQt(double *Qt);
	int SetRk(double *Rk, double f, ...);
	int SetZk(double *Zk, double f, ...);
*/
	//Zk,Hk,Ft;
};
