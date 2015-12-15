#include <iostream>
#include <cmath>

using namespace std;

void kvec(double* const k, const double* const x,const double mu);
void step(const double* const x, const double* const x5, double& max);
int main(){

	const double mu=0.012277471, tol=1e-5, q=0.5, tEnd=17.065216560157;
	const double a21=0.2, a31=3.0/40, a32=9.0/40, a41=44.0/45, a42=-56.0/15, a43=32.0/9, a51=19372.0/6561, a52=-25360.0/2187, a53=64448.0/6561, a54=-212.0/729, a61=9017.0/3168, a62=-355.0/33, a63=46732.0/5247, a64=49.0/176, a65=-5103.0/18656, a71=35.0/384, a72=0, a73=500.0/1113, a74=125.0/192, a75=-2187.0/6784, a76=11.0/84;
	const double b51=a71, b53=a73, b54=a74, b55=a75, b56=a76, b41=5179.0/57600, b43=7571.0/16695, b44=393.0/640, b45=-92097.0/339200, b46=187.0/2100, b47=1.0/40;	
	double t=0, dt=1e-3, x[4], x5[4], k1[4], k2[4], k3[4], k4[4], k5[4], k6[4], k7[4], max;	

	x[0]=0.994, x[2]=0; //x(0) bzw. y(0)
	x[1]=0, x[3]=-2.00158510637908; //dx(0)/dt bzw. dy(0)/dt
	
	cout << t << "\t" << dt << "\t" << x[0] << "\t" << x[2] << endl;	

	while(t<=tEnd){
		t+=dt;

		kvec(k1,x,mu);

		double tmp[4];
		for(int i=0; i<4; i++)
			tmp[i] = x[i]+dt*a21*k1[i];
 		kvec(k2,tmp,mu);

		for(int i=0; i<4; i++)
			tmp[i] = x[i]+dt*(a31*k1[i]+a32*k2[i]);
 		kvec(k3,tmp,mu);
		
		for(int i=0; i<4; i++)
			tmp[i] = x[i]+dt*(a41*k1[i]+a42*k2[i]+a43*k3[i]);
 		kvec(k4,tmp,mu);

		for(int i=0; i<4; i++)
			tmp[i] = x[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]);
 		kvec(k5,tmp,mu);

		for(int i=0; i<4; i++)
			tmp[i] = x[i]+dt*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]);
 		kvec(k6,tmp,mu);

		for(int i=0; i<4; i++)
			tmp[i] = x[i]+dt*(a71*k1[i]+a72*k2[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]);
 		kvec(k7,tmp,mu);

		for(int i=0; i<4; i++){
			x5[i] = x[i]+dt*(b51*k1[i]+b53*k3[i]+b54*k4[i]+b55*k5[i]+b56*k6[i]);			
			x[i]  += dt*(b41*k1[i]+b43*k3[i]+b44*k4[i]+b45*k5[i]+b46*k6[i]+b47*k7[i]);	
		}

		step(x,x5,max);

		//if (max>tol)
		dt *= q*pow((tol/max),0.2);

		cout << t << "\t" << dt << "\t" << x[0] << "\t" << x[2] << "\t" << endl;
	}	

	return 0;
}

void kvec(double* const k, const double* const x, const double mu){
	double r = sqrt((x[0]+mu)*(x[0]+mu)+x[2]*x[2]); 
	double s = sqrt((x[0]-1+mu)*(x[0]-1+mu)+x[2]*x[2]);	

	k[0] = x[1]; // dx/dt
	k[1] = x[0]+2*x[3]-(1-mu)*(x[0]+mu)/(r*r*r)-mu*(x[0]-1+mu)/(s*s*s); //d²x/dt²
	k[2] = x[3]; //dy/dt
	k[3] = x[2]-2*x[1]-(1-mu)*x[2]/(r*r*r)-mu*x[2]/(s*s*s); //d²y/dt²
}

void step(const double* const x, const double* const x5, double& max){
	double a;
	max=0;
	for(int i=0; i<4; i++){
			a = abs(x[i]-x5[i]);
			if (a>max)
				max = a;
	}
}
