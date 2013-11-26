#include <stdio.h>
#include <math.h>

void RK4(double t, double x[], double h, int n, int ic, 
	void(*F)(double t, double x[], double k[])){
	double K1[ic];
	double K2[ic];
	double K3[ic];
	double K4[ic];
	double R[ic];

	for(int i = 0; i < n; i++){
		F(t, x, K1);

		for(int j = 0; j < ic; j++)
			R[j] = x[j] + h/2*K1[j];
		F(t+h/2, R, K2);

		for(int j = 0; j < ic; j++)
			R[j] = x[j] + h/2*K2[j];
		F(t+h/2, R, K3);


		for(int j = 0; j < ic; j++)
			R[j] = x[j] + h*K3[j];
		F(t+h, R, K4);

		for(int j = 0; j < ic; j++)
			x[j] = x[j] + h/6*(K1[j]+2*(K2[j]+K3[j])+K4[j]);
		t = t + h;
		printf("%0.10f %0.10f\n", x[0], x[2]);
	}
}

void baseball(double t, double x[], double R[]){
	double k = .002*sqrt(x[1]*x[1]+x[3]*x[3]);
	double g = 32;
	R[0] = x[1];
	R[1] = -k*x[1];
	R[2] = x[3];
	R[3] = -k*x[3]-g;
}

void mercury(double t, double x[], double R[]){
	double G = 6.67384e-11;
	double M = 1.9891e30;
	double pow = sqrt(x[0]*x[0]+x[2]*x[2]);
	double comp = -M*G/(pow*pow*pow);
	R[0] = x[1];
	R[1] = x[0]*comp;
	R[2] = x[3];
	R[3] = x[2]*comp;
}

int main(){
	/*Initial conditions for basebal trajectory
	double t = 0.0;
	double h = 0.0647859;
	double v = 208.0;
	double theta = 43*M_PI/180;
	double x[4];
	x[0] = 0.0; x[1] = v*cos(theta);
	x[2] = 3.0; x[3] = v*sin(theta);
	RK4(t, x, h, 100, 4, baseball);
	*/
	//initial conditions for mercury
	double t = 0.0;
	double h = 76005.216;
	double x[4];
	x[0] = 4.6e10; x[1] = 0;
	x[2] = 0; x[3] = 5.898e4;
	RK4(t, x, h, 100, 4, mercury);
}
