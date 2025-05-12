#ifndef _WENO_H
#define _WENO_H

double weno3P(double *f)
{

	//assign value to v1, v2,...
	int k = 0;
	double v1 = *(f + k - 1);
	double v2 = *(f + k);
	double v3 = *(f + k + 1);

	//smoothness indicator
	double epsilon = 1.0e-16;
	double s1 = (v2 - v3)*(v2 - v3) + epsilon;
	double s2 = (v2 - v1)*(v2 - v1) + epsilon;

	//weights
	double a1 = 1.0e1/s1/s1;
	double a2 = 1.0e1/s2/s2;
	double tw1 = 1.0 / (a1 +a2);
	double w1 = tw1*a1;
	double w2 = tw1*a2;

	//return weighted average
	return  w1*(0.5*v2 + 0.5*v3)
		  + w2*(-0.5*v1 + 1.5*v2);
}

double weno3M(double *f)
{

	//assign value to v1, v2,...
	int k = 1;
	double v1 = *(f + k + 1);
	double v2 = *(f + k);
	double v3 = *(f + k - 1);

	//smoothness indicator
	double epsilon = 1.0e-16;
	double s1 = (v2 - v3)*(v2 - v3) + epsilon;
	double s2 = (v2 - v1)*(v2 - v1) + epsilon;

	//weights
	double a1 = 1.0e1/s1/s1;
	double a2 = 1.0e1/s2/s2;
	double tw1 = 1.0 / (a1 +a2);
	double w1 = tw1*a1;
	double w2 = tw1*a2;

	//return weighted average
	return  w1*(0.5*v2 + 0.5*v3)
		  + w2*(-0.5*v1 + 1.5*v2);
}


#endif