#include "disorder.h"

double *ini_torques(Ran * myran, double epsilon)
{
	double d = 1.0 / (Grid::mm - 1);
	double *disorder = new double[Grid::mm];
	for (int i = 0; i < Grid::mm; i++)
		disorder[i] = (-0.5 + i*d)*epsilon*2.0*PI;
	shuffle(disorder, Grid::mm, myran);
	return disorder;
}
