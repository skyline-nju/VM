#ifndef DISORDER_H
#define DISORDER_H
#include "grid.h"

double *ini_torques(Ran *myran, double epsilon);

template<class T>
void shuffle(T *a, int n, Ran *myran)
{
	for (int i = n - 1; i >= 0; i--)
	{
		int j = int(myran->doub() * (i + 1));
		if (j > i)
			j = i;
		else if (j < 0)
			j = 0;
		T tmp = a[i];
		a[i] = a[j];
		a[j] = tmp;
	}
}

#endif

