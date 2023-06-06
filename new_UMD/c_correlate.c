/*
 *Author : Kevin Jiguet-Covex
 *Date : June 6, 2023
 */


#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void compute_autocorrelation(const double *data, double *correlation, double *sigma, const int origin , const int length, const int firststep, const int lendata)//Calculates the correlation function and msd of a discrete signal. 
{

	double mean = 0;

	for (int i = firststep ; i < lendata; i = i + 1){

		mean = mean + data[i];

	}

	mean = mean/(lendata-firststep);//calculation of the mean of the data values

	for (int orig = firststep ; orig < (lendata-length) ; orig = orig + origin){//For each orig value, the autocorrelation function will be calculated

		double DeltaOrig=data[orig]-mean;

		for (int i=0 ; i<length ; i = i+1){//Calculation of the autocorrelation function, with a shift from O to length-1

			double element = DeltaOrig*(data[orig+i]-mean);
			correlation[i] = correlation[i] + element;
			sigma[i] = sigma[i] + pow(element,2);
		}
	}

	for (int i=0 ; i<length ; i=i+1){//Normalization by the length of the data (!! not by the mean !! This one is done in the python script)

		correlation[i] = correlation[i]/(lendata-length-firststep);
		sigma[i] = sigma[i]/(lendata-length-firststep);
    }

	for (int i=0 ; i<length ; i=i+1){

 		sigma[i] = sqrt(fabs(sigma[i]-pow(correlation[i],2)));
	}
}
