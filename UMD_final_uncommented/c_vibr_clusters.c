#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <stdbool.h>
#include <string.h>


void free_memory(double *tab){
	free(tab);
}

void vel_to_acc(double *Vel, int length, double timestep){

	for(int i=0 ; i<length-1 ; i++){
		Vel[i]=(Vel[i+1]-Vel[i])/timestep;
	}
	Vel[length-1]=0;
}


double* crossproduct(double *V1, double *V2){

	double *V3 = calloc(3,sizeof(double));

	V3[0] = V1[1]*V2[2]-V1[2]*V2[1];
	V3[1] = V1[2]*V2[0]-V1[0]*V2[2];
	V3[2] = V1[0]*V2[1]-V1[1]*V2[0];

	return V3;
}

double scalarproduct(double *V1, double *V2){

	double dot;

	dot = V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2];
	
	return dot;
}

double* compute_autocorrelation(double *pos, int nostep){//Compute the autocorrelation of the tab pos

	double* res;
	double normalization = nostep;
	int maxtau = nostep/2;
	
	res = calloc(maxtau,sizeof(double));
	double number;

	for(int i=0 ; i<maxtau ; i++){
		number = 0;
		for(int j=i ; j<nostep ; j++){
			number = number + pos[j]*pos[(j-i)];
		}
		res[i]=number/normalization;
		normalization--;
	}
	return res;
}

double*** relative(double *listPos, double *listVel, double *masses, int natom, int numsteps, int length){

	double *centralParameters;
	double *relativePos;
	double **relativeVel;
	double **PhiVel;
	double **ThetaVel;
	double **RadVel;
	double Direction[3];
	double norm;

	double phivel,thetavel,radvel;
	double tanvel[3];
	double vel[3];

	double ***velocityParameters = (double***)malloc(4*sizeof(double**));
	
	centralParameters = calloc(6,sizeof(double));
	relativePos = calloc(3*natom,sizeof(double));
	relativeVel = (double**)malloc(3*natom*sizeof(double*));
	PhiVel = (double**)malloc((natom-1)*sizeof(double*));
	ThetaVel = (double**)malloc((natom-1)*sizeof(double*));
	RadVel = (double**)malloc((natom-1)*sizeof(double*));

	for(int i=0 ; i<natom-1 ; i++){
		RadVel[i]=calloc(numsteps,sizeof(double));
		ThetaVel[i]=calloc(numsteps,sizeof(double));
		PhiVel[i]=calloc(numsteps,sizeof(double));
		
		relativeVel[3*i]=calloc(numsteps,sizeof(double));
		relativeVel[3*i+1]=calloc(numsteps,sizeof(double));
		relativeVel[3*i+2]=calloc(numsteps,sizeof(double));
	}

	relativeVel[3*(natom-1)]=calloc(numsteps,sizeof(double));
	relativeVel[3*(natom-1)+1]=calloc(numsteps,sizeof(double));
	relativeVel[3*(natom-1)+2]=calloc(numsteps,sizeof(double));

	velocityParameters[0] =	relativeVel;
	velocityParameters[1] = RadVel;
	velocityParameters[2] = ThetaVel;
	velocityParameters[3] = PhiVel;
	double normalization=0;

	for(int at=0 ; at<natom ; at++){
		normalization = normalization+masses[at];
	}
	for(int step=0 ; step<numsteps ; step++){
		
		for(int at=0 ; at<natom ; at++){
			centralParameters[0] = centralParameters[0]+listPos[3*numsteps*at+3*step]*masses[at];
			centralParameters[1] = centralParameters[1]+listPos[3*numsteps*at+3*step+1]*masses[at];
			centralParameters[2] = centralParameters[2]+listPos[3*numsteps*at+3*step+2]*masses[at];
			centralParameters[3] = centralParameters[3]+listVel[3*numsteps*at+3*step]*masses[at];
			centralParameters[4] = centralParameters[4]+listVel[3*numsteps*at+3*step+1]*masses[at];
			centralParameters[5] = centralParameters[5]+listVel[3*numsteps*at+3*step+2]*masses[at];
		}

		centralParameters[0] = centralParameters[0]/normalization;
		centralParameters[1] = centralParameters[1]/normalization;
		centralParameters[2] = centralParameters[2]/normalization;
		centralParameters[3] = centralParameters[3]/normalization;
		centralParameters[4] = centralParameters[4]/normalization;
		centralParameters[5] = centralParameters[5]/normalization;

		for(int at=0 ; at<natom ; at++){

	if(3*numsteps*at+3*step+2>=length){printf("ERROR : %d for %d !\n",3*numsteps*at+3*step+2,length);}

			relativePos[3*at] = listPos[3*numsteps*at+3*step] - centralParameters[0];
			relativePos[3*at+1] = listPos[3*numsteps*at+3*step+1] - centralParameters[1];
			relativePos[3*at+2] = listPos[3*numsteps*at+3*step+2] - centralParameters[2];
			relativeVel[3*at][step] = listVel[3*numsteps*at+3*step] - centralParameters[3];
			relativeVel[3*at+1][step] = listVel[3*numsteps*at+3*step+1] - centralParameters[4];
			relativeVel[3*at+2][step] = listVel[3*numsteps*at+3*step+2] - centralParameters[5];
		}
		
		for(int at=1 ; at<natom ; at++){
			Direction[0] = relativePos[3*at]-relativePos[0];
			Direction[1] = relativePos[3*at+1]-relativePos[1];
			Direction[2] = relativePos[3*at+2]-relativePos[2];

			norm = sqrt(pow(Direction[0],2)+pow(Direction[1],2)+pow(Direction[2],2));
			
			Direction[0] = Direction[0]/norm;
			Direction[1] = Direction[1]/norm;
			Direction[2] = Direction[2]/norm;			

			vel[0] = relativeVel[3*at][step];
			vel[1] = relativeVel[3*at+1][step];
			vel[2] = relativeVel[3*at+2][step];

			radvel = scalarproduct(vel,Direction);
			thetavel = scalarproduct(vel,crossproduct((double[]){1,0,0},Direction));
			phivel = scalarproduct(vel,crossproduct((double[]){0,1,0},Direction));
				
			RadVel[at-1][step] = radvel;
			ThetaVel[at-1][step] = thetavel;
			PhiVel[at-1][step] = phivel;
		}			
	
	}
	  
	return velocityParameters;
}

void correlation(double *PosMatrix, double *VelMatrix, double *masses, int natom, int numsteps, double *relativeVel_autocX,double *relativeVel_autocY,double *relativeVel_autocZ,double *RadVel_autoc,double *ThetaVel_autoc,double *PhiVel_autoc, int length, double timestep){

	//vel_to_acc(VelMatrix,3*numsteps*natom,timestep);

	double*** velocityParameters = relative(PosMatrix,VelMatrix,masses,natom,numsteps,length);

	int maxtau = numsteps/2;

	double *Vel_autocX;
	double *Vel_autocY;
	double *Vel_autocZ;
	double *Rad_autoc;
	double *Theta_autoc;
	double *Phi_autoc;

	for(int at=0 ; at<natom-1 ; at++){
		
		Vel_autocX = compute_autocorrelation(velocityParameters[0][3*at*0],numsteps);
		Vel_autocY = compute_autocorrelation(velocityParameters[0][3*at+1],numsteps);
		Vel_autocZ = compute_autocorrelation(velocityParameters[0][3*at+2],numsteps);
		Rad_autoc  = compute_autocorrelation(velocityParameters[1][at],numsteps);
		Theta_autoc= compute_autocorrelation(velocityParameters[2][at],numsteps);
		Phi_autoc  = compute_autocorrelation(velocityParameters[3][at],numsteps);				

		memcpy(&relativeVel_autocX[at*maxtau],Vel_autocX,maxtau*sizeof(double));
		memcpy(&relativeVel_autocY[at*maxtau],Vel_autocY,maxtau*sizeof(double));
		memcpy(&relativeVel_autocZ[at*maxtau],Vel_autocZ,maxtau*sizeof(double));
		memcpy(&RadVel_autoc[at*maxtau],Rad_autoc,maxtau*sizeof(double));
		memcpy(&ThetaVel_autoc[at*maxtau],Theta_autoc,maxtau*sizeof(double));
		memcpy(&PhiVel_autoc[at*maxtau],Phi_autoc,maxtau*sizeof(double));
	
	}

	Vel_autocX = compute_autocorrelation(velocityParameters[0][3*(natom-1)],numsteps);
	Vel_autocY = compute_autocorrelation(velocityParameters[0][3*(natom-1)+1],numsteps);
	Vel_autocZ = compute_autocorrelation(velocityParameters[0][3*(natom-1)+2],numsteps);
	memcpy(&relativeVel_autocX[(natom-1)*maxtau],Vel_autocX,maxtau*sizeof(double));
	memcpy(&relativeVel_autocY[(natom-1)*maxtau],Vel_autocY,maxtau*sizeof(double));
	memcpy(&relativeVel_autocZ[(natom-1)*maxtau],Vel_autocZ,maxtau*sizeof(double));	

	free(Vel_autocX);
	free(Vel_autocY);
	free(Vel_autocZ);
	free(Rad_autoc);
	free(Theta_autoc);
	free(Phi_autoc);

	free(velocityParameters[0][3*(natom-1)]);
	free(velocityParameters[0][3*(natom-1)+1]);
	free(velocityParameters[0][3*(natom-1)+2]);
	
	for(int i=0 ; i<natom-1 ; i++){
		free(velocityParameters[0][3*i]);
		free(velocityParameters[0][3*i+1]);
		free(velocityParameters[0][3*i+2]);
		free(velocityParameters[1][i]);
		free(velocityParameters[2][i]);
		free(velocityParameters[3][i]);
	}

	free(velocityParameters[0]);
	free(velocityParameters[1]);
	free(velocityParameters[2]);
	free(velocityParameters[3]);
}