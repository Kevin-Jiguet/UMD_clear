/*
 * Author: Kevin Jiguet-Covex
 * Date:   May 31, 2023
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <stdbool.h>

void free_memory(int *tab){
	free(tab);
}

int min(int a,int b){

	if(a>b){return b;}
	else{return a;}

}	

double min3(double a, double b, double c){

	double minimum = a;

	if(minimum>b){minimum=b;}
	if(minimum>c){minimum=c;}

	return minimum;

}


double* read_coord(char *MySnapshot, int nAtoms, int len){//reads the xcart coordinates of a snapshot and returns it in a tab

	int depth = 1;
	int atom = 0;
	int flagxcart = 0;
	double number = 0;
	double* Coordinates;
	int lineIndex=0;

	Coordinates = calloc(3*nAtoms,sizeof(double));
	
	for(int i=0 ; i<len ; i++){

		if(MySnapshot[i]=='\t'){
		

			if(flagxcart){
				Coordinates[3*atom+lineIndex-3]=42+0*number;
				flagxcart=0;
			}
			number=0;
			lineIndex++;
			depth=1;
		}

		else if(lineIndex==3||lineIndex==4||lineIndex==5){
			flagxcart = 1;
			if(isdigit(MySnapshot[i])){
				if(depth == 0){
					number=number*10;
					number=number+(MySnapshot[i]-'0');
				}
			
				else{
					depth++;
					number=number+(MySnapshot[i]-'0')/depth;
				}
			}
			

		}
		else if(MySnapshot[i]=='\n'){
			atom++;
			lineIndex=0;
			depth=1;
			number=0;
		}
	}
return Coordinates;
}

double* angles(int *polyhedras, int nPol, int M, int *MySnapshot, int CentMin, int CentMax, int OutMin, int OutMax, double acell0, double acell1, double acell2){

	double* AnglesList;
	int* OutList;
	int* CoordList;
	int at,ne1,ne2;	
	int xCent,yCent,zCent,Deltx1,Delty1,Deltz1,Deltx2,Delty2,Deltz2,valx1,valy1,valz1,valx2,valy2,valz2,deltzz,deltxx,deltyy,valxx,valyy,valzz,dist1,dist2,distij,angle;
	int centat=1;
	int nbounds=0 ;
	int ncent=0;
	int flagcoord=0;


	CoordList = calloc(M,sizeof(int));
	OutList = calloc(3*M,sizeof(int));
	AnglesList = calloc((OutMax-OutMin+1)*(OutMax-OutMin)/2+1+CentMax-CentMin+1,sizeof(double));
	AnglesList[0]=1;



	while(ncent<nPol){

		at = polyhedras[ncent];
		ncent++;
		if(at==-1){
			AnglesList[(int)AnglesList[0]]=-1;
			AnglesList[0]++;
			centat=1;
			nbounds = 0;
		}		
		
		else if(centat==1){
			centat=0;
			xCent = MySnapshot[3*at];
			yCent = MySnapshot[3*at+1];
			zCent = MySnapshot[3*at+2];
		}
		else{
			flagcoord=1;
			while(flagcoord){
				OutList[nbounds]=polyhedras[ncent];
				CoordList[3*nbounds]=MySnapshot[3*OutList[nbounds]];
				CoordList[3*nbounds+1]=MySnapshot[3*OutList[nbounds]+1];
				CoordList[3*nbounds+2]=MySnapshot[3*OutList[nbounds]+2];
				nbounds++;
				ncent++;				
				if(polyhedras[ncent]==-1){
					flagcoord=0;
				}
			}
			
			for(int i=0 ; i<nbounds ; i++){
								

				for(int j=i+1 ; j<nbounds ; j++){
					Deltx1 = CoordList[3*i]-xCent;
					Delty1 = CoordList[3*i+1]-yCent;
					Deltz1 = CoordList[3*i+2]-zCent;

					Deltx2 = CoordList[3*j]-xCent;
					Delty2 = CoordList[3*j+1]-yCent;
					Deltz2 = CoordList[3*j+2]-zCent;

					deltxx = CoordList[3*i]-CoordList[3*j]; 					
					deltyy = CoordList[3*i+1]-CoordList[3*j+1]; 					
					deltzz = CoordList[3*i+2]-CoordList[3*j+2];					
					
					valx1=min3(pow(Deltx1,2),pow(Deltx1-acell0,2),pow(Deltx1+acell0,2));													
					valy1=min3(pow(Delty1,2),pow(Delty1-acell1,2),pow(Delty1+acell1,2));
					valz1=min3(pow(Deltz1,2),pow(Deltz1-acell2,2),pow(Deltz1+acell2,2));

					dist1=valx1+valy1+valz1;

					valx2=min3(pow(Deltx2,2),pow(Deltx2-acell0,2),pow(Deltx2+acell0,2));
					valy2=min3(pow(Delty2,2),pow(Delty2-acell1,2),pow(Delty2+acell1,2));
					valz2=min3(pow(Deltz2,2),pow(Deltz2-acell2,2),pow(Deltz2+acell2,2));

					dist2=valx2+valy2+valz2;
				
					valxx = min3(pow(deltxx,2),pow(deltxx-acell0,2),pow(deltxx+acell0,2));
					valyy = min3(pow(deltyy,2),pow(deltyy-acell1,2),pow(deltyy+acell1,2));
					valzz = min3(pow(deltzz,2),pow(deltzz-acell2,2),pow(deltzz+acell2,2));

					distij = valxx+valyy+valzz;
					
					AnglesList[(int)AnglesList[0]] = acos((dist1+dist2-distij)/2);
					AnglesList[0]++;
				}
			}
		}
	}

	return AnglesList;
} 

int atom_to_indexB(int atom,int CentMin,int CentMax,int OutMin,int OutMax){

	int place;
	if(CentMax<=OutMax){
		if(atom<=CentMax){
			place = atom - CentMin;
		}
		else{
			place = atom - OutMin + CentMax - CentMin + 1;
		}
	}	
	else{
		if(atom<=OutMax){
			place = atom - OutMin;
		}
		else{
			place = atom - CentMin + OutMax - OutMin + 1;
		}
	}
}

int atom_to_indexAt(int atom,int CentMin,int CentMax,int OutMin,int OutMax){

	int indice;
	if(atom<=CentMax && atom>=CentMin){	
		indice = atom-CentMin;	
	}
	else{
		indice = atom-OutMin + CentMax - CentMin+1;	
	}

	return indice;
}

void reinit(int *AllAts,int CentMin, int CentMax,int OutMin,int OutMax){

	if(CentMax!=OutMax){
		for(int i=0 ; i<(CentMax-CentMin+1) ; i++){
			AllAts[i]=CentMin+i;
		}
		for(int i=0 ; i<(OutMax-OutMin+1) ; i++){
			AllAts[CentMax-CentMin+1+i]=OutMin+i;
		}
	}
	else{
		for(int i=0 ; i<(CentMax-CentMin+1) ; i++){
			AllAts[i]=CentMin+i;
		}
	}
}


void clusteringall(int *neighbors,const int *SnapshotBonds, const int *indBonds,int *AllAtoms, const int CentMin, const int CentMax, const int OutMin, const int OutMax, const int atom){

	int indiceAt, at, indiceCoord;
	indiceAt = atom_to_indexB(atom,CentMin,CentMax,OutMin,OutMax);
	for(int i = indBonds[indiceAt]; i<indBonds[indiceAt+1]; i++){
		
		at = SnapshotBonds[i];
		indiceCoord = atom_to_indexAt(at,CentMin,CentMax,OutMin,OutMax);
		
		if(AllAtoms[indiceCoord]!=-1){
			AllAtoms[indiceCoord]=-1;
			neighbors[neighbors[0]]=at;
			neighbors[0]++;
			clusteringall(neighbors,SnapshotBonds,indBonds,AllAtoms,CentMin,CentMax,OutMin,OutMax,at);
		}
	}
}

void clusteringrec(int *neighbors,const int *SnapshotBonds, const int *indBonds,int *AllAtoms, const int CentMin, const int CentMax, const int OutMin, const int OutMax, const int atom, const int r){

	if(r>0){

		int indiceAt, at, indiceCoord;
		indiceAt = atom_to_indexB(atom,CentMin,CentMax,OutMin,OutMax);

		for(int i = indBonds[indiceAt]; i<indBonds[indiceAt+1]; i++){
		
			at = SnapshotBonds[i];
			indiceCoord = atom_to_indexAt(at,CentMin,CentMax,OutMin,OutMax);
		
			if(AllAtoms[indiceCoord]!=-1 && (at<=OutMax && at>=OutMin)){
				AllAtoms[indiceCoord]=-1;
				neighbors[neighbors[0]]=at;
				neighbors[0]++;
				clusteringrec(neighbors,SnapshotBonds,indBonds,AllAtoms,CentMin,CentMax,OutMin,OutMax,at,r-1);
			}
		}
	}
}

int* fullclustering(const int *SnapshotBonds, const int *indBonds, const int CentMin, const int CentMax, const int OutMin, const int OutMax, const int M, const int r){
	
	int *AllAts;
	int *neighbors;
	int *Neighbors;
	int index = 0;
	int nAts,Min1,Min2,Max1,Max2,atom,len;		

	if(CentMax==OutMax){
		nAts = 1+CentMax-CentMin ;
	}
	else{
		nAts = 2+CentMax+OutMax-CentMin-OutMin ;	
	}

	if(r>0){
		len = min(nAts*(pow(M,r)+1)+1,nAts*nAts+1);
	}	
	else{
		len = 2*nAts+1;					
	}

	neighbors = calloc(len,sizeof(int));
	AllAts = calloc(nAts,sizeof(int));		

//Initialisation

	reinit(AllAts,CentMin,CentMax,OutMin,OutMax);	

	for(int i=0 ; i<len ; i++){
		neighbors[i]=-1;
	}


//Neighbors computation
	neighbors[0]=1;

	if(r>0){
		while(index<CentMax-CentMin+1){
			if(AllAts[index] != -1){
				atom=AllAts[index];
				AllAts[index]=-1;
				neighbors[neighbors[0]]=atom;
				neighbors[0]++;
				clusteringrec(neighbors,SnapshotBonds,indBonds,AllAts,CentMin,CentMax,OutMin,OutMax,atom,r);
			}		
			neighbors[0]++;
			index++;
			reinit(AllAts,CentMin,CentMax,OutMin,OutMax);	
		}
	}
	else{
		while(index<nAts){
			if(AllAts[index] != -1){
				atom=AllAts[index];
				AllAts[index]=-1;
				neighbors[neighbors[0]]=atom;
				neighbors[0]++;
				clusteringall(neighbors,SnapshotBonds,indBonds,AllAts,CentMin,CentMax,OutMin,OutMax,atom);
				neighbors[0]++;
			}					
			index++;		
		}
		while(neighbors[neighbors[0]]==-1){
			neighbors[0]--;
		}

		neighbors[0]++;
		
		free(AllAts);
		return neighbors;
	}
	
	Neighbors = calloc(neighbors[0]+1,sizeof(int));
		
	int ind=1;
	len=0;

	for(int i=1 ; i<(neighbors[0]+1) ; i++){
		if(neighbors[i]==-1){
			if(len<=1){
				ind--;
				Neighbors[ind]=-1;
				len=0;
				
			}
			else{
				Neighbors[ind]=-1;
				ind++;
				len=0;
			}
		}
		else{				
			Neighbors[ind]=neighbors[i];
			ind++;
			len++;
		}
	}

	while(Neighbors[ind]==-1){
		ind--;
	}

	Neighbors[0]=ind+1;

	free(neighbors);
	free(AllAts);

	return Neighbors;

}
	
