#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <stdbool.h>
#include <string.h>


void free_memory_int(int *tab){
	free(tab);
}


void free_memory_double(double *tab){
	free(tab);
}


int min(int a,int b){//returns the minimum of 2 numbers, except if one of them is negative (it will return the other).

	if(a>b||a<0){return b;}
	else{return a;}

}	

double min3(double a, double b, double c){//The minimum between 3 numbers.

	double minimum = a;

	if(minimum>b){minimum=b;}
	if(minimum>c){minimum=c;}

	return minimum;

}


double* angles(int *polyhedras, int nPol, int M, double *MySnapshot, int CentMin, int CentMax, int OutMin, int OutMax, double acell0, double acell1, double acell2){
	double* AnglesList = NULL;
	int* OutList = NULL;
	double* CoordList = NULL;
	int at,ne1,ne2;	
	double xCent,yCent,zCent,Deltx1,Delty1,Deltz1,Deltx2,Delty2,Deltz2,valx1,valy1,valz1,valx2,valy2,valz2,deltzz,deltxx,deltyy,valxx,valyy,valzz,dist1,dist2,distij,angle;
	int centat=1;
	int nbounds=0 ;
	int ncent=0;
	int indexPol=0;
	int flagcoord=0;

	CoordList = calloc(3*M,sizeof(double));
	OutList = calloc(M,sizeof(int));
	AnglesList = calloc(((OutMax-OutMin+1)*(OutMax-OutMin)/2+1+CentMax-CentMin+1),sizeof(double));

	if(CoordList == NULL || OutList == NULL || AnglesList == NULL){
		printf("Memory Allocation failure for CoordList, OutList or AnglesList in function 'angles'");
		return EXIT_FAILURE;
	}
	AnglesList[0]=1;
	while(ncent<=nPol){

		indexPol++;
		at = polyhedras[indexPol];
		if(at==-1){
			AnglesList[(int)AnglesList[0]]=-1;
			AnglesList[0]++;
			centat=1;
			ncent++;
			nbounds = 0;
		}		
		
		else if(centat==1){
			centat=0;
			//printf("central atom %d\n",at);
			xCent = MySnapshot[3*at];
			yCent = MySnapshot[3*at+1];
			zCent = MySnapshot[3*at+2];

		}
		else{
			flagcoord=1;
			while(flagcoord){
				at = polyhedras[indexPol];				
				OutList[nbounds]=at;
				//printf("bounds : %d on %d   ",nbounds,M);
				CoordList[3*nbounds]=MySnapshot[3*at];
				CoordList[3*nbounds+1]=MySnapshot[3*at+1];
				CoordList[3*nbounds+2]=MySnapshot[3*at+2];
				//printf("CoordLists : %f %f %f (%d) (at %d)\n",CoordList[3*nbounds],CoordList[3*nbounds+1],CoordList[3*nbounds+2],nbounds,OutList[nbounds]);
				nbounds++;
//				printf("atom %d, nbounds %d\n",at,nbounds);
				indexPol++;				
				if(polyhedras[indexPol]==-1){
					flagcoord=0;
					indexPol--;
				}

			}
			
			for(int i=0 ; i<nbounds ; i++){

				Deltx1 = CoordList[3*i]-xCent;
				Delty1 = CoordList[3*i+1]-yCent;
				Deltz1 = CoordList[3*i+2]-zCent;
				//printf("Delt1 = %f %f %f, i= %d\n",Deltx1,Delty1,Deltz1, i);							
				for(int j=i+1 ; j<nbounds ; j++){
					

					Deltx2 = CoordList[3*j]-xCent;
					Delty2 = CoordList[3*j+1]-yCent;
					Deltz2 = CoordList[3*j+2]-zCent;
//					printf("Delt2 = %f %f %f, j= %d\n",Deltx2,Delty2,Deltz2, j);


					deltxx = CoordList[3*i]-CoordList[3*j]; 					
					deltyy = CoordList[3*i+1]-CoordList[3*j+1]; 					
					deltzz = CoordList[3*i+2]-CoordList[3*j+2];					
//					printf("deltxx = %f %f %f\n",deltxx,deltyy,deltzz);

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
//					printf("ValXX = %f %f %f\n",valxx,valyy,valzz);

					distij = valxx+valyy+valzz;
//					printf("Valeurs : %f  %f  %f, cos = %f\n",distij,dist1,dist2,(dist1+dist2-distij)/(2*sqrt(dist1*dist2)));
					AnglesList[(int)AnglesList[0]] = acos((dist1+dist2-distij)/(2*sqrt(dist1*dist2)))/3.141592653589793*180;
					AnglesList[0]++;
				}
			}
		}
	}

	free(OutList);
	free(CoordList);
	return AnglesList;
} 

int atom_to_index(int atom,int CentMin,int CentMax,int OutMin,int OutMax){

	int place;
	if((CentMax-OutMax)*(CentMin-OutMin) > 0){
		if(CentMax<OutMax){
			if(atom<=CentMax){
				place = atom - CentMin;
			}
			else{
				place = atom - OutMin + CentMax - CentMin + 1;
			}
		}	
		else if (CentMax>OutMax){
			if(atom<=OutMax){
				place = atom - OutMin;
			}
			else{
				place = atom - CentMin + OutMax - OutMin + 1;
			}	
		}
	}
	else{
		place = atom - min(OutMin,CentMin);
	}
//	printf("place of %d = %d ; CM = %d ; Cm = %d ; OM = %d ; Om = %d",atom,place,CentMax,CentMin,OutMax,OutMin);
	return place;
}

void reinitAts(int *CentAts, int CentMin, int CentMax, int * OutAts, int OutMin, int OutMax){
	for(int i=0 ; i<CentMax-CentMin+1 ; i++){
		CentAts[i] = CentMin + i;
	}

	for(int i=0 ; i<OutMax-OutMin+1 ; i++){
		OutAts[i] = OutMin + i;
	}
}


void reinit(int *AllAts,int CentMin, int CentMax,int OutMin,int OutMax,int imax){


	for(int i = 0 ; i<imax ; i++){
		if((CentMax >= i && CentMin <=i)||(OutMax >= i && OutMin <= i)){
			AllAts[i] = i;
		} 
		else{
			AllAts[i]=-1;
		}
	}
}


void clusteringall(int *neighbors,const int *SnapshotBonds, const int *indBonds,int *AllAtoms, const int CentMin, const int CentMax, const int OutMin, const int OutMax, const int atom){

	int indiceAt, at, indiceCoord;
	indiceAt = atom_to_index(atom,CentMin,CentMax,OutMin,OutMax);
	for(int i = indBonds[indiceAt]; i<indBonds[indiceAt+1]; i++){
		at = SnapshotBonds[i];
		
		if(AllAtoms[at]!=-1){
			AllAtoms[at]=-1;
			neighbors[neighbors[0]]=at;
			neighbors[0]++;
			clusteringall(neighbors,SnapshotBonds,indBonds,AllAtoms,CentMin,CentMax,OutMin,OutMax,at);
		}
	}
}

void clusteringrec(int *neighbors,const int *SnapshotBonds, const int *indBonds,int *CentAts, const int CentMin, const int CentMax, int *OutAts, const int OutMin, const int OutMax, const int atom, const int r){

	if(r>0){

		int indiceAt, at, indiceCoord;
		indiceAt = atom_to_index(atom,CentMin,CentMax,OutMin,OutMax);

//		printf("\n The atom %d has %d links ; the outer are : \n",atom,indBonds[indiceAt+1]-indBonds[indiceAt]);
		
		for(int x=0 ; x<r ; x++){printf("\t");}
		for(int i = indBonds[indiceAt]; i<indBonds[indiceAt+1]; i++){
			at = SnapshotBonds[i];
			if(at<=OutMax && at>=OutMin){
//				printf("%d (root %d) ",at,atom);
				indiceCoord = at - OutMin ;
				if(OutAts[indiceCoord]!=-1){
					OutAts[indiceCoord]=-1;
					neighbors[neighbors[0]]=at;
					neighbors[0]++;
				}
				clusteringrec(neighbors,SnapshotBonds,indBonds,CentAts,CentMin,CentMax,OutAts,OutMin,OutMax,at,r-1);
				
			}
		}
	}
}

int* fullclustering(const int *SnapshotBonds, const int *indBonds, const int nAts, const int CentMin, const int CentMax, const int OutMin, const int OutMax, const int M, const int r){
	
	int *AllAts = NULL;
	int *neighbors = NULL;
	int *Neighbors = NULL;
	int *CentAts = NULL;
	int *OutAts = NULL;

	int maxClustSize = 0;
	int index = 0;
	int atom;
	int len;
	int imax;

	int meanClustSize = indBonds[nAts]/nAts;


	if(r>0){
		maxClustSize = min(pow(M,r)+1,OutMax-OutMin+2);
		len = indBonds[nAts]+2*nAts+1 ;
		CentAts = calloc(CentMax-CentMin+1,sizeof(int));
		OutAts = calloc(OutMax-OutMin+1,sizeof(int));

		if (CentAts == NULL||OutAts == NULL)
		{
		    printf("Memory Allocation failure for CentAts (length %d) or Outats (length %d)\n,",CentMax-CentMin+1,OutMax-OutMin+1);
		    return EXIT_FAILURE;
		}
		reinitAts(CentAts,CentMin,CentMax,OutAts,OutMin,OutMax);
	}	
	else{
		len = 2*nAts+1;					
		
		if(CentMax>OutMax){
			imax = CentMax;
		}
		else{
			imax = OutMax;
		}
	
		AllAts = calloc(imax,sizeof(int));
		if (AllAts == NULL)
		{
		    printf("Memory Allocation failure for AllAts : len is %d\n",len);
		    return EXIT_FAILURE;
		}
		reinit(AllAts,CentMin,CentMax,OutMin,OutMax,imax);	


	}

	neighbors = calloc(len,sizeof(int));
	if (neighbors == NULL)
	{
	    printf("Memory Allocation failure for neighbors : len is %d\n",len);
	    return EXIT_FAILURE;
	}


//Initialisation

	for(int i=0 ; i<len ; i++){
		neighbors[i]=-1;
	}
	neighbors[0]=1;



//Neighbors computation

	if(r>0){
		while(index<CentMax-CentMin+1){
			if(CentAts[index] != -1){
				atom=CentAts[index];
				CentAts[index]=-1;
				neighbors[neighbors[0]]=atom;
				neighbors[0]++;
				clusteringrec(neighbors,SnapshotBonds,indBonds,CentAts,CentMin,CentMax,OutAts,OutMin,OutMax,atom,r);
			}		
			neighbors[0]++;
			index++;
			reinitAts(CentAts,CentMin,CentMax,OutAts,OutMin,OutMax);	
			if(neighbors[0]>len-2*maxClustSize){
//				printf("about to reallocate ! %d to %d... ",len,len+2*maxClustSize+1);
				
				int *neighborstemp = (int *)realloc(neighbors,(len+2*maxClustSize+1)*sizeof(int));
				if(neighborstemp == NULL){
					printf("Memory allocation failure for neighbors \n");
					return EXIT_FAILURE;
				}
				neighbors = neighborstemp;
				len = len + 2*maxClustSize +1;
				for(int k=neighbors[0] ; k<len ; k++){
					neighbors[k] = -1;
				}													
			}
		}
	}
	else{	
		while(index<imax){
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
		free(OutAts);
		free(CentAts);
		return neighbors ;
	}
	
	Neighbors = calloc(neighbors[0]+2,sizeof(int));
	for(int i= 0 ; i<neighbors[0]+2 ; i++){
		Neighbors[i]=-1;
	}

        if (Neighbors == NULL)
        { 
              printf("Failure of Neighbors memory allocation\n");
              return EXIT_FAILURE;
        }
		
	int ind=1;
	len=0;

	for(int i=1 ; i<(neighbors[0]+1) ; i++){
		if(neighbors[i]==-1){
			if(len==1){
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

	ind = neighbors[0]+1;
	while(Neighbors[ind]==-1){
		ind--;
	}

	Neighbors[0]=ind+1;

	free(neighbors);
	free(AllAts);
	free(OutAts);
	free(CentAts);

	return Neighbors;
}
