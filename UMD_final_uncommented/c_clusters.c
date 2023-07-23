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


void reinit(int *Tab,int *Indexes,int n){
	int imax=0;
	int imin=0;
	for(int el = 0 ; el<n ; el++){
		imin = Indexes[2*el];
		imax = Indexes[2*el+1];
		for(int i = imin ; i<imax+1 ; i++){
			Tab[i] = i;
		}	
	}
}

bool is_in(int atom,int *Indexes, int n){
	bool result = false;
	for(int i=0 ; i<n ; i++){
		if((atom>=Indexes[2*i] && atom<=Indexes[2*i+1])){
			result = true;
		}
	}
	return result;
}

int num_atoms(int* Indexes, int n){
	int num = 0;
	for(int i=0 ; i<n ; i++){
		num = num + Indexes[2*i+1]-Indexes[2*i] + 1;
	}
	return num;
}

double min3(double a, double b, double c){//The minimum between 3 numbers.

	double minimum = a;

	if(minimum>b){minimum=b;}
	if(minimum>c){minimum=c;}

	return minimum;

}


double* angles(int *polyhedras, int nPol, int M, double *MySnapshot, int *CentIndexes, int nCent, int *OutIndexes, int nOut, double acell0, double acell1, double acell2){
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

	int numCent = num_atoms(CentIndexes,nCent);
	int numOut = num_atoms(OutIndexes,nOut);

	AnglesList = calloc(((M+1)*M/2*numCent+numCent),sizeof(double));

	if(CoordList == NULL || OutList == NULL || AnglesList == NULL){
		printf("Memory Allocation failure for CoordList (length %d), OutList (length %d) or AnglesList (length %d) in function 'angles'",3*M,M,(M+1)*M/2*numCent+numCent);
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

int atom_to_index(int atom,int *CentIndexes, int *OutIndexes, int nCent, int nOut){

	int place;



	return place;
}



void clusteringall(int *neighbors,const int *SnapshotBonds, const int *indBonds,int *AllAtoms, const int atom){

	int at;

	for(int i = indBonds[atom]; i<indBonds[atom+1]; i++){
		at = SnapshotBonds[i];
		if(AllAtoms[at]!=-1){
			AllAtoms[at]=-1;
			neighbors[neighbors[0]]=at;
			neighbors[0]++;
			clusteringall(neighbors,SnapshotBonds,indBonds,AllAtoms,at);
		}
	}
}

void clusteringrec(int *neighbors,const int *SnapshotBonds, const int *indBonds, int *OutAts, int *OutIndexes, const int nOut, const int atom, const int r){

	int at=0;
	if(r>0){
//		printf("\n The atom %d has %d links ; the outer are : \n",atom,indBonds[indiceAt+1]-indBonds[indiceAt]);
		
		for(int x=0 ; x<r ; x++){printf("\t");}
		for(int i = indBonds[atom]; i<indBonds[atom+1]; i++){
			at = SnapshotBonds[i];
			if(is_in(at,OutIndexes,nOut)){
//				printf("%d (root %d) ",at,atom);
				if(OutAts[at]!=-1){
					OutAts[at]=-1;
					neighbors[neighbors[0]]=at;
					neighbors[0]++;
				}
				clusteringrec(neighbors,SnapshotBonds,indBonds,OutAts,OutIndexes,nOut,at,r-1);
				
			}
		}
	}
}

int* fullclustering(const int *SnapshotBonds, const int *indBonds, const int natom, const int nAts, int *CentIndexes, int *OutIndexes, int nCent, int nOut, const int M, const int r){
	
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
	int prev_index=0;

	int numCent = num_atoms(CentIndexes,nCent);
	int numOut = num_atoms(OutIndexes,nCent);

	int meanClustSize = indBonds[nAts]/(nAts);


	if(r>0){
		maxClustSize = min(pow(M,r)+1,numOut+2);
		len = indBonds[nAts]+2*nAts+1 ;
		CentAts = calloc(natom,sizeof(int));
		OutAts = calloc(natom,sizeof(int));
	

		if (CentAts == NULL||OutAts == NULL)
		{
		    printf("Memory Allocation failure for CentAts (length %d) or Outats (length %d)\n,",natom,natom);
		    return EXIT_FAILURE;
		}

		for(int i=0 ; i<natom ; i++){
			CentAts[i]=-1;
			OutAts[i]=-1;
		}

		reinit(CentAts,CentIndexes,nCent);
		reinit(OutAts,OutIndexes,nOut);
	}	
	else{
		len = 2*nAts+2;					
		
		AllAts = calloc(natom,sizeof(int));
		if (AllAts == NULL)
		{
		    printf("Memory Allocation failure for AllAts : len is %d\n",len);
		    return EXIT_FAILURE;
		}
		for(int i=0 ; i<natom ; i++){
			AllAts[i]=-1;
		}
		reinit(AllAts,CentIndexes,nCent);
		reinit(AllAts,OutIndexes,nOut);
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
		while(index<natom){
			if(CentAts[index] != -1){
				atom=CentAts[index];
				CentAts[index]=-1;
				neighbors[neighbors[0]]=atom;
				neighbors[0]++;
				prev_index = neighbors[0];
				clusteringrec(neighbors,SnapshotBonds,indBonds,OutAts,OutIndexes,nOut,atom,r);
				if(neighbors[0]==prev_index){
					neighbors[0]--;
					neighbors[neighbors[0]]=-1;
				}
				else{
					neighbors[0]++;
				}

			}		


			index++;
			reinit(OutAts,OutIndexes,nOut);
			if(neighbors[0]>len-2*maxClustSize){
//				printf("about to reallocate ! %d to %d... ",len,len+2*maxClustSize+1);
				
				int *neighborstemp = (int *)realloc(neighbors,(len+2*maxClustSize+1)*sizeof(int));

				if(neighborstemp == NULL){printf("Memory allocation failure for neighbors \n");return EXIT_FAILURE;}

				neighbors = neighborstemp;
				len = len + 2*maxClustSize +1;
				for(int k=neighbors[0] ; k<len ; k++){
					neighbors[k] = -1;
				}													
			}
		}
	free(OutAts);
	free(CentAts);
	while(neighbors[neighbors[0]]==-1){
		neighbors[0]--;
	}
	return neighbors;
	}
	else{	
		while(index<natom){
			if(AllAts[index] != -1){
				atom=AllAts[index];
				AllAts[index]=-1;
				neighbors[neighbors[0]]=atom;
				neighbors[0]++;
				clusteringall(neighbors,SnapshotBonds,indBonds,AllAts,atom);
				neighbors[0]++;
			}					
			index++;
		}
		while(neighbors[neighbors[0]]==-1){
			neighbors[0]--;
		}
		
		free(AllAts);
		free(OutAts);
		free(CentAts);
		return neighbors;
	}
}
