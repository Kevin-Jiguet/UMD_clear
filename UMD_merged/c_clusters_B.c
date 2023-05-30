#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int min(int a,int b){

	if(a>b){return b;}
	else{return a;}

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
	
