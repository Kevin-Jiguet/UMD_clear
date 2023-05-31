#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int min(int a,int b){

	if(a>b){return b;}
	else{return a;}

}	

int atom_to_indexB(int atom,int CentMin,int CentMax,int OutMin,int OutMax){//Gives, in the BondIndexes tab, the index of the data pertaining to the atom.

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

int atom_to_indexAt(int atom,int CentMin,int CentMax,int OutMin,int OutMax){//Gives, in the AllAtoms tab, the index of the atom.

	int indice;
	if(atom<=CentMax && atom>=CentMin){	
		indice = atom-CentMin;	
	}
	else{
		indice = atom-OutMin + CentMax - CentMin+1;	
	}

	return indice;
}

void reinit(int *AllAts,int CentMin, int CentMax,int OutMin,int OutMax){//Reinitializes (refills) the AllAts list with all the atoms.

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


void clusteringall(int *neighbors,const int *SnapshotBonds, const int *indBonds,int *AllAtoms, const int CentMin, const int CentMax, const int OutMin, const int OutMax, const int atom){//fills the neighbors tab with the atoms of each cluster. Polymerization.

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

void clusteringrec(int *neighbors,const int *SnapshotBonds, const int *indBonds,int *AllAtoms, const int CentMin, const int CentMax, const int OutMin, const int OutMax, const int atom, const int r){//fills the neighbors tab with the atoms of each cluster. Coordination.

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

int* fullclustering(const int *SnapshotBonds, const int *indBonds, const int CentMin, const int CentMax, const int OutMin, const int OutMax, const int M, const int r){//returns the tab Neighbors, containing all the clusters, each separated by the value "-1".
	
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

//Memory allocation
	
	neighbors = calloc(len,sizeof(int));
	AllAts = calloc(nAts,sizeof(int));		

//Initialisation

	reinit(AllAts,CentMin,CentMax,OutMin,OutMax);	

	for(int i=0 ; i<len ; i++){
		neighbors[i]=-1;
	}


//Neighbors computation
	neighbors[0]=1;

	if(r>0){//Coordination
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
			reinit(AllAts,CentMin,CentMax,OutMin,OutMax);//Reinitialization because the same coordinating atom can be in different clusters
		}
	}
	else{//Polymerization
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
		}//No reinitialization because the same atom cannot be in two different polymers
		while(neighbors[neighbors[0]]==-1){
			neighbors[0]--;
		}

		neighbors[0]++;
		
		free(AllAts);
		return neighbors;
	}
	
	Neighbors = calloc(neighbors[0]+1,sizeof(int));//More condensed data
		
	int ind=1;
	len=0;
//Putting only the relevant data in the Neighbors tab
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
	
