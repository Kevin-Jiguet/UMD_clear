/*
 * Author: Kevin Jiguet-Covex
 * Date:   May 31, 2023
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <stdbool.h>

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

	int depth = 0;
	int atom = 0;
	int flagxcart = 0;
	double number = 0;
	double* Coordinates;
	int lineIndex=0;//index on the line of the number we are looking at
	Coordinates = calloc(3*nAtoms,sizeof(double));
	for(int i=0 ; i<len ; i++){

		if(MySnapshot[i]==' '){	
			if(flagxcart){
				Coordinates[3*atom+lineIndex-3]=number;
				flagxcart=0;
			}
			number=0;
			lineIndex++;
			depth=0;//The next digits will be in the integer part
		}

		else if(lineIndex==3||lineIndex==4||lineIndex==5){//If we are on the xcart coordinates
			flagxcart = 1;
			if(isdigit(MySnapshot[i])){
				if(depth == 0){
					number=number*10;
					number=number+(MySnapshot[i]-'0');
				}
			
				else{
					number=number+(double)(MySnapshot[i]-'0')/depth;
					depth=depth*10;

				}
			}
			else if(MySnapshot[i]=='.'){
				depth=10;//We're in the decimal part now
			}
			

		}
		else if(MySnapshot[i]=='\n'){
			atom++;//Next line, next atom ; we reinitialize the variables
			lineIndex=0;
			depth=0;
			number=0;
		}
	}
return Coordinates;
}

//Calculates the angles within the polyhedras if it contains at least 2 coordinated atoms
double* angles(const int *polyhedras, int nPol, int M, const double *MySnapshot, int CentMin, int CentMax, int OutMin, int OutMax, double acell0, double acell1, double acell2){
	
	double* AnglesList;
	int* OutList;
	double* CoordList;
	int at,ne1,ne2;	
	double xCent,yCent,zCent;
	double Deltx1,Delty1,Deltz1,Deltx2,Delty2,Deltz2,valx1,valy1,valz1,valx2,valy2,valz2,deltzz,deltxx,deltyy,valxx,valyy,valzz,dist1,dist2,distij,angle;
	double x1,y1,z1,x2,y2,z2;
	int centat=1;
	int nbounds=0 ;
	int ncent=0;
	int index=0;
	int flagcoord=0;

	CoordList = calloc(3*M,sizeof(double));
	OutList = calloc(M,sizeof(int));
	AnglesList = calloc((OutMax-OutMin+1)*(OutMax-OutMin)/2+1+CentMax-CentMin+1,sizeof(double));
	AnglesList[0]=1;



	while(ncent<nPol){
		
		index++;
		at = polyhedras[index];
		
		
		if(at==-1){
			centat = 1;//means that the next atom on the polyhedras list will be a central atom
			nbounds = 0;
			ncent++;
			AnglesList[(int)AnglesList[0]]=-1;//separation between clusters
			AnglesList[0]++;

		}		
		
		else if(centat==1){//we define the coordinates of the central atom
			centat = 0;//The next atom on the polyhedras list will be an adjacent atom
			xCent = MySnapshot[3*at];
			yCent = MySnapshot[3*at+1];
			zCent = MySnapshot[3*at+2];
		}

		else {
			flagcoord=1;
			while(flagcoord){//while we are looking at the same cluster
				at = polyhedras[index];
				if(at!=-1){
					CoordList[3*nbounds]=MySnapshot[3*at];//We fill the CoordList with the coordinates of outer atoms
					CoordList[3*nbounds+1]=MySnapshot[3*at+1];
					CoordList[3*nbounds+2]=MySnapshot[3*at+2];
					nbounds++;
					index++;
				}
				else{
					flagcoord=0;//Means that we're out of the cluster now
					index--;
				}
			}
			for(int i=0 ; i<nbounds ; i++){//We browse the CoordList to compute the angle for each pair of outer (adjacent) atoms

				x1 = CoordList[3*i];
				y1 = CoordList[3*i+1];
				z1 = CoordList[3*i+2];
								
				for(int j=i+1 ; j<nbounds ; j++){

					x2 = CoordList[3*j];
					y2 = CoordList[3*j+1];
					z2 = CoordList[3*j+2];

					Deltx1 = x1-xCent;
					Delty1 = y1-yCent;
					Deltz1 = z1-zCent;

					Deltx2 = x2-xCent;
					Delty2 = y2-yCent;
					Deltz2 = z2-zCent;

					deltxx = x1-x2; 					
					deltyy = y1-y2; 					
					deltzz = z1-z2;					
					
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

					
					if(dist1==0||dist2==0){
					AnglesList[(int)AnglesList[0]]=0.0;
					}
					else{
					AnglesList[(int)AnglesList[0]] = acos((dist1+dist2-distij)/(2*sqrt(dist1*dist2)))*180/3.141592653589793;
					}
					AnglesList[0]++;
				}
			}
		}
	}
	return AnglesList;
} int atom_to_indexB(int atom,int CentMin,int CentMax,int OutMin,int OutMax){//Gives, in the BondIndexes tab, the index of the data pertaining to the atom.

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
	
