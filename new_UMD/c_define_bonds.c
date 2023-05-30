#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <stdbool.h>

int* defineBonds(const char *lines,int len,int CentMin,int CentMax,int OutMin,int OutMax){
	
	int number=0;
	int nAts=0;
	int nMax=0;
	bool ligand = false ;
	bool newline = true ;
	int *BondsList;
	int *BondIndexes;

	BondsList = calloc(len,sizeof(int));
	BondIndexes = calloc((CentMax+OutMax-CentMin-OutMin+2+3),sizeof(int));
	
	BondIndexes[0]=1;
	BondIndexes[1]=0;
	BondsList[0]=1;
	
	

	for(int i=0 ; i<len ; i++){


		if(isdigit(lines[i])){number = number*10; number = number + (lines[i]-'0');}

		else if(lines[i] == '\t'){
			if(newline||ligand){
				if((number<=CentMax && number>=CentMin)||(number<=OutMax && number>=OutMin)){					
					if(newline){
						newline = false;
						ligand = true;
					}		
					else{	
						BondsList[BondsList[0]]=number; 
						BondsList[0]++;
						nAts++;
					}

	
				}	
				else{newline = false ; number=0;}
			}
			number=0;
		}
		
		else if(lines[i] == '\n'){

			if(nAts>nMax){nMax=nAts;}

			if(ligand){
			BondIndexes[BondIndexes[0]+1]=BondIndexes[BondIndexes[0]]+nAts;
			BondIndexes[0]++;
			}

			nAts=0;
			ligand = false;
			newline = true;
			number=0;
		}
	}

	
	BondIndexes[BondIndexes[0]+1]=nMax;
	BondIndexes[0]++;

	int *BList;

	BList = malloc((BondsList[0]+BondIndexes[0]+2)*sizeof(int));

	for(int i=1 ; i<BondsList[0]+1 ; i++){
		BList[i+1]=BondsList[i];
	}

	for(int i=1 ; i<BondIndexes[0]+1 ; i++){
		BList[BondsList[0]+1+i]=BondIndexes[i];
	}
	
	BList[0]=BondsList[0];
	BList[1]=BondIndexes[0];

	free(BondsList);
	free(BondIndexes);

	return BList ;
}
