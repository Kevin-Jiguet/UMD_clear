/*
 * Author: Kevin Jiguet-Covex
 * Date:   May 31, 2023
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <stdbool.h>

int* defineBonds(const char *lines,int len,int CentMin,int CentMax,int OutMin,int OutMax){//Returns BList, a tab listing all the coordinating atoms, concatenated with a the list of indexes that contains the information about which atome they pertain to.
											  //This only takes into account the atoms whose tag number is between CentMin and CentMax or between OutMin and Outmax.
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

		if(isdigit(lines[i])){number = number*10; number = number + (lines[i]-'0');}//Retrieving number (tag of atom) from string of char

		else if(lines[i] == '\t'){
			if(newline||ligand){
				if((number<=CentMax && number>=CentMin)||(number<=OutMax && number>=OutMin)){					
					if(newline){
						newline = false;
						ligand = true;
					}		
					else{	
						BondsList[BondsList[0]]=number; //Adding bound atom to the list
						BondsList[0]++;//Incrementing total number of bound atoms
						nAts++;//Counting atoms
					}

	
				}	
				else{newline = false ; number=0;}
			}
			number=0;
		}
		
		else if(lines[i] == '\n'){

			if(nAts>nMax){nMax=nAts;}

			if(ligand){
			BondIndexes[BondIndexes[0]+1]=BondIndexes[BondIndexes[0]]+nAts;//Incrementing by how many atoms are bound to the central atom 
			BondIndexes[0]++;//Incrementing the number of central atoms
			}
			//Reinitializing for the next central atom loop
			nAts=0;
			ligand = false;
			newline = true;
			number=0;
		}
	}

	
	BondIndexes[BondIndexes[0]+1]=nMax;
	BondIndexes[0]++;

	int *BList;//Tab that will contain every data in an unique row

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
