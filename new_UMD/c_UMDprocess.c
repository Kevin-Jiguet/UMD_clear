#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <stdbool.h>

//This script is intended to work in tandem with the python script umd_processes_fast.py
//It's main purpose is to read a snapshot and to extract the wanted information and only this

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

double* read_umd_values(char *MySnapshot, int nAtoms, int len, int X){//reads some specifics values of a snapshot and returns it in a tab

	int depth = 0;
	int atom = 0;
	int flagcoord = 0;//If we're on the right place on the line to read the values
	int flagexp = 0;//If we're reading an exponent
	int sign = 1;//Takes in account the fact that the value can be negative
	int signexp =1;//And the exponent too
	double number = 0;
	int exp = 0;
	double* Values;
	int lineIndex=0;//index on the line of the number we are looking at
	Values = calloc(3*nAtoms,sizeof(double));
	for(int i=0 ; i<len ; i++){

		if(MySnapshot[i]==' '){	
			if(flagcoord){
				Values[3*atom+lineIndex-X]=number*sign*pow(10,exp*signexp);
				flagcoord=0;
			}
			exp=0;
			flagexp=0;
			number=0;
			lineIndex++;
			sign = 1;
			signexp = 1;
			depth=0;//The next digits will be in the integer part
		}
		else if(MySnapshot[i]=='-'){
			if(flagexp){
				signexp=-1;
			}
			else{
				sign = -1;
			}
		}
		else if(lineIndex==X||lineIndex==X+1||lineIndex==X+2){//If we are on the vels coordinates
			flagcoord = 1;
			if(isdigit(MySnapshot[i])){
				if(flagexp){
					exp = exp*10;
					exp = exp + (MySnapshot[i]-'0');
				}
				else if(depth == 0){
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
			else if(MySnapshot[i]=='e'){
				exp=0;
				flagexp = 1;//From now on the digits are the exponent's
			}
			

		}
		else if(MySnapshot[i]=='\n'){
			atom++;
			lineIndex=0;
			exp=0;
			sign=1;
			signexp=1;
			depth=0;
			number=0;
			flagexp=0;
			flagcoord = 0;
		}
	}
return Values;
}


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
