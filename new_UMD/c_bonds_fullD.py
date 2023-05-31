#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double min(double a, double b, double c){
	double m=a;
	if(b<m){
		m=b;
	}
	if(c<m){
		m=c;
	}
return m;
}

int maxval(const int tab[], int L){

    int val;
    val=tab[0];
    for(int i=0 ; i<L ; i++){
        if(val<tab[i]){
            val=tab[i];
        }
    }
return val;
}

int* compute_fBonds(const double *MySnapshot, const double *BondTable, const int *CrystalTypes, const int nAtoms, const int ntypat, const double acell0, const double acell1, const double acell2, const int nCells){//returns the Bonds tab, containing the bound atoms. Two successive atoms are bound to the same central atom. The value "-1" is used to separate each central atom from the next one.

    int* DicoAtoms;
    int* nMesh;
    int* Mesh;
    int* MeshIndex;
    int* Bonding;
    int* nBonding;

    double x,y,z;
    int nX,nY,nZ,cellnum;

    double Lx = acell0/nCells;
    double Ly = acell1/nCells;
    double Lz = acell2/nCells;
    int ntotAts=0;
    int numCoord=0;
    int indexCoord=0;
    int index = 0;
    int cx,cy,cz;    
    int at,ne;
    double dx,dy,dz,valx,valy,valz,distij;


//Allocating memory
    DicoAtoms = calloc((3*nAtoms+1),sizeof(int));
    nMesh = calloc((nCells*nCells*nCells),sizeof(int));
    MeshIndex = calloc((nCells*nCells*nCells),sizeof(int));
    nBonding = calloc(nAtoms , sizeof(long int));

    for (int i=0 ; i<nAtoms ; i=i+1){//Determining to which specific sub-cell each atom belongs

        x = MySnapshot[3*i];
        y = MySnapshot[3*i+1];
        z = MySnapshot[3*i+2];

        nX = (int)floor(x/Lx);
        nY = (int)floor(y/Ly);
        nZ = (int)floor(z/Lz);

        DicoAtoms[3*i] = nX;
        DicoAtoms[3*i+1] = nY;
        DicoAtoms[3*i+2] = nZ;

        nMesh[nX*nCells*nCells+nY*nCells+nZ]++;

    }


    int M = maxval(nMesh,nCells*nCells*nCells);//maximum number of atoms in any cell ; used to determine the size of the tabs

//Allocating memory
    Mesh = calloc((nCells*nCells*nCells*M),sizeof(int));
    Bonding = calloc((nAtoms*M*27+1),sizeof(int));
    if(Mesh==NULL){fprintf(stderr,"Mem alloc failure");exit(1);}

    int Ats[M];//table which will contain the atoms in a central cell
    int Coord[M*27];//table which will contain the coordinating atoms in the 27 (3^3) surrounding cells
    
    Bonding[0]=M;

    for(int i = 1; i<nAtoms*M*27+1 ; i++){
	Bonding[i]=-1;
    }

    for(int iat = 0 ; iat<nAtoms ; iat++){

        nX = DicoAtoms[3*iat];
        nY = DicoAtoms[3*iat+1];
        nZ = DicoAtoms[3*iat+2];
    

        cellnum = nX*nCells*nCells+nY*nCells+nZ;
        index = MeshIndex[cellnum];
	if(cellnum>-1 && cellnum<nCells*nCells*nCells){Mesh[M*cellnum+index] = iat;}

	Mesh[M*cellnum+index] = iat;
	MeshIndex[cellnum]++;
    }
//Filling the Bonding tab
    for(int iat=0 ; iat<nAtoms ; iat++) {

        numCoord = 0;

        if(DicoAtoms[3*iat]!=-1){
            
	    nX = DicoAtoms[3*iat];
            nY = DicoAtoms[3*iat+1];
            nZ = DicoAtoms[3*iat+2];
            cellnum = nX*nCells*nCells+nY*nCells+nZ;
            index = MeshIndex[cellnum];
            
	    for(int iatom=0 ; iatom<index ; iatom++){
               
	        Ats[iatom]=Mesh[M*cellnum+iatom];
		DicoAtoms[3*Ats[iatom]]=-1;//Marking that we will have already calculated the bonds for this specific atom
            }
//Filling the Coord tab
            for(int i=-1 ; i<2 ; i++){
                    cx = (nX+i+nCells)%nCells;
                    for(int j=-1 ; j<2 ; j++){
                        cy = (nY+j+nCells)%nCells;
                        for(int k=-1 ; k<2 ; k++){
                            cz = (nZ+k+nCells)%nCells;
                            cellnum = cx*nCells*nCells+cy*nCells+cz;
                            indexCoord = MeshIndex[cellnum];
                         
			    for(int jatom = 0 ; jatom<indexCoord ; jatom++){
                                Coord[numCoord]=Mesh[M*cellnum+jatom];
                                numCoord++;
			
			       }
                            }
                        }
                    }


                 for(int iatom = 0 ; iatom<index ; iatom++){

                    at = Ats[iatom];
		    
                    for(int jatom = 0 ; jatom<numCoord ; jatom++){

                        ne = Coord[jatom];

                        if(at<ne){//Determining whether at and ne are bound together

                            dx = MySnapshot[at*3]-MySnapshot[ne*3];
                            dy = MySnapshot[at*3+1]-MySnapshot[ne*3+1];
                            dz = MySnapshot[at*3+2]-MySnapshot[ne*3+2];
                            valx = min(pow(dx,2),pow(dx-acell0,2),pow(dx+acell0,2));
                            valy = min(pow(dy,2),pow(dy-acell1,2),pow(dy+acell1,2));
                            valz = min(pow(dz,2),pow(dz-acell1,2),pow(dz+acell2,2));
                            distij = valx+valy+valz;

                            if(distij<BondTable[CrystalTypes[at]*ntypat+CrystalTypes[ne]]){//If they are, we fill the Bonding tab with the appropriate values
                                Bonding[1+at*M*27+nBonding[at]]=ne;
                                Bonding[1+ne*M*27+nBonding[ne]]=at;
				ntotAts = ntotAts+2;
                                nBonding[ne]++;
                                nBonding[at]++;

                                }
                            }
                        }
                    }
                }
            }


	int* Bonds;//Will contain the data of Bonding in a more condensed form
	
	Bonds = calloc(ntotAts+nAtoms+2,sizeof(int));
	Bonds[0]=1;
	for(int i=0 ; i<nAtoms ; i++){
		for(int j=0 ; j<M*27 ; j++){
			if(Bonding[i*M*27+j+1]!=-1){
				Bonds[Bonds[0]]=Bonding[i*27*M+j+1];
				Bonds[0]++;	
			}
		}
		Bonds[Bonds[0]]=-1;
		Bonds[0]++;;		
	}

	free(Bonding);
	free(Mesh);
	free(nMesh);
	free(MeshIndex);
	free(DicoAtoms);
	free(nBonding);

return Bonds;

}
