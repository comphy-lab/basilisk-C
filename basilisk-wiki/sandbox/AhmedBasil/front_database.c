#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_FRONTS 10
#define MAX_ELEMS 10000
#define MAX_POINTS 10000
#define aspmax 1.40
#define amin 0.6
#define amax 0.8

/**
This file maintains database of Front.
 Include functions
		to  add or delete element to the front.
*/

struct Front{
	int *PConnectNext;
	int *PConnectPrev;
	double **PCoords;
	int NP,NPEmpty;
	int FirstPoint;
	int FirstEmptyPoint;

	int *EConnectNext;
	int *EConnectPrev;
	int **ENbr;
	int **ECorner;
	int NE,NEEmpty;
	int FirstElement;
	int FirstEmptyElem;
};

/**
Create an empty front;
*/
void create_empty_front(struct Front *front){
	struct Front *f = (struct Front *)malloc(sizeof(struct Front *));

	f->PConnectNext = (int *)malloc((MAX_POINTS + 1) * sizeof(int));
	f->PConnectPrev = (int *)malloc((MAX_POINTS + 1) * sizeof(int));
	f->PConnectNext = f->PConnectNext + 1;
	f->PCoords = (double **)malloc(MAX_POINTS * sizeof(double *));
	for(int l=0; l<MAX_POINTS; ++l){
		f->PCoords[l] = (double *)malloc(3 * sizeof(double));
	}

	f->EConnectNext = (int *)malloc((MAX_ELEMS + 1) * sizeof(int));
	f->EConnectPrev = (int *)malloc((MAX_ELEMS + 1) * sizeof(int));
	f->EConnectNext = f->EConnectNext + 1;
	f->ENbr = (int **)malloc(MAX_ELEMS * sizeof(int *));
	f->ECorner = (int **)malloc(MAX_ELEMS * sizeof(int *));
	for(int l=0; l<MAX_ELEMS; ++l){
		f->ENbr[l] = (int *)malloc(3 * sizeof(int));
		f->ECorner[l] = (int *)malloc(3 * sizeof(int));
	}

	//Set  default values
	f->NP = 0;
	f->NPEmpty = MAX_POINTS;
	f->NE = 0;
	f->NEEmpty = MAX_ELEMS;

	//Setting default connectivity(Cyclic connection)
	int *connectf,*connectb;

	connectf = f->PConnectNext;
	connectb = f->PConnectPrev;
	for (int l=0; l<=MAX_POINTS ; ++l){
		connectf[l-1] = l;
		connectb[l] = l-1;
	}
	f->FirstPoint = -1;
	f->FirstEmptyPoint = 0;

	connectf = f->EConnectNext;
	connectb = f->EConnectPrev;	
	for (int l=0; l<=MAX_ELEMS ; ++l){
		connectf[l-1] = l;
		connectb[l] = l-1;
	}
	f->FirstElement = -1;
	f->FirstEmptyElem = 0;

	front = f;	
}
 
/**
List containing all fronts.
*/
struct FrontList{
	struct Front **FArray;
	int NF;
};

void initialize_front_list(struct FrontList *front_list){
	struct FrontList *fl = (struct FrontList *)malloc(sizeof(struct FrontList *));
	fl->NF = 0;
	fl->FArray	= (struct Front **)malloc(MAX_FRONTS*sizeof(struct Front *)); 
	
	front_list = fl; 
}

/**
Push a front to front_list.
*/
int push_to_front_list(struct FrontList *front_list,struct Front *front){
	if (front_list->NF > MAX_FRONTS){
		//ERROR: cannot add one more front
		return(0); 
	}
	create_empty_front(front);
	front_list->FArray[front_list->NF]	=	front;
	front_list->NF += 1;		
	return(1);
}

/**
Functions related to front.
*/

/**
Objects (Points or Elements) are kept as a linked array.
One object is connected to two neighbours  (Next and Prev).
Total number of used and unused objects is constant.
First Object(FO) and First Empty Object (FEO) are two adjacent
Objects in the list. All the objects from FO connected in the Prev(<-) 
direction constitute the used objects.
*/
/*
Example:
	5 <- 2 <- 1 <- 8 <- 4 <- || -> 3 -> 6 ->7
	                    FO        FEO

In the above example objects 5,2,1,8,4 are used objects. 
Adding an object to used list is straight  forward. Just
retag FO and FEO. See below, 
	5 <- 2 <- 1 <- 8 <- 4 <- 3 || -> 6 ->7
	                        FO      FEO

Let's say we want to delete 1(Delete Object,DO);
	5 <- 2 <- 1 <- 8 <- 4 <- || -> 3 -> 6 ->7
	          DO       FO         FEO
we will replace DO between FO and FEO and retag FEO
	5 <- 2 <- 8 <- 4 <- || -> 1 -> 3 -> 6 ->7
	               FO         FEO

*/
int delete_obj(int *ConnectNext, int *ConnectPrev, int *OldObj, int *FirstObj,
                int *FirstEmptyObj, int *NumObj, int *NumEmptyObj){
	if (*NumObj == 0){
		//ERROR: cannot delete from an empty array of objects.
		return(0);
	}

	ConnectNext[ConnectPrev[*OldObj]] = ConnectNext[*OldObj]; 
	ConnectPrev[ConnectNext[*OldObj]] = ConnectPrev[*OldObj]; 
	ConnectNext[*OldObj] = *FirstEmptyObj;
	ConnectPrev[*OldObj] = *FirstObj;
	ConnectNext[*FirstObj] = *OldObj;
	ConnectPrev[*FirstEmptyObj] = *OldObj;

	*FirstEmptyObj = *OldObj;
	*NumObj -= 1;
	*NumEmptyObj += 1;

	return(1);
}

int add_obj(int *ConnectNext, int *ConnectPrev, int *NewObj, int *FirstObj,
             int *FirstEmptyObj, int *NumObj, int *NumEmptyObj){
	if (*NumEmptyObj == 0){
		//ERROR: cannot add to a completely-filled array.
		return(0);
	}

	*NewObj = *FirstEmptyObj;
	*FirstObj = *FirstEmptyObj;
	*FirstEmptyObj = ConnectNext[*FirstEmptyObj];
	*NumObj += 1;
	*NumEmptyObj -= 1;

	return(1);
}

/**
Function nflocal().
*/  
	
int nflocal(int **Corner, int **Nbr,
            int *m, int *n0, int *n1, int *n2,
            int *mc, int *nc0, int *nc1, int *nc2, int *p){

	*n1 = (*n0 + 1)%3;
	*n2 = (*n0 + 2)%3;
	p[0] = Corner[*m][*n0];  
	p[1] = Corner[*m][*n1];  
	p[2] = Corner[*m][*n2];  

	*mc = Nbr[*m][*n0];
	for(int i=0 ; i<3 ; ++i){
		if(Corner[*mc][i] == p[0]) {
			*nc0= i;
		}
	}
	*nc1 = (*nc0 + 1)%3;
	*nc2 = (*nc0 + 2)%3;
	p[3] = Corner[*mc][*nc1];

	int vertex;
	for (int i=0 ; i<3; ++i) {
		vertex = Corner[Nbr[*m][*n1]][i];
		if( vertex != p[1] && vertex != p[2])		p[4] = vertex; 
		vertex = Corner[Nbr[*m][*n2]][i];
		if( vertex != p[2] && vertex != p[0])		p[5] = vertex; 	
	}

	return(1);
}

/**
Delete an element (And also one of its neighbour).
*/
int delete_front_element(struct Front *fr, int elem_id, int vertex_id){
	int m, n0, n1, n2, 
	    mc, nc0, nc1, nc2,
	    m1, m2, mc0, mc1, p[8],
	    mtemp, ntemp,
	    vertex,
	    **Corner, **Nbr;
	double **Coords;

	m = elem_id;
	n0 = vertex_id;

	Corner = fr->ECorner;
	Nbr = fr->ENbr;
	Coords = fr->PCoords;

	nflocal(fr->ECorner, fr->ENbr, &m, &n0, &n1, &n2,
	         &mc, &nc0, &nc1, &nc2, p);
	
	m1 = Nbr[m][n1];
	m2 = Nbr[m][n2];
	mc0 = Nbr[mc][nc0];
	mc1 = Nbr[mc][nc1];

	for(int i=0; i<3; ++i){
		vertex = Corner[mc0][i];
		if(vertex != p[0] && vertex != p[3]) p[6] = vertex;
		vertex = Corner[mc1][i];
		if(vertex != p[3] && vertex != p[1]) p[7] = vertex;
	}
	
	//Conditions to abandon deleting.
	for(int i=0; i<3; ++i){
		for (int j=0; j<3; ++j){
			if( Nbr[m1][i] == Nbr[m2][j] && Nbr[m1][i] != m) {
				//i.e if m1 and m2 have a common neighbour other than m
				return(0); 
			}
			if( Nbr[mc0][i] == Nbr[mc1][j] && Nbr[mc0][i] != mc) {
				//i.e if mc0 and mc1 have a common neighbour other than mc
				return(0);
			}
		}
	}
	
	for(int i=0; i<3; ++i){
		Coords[p[0]][i] = 0.50 *( Coords[p[0]][i] + Coords[p[1]][i])-
		                  0.125 *( Coords[p[2]][i] + Coords[p[3]][i])+
		                  0.0625 *( Coords[p[4]][i]+ Coords[p[5]][i] + Coords[p[6]][i] + Coords[p[7]][i]);
		
	}

	//Delete Elements m,mc and Point p[1]
	delete_obj(fr->EConnectNext, fr->EConnectPrev, &m, &(fr->FirstElement),
	           &(fr->FirstEmptyElem), &(fr->NE), &(fr->NEEmpty) );
	delete_obj(fr->EConnectNext, fr->EConnectPrev, &mc, &(fr->FirstElement),
	           &(fr->FirstEmptyElem), &(fr->NE), &(fr->NEEmpty) );
	delete_obj(fr->PConnectNext, fr->PConnectPrev, p+1, &(fr->FirstPoint),
	           &(fr->FirstEmptyPoint), &(fr->NP), &(fr->NPEmpty) );

	//Resetting neighbour elements around deleed elements
	for(int i=0; i<3; ++i){
		if(Nbr[m1][i] == m) Nbr[m1][i] = m2;
		if(Nbr[m2][i] == m) Nbr[m2][i] = m1;
		if(Nbr[mc0][i] == mc) Nbr[mc0][i] = mc1;
		if(Nbr[mc1][i] == mc) Nbr[mc1][i] = mc0;
	}
	
	//Reset corner of all elements having P[1] as a corner from P[1] to P[0]
	mtemp = m1;
	while(1){
		for(int i=0; i<3; ++i){
			if(Corner[mtemp][i] == p[1]) {
				Corner[mtemp][i] = p[0];
				ntemp = i;
			}
		}
		if(Nbr[mtemp][ntemp] == mc) break;
		mtemp = Nbr[mtemp][ntemp]; 
	}
	
		
	return(1);	//Delete successful		
}

/**
Add an element.
*/
int add_front_element(struct Front *fr, int elem_id, int vertex_id){
	int m, n0, n1, n2, 
	    mc, nc0, nc1, nc2,
	    m1, m2, mc0, mc1, p[8],
	    newpt, newm, newmc,
	    vertex,
	    **Corner, **Nbr;
  double **Coords;

	m = elem_id;
	n0 = vertex_id;
	
	Corner = fr->ECorner;
	Nbr = fr->ENbr;
	Coords = fr->PCoords;

	nflocal(fr->ECorner, fr->ENbr, &m, &n0, &n1, &n2,
	         &mc, &nc0, &nc1, &nc2, p);
	
	m1 = Nbr[m][n1];
	m2 = Nbr[m][n2];
	mc0 = Nbr[mc][nc0];
	mc1 = Nbr[mc][nc1];

	for(int i=0; i<3; ++i){
		vertex = Corner[mc0][i];
		if(vertex != p[0] && vertex != p[3]) p[6] = vertex;
		vertex = Corner[mc1][i];
		if(vertex != p[3] && vertex != p[1]) p[7] = vertex;
	}
		

	//Create Elements newm, newmc and Point newpt
	add_obj(fr->EConnectNext, fr->EConnectPrev, &newm, &(fr->FirstElement),
	           &(fr->FirstEmptyElem), &(fr->NE), &(fr->NEEmpty) );
	add_obj(fr->EConnectNext, fr->EConnectPrev, &newmc, &(fr->FirstElement),
	           &(fr->FirstEmptyElem), &(fr->NE), &(fr->NEEmpty) );
	add_obj(fr->PConnectNext, fr->PConnectPrev, &newpt, &(fr->FirstPoint),
	           &(fr->FirstEmptyPoint), &(fr->NP), &(fr->NPEmpty) );

	//Coordinates of new point
	for(int i=0; i<3; ++i){
		Coords[newpt][i] = 0.50 *( Coords[p[0]][i] + Coords[p[1]][i])-
		                   0.125 *( Coords[p[2]][i] + Coords[p[3]][i])+
		                   0.0625 *( Coords[p[4]][i]+ Coords[p[5]][i] + Coords[p[6]][i] + Coords[p[7]][i]);
		
	}

	//Resetting Corners and Nbrs
	Corner[m][n0] = newpt;
	Corner[mc][nc0] = newpt;
	
	Corner[newm][0] = newpt;
	Corner[newm][1] = p[2];
	Corner[newm][2] = p[0];

	Corner[newmc][0] = newpt;
	Corner[newmc][1] = p[0];
	Corner[newmc][2] = p[3];
	
	Nbr[m][n2] = newm;
	Nbr[mc][nc0] = newmc;

	Nbr[newm][0] = m;
	Nbr[newm][1] = m2;
	Nbr[newm][2] = newmc;
		
	Nbr[newmc][0] = newm;
	Nbr[newmc][1] = mc0;
	Nbr[newmc][2] = mc;

	for (int i=0; i<3; ++i){
		if(Nbr[m2][i] == m)    Nbr[m2][i] = newm;
		if(Nbr[mc0][i] == mc)    Nbr[mc0][i] = newmc;
	}

	return(1);		
}
/**
Finding distance between two points in space
*/
double distance(double *CoordP1, double *CoordP2){
	double d = (CoordP1[0] - CoordP2[0]) * (CoordP1[0] - CoordP2[0]) +
	           (CoordP1[1] - CoordP2[1]) * (CoordP1[1] - CoordP2[1]) +
	           (CoordP1[2] - CoordP2[2]) * (CoordP1[2] - CoordP2[2]) ;
	return (sqrt(d));
}
/**
Finding Area of a triangular element
*/

double element_area(double **Coords, int p0, int p1, int p2){
	double area, x0, y0, z0, x1, y1, z1;
	x0 = Coords[p1][0] - Coords[p0][0]; 		
	y0 = Coords[p1][1] - Coords[p0][1]; 		
	z0 = Coords[p1][2] - Coords[p0][2]; 		
	x1 = Coords[p2][0] - Coords[p0][0]; 		
	y1 = Coords[p2][1] - Coords[p0][1]; 		
	z1 = Coords[p2][2] - Coords[p0][2]; 		

	//Area  = 0.5 | u x v |	
	area = 0.50 * sqrt((y0*z1 - y1*z0)*(y0*z1 - y1*z0) + 
	                   (x1*z0 - x0*z1)*(x1*z0 - x0*z1) +
	                   (x0*y1 - x1*y0)*(x0*y1 - x1*y0) );
	return (area);
}
/**
Finding aspect ratio of a triangle. 
Aspect ratio is 1.0 for an equilateral triangle.
Aspect ratio >= 1.
For quality mesh, aspect ratio is limited by aspmax.
*/

double aspect_ratio(int elem, int **Corner, double **Coords){
	int p0, p1, p2;
	double AR, s0, s1, s2, s, area;
	p0 = Corner[elem][0];
	p1 = Corner[elem][1];
	p2 = Corner[elem][2];
	s0 = distance(Coords[p0], Coords[p1]);
	s1 = distance(Coords[p1], Coords[p2]);
	s2 = distance(Coords[p2], Coords[p0]);
	s = (s0 + s1 + s2)/3.0;
	area = element_area(Coords, p0, p1, p2);
	AR = 0.25*sqrt(3)*s*s/area;
	return(AR);
}

/**
Finding the quality of a grid.
Following function returns 0 if all the three 
quality standards are satisfied. Returns 1 otherwise
*/
int front_quality(struct Front *fr){
	int *ConnectPrev, **Corner,  NE, FirstElement,
	    elem, ielem;
	double **Coords, area, ar;

	ConnectPrev = fr->EConnectPrev;
	Coords = fr->PCoords;
	Corner = fr->ECorner;
	NE = fr->NE;
	FirstElement = fr->FirstElement;	

	double aamin = amin*amin*1.7/4.0;
	double aamax = amax*amax*1.7/4.0;

	ielem = 0;
	elem = FirstElement;
	while(ielem < NE) {	
		area = element_area(Coords, Corner[elem][0], Corner[elem][1], Corner[elem][2]);
		if(area < aamin) return(1);
		if(area > aamax) return(1);
		ar = aspect_ratio(elem, Corner, Coords);
		if(ar > aspmax) return(1);
		
		ielem++;
		elem = ConnectPrev[elem];
	}
	
	return(0);
}
/**
Checking the quality of each element and delete/add 
elements.
*/

void regrid_front(struct Front *fr){
	int *ConnectPrev, *ConnectNext, **Corner, **Nbr,
			*NE, *FirstElement,
	    elem, ielem, p0, p1, p2, n0, n1, n2, ntemp, irough;
	double **Coords,
	       s0, s1, s2, stemp;

	ConnectPrev = fr->EConnectPrev;
	Coords = fr->PCoords;
	Corner = fr->ECorner;
	NE = &(fr->NE);
	FirstElement = &(fr->FirstElement);	

	for(int ipass = 0; ipass<6; ++ipass){

		//Adding elements if needed
		ielem = 0;
		elem = *FirstElement;
		while(ielem < *NE) {

			p0 = Corner[elem][0];
			p1 = Corner[elem][1];
			p2 = Corner[elem][2];
			s0 = distance(Coords[p0], Coords[p1]);
			s1 = distance(Coords[p1], Coords[p2]);
			s2 = distance(Coords[p2], Coords[p0]);
			n0 = 0; n1 = 1; n2 = 2;		

			//Sorting sidelengths in increasing order s0 <= s1 <= s2
			if(s0 > s1){
				stemp = s0;	s0 = s1; s1 = stemp;
				ntemp = n0;	n0 = n1; n1 = ntemp;
			}
			if(s1 > s2){
				stemp = s1;	s1 = s2; s2 = stemp;
				ntemp = n1;	n1 = n2; n2 = ntemp;
			}
			if(s0 > s1){
				stemp = s0;	s0 = s1; s1 = stemp;
				ntemp = n0;	n0 = n1; n1 = ntemp;
			}
				
		  if( (s2 > amax) || 
					((s0 > amin ) && (aspect_ratio(elem, Corner, Coords) > aspmax) ) ) {
					add_front_element(fr, elem, n2);
					//if an element added, start the loop from the beginning
					ielem = 0;
					elem = *FirstElement;
					continue;					
			}

			ielem++;
			elem = ConnectPrev[elem];
		}

		//Deleting elements if needed
		ielem = 0;
		elem = *FirstElement;
		while(ielem < *NE) {

			p0 = Corner[elem][0];
			p1 = Corner[elem][1];
			p2 = Corner[elem][2];
			s0 = distance(Coords[p0], Coords[p1]);
			s1 = distance(Coords[p1], Coords[p2]);
			s2 = distance(Coords[p2], Coords[p0]);
			n0 = 0; n1 = 1; n2 = 2;		

			if(s0 > s1){
				stemp = s0;	s0 = s1; s1 = stemp;
				ntemp = n0;	n0 = n1; n1 = ntemp;
			}
			if(s1 > s2){
				stemp = s1;	s1 = s2; s2 = stemp;
				ntemp = n1;	n1 = n2; n2 = ntemp;
			}
			if(s0 > s1){
				stemp = s0;	s0 = s1; s1 = stemp;
				ntemp = n0;	n0 = n1; n1 = ntemp;
			}
	
		  if(s0 < amin) {
					if (delete_front_element(fr, elem, n2)) {
						//if delete was succesful, start the loop from the beginning
						ielem = 0;
						elem = *FirstElement;
						continue;					
					}
			}

			ielem++;
			elem = ConnectPrev[elem];
		}
		
		if(front_quality(fr) == 0) {
			break;
		}
	}
}