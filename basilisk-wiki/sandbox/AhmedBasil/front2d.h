#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <dirent.h>
#include <errno.h>

#include <string.h>

#include "draw.h"

#include "front_utils.h"
#include "vofi.h"

#define MAX_FRONTS 10
#define MAX_ELEMS 10000
#define MAX_POINTS 10000
//#define amin 0.005
//#define amax 0.010
#define amin 0.0031250
#define amax 0.0046875


#ifndef LEVEL
# define LEVEL 8
#endif

scalar c_refer[];
/**
This file maintains database of Front.
 Include functions
		to  add or delete element to the front.
*/

typedef struct{
	int *PConnectNext;
	int *PConnectPrev;
	double **PCoords;
	int *NP, *NPEmpty;
	int *FirstPoint;
	int *FirstEmptyPoint;

	int *EConnectNext;
	int *EConnectPrev;
	int **ENbr;
	int **ECorner;
	int *NE,*NEEmpty;
	int *FirstElement;
	int *FirstEmptyElement;
}Front;

typedef struct{
	Front *fr;
	int ndiv; 
	double delta;
	int **numpoints;			//number of marker points occcupied by each cell at level MAX_LEVEL
	int **iselems;			  //=1 if cell contains a marker point ; should be bool** rather than int **
	int **numbuffer;		  //used to store temporary variables	
	int **numpoints_cum;	//cumulative number of points 
	int **elem_id;	  		//element id of any element in the cell (if any) 
	int *point_id;				//marker points rearranged
}Front2vof;


Front *fr_circle;
Front2vof *f2v_circle;

/**
foreach_frontpoint(fr) is macro that iterates through each marker points of the front*/
@def foreach_frontpoint(fr) 
	{
		int *PrevPoint, NP, FirstPoint, frontpoint, ipt;
		double **Coords;

		PrevPoint  = fr->PConnectPrev;
		Coords     = fr->PCoords;
		NP         = *(fr->NP);
		FirstPoint = *(fr->FirstPoint);	
		NOT_UNUSED(Coords);

		frontpoint = FirstPoint;
		ipt = 0;
		while(ipt < NP){
@
@def end_foreach_frontpoint()
			ipt++;
			frontpoint = PrevPoint[frontpoint];
		}
	}
@

/**
foreach_frontelement(fr) is macro that iterates through each marker elements of the front*/
@def foreach_frontelement(fr) 
	{
		int *PrevElement, **Corner, **Nbr, NE, FirstElement, frontelement, ielem;
		double **Coords;

		PrevElement  = fr->EConnectPrev;
		Coords       = fr->PCoords; 
		Corner       = fr->ECorner; 
		Nbr          = fr->ENbr;
		NE           = *(fr->NE);
		FirstElement = *(fr->FirstElement);	
		NOT_UNUSED(Coords);
		NOT_UNUSED(Corner);
		NOT_UNUSED(Nbr);

		frontelement = FirstElement;
		ielem = 0;
		while(ielem < NE){
@
@def end_foreach_frontelement()
			ielem++;
			frontelement = PrevElement[frontelement];
		}
	}
@

/**
Create an empty front;
*/
Front *create_empty_front(){
		
	Front *f = (Front *)malloc(sizeof(Front));

	f->PConnectNext = (int *)malloc((MAX_POINTS + 1) * sizeof(int));
	f->PConnectPrev = (int *)malloc((MAX_POINTS + 1) * sizeof(int));
	f->PConnectNext = f->PConnectNext + 1;
	f->PCoords = (double **)malloc(MAX_POINTS * sizeof(double *));
	for(int l=0; l<MAX_POINTS; ++l){
		f->PCoords[l] = (double *)malloc(2 * sizeof(double));
	}

	f->EConnectNext = (int *)malloc((MAX_ELEMS + 1) * sizeof(int));
	f->EConnectPrev = (int *)malloc((MAX_ELEMS + 1) * sizeof(int));
	f->EConnectNext = f->EConnectNext + 1;
	f->ENbr = (int **)malloc(MAX_ELEMS * sizeof(int *));
	f->ECorner = (int **)malloc(MAX_ELEMS * sizeof(int *));
	for(int l=0; l<MAX_ELEMS; ++l){
		f->ENbr[l] = (int *)malloc(2 * sizeof(int));
		f->ECorner[l] = (int *)malloc(2 * sizeof(int));
	}

	f->NP = (int *)malloc(sizeof(int));
	f->NPEmpty = (int *)malloc(sizeof(int));
	f->FirstPoint = (int *)malloc(sizeof(int));
	f->FirstEmptyPoint = (int *)malloc(sizeof(int));
	f->NE = (int *)malloc(sizeof(int));
	f->NEEmpty = (int *)malloc(sizeof(int));
	f->FirstElement = (int *)malloc(sizeof(int));
	f->FirstEmptyElement = (int *)malloc(sizeof(int));

	//Set  default values
 	
	*(f->NP) = 0;
	*(f->NPEmpty) = MAX_POINTS;
	*(f->NE) = 0;
	*(f->NEEmpty) = MAX_ELEMS;

	//Setting default connectivity(Cyclic connection)
	int *connectf,*connectb;

	connectf = f->PConnectNext;
	connectb = f->PConnectPrev;
	for (int l=0; l<=MAX_POINTS ; ++l){
		connectf[l-1] = l;
		connectb[l] = l-1;
	}
	*(f->FirstPoint) = -1;
	*(f->FirstEmptyPoint) = 0;

	connectf = f->EConnectNext;
	connectb = f->EConnectPrev;	
	for (int l=0; l<=MAX_ELEMS ; ++l){
		connectf[l-1] = l;
		connectb[l] = l-1;
	}
	*(f->FirstElement) = -1;
	*(f->FirstEmptyElement) = 0;

	return(f);
}

/**
Free memory*/
void free_front(Front *fr){
	if(fr == NULL) return;

	fr->PConnectNext = fr->PConnectNext - 1;
	free(fr->PConnectNext);	
	free(fr->PConnectPrev);	
	for(int l=0; l<MAX_POINTS; ++l){
		free(fr->PCoords[l]);
	}
	free(fr->PCoords);
	fr->EConnectNext = fr->EConnectNext - 1;
	free(fr->EConnectNext);	
	free(fr->EConnectPrev);	
	for(int l=0; l<MAX_ELEMS; ++l){
		free(fr->ENbr[l]);
		free(fr->ECorner[l]);
	}
	free(fr->ENbr);
	free(fr->ECorner);

	free(fr->NP);
	free(fr->NPEmpty);
	free(fr->FirstPoint);
	free(fr->FirstEmptyPoint);
	free(fr->NE);
	free(fr->NEEmpty);
	free(fr->FirstElement);
	free(fr->FirstEmptyElement);
	
	free(fr);
	return;
}
 
/**
List containing all fronts.
*/
typedef struct{
	Front **FArray;
	int NF;
}FrontList;

FrontList *initialize_front_list(){
	FrontList *fl = (FrontList *)malloc(sizeof(FrontList *));
	fl->NF = 0;
	fl->FArray	= (Front **)malloc(MAX_FRONTS*sizeof(Front *)); 
	
	return(fl); 
}

/**
Push a front to front_list.
*/
int push_to_front_list(FrontList *front_list, Front *front){
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
Delete an element.
*/
int delete_front_element(Front *fr, int elem_id){
	int **Corner, **Nbr, p[4], n0, n1;
	double **Coords;

	Corner = fr->ECorner;
	Nbr = fr->ENbr;
	Coords = fr->PCoords;

	n0 = Nbr[elem_id][0];
	n1 = Nbr[elem_id][1];
	p[0] = Corner[elem_id][0];	
	p[1] = Corner[elem_id][1];	
	//p[2] = Corner[n0][0];	
	//p[3] = Corner[n1][1];	
	
	//Conditions to abandon deleting.
	//I don't know if any exists.
	
	for(int i=0; i<2; ++i){
		Coords[p[0]][i] = 0.50 *( Coords[p[0]][i] + Coords[p[1]][i]);
		//Need to redefine it.
		
	}

	//Delete Element elem_id and p[1]
	delete_obj(fr->EConnectNext, fr->EConnectPrev, &elem_id, fr->FirstElement,
	           fr->FirstEmptyElement, fr->NE, fr->NEEmpty );
	delete_obj(fr->PConnectNext, fr->PConnectPrev, p+1, fr->FirstPoint,
	           fr->FirstEmptyPoint, fr->NP, fr->NPEmpty );

	//Resetting neighbour elements around deleed elements
	Nbr[n0][1] = n1;
	Nbr[n1][0] = n0;
	
	//Reset corner of all elements having P[1] as a corner from P[1] to P[0]
	//Corner[n0][1] = p[0];
	Corner[n1][0] = p[0];
		
	return(1);	//Delete successful		
}

/**
Add an element.
*/
int add_front_element(Front *fr, int elem_id){
	int **Corner, **Nbr, p[4], n1,
	    newpt, newelem;
	double **Coords;

	Corner = fr->ECorner;
	Nbr = fr->ENbr;
	Coords = fr->PCoords;

	//n0 = Nbr[elem_id][0];
	n1 = Nbr[elem_id][1];

	p[0] = Corner[elem_id][0];	
	p[1] = Corner[elem_id][1];	
	//p[2] = Corner[n0][0];	
	//p[3] = Corner[n1][1];	
	newpt = 0;	  //to avoid warning
	newelem = 0;	// do

	//Create Element newelem, and Point newpt
	add_obj(fr->EConnectNext, fr->EConnectPrev, &newelem, fr->FirstElement,
	           fr->FirstEmptyElement, fr->NE, fr->NEEmpty );
	add_obj(fr->PConnectNext, fr->PConnectPrev, &newpt, fr->FirstPoint,
	           fr->FirstEmptyPoint, fr->NP, fr->NPEmpty );

	//Coordinates of new point
	for(int i=0; i<2; ++i){
		Coords[newpt][i] = 0.50 *( Coords[p[0]][i] + Coords[p[1]][i]);
		//Need to redefine it.
		
	}

	//Resetting Corners and Nbrs
	Nbr[elem_id][1] = newelem;
	Nbr[newelem][0] = elem_id;	
	Nbr[newelem][1] = n1;	
	Nbr[n1][0] = newelem;

	Corner[elem_id][1] = newpt;	
	Corner[newelem][0] = newpt;
	Corner[newelem][1] = p[1];

	return(1);		
}
/**
elem_length
*/
double elem_length(double **Coords, int **Corner, int elem){
	double *p0 = Coords[Corner[elem][0]];
	double *p1 = Coords[Corner[elem][1]];
	double d = (p1[0] - p0[0]) * (p1[0] - p0[0]) +
	           (p1[1] - p0[1]) * (p1[1] - p0[1]);
	return (sqrt(d));
}

/**
Finding the quality of a grid.
Following function returns 0 if all the three 
quality standards are satisfied. Returns 1 otherwise
*/
int front_quality(Front *fr){

	double elength;
	foreach_frontelement(fr) {
		elength = elem_length(Coords, Corner, frontelement);	
		if(elength < amin) return(1);
		if(elength > amax) return(1);
	}
	
	return(0);
}

/**
Checking the quality of each element and delete/add 
elements.
*/

void regrid_front(Front *fr){
	int *ConnectPrev, **Corner,
			*NE, *FirstElement,
	    elem, ielem;
	double **Coords, elength;

	ConnectPrev = fr->EConnectPrev;
	Coords = fr->PCoords;
	Corner = fr->ECorner;
	NE = fr->NE;
	FirstElement = fr->FirstElement;	

	for(int ipass = 0; ipass<6; ++ipass){

		//Adding elements if needed
		ielem = 0;
		elem = *FirstElement;
		while(ielem < *NE) {

			elength = elem_length(Coords, Corner, elem);

		  if( elength > amax ) { 
				add_front_element(fr, elem);
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

			elength = elem_length(Coords, Corner, elem);

		  if( elength < amin ) { 
				delete_front_element(fr, elem);
						//After deleting, start the loop from the beginning
				ielem = 0;
				elem = *FirstElement;
				continue;					
				
			}

			ielem++;
			elem = ConnectPrev[elem];
		}
		
		if(front_quality(fr) == 0) {
			break;
		}

	}
}

/**
Write the details of front into a file.
Points and Elements will be renamed such that point_id
and element_id falls between 0 to NP-1 and 0 to NE-1 respectively.
For element without a neighbour (elements at boundary), neighbour element is designated by -1.
*/
int write_front(Front *fr){

	int point_redirect[MAX_POINTS];
	int elem_red[MAX_ELEMS+1];
	int *elem_redirect = elem_red+1;
	elem_redirect[-1] = -1;

	FILE * fp;
	fp = fopen("data_front.txt","w");
	fprintf(fp, "NUM_POINTS\n{\n%d\n}\n", *(fr->NP)); 
	fprintf(fp, "NUM_ELEMENTS\n{\n%d\n}\n", *(fr->NE)); 
	//writing point coordinates
	fprintf(fp, "POINTS\n{\n");
		foreach_frontpoint(fr) {
			point_redirect[frontpoint] = ipt;
			fprintf(fp,"%f %f\n", Coords[frontpoint][0], Coords[frontpoint][1]);
		}
	fprintf(fp,"}\n");

	//writing element details
	fprintf(fp, "ELEMENTS\n{\n");
		foreach_frontelement(fr) {
			elem_redirect[frontelement] = ielem;
		}

		fprintf(fp, "CORNER\n{\n");
		foreach_frontelement(fr) {
			fprintf(fp,"%d %d\n", point_redirect[Corner[frontelement][0]], point_redirect[Corner[frontelement][1]] );
		}
		fprintf(fp,"}\n");

		fprintf(fp, "NEIGHBOUR\n{\n");
		foreach_frontelement(fr) {
			fprintf(fp,"%d %d\n", elem_redirect[Nbr[frontelement][0]], elem_redirect[Nbr[frontelement][1]] );
		}
		fprintf(fp,"}\n");

	fprintf(fp,"}\n");

	fclose(fp);

	return 0;	
}

/**
Write the details of front into a file.
Written in 2 different files data_front_xyz.dat and data_front_connectivity.dat.
*/
int write_front_python(double t, Front *fr){

	int point_redirect[MAX_POINTS];

	FILE * fp1;
	FILE * fp2;
	char filename1[100];
	char filename2[100];
	sprintf(filename1,"pythondirectory/front_xyz_t%.6f.dat", t);
	sprintf(filename2,"pythondirectory/front_connectivity_t%.6f.dat", t);
	fp1 = fopen(filename1, "w");
	fp2 = fopen(filename2, "w");
	//writing point coordinates
		foreach_frontpoint(fr){
			point_redirect[frontpoint] = ipt;
			fprintf(fp1,"%f\n%f\n", Coords[frontpoint][0], Coords[frontpoint][1]);
		}

	//writing connectivity

		foreach_frontelement(fr){
			fprintf(fp2,"%d\n%d\n", point_redirect[Corner[frontelement][0]], point_redirect[Corner[frontelement][1]]);
		}
	fclose(fp1);	
	fclose(fp2);	

	return 0;
}
/**
print Front details on stdout
*/
int write_front_stdout(Front *fr){

	fprintf(stdout, "\n{FRONT}\n"); 
	fprintf(stdout, "NUM_POINTS\n{\n%d\n}\n", *(fr->NP)); 
	fprintf(stdout, "NUM_ELEMENTS\n{\n%d\n}\n", *(fr->NE)); 
	//writing point coordinates
	fprintf(stdout, "POINTS\n{\n");
		foreach_frontpoint(fr){
			fprintf(stdout,"%d (%f, %f)\n", frontpoint, Coords[frontpoint][0], Coords[frontpoint][1] );
		}
	fprintf(stdout,"}\n");

	//writing element details
	fprintf(stdout, "ELEMENTS\n{\n");
		foreach_frontelement(fr){
			fprintf(stdout,"[%d] C(%d %d)",frontelement, Corner[frontelement][0], Corner[frontelement][1]);
			fprintf(stdout," N(%d %d)\n", Nbr[frontelement][0], Nbr[frontelement][1]);
		}

	fprintf(stdout,"}\n");

	return 0;

}
/**
Read the details of front from a file.
*/

int read_front(Front *fr){

	FILE *fp;

	fp = fopen("data_front.txt","r");

	if(fp == NULL){
		printf("\nCouldn't open specified file");
		return(-1);
	}

	char sbuf[20];

	int **Corner, **Nbr, NP, NE, n0, n1;
	double **Coords, x, y;

	Coords = fr->PCoords;
	Corner = fr->ECorner;
	Nbr = fr->ENbr;

	//read NP
	while(fscanf(fp,"%s",sbuf) != EOF){
		if(strcmp(sbuf,"NUM_POINTS") == 0){
			fscanf(fp,"%s",sbuf);
			if(strcmp(sbuf,"{"))	{
				printf("\nCouldn't find an opening brace");
				return (-1);
			}
			fscanf(fp,"%d",&NP);
			break; 
		}
	}
	while(fscanf(fp,"%s",sbuf) != EOF){
		if(strcmp(sbuf,"NUM_ELEMENTS") == 0){
			fscanf(fp,"%s",sbuf);
			if(strcmp(sbuf,"{"))	{
				printf("\nCouldn't find an opening brace");
				return (-1);
			}
			fscanf(fp,"%d",&NE);
			break; 
		}
	}

	//read coordinates
	if(NP<0 || NP>MAX_POINTS) {
		printf("\nNP not within limit");
		return(-1);
	}
	*(fr->NP) = NP;
	*(fr->FirstPoint) = NP - 1;
	*(fr->FirstEmptyPoint) = NP;

	while(fscanf(fp,"%s",sbuf) != EOF){
		if(strcmp(sbuf,"POINTS") == 0){
			fscanf(fp,"%s",sbuf);
			if(strcmp(sbuf,"{"))	{
				printf("\nCouldn't find an opening brace");
				return (-1);
			}
			break; 
		}
	}

	for( int ip = 0; ip<NP; ++ip){
		fscanf(fp,"%lf",&x);
		fscanf(fp,"%lf",&y);
		Coords[ip][0] = x;
		Coords[ip][1] = y;
	}
	//read corners and neighbours
	if(NE<0 || NE>MAX_ELEMS) {
		printf("\nNE not within limit");
		return(-1);
	}
	*(fr->NE) = NE;
	*(fr->FirstElement) = NE - 1;
	*(fr->FirstEmptyElement) = NE;

	while(fscanf(fp,"%s",sbuf) != EOF){
		if(strcmp(sbuf,"CORNER") == 0){
			fscanf(fp,"%s",sbuf);
			if(strcmp(sbuf,"{"))	{
				printf("\nCouldn't find an opening brace");
				return (-1);
			}
			break; 
		}
	}

	for( int ip = 0; ip<NE; ++ip){
		fscanf(fp,"%d",&n0);
		fscanf(fp,"%d",&n1);
		Corner[ip][0] = n0;
		Corner[ip][1] = n1;
	}

	while(fscanf(fp,"%s",sbuf) != EOF){
		if(strcmp(sbuf,"NEIGHBOUR") == 0){
			fscanf(fp,"%s",sbuf);
			if(strcmp(sbuf,"{"))	{
				printf("\nCouldn't find an opening brace");
				return (-1);
			}
			break; 
		}
	}

	for( int ip = 0; ip<NE; ++ip){
		fscanf(fp,"%d",&n0);
		fscanf(fp,"%d",&n1);
		Nbr[ip][0] = n0;
		Nbr[ip][1] = n1;
	}
	

	fclose(fp);	

	return(0);

}

/*read and write in binary mode*/
bool dump_front(struct Dump p) {
	
  FILE * fp = p.fp;
  char def[] = "dump_front", * file = p.file ? p.file : p.fp ? NULL : def;

  char * name = NULL;
  if (file) {
    name = (char *) malloc (strlen(file) + 2);
    strcpy (name, file);
    if ((fp = fopen (name, "w")) == NULL) {
      perror (name);
      exit (1);
    }
		else
			fprintf(stdout,"\nwriting to %s", name);
  }
	int point_redirect[MAX_POINTS];
	int elem_red[MAX_ELEMS+1];
	int *elem_redirect = elem_red+1;
	elem_redirect[-1] = -1;

	double x, y;// coordinates;
	int c1, c2, n1, n2; // corner points and neighbour elems
	Front *fr = fr_circle;

	//write NP and NE
	int NP = *(fr->NP);
	int NE = *(fr->NE);
	fwrite(&NP, sizeof(int), 1, fp); 
	fwrite(&NE, sizeof(int), 1, fp);
	//fprintf(stdout,"\nNp NE %d %d", NP, NE); 
	//writing point coordinates
	foreach_frontpoint(fr) {
		//redirect pts so that point_id lies bw 0 & NP-1
		point_redirect[frontpoint] = ipt;
		x = Coords[frontpoint][0];
		y = Coords[frontpoint][1];
		fwrite(&x, sizeof(double), 1, fp); 
		fwrite(&y, sizeof(double), 1, fp); 
		//fprintf(stdout,"\nx y %f %f", x, y); 
	}

	//redirect elems so that element_id lies bw 0 & NE-1
	foreach_frontelement(fr) {
		elem_redirect[frontelement] = ielem;
	}

	//write corners
	foreach_frontelement(fr) {
		c1 = point_redirect[Corner[frontelement][0]];
		c2 = point_redirect[Corner[frontelement][1]];
		fwrite(&c1, sizeof(int), 1, fp); 
		fwrite(&c2, sizeof(int), 1, fp); 
	}

	//write Nbrs
	foreach_frontelement(fr) {
		n1 = elem_redirect[Nbr[frontelement][0]];
		n2 = elem_redirect[Nbr[frontelement][1]];
		fwrite(&n1, sizeof(int), 1, fp); 
		fwrite(&n2, sizeof(int), 1, fp); 
	}

  if (file) {
    fclose (fp);
    free (name);
  }

	return 0;	
}
/*
struct Restore_front {
	bool status;
	double *xs;
	double *xe;
	double *ys;
	double *ye;
	int ne;
}*/

bool restore_front(struct Dump p) {

  FILE * fp = p.fp;
  char def[] = "dump_front", * file = p.file ? p.file : p.fp ? NULL : def;

  char * name = NULL;
  if (file) {
    name = (char *) malloc (strlen(file) + 2);
    strcpy (name, file);
    if ((fp = fopen (name, "r")) == NULL) {
      //perror (name);
      //exit (1);
			fprintf(stdout,"\nRestore fr:failed\n");
			return(0);
    }
  }

	//bool alloc;
	if(fr_circle == NULL) {
		fr_circle = create_empty_front();
		//alloc = 1;
	}
	Front *fr = fr_circle;

	double x,y;
	int a,b;
	double **Coords = fr->PCoords;
	int **Nbr = fr->ENbr;
	int **Corner = fr->ECorner;

	//read NP and NE
	if(fread(&a, sizeof(int), 1, fp) < 1) return(0); 
	if(fread(&b, sizeof(int), 1, fp) < 1) return(0); 

	*(fr->NP) = a;
	*(fr->FirstPoint) = a - 1;
	*(fr->FirstEmptyPoint) = a;
	*(fr->NE) = b;
	*(fr->FirstElement) = b - 1;
	*(fr->FirstEmptyPoint) = b;

	fprintf(ferr, "\nNP %d,NE %d", a, b);

	//writing point coordinates
	for(int i=0; i<*(fr->NP); ++i) {
		if(fread(Coords[i], sizeof(double), 2, fp) < 2) 
			return 0;
		fprintf(ferr, "< %f, %f>", Coords[i][0], Coords[i][1]);
	}
		

	//read corners
	for(int i=0; i<*(fr->NE); ++i) {
		if(fread(Corner[i], sizeof(int), 2, fp) < 2)
			return 0;
	}
	//read nbrs
	for(int i=0; i<*(fr->NE); ++i) {
		if(fread(Nbr[i], sizeof(int), 2, fp) < 2)
			return 0;
	}

  if (file) {
    fclose (fp);
    free (name);
  }

	//if(alloc)
		//free_front(fr);
		
	fprintf(stdout,"\nRestore fr:success");
	return 1;	
}

bool draw_front(struct Dump p) {

  FILE * fp = p.fp;
  char def[] = "dump_front", * file = p.file ? p.file : p.fp ? NULL : def;

  char * name = NULL;
  if (file) {
    name = (char *) malloc (strlen(file) + 2);
    strcpy (name, file);
    if ((fp = fopen (name, "r")) == NULL) {
      //perror (name);
      //exit (1);
			fprintf(stdout,"\nRestore fr:failed\n");
			return(0);
    }
  }

	Front *fr = create_empty_front();

	int a,b;
	double **Coords = fr->PCoords;
	int **Nbr = fr->ENbr;
	int **Corner = fr->ECorner;

	//read NP and NE
	if(fread(&a, sizeof(int), 1, fp) < 1) return(0); 
	if(fread(&b, sizeof(int), 1, fp) < 1) return(0); 

	*(fr->NP) = a;
	*(fr->FirstPoint) = a - 1;
	*(fr->FirstEmptyPoint) = a;
	*(fr->NE) = b;
	*(fr->FirstElement) = b - 1;
	*(fr->FirstEmptyPoint) = b;

	fprintf(ferr, "\nFront %d points,%d elems", a, b);

	//writing point coordinates
	for(int i=0; i<*(fr->NP); ++i) {
		if(fread(Coords[i], sizeof(double), 2, fp) < 2) 
			return 0;
		//fprintf(ferr, "< %f, %f>", Coords[i][0], Coords[i][1]);
	}
		

	//read corners
	for(int i=0; i<*(fr->NE); ++i) {
		if(fread(Corner[i], sizeof(int), 2, fp) < 2)
			return 0;
	}
	//read nbrs
	for(int i=0; i<*(fr->NE); ++i) {
		if(fread(Nbr[i], sizeof(int), 2, fp) < 2)
			return 0;
	}

  if (file) {
    fclose (fp);
    free (name);
  }

	fprintf(ferr,"\nRestore fr:success");
	//draw here
  bview * view = draw();
	
#if dimension == 2
	foreach_frontpoint(fr) {
		glBegin (GL_POINTS);
			glvertex2d(view, Coords[frontpoint][0], Coords[frontpoint][1]);
		glEnd();
		view->ni++;
	}
	foreach_frontelement(fr) {
		glBegin (GL_LINES);
			glvertex2d(view, Coords[Corner[frontelement][0]][0], Coords[Corner[frontelement][0]][1]);
			glvertex2d(view, Coords[Corner[frontelement][1]][0], Coords[Corner[frontelement][1]][1]);
		glEnd();
		view->ni++;
	}
#elif dimension == 3
	assert(false);	//not yet implemented. Is straghtforward
#endif

	free_front(fr);
		
	return 1;	
}

/**
circle front
*/

int add_circle(Front *fr, double xc, double yc,
                double radin, int np){

	int **Corner, **Nbr;
	double **Coords,  _PI, theta, dtheta;

	if(np > MAX_POINTS || np < 10) {
		printf("np should be 10 and MAX_POINTS");
		return(-1);
	}

	Coords = fr->PCoords;
	Corner = fr->ECorner;
	Nbr    = fr->ENbr;

	_PI = 4.0 * atan(1.0);
	dtheta = 2.0*_PI/((double) np);

	*(fr->NP) = np;
	*(fr->FirstPoint) = np-1;
	*(fr->FirstEmptyPoint) = np;
	*(fr->NE) = np;
	*(fr->FirstElement) = np-1;
	*(fr->FirstEmptyElement) = np;

	theta = 0.0;
	for(int pid=0; pid<np; ++pid, theta += dtheta){
		Coords[pid][0] = xc + radin*cos(theta);
		Coords[pid][1] = yc + radin*sin(theta);

		Corner[pid][0] = pid;
		Corner[pid][1] = pid+1;

		Nbr[pid][0] = pid-1;
		Nbr[pid][1] = pid+1;
	}

	Corner[np-1][1] = 0;
	Nbr[0][0] = np-1;
	Nbr[np-1][1] = 0;
	
	return(0);

}
/**
Routines related to front2vof*/

Front2vof *create_empty_front2vof(){

	Front2vof *f2v = (Front2vof *)malloc(sizeof(Front2vof));
	f2v->ndiv      = 1 << LEVEL;
	f2v->delta     = L0/(f2v->ndiv);

	f2v->numpoints      = allocate_int_2d_array(f2v->ndiv, f2v->ndiv, GHOSTS);
	f2v->iselems        = allocate_int_2d_array(f2v->ndiv, f2v->ndiv, GHOSTS);
	f2v->numbuffer      = allocate_int_2d_array(f2v->ndiv, f2v->ndiv, GHOSTS);
	f2v->numpoints_cum  = allocate_int_2d_array(f2v->ndiv, f2v->ndiv, GHOSTS);
	f2v->elem_id        = allocate_int_2d_array(f2v->ndiv, f2v->ndiv, GHOSTS);
	f2v->point_id       = (int *)malloc(MAX_POINTS*sizeof(int));
	
	return f2v; 
}

void free_front2vof(Front2vof *f2v){
	deallocate_int_2d_array(f2v->numpoints, GHOSTS);
	deallocate_int_2d_array(f2v->iselems, GHOSTS);
	deallocate_int_2d_array(f2v->numbuffer, GHOSTS);
	deallocate_int_2d_array(f2v->numpoints_cum, GHOSTS);
	deallocate_int_2d_array(f2v->elem_id, GHOSTS);
	free(f2v->point_id);
	free(f2v);
}

void reorder_front(Front2vof *f2v, Front *fr) {
	int **numpoints = f2v->numpoints;
	int **iselems = f2v->iselems;
	int **numbuffer = f2v->numbuffer;
	int **numpoints_cum = f2v->numpoints_cum;
	int **elem_id = f2v->elem_id;
	int *point_id = f2v->point_id;
	int ndiv = f2v->ndiv;
	int cumulative = 0;
	int I_, J_;
	double delta = f2v->delta;

	//reorder points
	for(int i = -GHOSTS; i<ndiv + GHOSTS; ++i) {
		for(int j = -GHOSTS; j<ndiv + GHOSTS; ++j) {
			numpoints[i][j] = 0;
			numbuffer[i][j] = 0;
		}
	}
	foreach_frontpoint(fr){
		I_ =  (Coords[frontpoint][0]-X0)/delta;
		J_ =  (Coords[frontpoint][1]-Y0)/delta;
		numpoints[I_][J_]++; 
	}
	for(int i = -GHOSTS; i<ndiv + GHOSTS; ++i) {
		for(int j = -GHOSTS; j<ndiv + GHOSTS; ++j) {
			numpoints_cum[i][j] = cumulative;
			cumulative += numpoints[i][j];
		}
	}
	foreach_frontpoint(fr){
		I_ =  (Coords[frontpoint][0]-X0)/delta;
		J_ =  (Coords[frontpoint][1]-Y0)/delta;
		point_id[numpoints_cum[I_][J_] + numbuffer[I_][J_]] = frontpoint; 
		numbuffer[I_][J_]++;
	}
	//reorder elems
	for(int i = -GHOSTS; i<ndiv + GHOSTS; ++i) {
		for(int j = -GHOSTS; j<ndiv + GHOSTS; ++j) {
			iselems[i][j] = 0;
		}
	}
	int I_2, J_2, I_3, J_3;
	double a,b,c; //constants of the line ax+by = c , passing through the endpoints of the element 
	double xipt, yipt; //intercepts
	foreach_frontelement(fr){
		I_  =  (Coords[Corner[frontelement][0]][0] - X0)/delta;
		J_  =  (Coords[Corner[frontelement][0]][1] - Y0)/delta;
		I_2 =  (Coords[Corner[frontelement][1]][0] - X0)/delta;
		J_2 =  (Coords[Corner[frontelement][1]][1] - Y0)/delta;			
		iselems[I_][J_] = 1; 
		elem_id[I_][J_] = frontelement; 
	
		if( (I_ != I_2) || (J_ != J_2) ){
			iselems[I_2][J_2] = 1; 
			elem_id[I_2][J_2] = frontelement; 
		
			//   a*x    +     b*y   =            c 
			//(y0-y1)*x + (x1-x0)*y = (y0-y1)*x0 + (x1-x0)*y0	
			a = (Coords[Corner[frontelement][0]][1] -  Coords[Corner[frontelement][1]][1]);
			b = (Coords[Corner[frontelement][1]][0] -  Coords[Corner[frontelement][0]][0]);
			c = a*Coords[Corner[frontelement][0]][0] + b*Coords[Corner[frontelement][0]][1];

			for(int i = MIN(I_, I_2)+1; i <= MAX(I_, I_2); ++i) {
				yipt = (c - a*(i*delta))/b;
				J_3 = (yipt - Y0)/delta;
				iselems[i][J_3] = 1; 
				elem_id[i][J_3] = frontelement; 
			}
			for(int j = MIN(J_, J_2)+1; j <= MAX(J_, J_2); ++j) {
				xipt = (c - b*(j*delta))/a;
				I_3 = (xipt - X0)/delta;
				iselems[I_3][j] = 1; 
				elem_id[I_3][j] = frontelement; 
			}
		} 
	}
}

void front_to_vof(Front2vof *f2v, Front *fr) {
	reorder_front(f2v_circle, fr_circle);
	foreach() 
		c_refer[] = 0.0;

	int **numpoints_cum = f2v->numpoints_cum;
	int **numpoints = f2v->numpoints;
	int **elem_id = f2v->elem_id;
	int **iselems = f2v->iselems;
	int *point_id = f2v->point_id;
	int ndiv = f2v->ndiv;
	int *list_ = (int *)malloc(20*sizeof(int));	
									//Collection of points/elems in the stencil
	int nlist_;			//number of points/elems in the list
	int I_0, I_, J_0, J_ , I_2, J_2, dLevel;
	int elist_;
	double **coords = fr->PCoords;
	int **nbr = fr->ENbr;
	int **corner = fr->ECorner;
	double **A; A = (double **)malloc(2*sizeof(double *));
	A[0] = (double *)malloc(2*sizeof(double));
	A[1] = (double *)malloc(2*sizeof(double));
	double xc, yc, R, fh;
	double xy[2];

	foreach() {
		if (level > LEVEL) {
			fprintf(stdout, "ERROR ");
			exit(-1);
		}
	
		dLevel = 1 << (LEVEL - level);
		I_0 = point.i - GHOSTS;
		J_0 = point.j - GHOSTS;
		I_  = I_0*dLevel;
		J_  = J_0*dLevel;
		
		nlist_ = 0;
		for(int i = I_; i< I_+dLevel; ++i) {
			for(int j = J_; j< J_+dLevel; ++j){
				if(iselems[i][j]) {
					nlist_ = 1;
					elist_ = elem_id[i][j];
					break;
				}
			}
		}
		if(!nlist_) {
				//c_refer[] = (double) (c_refer[] > 0.5)	;
				c_refer[] = (double) (c_refer[] > 0.5)	;
		}
		else {
		//fprintf(stdout, "\n-- <%d|%d,%d>", level, I_0, J_0);
		//fprintf(stdout, "<<%d>>", elist_);
			/*
			nlist_ = 0;
			for(int i = max(I_-dLevel,-GHOSTS); i< min(I_+2*dLevel,ndiv+GHOSTS); ++i) {
				for(int j = max(J_-dLevel,-GHOSTS); j< min(J_+2*dLevel,ndiv+GHOSTS); ++j){
					for(int k=0; k<numpoints[i][j]; ++k) {
					 list_[nlist_] = point_id[numpoints_cum[i][j]+k]; 
					 nlist_++;
					}
				}
			}

		  //Implement Here: Spherical or Planar method to find volume of frac	
			if (circle_fit(coords, list_, nlist_, A, &xc, &yc, &R) ) {
				xy[0] = x -Delta/2.0;
				xy[1] = y -Delta/2.0;
				vofi_xc = 0.50;
				vofi_yc = 1.00;
				vofi_rc = 0.25;
				fh = Get_fh(circle, xy, Delta, 2, 0);
				c_refer[] = Get_cc(circle, xy, Delta, fh, 2);
				//fprintf(stdout, "|cref %f, c %f ", c_refer[],f[] );
			}
			else {
			*/
			
				list_[0] = elist_;	I_2 = (coords[corner[elist_][0]][0]-X0)/Delta; J_2 = (coords[corner[elist_][0]][1]-Y0)/Delta;				
				nlist_ = 1;
				
				//left direction
				while((I_2 == I_0) && (J_2 == J_0)) {
					elist_ = nbr[elist_][0];
					list_[nlist_] = elist_;
					nlist_++;
					I_2 = (coords[corner[elist_][0]][0]-X0)/Delta; J_2 = (coords[corner[elist_][0]][1]-Y0)/Delta;
					//fprintf(stdout,"**<%d|(%d,%d)] ",elist_, I_2, J_2);
				}
				
				//right direction
				elist_ = list_[0]; I_2 = (coords[corner[elist_][1]][0]-X0)/Delta; J_2 = (coords[corner[elist_][1]][1]-Y0)/Delta;
				while((I_2 == I_0) && (J_2 == J_0)) {
					elist_ = nbr[elist_][1];
					list_[nlist_] = elist_;
					nlist_++;
					I_2 = (coords[corner[elist_][1]][0]-X0)/Delta; J_2 = (coords[corner[elist_][1]][1]-Y0)/Delta;
					//fprintf(stdout,"**<%d|(%d,%d)] ",elist_, I_2, J_2);
				}
		
				xy[0] = x -Delta/2.0;
				xy[1] = y -Delta/2.0;
				c_refer[] = 1.0 - piecewise_f2v(list_, nlist_, corner, coords, xy, Delta);
		
				//fprintf(stdout," <%f,%f>",c_refer[],f[]);	
				
				//c_refer[] = f[];//Get_cc(circle, xy, Delta, fh, 2);
			//}
		}
	}
	free(A[0]); 
	free(A[1]); 
	free(A);
	free(list_); 
}

event front_default (i = 0) {
#if TREE
  //for (scalar c in interfaces)
    c_refer.refine = c_refer.prolongation = fraction_refine;
#endif
}

event init (i = 0) {
  fraction (c_refer, sq(x - 0.5) + sq(y -1.0) - sq(0.25));
	/*adding circle to front*/
	fr_circle = create_empty_front();
	add_circle(fr_circle, 0.5, 1.0, 0.250, 400);
	dump_front();
	restore_front();
	f2v_circle = create_empty_front2vof();
	//front_to_vof(f2v_circle, fr_circle);
	
	/**
	Directory to write dump files*/
	DIR* dir = opendir("dumpdirectory");
	if (dir) {
		//fprintf(stdout, "\nDirectory exists and accesible");
    closedir(dir);
	} else if (ENOENT == errno) {
		fprintf(stdout, "\nDirectory doesn't exist; Exiting..");
		exit(-1);	
	} else {
		fprintf(stdout, "\nDirectory 'dumpdirectory' exists But can't open. Exiting..");
		exit(-1);	
	}
	/**
	Directory to write front_data files*/
	dir = opendir("pythondirectory");
	if (dir) {
		//fprintf(stdout, "\nDirectory exists and accesible");
    closedir(dir);
	} else if (ENOENT == errno) {
		fprintf(stdout, "\nDirectory 'pythondirectory' doesn't exist; Exiting..");
		exit(-1);	
	} else {
		fprintf(stdout, "\nDirectory 'pythondirectory' exists But can't open. Exiting..");
		exit(-1);	
	}
}

/**
Advancing front*/

event advance_front(i++) {

	Front *fr = fr_circle;
	int     is, js;
	double  pcoord[2], refcoord[2], wgt;

	double xmarker[2];
	double umarker, vmarker;
	double ustencil[5][5], vstencil[5][5];

	//Point point;
	
	foreach_frontpoint(fr) {
		
		xmarker[0] = Coords[frontpoint][0];
		xmarker[1] = Coords[frontpoint][1];
	
		//Locate the cell that contain the marker point	
		Point point = locate(xmarker[0],xmarker[1]);

		is = (xmarker[0] > x) ;
		js = (xmarker[1] > y) ;
		
		for (int ii = 0; ii<5; ++ii) {
			for (int jj = 0; jj<5; ++jj) {
				ustencil[ii][jj] = u.x[ii-2,jj-2];
				vstencil[ii][jj] = u.y[ii-2,jj-2];
			}	
		}
		
		umarker = 0.0;
		vmarker = 0.0;
		refcoord[0] = (xmarker[0] - x)/Delta;
		refcoord[1] = (xmarker[1] - y)/Delta;
	
		//Find velocity at the location of marker from grid velocity
		for(int i=0; i<4; ++i) {
			for(int j=0; j<4; ++j) {

				pcoord[0] = (double) is+i-2;		
				pcoord[1] = (double) js+j-2;		
				
				wgt = w_ijk(pcoord, refcoord);
 				umarker += wgt* ustencil[is+i][js+j];
 				vmarker += wgt* vstencil[is+i][js+j];
			}
		}
	
		//Update location of marker point	
		Coords[frontpoint][0] +=  umarker*dt;	
		Coords[frontpoint][1] +=  vmarker*dt;	

	}
	
	fprintf(stdout, "\n	EVENT(%f) :Advance Front ---", t);
}

event event_front_to_vof(i++) {
	//front_to_vof(f2v_circle, fr_circle);
	fprintf(stdout, "\n	EVENT(%f) :Front to Vof ---", t);
}


event free_front_memory(t=end) {
	if(fr_circle) free_front(fr_circle);
	if(f2v_circle) free_front2vof(f2v_circle);
}

event regrid(i+=10){
	regrid_front(fr_circle);
	//write_front_stdout(fr_circle);
}

/**
Output instance_files reqd for analysis*/
event write_outputfiles(t+= 0.1){
	write_front_python(t, fr_circle);
	char Dfile[40],DFfile[40]; 
	sprintf(Dfile, "dumpdirectory/dump_t%.1f", t);
	sprintf(DFfile, "dumpdirectory/dumpf_t%.1f", t);
	dump(Dfile);
	dump_front(DFfile);
}
