/**
# Rain term in saint venant

## Rain source term
The rain is simply added to the mass equation of Saint-Venant
*/

// Intensity of the rainfall
scalar rain[];
double maxrain;

event updaterain( i++ ){
	foreach(){
	
		h[] += dt*rain[];
	}
	boundary({h});
}


/** 
## Read_Lameeau function
This function reads the "lame d'eau" given by
Meteo-france. THe arguments of the routine are : 

* name : the name of the .asc files 

* start : the starting hour. format hhmm or hh

* end : the ending hour : format hhmm or hh

* step : the step between two files in minutes 

* t : the time
 
* mult : a multiplicator of the readen file 

* linear : (true/false) activating the bilinear interpolation between two points 
* timelin : (true/flase) activating the linear interpolation between two times

 */
#if TREE //input.h does not work without quadtree, so does read_lameau
#include "input.h"


static void coarsen_maxrain (Point point, scalar s){
  s[] = max(max(fine(s,0,0),fine(s,1,0)),max(fine(s,0,1),fine(s,1,1)));
}


event defaults( i = 0 ){
 rain.refine = rain.prolongation = refine_injection;
 rain.coarsen = coarsen_maxrain;
}

struct InputRain {
  char * name;
  int start;
  int end;
  int step;
  double t;
  double mult;
  bool linear;
  bool timelin;
};

scalar rain1[], rain2[];
void read_lameeau(struct InputRain p){
	// Default value
  	int time;
  	if( p.t == 0 )
  		time = (int)t;
  	else
  		time = (int)p.t;
  	if( p.mult == 0 )
  		p.mult = 1;
  	if( p.step == 0 )
  		p.step = 5*60;
  

  	int hstart = 0 , mstart = 0; 
  	if ( p.start <= 24 ){   // just hour
  		hstart = p.start;
  	}
 	else if ( p.start < 2400 && p.start > 100 ) {
 	 	hstart = p.start/100;
 	 	mstart = p.start%100;
 	}
 	else 
 		fprintf(stderr,"#warning : wrong hour format in read_lameeau routine\n#use hhmm format or hh format only\n");
  
 	time = time + hstart*3600 + mstart*60 + 5 ;
  	int timea = time + hstart*3600 + mstart*60 + 5;
  
  	int sec = time%60;
  	time = (time - sec)/60;
 	int min = time%60 ;
 	time = (time - min)/60;
  	int hour = time%24 ;
  
  	if( min%5 != 0 )
  		min = min - min%5;

	char str[100], str1[100];
	strcpy(str,p.name);
	if(hour < 10)  //add 0 to the right place
	  strcat(str,"0");
	sprintf(str1, "%d", hour);
	strcat(str,str1);
	if(min < 10) // again
	  strcat(str,"0");
	sprintf(str1, "%d", min); // here we read the next lameeau
	strcat(str,str1);
	strcat(str,".asc");
	
  // Reading files
	FILE * fpo = fopen(str,"r");
	if (fpo == NULL)
		fprintf(stderr,"# skiping file that does not exist : %s\n",str);
	else{
		input_grd(rain1, fp = fpo, linear = p.linear );
		foreach()
			rain1[] *= p.mult*1e-3/(5*60);
		if( p.timelin ) {
			strcpy(str,p.name);
			if( min + 5 >= 60 )
				min -= 60, hour++;
			if(hour < 10)  //add 0 to the right place
				strcat(str,"0");
	
			sprintf(str1, "%d", hour);
			strcat(str,str1);
			if(min + 5 < 10) // again
				strcat(str,"0");
			sprintf(str1, "%d", min + 5);
			strcat(str,str1);
			strcat(str,".asc");


			fpo = fopen(str,"r");
			if (fpo == NULL)
				fprintf(stderr,"# skiping file that does not exist : %s\n",str);
			else{
				input_grd(rain2, fp = fpo, linear = p.linear );
				foreach()
					rain2[] *= p.mult*1e-3/(5*60);
			}
		}
	}
	 
	 if( p.timelin ) {
	 	 int tt = timea%(5*60);
	 	 foreach()
	 	 	rain[] = rain1[] + tt*(rain2[]-rain1[])/(5*60); // Linear interpolation in time
	 }
	 else{
	 	 foreach()
	 	 	rain[] =rain1[];
	 }
}

#endif

/**
## Link to the homepage
* [Homepage](Readme)
*/

