/** Page containing all diagnostics functions and events of the Krabbendijke case. Note that the imported output_silces.h can be found [here](http://basilisk.fr/sandbox/vheusinkveld/myfunctions/output/sliceOf3D.h) */

#include "utils.h"
#if dimension == 3
	#include "lambda2.h"
#endif
#include "output_slices.h"

/** Define variables and structures to do: diagnostics, ouput data, output movies. */

struct sDiag dia; 			// Diagnostics

struct sDiag {
	double Ekin;			// Total kinetic energy
	double EkinOld;			// Track changes in kin energy 
	double WdoneOld;		// Track changes in work done 
	double rotVol;			// Diagnosed rotor volum
	double bE;			// Buoyancy energy
	double bEold;			// Track changes in buoyancy energy
	double diss;			// Diagnose dissipation
};

struct sEquiDiag {
    int level;			        // Level for which diagnostics should be done
    int ii; 			        // Keep track of how many additions are done
    double dtDiag;
    double startDiag;
    double endDiag;
    double dtOutput;
};

// Create variables to make a folder structure
struct sOutput {
    double dtDiag;
    double dtVisual;
    double dtSlices;
    double dtProfile;
    double startAve;
    double dtAve;
    char main_dir[12];
    char dir[30];
    char dir_profiles[60];
    char dir_slices[60];
    char dir_equifields[60];
    char dir_strvel[60];
    char dir_refvel[60];
    char dir_diffbins[60];
    char dir_dts[60];
    int sim_i;
};
	
struct sbViewSettings {
    double phi; 			// Phi for 3D bview movie
    double theta;			// Theta for 3D bview movie
    double sphi; 			// Polar angle for sliced image RED
    double stheta;			// Azimuthal angle for sliced image RED
};

/** Initialize structures */
struct sOutput out = {.dtDiag = 1., .dtVisual=1., .dtSlices=5., .dtProfile=60., .main_dir="results", .sim_i=0};

struct sEquiDiag ediag = {.level = 5, .ii = 0, .startDiag = 0., .dtDiag = 0., .dtOutput = 0.};

struct sbViewSettings bvsets = {.phi=0., .theta=0., .sphi=0., .stheta=0.};

event init(i = 0){
	bvsets.phi = 0.;
	bvsets.theta = -M_PI/6.;
	bvsets.sphi = 0.;
	bvsets.stheta = 0.;
}

/** Diagnosing: kinetic energy, diagnosed rotor volume, buoyancy energy, ammount of cells used.*/
event diagnostics (t+=out.dtDiag){
	int n = 0.;
	scalar ekin[]; 		// Kinetic energy field 
	double tempVol = 0.;    // Temp volume 
	double tempEkin = 0.;   // Temp kinetic energy
	double tempDiss = 0.;   // Temp dissipation 
	double maxVel = 0.;     // Maximum velocity in fan
	double bEnergy = 0.;    // Buoyant energy
		
	/** Loop over cells to get diagnostics */ 
	foreach(reduction(+:n) reduction(+:tempVol) reduction(+:tempEkin) 
		reduction(max:maxVel) reduction(+:bEnergy) reduction(+:tempDiss)) {

		tempVol += dv()*fan[];
		if(y + Delta/2. <= rot.y0){
		     bEnergy += dv()*y*(b[] - STRAT(y));
		}
		foreach_dimension() {
			ekin[] += sq(u.x[]);
		}
		maxVel = max(maxVel, sq(ekin[]));
		ekin[] *= 0.5*rho[]*dv();	
		tempEkin += ekin[];
		n++;
	}
	/** Assign values to respective global sturcture vars */ 
	dia.diss = 1.*tempDiss;
	dia.bE = 1.*bEnergy;
	dia.rotVol = 1.*tempVol;
	dia.Ekin = 1.*tempEkin;
	
	if (pid() == 0){
	/** Write away simulation data and case setup for main thread */ 
	char nameOut[90];
	char nameCase[90];
     	snprintf(nameOut, 90, "./%s/output", out.dir);
     	snprintf(nameCase, 90, "./%s/case", out.dir);
	static FILE * fpout = fopen(nameOut, "w");
	static FILE * fpca = fopen(nameCase, "w");

	if(t==0.){
		fprintf(fpca,"L0\tinversion\thubU\tTref\tLambda\txr\tyr\tzr\ttheta\tphit\tr\tW\tP\tcu\trampT\tmaxlvl\tminlvl\teps\n");
		fprintf(fpca, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%g\n", 
				L0,TREF/gCONST*STRAT(rot.y0), WIND(rot.y0),TREF, Lambda, rot.x0, rot.y0, rot.z0, rot.theta, rot.phit, rot.R, rot.W, rot.P, rot.cu, rot.rampT, maxlevel, minlevel, eps);
		
	        fprintf(stderr,"n\tred\tEkin\tWork\tbE\n");
		fprintf(fpout,"i\tt\tn\tred\tEkin\tWork\tbE\n");
	}
	fprintf(fpout, "%d\t%g\t%d\t%g\t%g\t%g\t%g\n",
		i,t,n,(double)((1<<(maxlevel*3))/n),dia.Ekin,rot.Work, dia.bE);
	
	fprintf(stderr, "n=%d\t%g\t%g\t%g\t%g\n",n,(double)((1<<(maxlevel*dimension))/n),dia.Ekin,rot.Work,dia.bE);
	
	fflush(fpout);
	fflush(fpca);	

	}

	dia.EkinOld = 1.*dia.Ekin;
	dia.WdoneOld = 1.*rot.Work;
	dia.bEold = 1.*dia.bE;
}

#if dimension == 3
/** diagnose velocity of the passing wind machine*/
event refvelocity(t+=1) {
    if(pid() == 0) {
    char nameRefvel[90];
    snprintf(nameRefvel, 90, "%st=%05g", out.dir_refvel, t);
    FILE * fpstr = fopen(nameRefvel, "w");
    fprintf(fpstr, "x,v,vx,vy,vz\n");
        
    double length = 300.;
    int ntot = 300.;
    double vels1[ntot];	

    double xf0 = rot.x0 - length/2*cos(M_PI/6);
    double yf0 = 3;
    double zf0 = rot.z0 + length/2*sin(M_PI/6);

    for(int n = 0; n < ntot; n++) {
	double dist = length*n/ntot;
	double xx = xf0 + dist*cos(-M_PI/6); 
	double yy = yf0; 
	double zz = zf0 + dist*sin(-M_PI/6); 

	double valx1 = interpolate(u.x, xx, yy, zz);
	double valy1 = interpolate(u.y, xx, yy, zz);
	double valz1 = interpolate(u.z, xx, yy, zz);
         
        vels1[n] = sqrt(sq(valx1) + sq(valy1) + sq(valz1));
	fprintf(fpstr, "%g,%g,%g,%g,%g,%g\n", xx, zz, vels1[n], valx1, valy1, valz1);     
    }  
    fclose(fpstr); 
    }
}
/** 'DTS' measurments*/
event dts_meas(t += 1) {
    if(pid() == 0) {
	char nameDtshorS[90];
    	snprintf(nameDtshorS, 90, "%shorS_t=%05g", out.dir_dts, t);
        FILE * fpstrhorS = fopen(nameDtshorS, "w");
        fprintf(fpstrhorS, "x,y,z,b\n");
        
	double lengthhorS = 600.;
        int ntothorS = 600;
	
	double xf0S = rot.x0 - lengthhorS/2*cos(M_PI/6);
	double yf0S = 1;
	double zf0S = rot.z0 + lengthhorS/2*sin(M_PI/6);

	for(int n = 0; n <= ntothorS; n++) {
	    double dist = lengthhorS*n/ntothorS;
	    double xx = xf0S + dist*cos(-M_PI/6); 
	    double yy = yf0S; 
	    double zz = zf0S + dist*sin(-M_PI/6); 
       
	    double valb = interpolate(b, xx, yy, zz);
          
	    fprintf(fpstrhorS, "%g,%g,%g,%g\n", xx, yy, zz, valb);     
	}  
 	fclose(fpstrhorS); 

	char nameDtshorL[90];
    	snprintf(nameDtshorL, 90, "%shorL_t=%05g", out.dir_dts, t);
        FILE * fpstrhorL = fopen(nameDtshorL, "w");
        fprintf(fpstrhorL, "x,y,z,b\n");
        
	double lengthhorL = 600.;
        int ntothorL = 600;
	
	double xf0L = rot.x0 - lengthhorL/2*cos(120*M_PI/180);
	double yf0L = 1;
	double zf0L = rot.z0 + lengthhorL/2*sin(120*M_PI/180);

	for(int n = 0; n <= ntothorL; n++) {
	    double dist = lengthhorL*n/ntothorL;

	    double xx = xf0L + dist*cos(-120*M_PI/180); 
	    double yy = yf0L; 
	    double zz = zf0L + dist*sin(-120*M_PI/180); 
       
	    double valb = interpolate(b, xx, yy, zz);
          
	    fprintf(fpstrhorL, "%g,%g,%g,%g\n", xx, yy, zz, valb);     
	}  
 	fclose(fpstrhorL); 

	char nameDtsver[90];
    	snprintf(nameDtsver, 90, "%sverS_t=%05g", out.dir_dts, t);
        FILE * fpstrver = fopen(nameDtsver, "w");
        fprintf(fpstrver, "x,y,z,b\n");
        
	double lengthver = 20.;
        int ntotver = 160;

	for(int n = 0; n <= ntotver; n++) {
	    double dist = lengthver*n/ntotver;

	    double xx = rot.x0 + 30*cos(-M_PI/6); 
	    double yy = dist; 
	    double zz = rot.z0 + 30*sin(-M_PI/6); 
       
	    double valb = interpolate(b, xx, yy, zz);
          
	    fprintf(fpstrver, "%g,%g,%g,%g\n", xx, yy, zz, valb);     
	}  
 	fclose(fpstrver); 
    }
    
}


#endif

/** Ouputting movies in 2- or 3-D*/

#if dimension == 2
event movies(t += 0.5) {
    vertex scalar omega[]; 	// Vorticity
    scalar lev[];	 	// Grid depth
    scalar ekinRho[]; 		// Kinetic energy

    foreach() {
        omega[] = ((u.y[1,0] - u.y[-1,0]) - (u.x[0,1] - u.x[0,-1]))/(2*Delta); // Curl(u) 
        ekinRho[] = 0.5*rho[]*(sq(u.x[]) + sq(u.y[]));
        lev[] = level;
    }

    boundary ({b, lev, omega, ekinRho});
    output_ppm (b, file = "ppm2mp4 ./results/buoyancy.mp4", n = 1<<maxlevel, linear = true, max=STRAT(L0), min=STRAT(0.));
    output_ppm (ekinRho, file = "ppm2mp4 ./results/ekin.mp4", n = 1<<maxlevel, min = 0, max = 1.*sq(rot.cu));
    output_ppm (omega, file = "ppm2mp4 ./results/vort.mp4", n = 1<<maxlevel, linear = true); 
    output_ppm (lev, file = "ppm2mp4 ./results/grid_depth.mp4", n = 1<<maxlevel, min = minlevel, max = maxlevel);
}
#elif dimension == 3
event movies(t += out.dtVisual) {
    scalar l2[];

    lambda2(u,l2);
    boundary({l2});

    view(fov=25, tx = 0., ty = 0., phi=M_PI/2., theta=M_PI/2., width = 800, height = 800);
    
    translate(-rot.x0,-rot.y0,-rot.z0) {
        box(notics=false);
        isosurface("l2", v=-0.02, color="b", min=STRAT(0.9*roughY0h), max=STRAT(1.1*rot.y0));
	draw_vof("fan", fc = {1,0,0});
    }
   
    translate(-rot.x0,0,-rot.z0){
        squares("u.y", n = {0,1,0}, alpha=2, min=-1, max=1);
    }

    /** Save file with certain fps*/
    char nameVid4[90];
    snprintf(nameVid4, 90, "ppm2mp4 -r %g ./%s/visual_3d_ux.mp4", 10., out.dir);
    save(nameVid4);
    clear();

    view(fov=25, tx = 0., ty = 0., phi=M_PI/2., theta=M_PI/2., width = 800, height = 800);
    
    translate(-rot.x0,-rot.y0,-rot.z0) {
        box(notics=false);
        isosurface("l2", v=-0.02, color="b", min=STRAT(0.9*roughY0h), max=STRAT(1.1*rot.y0));
	draw_vof("fan", fc = {1,0,0});
    }
   
    translate(-rot.x0,0,-rot.z0){
        squares("u.y", n = {0,1,0}, alpha=2, min=-.05, max=.05);
    }

    /** Save file with certain fps*/
    char nameVid3[90];
    snprintf(nameVid3, 90, "ppm2mp4 -r %g ./%s/visual_3d_uy.mp4", 10., out.dir);
    save(nameVid3);
    clear();

    view(fov=25, tx = 0., ty = 0., phi=M_PI/2., theta=M_PI/2., width = 800, height = 800);
    
    translate(-rot.x0,-rot.y0,-rot.z0) {
        box(notics=false);
        isosurface("l2", v=-0.02, color="b", min=STRAT(roughY0h), max=STRAT(2.*rot.y0));
	draw_vof("fan", fc = {1,0,0});
    }
   
    translate(-rot.x0,0,-rot.z0){
        squares("b", n = {0,1,0}, alpha=2, min=0.5*STRAT(2), max=1.5*STRAT(2));
    }

    /** Save file with certain fps*/
    char nameVid1[90];
    snprintf(nameVid1, 90, "ppm2mp4 -r %g ./%s/visual_3d_b2.mp4", 10., out.dir);
    save(nameVid1);
    clear();

    view(fov=25, tx = 0., ty = 0., phi=bvsets.phi, theta=bvsets.theta, width = 1200, height = 1200);
    
    translate(-rot.x0,-rot.y0,-rot.z0) {
        box(notics=false);
        isosurface("l2", v=-0.02, color="b", min=STRAT(0.), max=STRAT(1.5*rot.y0));
	draw_vof("fan", fc = {1,0,0});
    }
    translate(-rot.z0,-rot.y0, -L0){
      	squares("u.x", n = {0,0,1}, alpha=rot.z0, min=-1, max=1);
        cells(n = {0,0,1}, alpha = rot.z0);
    }
    
    translate(0.,-rot.y0,-rot.z0){
        squares("b", n = {1,0,0}, alpha=rot.x0, min=STRAT(0.), max=STRAT(1.5*rot.y0));
    }

    /** Save file with certain fps*/
    char nameVid2[90];
    snprintf(nameVid2, 90, "ppm2mp4 -r %g ./%s/visual_3d_b.mp4", 10., out.dir);
    save(nameVid2);
    clear();
}

/** Take relevant field slices and write away */
event slices(t=out.dtSlices; t+=out.dtSlices) {
    char nameSlice[90];
    coord slice = {1., 0., 1.};
    int res = L0/2;

    for(double yTemp = 0.5; yTemp<=1; yTemp+=0.5) {
	slice.y = yTemp/L0;
   
    	snprintf(nameSlice, 90, "%st=%05gy=%03g", out.dir_slices, t, yTemp);
    	FILE * fpsli = fopen(nameSlice, "w");
    	output_slice(list = (scalar *){b}, fp = fpsli, n = res, linear = true, plane=slice);
    	fclose(fpsli);
    }

    for(double yTemp = 2; yTemp<=4.; yTemp+=2.) {
	slice.y = yTemp/L0;
   
    	snprintf(nameSlice, 90, "%st=%05gy=%03g", out.dir_slices, t, yTemp);
    	FILE * fpsli = fopen(nameSlice, "w");
    	output_slice(list = (scalar *){b}, fp = fpsli, n = res, linear = true, plane=slice);
    	fclose(fpsli);
    }


}
#endif

/** Usefull functions */ 

/** Checks if required folders exists, if not they get created. */
void sim_dir_create(){
   
    sprintf(out.dir, "./%s/%s%02d", out.main_dir, sim_ID, out.sim_i);
    sprintf(out.dir_profiles, "%s/profiles/", out.dir);
    sprintf(out.dir_slices, "%s/slices/", out.dir);
    sprintf(out.dir_equifields, "%s/equifields/", out.dir);
    sprintf(out.dir_strvel, "%s/strvel/", out.dir);
    sprintf(out.dir_refvel, "%s/refvel/", out.dir);
    sprintf(out.dir_diffbins, "%s/diffbins/", out.dir);
    sprintf(out.dir_dts, "%s/dts/", out.dir);

    
    if (pid() == 0){
    struct stat st = {0};
    if (stat(out.main_dir, &st) == -1) {
        mkdir(out.main_dir, 0777);
    }
   
    if (stat(out.dir, &st) == -1) {
        mkdir(out.dir, 0777);
    }
    if (stat(out.dir_slices, &st) == -1) {
        mkdir(out.dir_slices, 0777);
    }
    if (stat(out.dir_profiles, &st) == -1) {
        mkdir(out.dir_profiles, 0777);
    }  
    if (stat(out.dir_equifields, &st) == -1) {
        mkdir(out.dir_equifields, 0777);
    }  
    if (stat(out.dir_strvel, &st) == -1) {
	mkdir(out.dir_strvel, 0777);
    }
    if (stat(out.dir_refvel, &st) == -1) {
	mkdir(out.dir_refvel, 0777);
    }
    if (stat(out.dir_diffbins, &st) == -1) {
	mkdir(out.dir_diffbins, 0777);
    }
    if (stat(out.dir_dts, &st) == -1) {
	mkdir(out.dir_dts, 0777);
    }

    }
}


