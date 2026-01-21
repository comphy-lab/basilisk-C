event materialloss( i += 500) {

	foreach(reduction(+:totalphi1)  )  {

		#ifdef AXI
			totalphi1 += 2.*pi*y*(1.0-phi1[])*sq(Delta);
  
		#else
			totalphi1 += (1.0-phi1[])*pow(Delta, dimension);
		#endif
	}

   fp1 = fopen("mass.txt","a");
   fp2 = fopen("masspercent.txt","a");

   fprintf(fp1,"%g  %.16e  \n",t,totalphi1);
   fprintf(fp2,"%g  %.16e  \n",t,totalphi1/totalphi10);
   fclose(fp1); 
   fclose(fp2);

  assert (totalphi1 > 1.e-14);

  totalphi1 = 0.;

}  

event interface_ice(i += 50) {
    double  xmin=50000.,xmax=0.,ymin=50000.,ymax=0.;
	double kl,xminf=0.,xmaxf=0.,yminf=0.,ymaxf=0.;
	
	
		foreach(serial)
	{   
		if((phi1[]-0.5)*( phi1[1,0] - 0.5)<=0.) {
		  kl = (phi1[1,0] - phi1[]) / Delta;
		xmaxf = (0.5 - phi1[]) / kl + x;
		if(xmaxf > xmax) xmax = xmaxf;
		}
		
		if((phi1[]-0.5)*( phi1[1,0] - 0.5)<=0.) {
		  kl = (phi1[1,0] - phi1[]) / Delta;
		  xminf = (0.5 - phi1[]) / kl + x;
		if(xminf < xmin) xmin = xminf;
		}
		
		if((phi1[]-0.5)*( phi1[0,1] - 0.5)<=0.) {
		  kl = (phi1[0,1]  - phi1[]) / Delta;
		  ymaxf = (0.5 - phi1[]) / kl + y;
		if(ymaxf > ymax) ymax = ymaxf;
		}
		
		if((phi1[]-0.5)*( phi1[0,1]  - 0.5)<=0.) {
		  kl = (phi1[0,1]  - phi1[]) / Delta;
		  yminf = (0.5 - phi1[]) / kl + y;
		if(yminf < ymin) ymin = yminf;
		}	

		
	}
		if(xmin>1000. ) xmin=0.0;
		if(ymin>1000. ) ymin=0.0;

	fp1 = fopen("0interface.plt", "a+");
	if (t<0.0000001) {fprintf(fp1, "variables = t,xmax, xmin,r1,ymax,ymin,r2  \n");}
	fprintf(fp1,"%e, %e, %e, %e , %e, %e, %e \n",t,xmax, xmin,fabs(xmin-xmax)/2.0,ymax,ymin,fabs(ymax-ymin)/2.0);
	fclose(fp1);
}

event writingFiles ( t = dtoutput ;t += dtoutput ;t<=endoutput) {
  dump (file = "dump");
  char nameOut[80];
  sprintf (nameOut, "snapshot-%5.4f", t);
  scalar *my_list={u,p,phi1,T};
  dump (file = nameOut,list=my_list);
}

void wirte_ASCII_of_str(char * str, FILE * file)
{
int value = 0;
while ((*str) != '\0')
{
value = (int)*str;
fwrite(&value, sizeof(int), 1, file);
str++;
}
//char null_char[] = "";
//value = (int)*null_char;
value = 0;
fwrite(&value, sizeof(int), 1, file);
}