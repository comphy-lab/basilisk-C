#ifndef OUTPUT_PVD_H
#define OUTPUT_PVD_H
/** Example Usage:
 	if(pid()==0) {
		bool firstTimeWritten = false;
		char pvd_name[]="pvd_name.pvd";	  
		fp = fopen(pvd_name, "r+");
		if( (i == 0) ||  (fp == NULL) ) {
			fp = fopen(pvd_name,"w");
			firstTimeWritten = true;
		}	 
		output_pvd(vtk_name, t, fp, firstTimeWritten);
		fclose(fp);
	}
 */
void output_pvd(char* name, double t, FILE* fp, bool firstTimeWritten){	
	char head[] ="<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n\t<Collection>\n";
	char tail[] = "\t</Collection>\n</VTKFile>\n";
	
	if (firstTimeWritten == true) // Write Head
		fprintf(fp,"%s", head);
	else // Overwrite old Tail with new Timestep and Tail
		fseek(fp, -strlen(tail), SEEK_END);
	
	fprintf(fp,"\t\t<DataSet timestep=\"%g\" group=\"\" part=\"0\" file=\"%s\"/>\n",t, name);	
	fprintf(fp,"%s",tail);
}
#endif