#ifndef OUTPUT_PVD_H
#define OUTPUT_PVD_H

/**
see [output_htg.h](output_htg.h) on how to use this code
*/

#if _MPI
void output_pvd_mpiio(char* name, double t, MPI_File fp, bool firstTimeWritten){	
	char head[] ="<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n\t<Collection>\n";
	char tail[] = "\t</Collection>\n</VTKFile>\n";
	
	if (firstTimeWritten == true) // Write Head
		MPI_File_write(fp,&head, strlen(head), MPI_CHAR, MPI_STATUS_IGNORE); 
	else
		MPI_File_seek(fp, -strlen(tail), MPI_SEEK_END);
	
	char buffer[100];
	snprintf(buffer,sizeof(buffer),"\t\t<DataSet timestep=\"%g\" group=\"\" part=\"0\" file=\"%s\"/>\n",t, name);
	
	MPI_File_write(fp, &buffer, strlen(buffer), MPI_CHAR, MPI_STATUS_IGNORE); 
	MPI_File_write(fp, &tail, 	strlen(tail), 	MPI_CHAR, MPI_STATUS_IGNORE); 
}
#endif

void output_pvd(char* name, double t, FILE* fp, bool firstTimeWritten){	
	char head[] ="<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n\t<Collection>\n";
	char tail[] = "\t</Collection>\n</VTKFile>\n";
	
	if (firstTimeWritten == true) // Write Head
		fprintf(fp,"%s", head);
	else
		fseek(fp, -strlen(tail), SEEK_END);
	
	fprintf(fp,"\t\t<DataSet timestep=\"%g\" group=\"\" part=\"0\" file=\"%s\"/>\n",t, name);	
	fprintf(fp,"%s",tail);
}

#endif