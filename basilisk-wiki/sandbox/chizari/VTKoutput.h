/**
The purpose is to export a vtk ascii file for using in Paraview. There are two main functions, one for the contour plot "output_paraview_CC" and the other for the interface "output_paraview_IF".

These functions may be used as follow,

	char name[500];
	sprintf(name, "%s-%.4f-CC.vtk", "out", t);
	output_paraview_CC(name, t, 90.0, {f, u.x, u.y});
        sprintf(name, "%s-%.4f-F1.vtk", "out", t);
        output_paraview_IF(name, t, 90.0, f);
*/

/**
Creating the VTK format for the contour plot will be done by the following function, note that the "rotationangle" variable is for rotating the whole area CCW.
*/

void output_paraview_CC(char *name, double time, double rotationangle, scalar *outlist)
{
	scalar *list = dump_list(outlist);
	;
	int i, cellnumber;
	//	int varNO = list_len(list);
	double xxx[4], yyy[4];
	const double angle = rotationangle * R_PI / 180.0;
	FILE *fp;
	scalar vartmp[];
	;
	cellnumber = 0;
	foreach ()
		cellnumber++;
	;
	fp = fopen(name, "w");
	fprintf(fp, "# vtk DataFile Version 2.0\r\n");
	fprintf(fp, "Unstructured Grid\r\n");
	fprintf(fp, "ASCII\r\n");
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\r\n");
	fprintf(fp, "POINTS %d float\r\n", cellnumber * 4);
	foreach_leaf()
	{
		xxx[0] = x - 0.50 * Delta;
		xxx[1] = x + 0.50 * Delta;
		xxx[2] = x + 0.50 * Delta;
		xxx[3] = x - 0.50 * Delta;
		yyy[0] = y - 0.50 * Delta;
		yyy[1] = y - 0.50 * Delta;
		yyy[2] = y + 0.50 * Delta;
		yyy[3] = y + 0.50 * Delta;
		fprintf(fp, "%.9f %.9f 0.0\r\n", xxx[0] * cos(angle) - yyy[0] * sin(angle), yyy[0] * cos(angle) + xxx[0] * sin(angle));
		fprintf(fp, "%.9f %.9f 0.0\r\n", xxx[1] * cos(angle) - yyy[1] * sin(angle), yyy[1] * cos(angle) + xxx[1] * sin(angle));
		fprintf(fp, "%.9f %.9f 0.0\r\n", xxx[2] * cos(angle) - yyy[2] * sin(angle), yyy[2] * cos(angle) + xxx[2] * sin(angle));
		fprintf(fp, "%.9f %.9f 0.0\r\n", xxx[3] * cos(angle) - yyy[3] * sin(angle), yyy[3] * cos(angle) + xxx[3] * sin(angle));
	}
	fprintf(fp, "CELLS %d %d\r\n", cellnumber, cellnumber + cellnumber * 4);
	for (i = 0; i < cellnumber; i++)
		fprintf(fp, "%d %d %d %d %d\r\n", 4, i * 4 + 0, i * 4 + 1, i * 4 + 2, i * 4 + 3);
	fprintf(fp, "CELL_TYPES %d\r\n", cellnumber);
	for (i = 0; i < cellnumber; i++)
		fprintf(fp, "7\r\n");
	fprintf(fp, "CELL_DATA %d\r\n", cellnumber);
	for (scalar s in list)
	{
		if (strcmp(s.name, "cm"))
		{
			fprintf(fp, "SCALARS %s float\r\n", s.name);
			fprintf(fp, "LOOKUP_TABLE default\r\n");
			foreach_leaf()
				fprintf(fp, "%f\r\n", s[]);
		}
	};
	fclose(fp);
	return;
}

/**
The following function will find the intersection points for 2D.
*/

int findintersectionpoints(double nx, double ny, double alpha, double xc, double yc, double delta, double xi[2], double yi[2], char typei[2])
{
	int ppp = 0;
	double xtmp[2], ytmp[2], underflow = 1.0e-6;
	if (fabs(nx) < underflow)
	{
		ytmp[0] = (alpha - nx * (-0.50)) / ny;
		ytmp[1] = (alpha - nx * (+0.50)) / ny;
		xi[ppp] = xc + (-0.50) * delta;
		yi[ppp] = yc + (ytmp[0]) * delta;
		typei[ppp] = 'l';
		(ppp)++;
		xi[ppp] = xc + (+0.50) * delta;
		yi[ppp] = yc + (ytmp[1]) * delta;
		typei[ppp] = 'r';
		(ppp)++;
	}
	else if (fabs(ny) < underflow)
	{
		xtmp[0] = (alpha - ny * (-0.50)) / nx;
		xtmp[1] = (alpha - ny * (+0.50)) / nx;
		xi[ppp] = xc + (xtmp[0]) * delta;
		yi[ppp] = yc + (-0.50) * delta;
		typei[ppp] = 'b';
		(ppp)++;
		xi[ppp] = xc + (xtmp[1]) * delta;
		yi[ppp] = yc + (+0.50) * delta;
		typei[ppp] = 't';
		(ppp)++;
	}
	else
	{
		xtmp[0] = (alpha - ny * (-0.50)) / nx;
		xtmp[1] = (alpha - ny * (+0.50)) / nx;
		ytmp[0] = (alpha - nx * (-0.50)) / ny;
		ytmp[1] = (alpha - nx * (+0.50)) / ny;

		if (-0.50 <= ytmp[0] && ytmp[0] <= +0.50)
		{
			xi[ppp] = xc + (-0.50) * delta;
			yi[ppp] = yc + (ytmp[0]) * delta;
			typei[ppp] = 'l';
			(ppp)++;
		}
		if (-0.50 <= xtmp[0] && xtmp[0] <= +0.50)
		{
			xi[ppp] = xc + (xtmp[0]) * delta;
			yi[ppp] = yc + (-0.50) * delta;
			typei[ppp] = 'b';
			(ppp)++;
		}
		if (-0.50 <= ytmp[1] && ytmp[1] <= +0.50)
		{
			xi[ppp] = xc + (+0.50) * delta;
			yi[ppp] = yc + (ytmp[1]) * delta;
			typei[ppp] = 'r';
			(ppp)++;
		}
		if (-0.50 <= xtmp[1] && xtmp[1] <= +0.50)
		{
			xi[ppp] = xc + (xtmp[1]) * delta;
			yi[ppp] = yc + (+0.50) * delta;
			typei[ppp] = 't';
			(ppp)++;
		}
	}
	return true;
}

/** And the VTK format of the interface will be, */

void output_paraview_IF(char *name, double time, double rotationangle, scalar intrfc)
{
	int interfacepoints, cellnumber, i;
	char typei[2], *typeintersect;
	double xi[2], yi[2], *xintersect, *yintersect;
	FILE *fp;
	const double angle = rotationangle * R_PI / 180.0;
	scalar alpha[];
	vector n[];
	;
	cellnumber = 0;
	foreach ()
		cellnumber++;
	xintersect = (double *)calloc(sizeof(double), cellnumber);
	yintersect = (double *)calloc(sizeof(double), cellnumber);
	typeintersect = (char *)calloc(sizeof(char), cellnumber);
	;
	reconstruction(intrfc, n, alpha);
	interfacepoints = 0;
	foreach_leaf()
	{
		if (intrfc[] > R_VOFLIMIT && intrfc[] < 1.0 - R_VOFLIMIT)
		{
			findintersectionpoints(n.x[], n.y[], alpha[], x, y, Delta, xi, yi, typei);
			xintersect[interfacepoints] = xi[0] * cos(angle) - yi[0] * sin(angle);
			yintersect[interfacepoints] = yi[0] * cos(angle) + xi[0] * sin(angle);
			typeintersect[interfacepoints] = typei[0];
			interfacepoints++;
			xintersect[interfacepoints] = xi[1] * cos(angle) - yi[1] * sin(angle);
			yintersect[interfacepoints] = yi[1] * cos(angle) + xi[1] * sin(angle);
			typeintersect[interfacepoints] = typei[1];
			interfacepoints++;
		}
	};
	fp = fopen(name, "w");
	fprintf(fp, "# vtk DataFile Version 2.0\r\n");
	fprintf(fp, "Unstructured Grid\r\n");
	fprintf(fp, "ASCII\r\n");
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\r\n");
	switch (interfacepoints)
	{
	case 0:
	{
		fprintf(fp, "POINTS %d float\r\n", 2);
		foreach_leaf()
		{
			fprintf(fp, "%.9f %.9f 0.0\r\n", x, y);
			fprintf(fp, "%.9f %.9f 0.0\r\n", x, y);
			break;
		}
		fprintf(fp, "CELLS %d %d\r\n", 1, 3);
		fprintf(fp, "%d %d %d\r\n", 2, 0, 1);
		fprintf(fp, "CELL_TYPES %d\r\n", 1);
		fprintf(fp, "3\r\n");
		break;
	}
	default:
	{
		fprintf(fp, "POINTS %d float\r\n", interfacepoints);
		for (i = 0; i < interfacepoints; i++)
			fprintf(fp, "%.9f %.9f 0.0\r\n", xintersect[i], yintersect[i]);
		fprintf(fp, "CELLS %d %d\r\n", interfacepoints / 2, interfacepoints + interfacepoints / 2);
		for (i = 0; i < interfacepoints / 2; i++)
			fprintf(fp, "%d %d %d\r\n", 2, i * 2 + 0, i * 2 + 1);
		fprintf(fp, "CELL_TYPES %d\r\n", interfacepoints / 2);
		for (i = 0; i < interfacepoints / 2; i++)
			fprintf(fp, "3\r\n");
		break;
	}
	};
	free(xintersect);
	free(yintersect);
	fclose(fp);
	return;
}