/**
This header file is calculating isolines of a scalar variable. The main function is "isoline",

int isoline (scalar varB, const double Value, int *loopnumber, int *loopsizes, double ***xyLoops, int *loopleafsizes, int **loopleafindex)

* varB: The scalar variable which the isoline will be calculated from,
* Value: The value that the function will follow to construct the isoline,
* loopnumber: The number of loops of isolines that we may have in the domain, maximum is 100,
* loopsizes[l]: Size of each loop "l" of isolines
* xyLoops[l][p][0]: The points "p" coordinates of isolines for each loop "l", note that the coordinates are the last index, [0] for the x and [1] for the y coordinate,

The last two parameters need some explainations, assuming,

*/

/*
i = 0;
foreach_leaf()
{
    // do some calculation on a cell
    ...
    // finish the calculations
    i++;
}
*/

/**
The "i" in the above example is showing the "leaf-index". Now, the last two parameters on the isoline function are connected to the leaf index as follow,

* loopleafsizes[l]: This is the loop size of the leaf index, note that this might be different from the previous size of the loop since the isoline may have more points than cells,
* loopleafindex[l][j]: This is the leaf-indices of the isoline,

The isoline-cells may be accessed as,
*/

/*
    i = 0;
    foreach_leaf()
    {
        for (l = 0; l < loopnumber; l++)
        {
            for (j = 0; j < loopleafsizes[l]; j++)
            {
                if (i == loopleafindex[l][j])
                {
                    // iso-cell values will be accessed here, such as f[], u.x[], ...
                }
            }
        }
        i++;
    }
*/

/**
A sample of using the isoline is,

*/

/*
		const int LoopNO = 100;
		cellnumber = 0;
		foreach ()
			cellnumber++;
		;
		// allocating variables
		;
		loopsizes = (int *)calloc(LoopNO, sizeof(int));
		xyLoops = (double ***)calloc(LoopNO, sizeof(double **));
		loopleafsizes = (int *)calloc(LoopNO, sizeof(int));
		loopleafindex = (int **)calloc(LoopNO, sizeof(int *));
		for (i = 0; i < LoopNO; i++)
		{
			loopleafindex[i] = (int *)calloc(cellnumber, sizeof(int));
			xyLoops[i] = (double **)calloc(cellnumber, sizeof(double *));
			for (c = 0; c < cellnumber; c++)
				xyLoops[i][c] = (double *)calloc(2, sizeof(double));
		}
		;
		// find isoline points in loop and loopleaf
		;
		isoline (f, 0.50, &loopnumber, loopsizes, xyLoops, loopleafsizes, loopleafindex);
                ;
                // ... data analysis ...
                ;
                // free the memory
                ;
		free (loopsizes);
		free (loopleafsizes);
		for (i = 0; i < LoopNO; i++)
		{
			free (loopleafindex[i]);
			for (c = 0; c < cellnumber; c++)
				free (xyLoops[i][c]);
			free (xyLoops[i]);
			for (c = 0; c < loopnumber; c++)
				free (vars[i][c]);
		}
		free (loopleafindex);
		free (xyLoops);
*/

#define ROUNDOFDMINCOEFFICIENT                          0.0010

int findnextconnection (double ***xyCIsPn, int noCeIsP, char **flag, char **type, double ***xyLoops, double xcurrent, double ycurrent, int loop, int *pointloop, double dmin, int *leafindex, int **loopleafindex)
{
    int l, p, c, e1, iflag;
    ;
    iflag = 0;
    l = loop;
    p = *pointloop;
    for (c = 0; c < noCeIsP; c++)
    {
        if (flag[c][0] == 'n' && flag[c][1] == 'n')
        {
            for (e1 = 0; e1 < 2; e1++)
            {
                if (fabs(xyCIsPn[c][e1][0] - xcurrent) < ROUNDOFDMINCOEFFICIENT * dmin && fabs(xyCIsPn[c][e1][1] - ycurrent) < ROUNDOFDMINCOEFFICIENT * dmin)
                {
                    p++;
                    xyLoops[l][p][0] = xyCIsPn[c][1 - e1][0];
                    xyLoops[l][p][1] = xyCIsPn[c][1 - e1][1];
                    flag[c][e1] = 'y';
                    flag[c][1 - e1] = 'y';
                    loopleafindex[l][p] = leafindex[c];
                    ;
                    e1 = 2;
                    c = noCeIsP;
                    iflag = 1;
                }
            }
        }
    }
    *pointloop = p;
    return iflag;
}

int leafindexremoveduplicates (int loopnumber, int *loopsizes, int *loopleafsizes, int **loopleafindex)
{
    int i, c, itmp, ctmp;
    ;
    for (i = 0; i < loopnumber; i++)
    {
        ctmp = 0;
        itmp = loopleafindex[i][ctmp];
        for (c = 0; c < loopsizes[i]; c++)
        {
            if (itmp == loopleafindex[i][c])
                ;
            else
            {
                ctmp++;
                loopleafindex[i][ctmp] = loopleafindex[i][c];
                itmp = loopleafindex[i][c];
            }
        }
        loopleafsizes[i] = ++ctmp;
    }
    return true;
}

double findlocation (double S0, double S1, double V0, double V1, double Vc)
{
    return S0 + (S1 - S0) / (V1 - V0) * (Vc - V0);
}

int isolineonecell (double *xyDvars, const double Value, double **xyCIsPn)
{
    int i = 0;
    if ( (xyDvars[4] < Value && Value <= xyDvars[5]) || (xyDvars[4] > Value && Value >= xyDvars[5]) )
    {
        xyCIsPn[i][0] = findlocation (xyDvars[0] + 0.50 * xyDvars[2], xyDvars[0] - 0.50 * xyDvars[2], xyDvars[4], xyDvars[5], Value);
        xyCIsPn[i][1] = xyDvars[1] + 0.50 * xyDvars[2];
        i++;
    }
    if ( (xyDvars[5] < Value && Value <= xyDvars[6]) || (xyDvars[5] > Value && Value >= xyDvars[6]) )
    {
        xyCIsPn[i][0] = xyDvars[0] - 0.50 * xyDvars[2];
        xyCIsPn[i][1] = findlocation (xyDvars[1] + 0.50 * xyDvars[2], xyDvars[1] - 0.50 * xyDvars[2], xyDvars[5], xyDvars[6], Value);
        i++;
    }
    if ( (xyDvars[6] < Value && Value <= xyDvars[7]) || (xyDvars[6] > Value && Value >= xyDvars[7]) )
    {
        xyCIsPn[i][0] = findlocation (xyDvars[0] - 0.50 * xyDvars[2], xyDvars[0] + 0.50 * xyDvars[2], xyDvars[6], xyDvars[7], Value);
        xyCIsPn[i][1] = xyDvars[1] - 0.50 * xyDvars[2];
        i++;
    }
    if ( (xyDvars[7] < Value && Value <= xyDvars[4]) || (xyDvars[7] > Value && Value >= xyDvars[4]) )
    {
        xyCIsPn[i][0] = xyDvars[0] + 0.50 * xyDvars[2];
        xyCIsPn[i][1] = findlocation (xyDvars[1] - 0.50 * xyDvars[2], xyDvars[1] + 0.50 * xyDvars[2], xyDvars[7], xyDvars[4], Value);
        i++;
    }
    return i;
}

int isolineallcells (double **xyDvars, int cellno, const double Value, double ***xyCIsPn, double dmin, int *leafindex)
{
    int c, ciso, pcell;
    ciso = 0;
    for (c = 0; c < cellno; c++)
    {
        pcell = isolineonecell (xyDvars[c], Value, xyCIsPn[ciso]);
        leafindex[ciso] = leafindex[c];
        if (pcell == 2)
        {
            if (fabs(xyCIsPn[ciso][0][0] - xyCIsPn[ciso][1][0]) >= ROUNDOFDMINCOEFFICIENT * dmin || fabs(xyCIsPn[ciso][0][1] - xyCIsPn[ciso][1][1]) >= ROUNDOFDMINCOEFFICIENT * dmin)
                ciso++;
        }
    }
    return ciso;
}

int findminmax (double *var, int NO, double *min, double *max)
{
    int i;
    *min = var[0];
    for (i = 1; i < NO; i++)
        if (*min > var[i])
            *min = var[i];
    *max = var[0];
    for (i = 1; i < NO; i++)
        if (*max < var[i])
            *max = var[i];
    return true;
}

int findcellpoints (scalar varB, double xc, double yc, double DL, double *varI)
{
    int i;
    double xxx[5], yyy[5];
    xxx[0] = xc;
    xxx[1] = xc + DL * 0.50;
    xxx[2] = xc - DL * 0.50;
    xxx[3] = xc - DL * 0.50;
    xxx[4] = xc + DL * 0.50;
    yyy[0] = yc;
    yyy[1] = yc + DL * 0.50;
    yyy[2] = yc + DL * 0.50;
    yyy[3] = yc - DL * 0.50;
    yyy[4] = yc - DL * 0.50;
    for (i = 0; i < 5; i++)
    {
        varI[i] = interpolate (varB, xxx[i], yyy[i]);
        if (varI[i] == nodata)
            return false;
    }
    return true;
}

int separatemidendpoints (double ***xyCIsPn, int noCeIsP, char **flag, char **type, double dmin)
{
    int i, j, e1, e2;
    ;
    for (i = 0; i < noCeIsP; i++)
    {
        for (e1 = 0; e1 < 2; e1++)
        {
            for (j = i + 1; j < noCeIsP; j++)
            {
                for (e2 = 0; e2 < 2; e2++)
                {
                    if (fabs(xyCIsPn[i][e1][0] - xyCIsPn[j][e2][0]) < ROUNDOFDMINCOEFFICIENT * dmin && fabs(xyCIsPn[i][e1][1] - xyCIsPn[j][e2][1]) < ROUNDOFDMINCOEFFICIENT * dmin)
                    {
                        type[i][e1] = 'm';
                        type[j][e2] = 'm';
                        ;
                        e2 = 2;
                        j = noCeIsP;
                    }
                }
            }
        }
    }
    return true;
}

int connectisolinepoints (double ***xyCIsPn, int noCeIsP, char **flag, char **type, double ***xyLoops, int *loopsizes, double dmin, int *leafindex, int **loopleafindex)
{
    int c, l, p, e1;
    char cflag;
    ;
    separatemidendpoints (xyCIsPn, noCeIsP, flag, type, dmin);
    l = -1;
    for (c = 0; c < noCeIsP; c++)
    {
        p = 0;
        cflag = 'n';
        if (flag[c][0] == 'n' && flag[c][1] == 'n')
        {
            for (e1 = 0; e1 < 2; e1++)
            {
                if (type[c][e1] == 'e' && type[c][1 - e1] == 'm')
                {
                    l++;
                    ;
                    xyLoops[l][p][0] = xyCIsPn[c][e1][0];
                    xyLoops[l][p][1] = xyCIsPn[c][e1][1];
                    flag[c][e1] = 'y';
                    loopleafindex[l][p] = leafindex[c];
                    ;
                    p++;
                    xyLoops[l][p][0] = xyCIsPn[c][1 - e1][0];
                    xyLoops[l][p][1] = xyCIsPn[c][1 - e1][1];
                    flag[c][1 - e1] = 'y';
                    loopleafindex[l][p] = leafindex[c];
                    ;
                    cflag = 'y';
                    e1 = 2;
                }
            }
        }
        switch (cflag)
        {
        case 'y':
        {
            while (findnextconnection (xyCIsPn, noCeIsP, flag, type, xyLoops, xyLoops[l][p][0], xyLoops[l][p][1], l, &p, dmin, leafindex, loopleafindex) == 1)
                ;
            switch (p)
            {
            case 0:
            {
                l--;
                break;
            }
            default:
            {
                loopsizes[l] = p;
                break;
            }
            }
            break;
        }
        }
    }
    return l + 1;
}

int checkisolineonecell (double *xyDvars, const double Value)
{
    int i = 0;
    if ( (xyDvars[4] < Value && Value <= xyDvars[5]) || (xyDvars[4] > Value && Value >= xyDvars[5]) )
        i++;
    if ( (xyDvars[5] < Value && Value <= xyDvars[6]) || (xyDvars[5] > Value && Value >= xyDvars[6]) )
        i++;
    if ( (xyDvars[6] < Value && Value <= xyDvars[7]) || (xyDvars[6] > Value && Value >= xyDvars[7]) )
        i++;
    if ( (xyDvars[7] < Value && Value <= xyDvars[4]) || (xyDvars[7] > Value && Value >= xyDvars[4]) )
        i++;
    return i;
}

int splitacelltofour (scalar varB, const double dmin, const double Value, double **xyDvars, int c, int *cnplus, int *leafindex)
{
    int i, ctmp, pctmp[4], cs, cnplustmp;
    double vtmp[5], xyDtmp[4][8];
    ;
    ctmp = 0;
    xyDtmp[ctmp][0] = xyDvars[c][0] + 0.250 * xyDvars[c][2];
    xyDtmp[ctmp][1] = xyDvars[c][1] + 0.250 * xyDvars[c][2];
    xyDtmp[ctmp][2] = 0.50 * xyDvars[c][2];
    findcellpoints (varB, xyDtmp[ctmp][0], xyDtmp[ctmp][1], xyDtmp[ctmp][2], vtmp);
    for (i = 0; i < 5; i++)
        xyDtmp[ctmp][i + 3] = vtmp[i];
    pctmp[ctmp] = checkisolineonecell (xyDtmp[ctmp], Value);
    ;
    ctmp = 1;
    xyDtmp[ctmp][0] = xyDvars[c][0] - 0.250 * xyDvars[c][2];
    xyDtmp[ctmp][1] = xyDvars[c][1] + 0.250 * xyDvars[c][2];
    xyDtmp[ctmp][2] = 0.50 * xyDvars[c][2];
    findcellpoints (varB, xyDtmp[ctmp][0], xyDtmp[ctmp][1], xyDtmp[ctmp][2], vtmp);
    for (i = 0; i < 5; i++)
        xyDtmp[ctmp][i + 3] = vtmp[i];
    pctmp[ctmp] = checkisolineonecell (xyDtmp[ctmp], Value);
    ;
    ctmp = 2;
    xyDtmp[ctmp][0] = xyDvars[c][0] - 0.250 * xyDvars[c][2];
    xyDtmp[ctmp][1] = xyDvars[c][1] - 0.250 * xyDvars[c][2];
    xyDtmp[ctmp][2] = 0.50 * xyDvars[c][2];
    findcellpoints (varB, xyDtmp[ctmp][0], xyDtmp[ctmp][1], xyDtmp[ctmp][2], vtmp);
    for (i = 0; i < 5; i++)
        xyDtmp[ctmp][i + 3] = vtmp[i];
    pctmp[ctmp] = checkisolineonecell (xyDtmp[ctmp], Value);
    ;
    ctmp = 3;
    xyDtmp[ctmp][0] = xyDvars[c][0] + 0.250 * xyDvars[c][2];
    xyDtmp[ctmp][1] = xyDvars[c][1] - 0.250 * xyDvars[c][2];
    xyDtmp[ctmp][2] = 0.50 * xyDvars[c][2];
    findcellpoints (varB, xyDtmp[ctmp][0], xyDtmp[ctmp][1], xyDtmp[ctmp][2], vtmp);
    for (i = 0; i < 5; i++)
        xyDtmp[ctmp][i + 3] = vtmp[i];
    pctmp[ctmp] = checkisolineonecell (xyDtmp[ctmp], Value);
    ;
    cs = 0;
    cnplustmp = *cnplus;
    for (ctmp = 0; ctmp < 4; ctmp++)
    {
        if (pctmp[ctmp] == 2 && cs == 0)
        {
            for (i = 0; i < 8; i++)
                xyDvars[c][i] = xyDtmp[ctmp][i];
            cs++;
            if (xyDvars[c][2] > 1.50 * dmin)
                splitacelltofour (varB, dmin, Value, xyDvars, c, cnplus, leafindex);
        }
        else if (pctmp[ctmp] == 2)
        {
            cnplustmp = *cnplus;
            for (i = 0; i < 8; i++)
                xyDvars[cnplustmp][i] = xyDtmp[ctmp][i];
            leafindex[cnplustmp] = leafindex[c];
            (*cnplus)++;
            if (xyDvars[cnplustmp][2] > 1.50 * dmin)
                splitacelltofour (varB, dmin, Value, xyDvars, cnplustmp, cnplus, leafindex);
        }
    }
    ;
/*    ctmp = *cnplus;
    xyDvars[ctmp][0] = xyDvars[c][0] - 0.250 * xyDvars[c][2];
    xyDvars[ctmp][1] = xyDvars[c][1] + 0.250 * xyDvars[c][2];
    xyDvars[ctmp][2] = 0.50 * xyDvars[c][2];
    findcellpoints (varB, xyDvars[ctmp][0], xyDvars[ctmp][1], xyDvars[ctmp][2], vtmp);
    for (i = 0; i < 5; i++)
        xyDvars[ctmp][i + 3] = vtmp[i];
    (*cnplus)++;
    if (xyDvars[ctmp][2] > 1.50 * dmin)
        splitacelltofour (varB, dmin, Value, xyDvars, ctmp, cnplus);
    ;
    ctmp = *cnplus;
    xyDvars[ctmp][0] = xyDvars[c][0] - 0.250 * xyDvars[c][2];
    xyDvars[ctmp][1] = xyDvars[c][1] - 0.250 * xyDvars[c][2];
    xyDvars[ctmp][2] = 0.50 * xyDvars[c][2];
    findcellpoints (varB, xyDvars[ctmp][0], xyDvars[ctmp][1], xyDvars[ctmp][2], vtmp);
    for (i = 0; i < 5; i++)
        xyDvars[ctmp][i + 3] = vtmp[i];
    (*cnplus)++;
    if (xyDvars[ctmp][2] > 1.50 * dmin)
        splitacelltofour (varB, dmin, Value, xyDvars, ctmp, cnplus);
    ;
    ctmp = *cnplus;
    xyDvars[ctmp][0] = xyDvars[c][0] + 0.250 * xyDvars[c][2];
    xyDvars[ctmp][1] = xyDvars[c][1] - 0.250 * xyDvars[c][2];
    xyDvars[ctmp][2] = 0.50 * xyDvars[c][2];
    findcellpoints (varB, xyDvars[ctmp][0], xyDvars[ctmp][1], xyDvars[ctmp][2], vtmp);
    for (i = 0; i < 5; i++)
        xyDvars[ctmp][i + 3] = vtmp[i];
    (*cnplus)++;
    if (xyDvars[ctmp][2] > 1.50 * dmin)
        splitacelltofour (varB, dmin, Value, xyDvars, ctmp, cnplus);
    ;
    ctmp = *cnplus;
    xyDvars[c][0] = xyDvars[c][0] + 0.250 * xyDvars[c][2];
    xyDvars[c][1] = xyDvars[c][1] + 0.250 * xyDvars[c][2];
    xyDvars[c][2] = 0.50 * xyDvars[c][2];
    findcellpoints (varB, xyDvars[c][0], xyDvars[c][1], xyDvars[c][2], vtmp);
    for (i = 0; i < 5; i++)
        xyDvars[c][i + 3] = vtmp[i];
    if (xyDvars[c][2] > 1.50 * dmin)
        splitacelltofour (varB, dmin, Value, xyDvars, c, cnplus);
*/    ;
    return true;
}

int convertallcellsizestosmallest (scalar varB, const double Value, double **xyDvars, int cellno, int *cnplus, double dmin, int *leafindex)
{
    int c;
    *cnplus = cellno;
    for (c = 0; c < cellno; c++)
        if (1.50 * dmin < xyDvars[c][2])
            splitacelltofour (varB, dmin, Value, xyDvars, c, cnplus, leafindex);
    return true;
}

int isoline (scalar varB, const double Value, int *loopnumber, int *loopsizes, double ***xyLoops, int *loopleafsizes, int **loopleafindex)
{
    int i, c, l, cellno, cellnototal, cnplus, noCeIsP, *leafindex;
    char **flag, **type;
    double **xyDvars, varI[5], varImin, varImax, ***xyCIsPn, dmin;
    ;
    cellnototal = 0;
    foreach()
        cellnototal++;
    ;
    leafindex = (int *)calloc(cellnototal, sizeof(int));
    xyDvars = (double **)calloc(cellnototal, sizeof(double *));
    xyCIsPn = (double ***)calloc(cellnototal, sizeof(double **));
    for (i = 0; i < cellnototal; i++)
    {
        xyDvars[i] = (double *)calloc(8, sizeof(double));
        xyCIsPn[i] = (double **)calloc(4, sizeof(double *));
        for (c = 0; c < 4; c++)
            xyCIsPn[i][c] = (double *)calloc(2, sizeof(double));
    }
    ;
    // find cells in the interested area
    ;
    c = 0;
    l = 0;
    foreach_leaf()
    {
        if (findcellpoints (varB, x, y, Delta, varI))
        {
            findminmax (varI, 5, &varImin, &varImax);
            if (varImin <= Value && Value <= varImax)
            {
                xyDvars[c][0] = x;
                xyDvars[c][1] = y;
                xyDvars[c][2] = Delta;
                xyDvars[c][3] = varB[];
                for (i = 1; i < 5; i++)
                    xyDvars[c][i + 3] = varI[i];
                leafindex[c] = l;
                c++;
            }
        }
        l++;
    }
    cellno = c;
    dmin = xyDvars[0][2];
    for (c = 1; c < cellno; c++)
        if (dmin > xyDvars[c][2])
            dmin = xyDvars[c][2];
    ;
    // convert all cells to same size
    ;
    cnplus = cellno;
    convertallcellsizestosmallest (varB, Value, xyDvars, cellno, &cnplus, dmin, leafindex);
    printf ("old cells: %d, new cells: %d\r\n", cellno, cnplus);
    cellno = cnplus;
    ;
    // isoline calculation for each cell
    ;
    noCeIsP = isolineallcells (xyDvars, cellno, Value, xyCIsPn, dmin, leafindex);
    printf ("total iso-line cells: %d\r\n", noCeIsP);
    ;
    // connect isoline points
    ;
    flag = (char **)calloc(noCeIsP, sizeof(char *));
    type = (char **)calloc(noCeIsP, sizeof(char *));
    for (i = 0; i < noCeIsP; i++)
    {
        flag[i] = (char *)calloc(2, sizeof(char));
        type[i] = (char *)calloc(2, sizeof(char));
    }
    ;
    for (i = 0; i < noCeIsP; i++)
    {
        flag[i][0] = 'n';
        flag[i][1] = 'n';
        type[i][0] = 'e';
        type[i][1] = 'e';
    }
    *loopnumber = connectisolinepoints (xyCIsPn, noCeIsP, flag, type, xyLoops, loopsizes, dmin, leafindex, loopleafindex);
    printf ("number of loops: %d\r\n", *loopnumber);
    for (i = 0; i < *loopnumber; i++)
        printf ("loop points: %d\r\n", loopsizes[i]);
	leafindexremoveduplicates (*loopnumber, loopsizes, loopleafsizes, loopleafindex);
    ;
    for (i = 0; i < cellnototal; i++)
        free(xyDvars[i]);
    free(xyDvars);
    ;
    for (i = 0; i < cellnototal; i++)
    {
        for (c = 0; c < 4; c++)
            free(xyCIsPn[i][c]);
        free(xyCIsPn[i]);
    }
    free(xyCIsPn);
    ;
    for (i = 0; i < noCeIsP; i++)
    {
        free (flag[i]);
        free (type[i]);
    }
    free (flag);
    free (type);
    ;
    free (leafindex);
    ;
    return noCeIsP;
}
