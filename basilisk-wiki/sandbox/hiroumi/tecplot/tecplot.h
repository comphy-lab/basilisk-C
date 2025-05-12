/**
This header file is included to output scalars at cell center in Tecplot binary format. The coordinates (i.e. x, y, z) are given at vertex.
*/

#include <TECIO.h>
struct OutputTec{
  scalar * tec_cc;
  char varname[400];
  char extname[400];
#if CANTERA
  scalar * list_SP;
#endif
};

#if CANTERA
extern char **sname;
#endif

void output_tec(struct OutputTec tec,int i){
   INTEGER4 nNodes,nCells,
            idebug,itec,dIsDouble,vIsDouble,zoneType,
            strandID,parentZn,isBlock,
            nFConns, fNMode, shrConn, fileType, nFaces,
            connectivityCount, iCellMax, jCellMax, kCellMax;
   INTEGER4 fileFormat;
   vertex scalar idv = new scalar; scalar idc = new scalar;
  
//Give ID number to every vertex. 
   int id=0;
   foreach_vertex() {
     #pragma omp critical (node)
     {
     id+=1;
     idv[]= id;
     }
   }
   nNodes=id;

//Give ID number to every cell center. 
   id =0;
   foreach() {
     #pragma omp critical (center)
     {
     idc[]= id;
     id+=1;
     }
   }
   nCells=id;

   int numcc =list_len(tec.tec_cc), numspc=0;
#if CANTERA
   numspc=list_len(tec.list_SP);
#endif
// "dimension" corresponds to the nubmer of the coordinate axes, numcc is the number of variables at cell center, and numspc is that of chemical species if defined with Cantera library.
   int valueLocation[dimension + numcc + numspc];
   for(int k=0;k<dimension;k++)valueLocation[k]=1; //x,y (,z) are given at vertexes.
   for(int k=dimension;k<dimension+numcc+numspc;k++)valueLocation[k]=0; //Scalars are given at cell centers.
   int *connectivity;
   double *xn,*yn,*sn,*sc;
   #if dimension ==3
   double *zn;
   #endif

   fileFormat = 0;
   fileType   = 0;
   idebug     = 0;
   vIsDouble  = 1;

   #if dimension==2
   zoneType  = 3;      /* FEQuadrilateral */
   #elif dimension==3
   zoneType  = 5;      /* FEBrick */
   #endif
   nFaces    = 0;
   iCellMax  = 0;
   jCellMax  = 0;
   kCellMax  = 0;
   strandID  = 1;     /* StaticZone */
   parentZn  = 0;      /* No Parent */
   isBlock   = 1;      /* Block */
   nFConns   = 0;
   fNMode    = 0;
   dIsDouble = 1;
   shrConn   = 0;
#if CANTERA
   for (int k=0;k<numspc;k++){
       strcat(tec.varname," ");
       strcat(tec.varname,sname[k]);
   }
#endif

   strcat(tec.varname,tec.extname); // Add the name of extra variables.
   char filename[40];
  
//The name of output file is tecdata_step + "number of step" _zone + "pid" + .plt.
   sprintf(filename,"tecdata_step%05d_zone%04d.plt",i,pid());

   itec = TECINI142((char*)"BASHILISK DATASET",
                    tec.varname,
                    filename,
                    (char*)".",
                    &fileFormat,
                    &fileType,
                    &idebug,
                    &vIsDouble);

   itec = TECZNE142((char*)"Primary Zone",
                    &zoneType,
                    &nNodes,
                    &nCells,
                    &nFaces,
                    &iCellMax,
                    &jCellMax,
                    &kCellMax,
                    &t,
                    &strandID,
                    &parentZn,
                    &isBlock,
                    &nFConns,
                    &fNMode,
                    0,              /* TotalNumFaceNodes */
                    0,              /* NumConnectedBoundaryFaces */
                    0,              /* TotalNumBoundaryConnections */
                    0,//NULL,           /* PassiveVarList */
                    valueLocation,  /* ValueLocation = Nodal */
                    0,//NULL,           /* SharVarFromZone */
                    &shrConn);

   xn  = (double*)malloc(nNodes * sizeof(double));
   yn  = (double*)malloc(nNodes * sizeof(double));
   #if dimension>2
   zn  = (double*)malloc(nNodes * sizeof(double));
   #endif
   foreach_vertex() {
      int idvv = (int)idv[]-1;
      xn[idvv]=x;
      yn[idvv]=y;
      #if dimension>2
      zn[idvv]=z;
      #endif
   }

   itec = TECDAT142(&nNodes, xn, &dIsDouble);
   itec = TECDAT142(&nNodes, yn, &dIsDouble);
   #if dimension>2
   itec = TECDAT142(&nNodes, zn, &dIsDouble);
   #endif
   free(xn);
   free(yn);
   #if dimension>2
   free(zn);
   #endif

   sc  = (double*)malloc(nCells * sizeof(double));
   for (scalar s in tec.tec_cc){
     foreach() {
        int idcc = (int)idc[];
        sc[idcc]=s[];
     }
     itec = TECDAT142(&nCells, sc, &dIsDouble);
   }
   free(sc);

   #if dimension==2
   connectivityCount = 4 * nCells;
   #elif dimension==3
   connectivityCount = 8 * nCells;
   #endif
   connectivity = (INTEGER4*)malloc(connectivityCount * sizeof(INTEGER4));

   //Specify the connectivity between the IDs of cell center and vertex.
   foreach(){
     INTEGER4 itecc=(INTEGER4)idc[];
     #if dimension==2
     connectivity[itecc*4]   =(INTEGER4)idv[0,0];
     connectivity[itecc*4+1] =(INTEGER4)idv[1,0];
     connectivity[itecc*4+2] =(INTEGER4)idv[1,1];
     connectivity[itecc*4+3] =(INTEGER4)idv[0,1];
     #elif dimension==3
     connectivity[itecc*8]   =(INTEGER4)idv[0,0,0];
     connectivity[itecc*8+1] =(INTEGER4)idv[1,0,0];
     connectivity[itecc*8+2] =(INTEGER4)idv[1,1,0];
     connectivity[itecc*8+3] =(INTEGER4)idv[0,1,0];
     connectivity[itecc*8+4] =(INTEGER4)idv[0,0,1];
     connectivity[itecc*8+5] =(INTEGER4)idv[1,0,1];
     connectivity[itecc*8+6] =(INTEGER4)idv[1,1,1];
     connectivity[itecc*8+7] =(INTEGER4)idv[0,1,1];
     #endif
   }

   itec = TECNODE142(&connectivityCount, connectivity);
   free(connectivity);
   delete({idv});
   delete({idc});
   itec = TECEND142();
}