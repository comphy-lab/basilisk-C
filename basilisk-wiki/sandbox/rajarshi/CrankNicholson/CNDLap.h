/**
# Crank Nicholson RHS calculator

   This Header file has a function which is used to solve the Crank Nicholson Scheme for a time evolving equation with a Laplacian
   The function is used in addition with the Poisson-Helmholtz solver. 
*/


void DLap(scalar * DataP, scalar * CLaplacianP, double Diff, double dt) {

   scalar Data = DataP[0], CLaplacian = CLaplacianP[0];
   face vector GradP[], GradAvg[];

   foreach_face()
        GradP.x[] = (Data[-2,0]-27*Data[-1,0]+27*Data[]-Data[1,0])/(24*Delta);
   boundary((scalar *){GradP});

   foreach_face()
      GradAvg.x[] = (-17*GradP.x[0,-2] + 308*GradP.x[0,-1] + 5178*GradP.x[0,0] + 308*GradP.x[0,1] - 17*GradP.x[0,2])/5760;
   
   foreach()
      CLaplacian[] = -( (GradAvg.x[1,0]-GradAvg.x[0,0]+GradAvg.y[0,1]-GradAvg.y[0,0])/Delta + (2./(Diff*dt))*Data[] );

   boundary({CLaplacian});
}