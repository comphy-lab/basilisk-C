/**
# Second order structure function
The functions within this file aim to calculate a second-order longitudional structure function from a vector field. We can define a distance vector $\mathbf{l}$ separating two points in space $(\mathbf{x}_1,\ \mathbf{x}_1+\mathbf{l})$. Generally speaking, the vectors $\mathbf{u}$ at these points are different. The difference in the magnitude of their projections on the $\mathbf{l}$ vector may have some interesting properties. We define:  
   
   $$\delta v_{\|}(\mathbf{x}_1,\mathbf{l})= \left( \mathbf{v}(\mathbf{x}_1)-\mathbf{v}(\mathbf{x}_1)\right) \cdot \frac{\mathbf{l}}{l},$$

where $l$ is $\|\mathbf{l}\|$. By definition, from a statistical perspective, within a homogeneous and isotropic region, one would expect $\|\delta v_{\|}(\mathbf{x_1,l})\|$ only to be dependend on $l$, i.e. the length of the separation vector. Therefore we can define a second order statistic related to a vector field ($\mathbf{u}$) according to:

$$S_2(l)=\langle \delta v_{\|}(\mathbf{x}_1,\mathbf{l})^2 \rangle$$ 

where the variable $\langle x \rangle$ represents a region-averaged value of a dummy variable $x$. The code belows aims to evaluate some discretized version of $S_2(l)$ from a vector field on a three-dimensional-tree grid.

There are many sub-functions so that maybe later, it may be modified/extended to e.g. 2D and 1D.

First, a function is defined that selects `j` cells for analysis and stores the relevant data in an array.
*/

int obtain_a_cache(int j,double d[][6],int lev[],vector u){
#if _MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  srand(time(NULL));// In order to obtain a smooth-ish pair distances distribution we add a stogastic component;
  int index[j];
  int n=0;
  foreach()
    n++;
  double meaninc=(double)(n)/(double)(j); //The mean increment to obtain approx j cells.
  index[0]=0;
  int g=1;
  while (g<j){
    double randsteprel =((sign(noise())*pow((rand()/(double)RAND_MAX),0.2))+1.0);
    int inc = (int)((1.0*randsteprel*meaninc)+1.0);
    index[g]=index[g-1]+inc;
    g++;
  }
  n=0;
  int jj=0;
  foreach(){
    if(n==index[jj]){
      d[jj][0]=x; d[jj][1]=y; d[jj][2]=z;
      d[jj][3]=u.x[]; d[jj][4]=u.y[]; d[jj][5]=u.z[];
      lev[jj]=level;
      jj++;
    }
    n++;
  }
  printf("there are %d in cache and %d requested\n",jj,j);
  return jj; //Return jj (<=j) the number of cells in the cache
}
/**
Since this function is not necessarily cheap, each thread does part of the calculations. This function sets the variables that determine what point combinations need to be assessed by the specific thread.
*/
void set_g_and_gg(int ggg[2], int j){
  int total = j*(j-1)/2;
  int start = (pid()*total/npe())+1;
  int m = 0;
  for (int g=0;g<j;g++){
    for (int gg=g+1;gg<j;gg++){
      m++;
      if (m==start){
	ggg[0]=g;
	ggg[1]=gg;
        printf("%d %d %d %d %d %d\n",ggg[0],ggg[1],total,pid(),start,j);
	return;
      }
    }
  }
}
/**
This function actually calculates the the structure function for each thread.
*/
int structure2(double cach[][6],int lev[],int j,double L,int Len,double dvarr[][4]){
  int g=0;
  int gg=1;
#if _MPI
  int ggg[2];
  set_g_and_gg(ggg,j);
  g=ggg[0];
  gg=ggg[1];
#endif
  
  int iter=1;
  int todo=((j*(j-1))/2)/npe();
  double Del = L/(double)Len;
  double l;
  while (iter<=todo){
    l=pow(sq(cach[g][0]-cach[gg][0])+
	  sq(cach[g][1]-cach[gg][1])+
	  sq(cach[g][2]-cach[gg][2]),0.5);
    if (l>0.&&l<L0){ //distance within range
      double weight = (double)(1<<(3*(2*depth()-lev[gg]-lev[g]))); //weigh with cells' volumes
      double dv = (sq((cach[g][3]-cach[gg][3])*(cach[g][0]-cach[gg][0]))+
		   sq((cach[g][4]-cach[gg][4])*(cach[g][1]-cach[gg][1]))+
		   sq((cach[g][5]-cach[gg][5])*(cach[g][2]-cach[gg][2])))/(sq(l));
      int index = (int)((l/Del)+0.5);
      if  (index<Len){
	dvarr[index][1]+=dv*weight;
	dvarr[index][2]+=weight;
	dvarr[index][3]+=1.;
      }
    }
    gg++; 
    if (gg==j){
      g++;
      gg=g+1;
    }
    iter++;
  }
  return 1;
}
/**
This function can be used to print the result to a file
*/
void printdvarr(FILE * fp, double dvarr[][4], int Len){
  int it = 0;
  while (it<Len){
    for (int gg=0;gg<4;gg++)
      fprintf(fp,"%g\t",dvarr[it][gg]);
    fprintf(fp,"\n");
    it++;
  }
}
/**
This function comminucates the points that each thread found in the `obtain_a_cache` function. Note that this requires all-to-all communication, a feature that is warranted by the non-local nature of the statistic.
*/
#if _MPI
int combine_MPI_caches(double cach[][6],int lev[],int j){
  int tot=0;
  MPI_Allreduce(&j,&tot,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  double cachsend[j*6];
  double * cachrecv= (double *) malloc(sizeof(double)* tot * 6);
  int  levrecv[tot];
  for (int g=0;g<j;g++){
    for (int gg=0;gg<6;gg++)
      cachsend[6*g+gg]=cach[g][gg];
  }
  MPI_Allgather(&cachsend[0],6*j,MPI_DOUBLE,cachrecv,6*j,MPI_DOUBLE,MPI_COMM_WORLD);
  MPI_Allgather(&lev[0],j,MPI_INT,&levrecv[0],j,MPI_INT,MPI_COMM_WORLD);
  for (int g=0;g<tot;g++){
    lev[g]=levrecv[g];
    for (int gg=0;gg<6;gg++)
      cach[g][gg]=cachrecv[6*g+gg];
  }
  return tot;
  
}
#endif
/**
## User interface function 

The function requires a file pointer to write the result, the vector field, the number of cells that will be used ($n$, so that there are $\mathcal{O}(n^2)$ queries), the maximum separation distance, and the number of bins for the discretized represenation of $l$.    
*/

int long2structure(FILE * fp,vector u,int n,double L,int Len){
  double Del = L/(double)Len;
  double dvarr[Len][4];
  int it=0;
  while (it<Len){
    for (int gg=0;gg<4;gg++)
      dvarr[it][gg]=0.;
    it++;
  }
  
  double cach[n][6];
  int lev[n];
  int locpnts = ((n/npe())+0.5);
  int j=obtain_a_cache(locpnts,cach,lev,u);
#if _MPI
  int d=j;
  j=combine_MPI_caches(cach,lev,d);
#endif
  structure2(cach,lev,j,L,Len,dvarr);
  it=0;
#if _MPI
  double dvarrrecv[Len][4]; // An array to store the reduced version of *dvarr*, according to MPI_SUM
  MPI_Reduce(&dvarr[0][0],&dvarrrecv[0][0],Len*4,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if (pid()==0){
    while (it<Len){
      for (int gg=0;gg<4;gg++)
	dvarr[it][gg]=dvarrrecv[it][gg];
      it++;
    }
    it=0;
#endif
    double aa;
    while (it<Len){
      if (dvarr[it][2]>0.){
	aa = dvarr[it][1]/dvarr[it][2];
	dvarr[it][1] = aa;
      }
      dvarr[it][0]=it*Del;
      it++;
    }
    printdvarr(fp,dvarr,Len);
    
#if _MPI
  }
#endif
  return 1; 
}
/**
## test
* [Structure function of linearly varying vector fields](test_sf.c)

## Usage

* [Les of isotropic turbulence](isotropicLES.c)
*/