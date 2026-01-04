/**
#Thinning of a scalar field

##Restriction/prolongation for booleans

*/
int layers = 1;
scalar skeleton[];

static inline void myrestrict(Point point, scalar s){
  double sum = 0.;
  foreach_child()
    if(s[])sum=1.;
  s[] = sum;
}

static inline void myprolongation(Point point, scalar s){
  double sum = s[];
  foreach_child()
    s[] = sum;
}

struct thinningStruct{
  scalar c;
  int thickness; // optional
};
/**
##2D

Method taken from [Guo and Hall, 1989](#guo1989parallel). 
*/

void thinningIteration1(scalar c, bool * modif){

  scalar cp[];

  c.restriction = cp.restriction = myrestrict;
  c.prolongation = cp.prolongation = myprolongation;

  foreach(){
    cp[] = c[];
  }
  boundary({cp});
  restriction({cp});

  foreach(){
    int p0 = (int) cp[ -1, -1];
    int p1 = (int) cp[ -1,  0];
    int p2 = (int) cp[ -1,  1];
    int p3 = (int) cp[  0,  1];
    int p4 = (int) cp[  1,  1];
    int p5 = (int) cp[  1,  0];
    int p6 = (int) cp[  1, -1];
    int p7 = (int) cp[  0, -1];
    int C = ((1-p1) & ( p2 | p3)) +
    ((1-p3) & ( p4 | p5)) +
    ((1-p5) & ( p6 | p7)) +
    ((1-p7) & ( p0 | p1));
    if(C == 1 && c[]) {
                                        /// calculate N
      int N1 = (p0 | p1) +
      (p2 | p3) +
      (p4 | p5) +
      (p6 | p7);
      int N2 = (p1 | p2) +
      (p3 | p4) +
      (p5 | p6) +
      (p7 | p0);
      int N = min(N1,N2);
      if ((N == 2) || (N == 3)) {
                                                /// calculate criteria 3
        int c3 = ( p1 | p2 | (1-p4)) & p3;
        if(c3 == 0) {
          c[] = 0; 
          * modif = 1;
        }
      }
    }
  }
  boundary({c});
  restriction({c});
}

void thinningIteration2(scalar c, bool * modif){
  scalar cp[];

  c.restriction = cp.restriction = myrestrict;
  c.prolongation = cp.prolongation = myprolongation;

  foreach(){
    cp[] = c[];
  }
  boundary({cp});
  restriction({cp});

  foreach(){
    int p0 = (int) cp[ -1, -1];
    int p1 = (int) cp[ -1,  0];
    int p2 = (int) cp[ -1,  1];
    int p3 = (int) cp[  0,  1];
    int p4 = (int) cp[  1,  1];
    int p5 = (int) cp[  1,  0];
    int p6 = (int) cp[  1, -1];
    int p7 = (int) cp[  0, -1];
    int C = ((1-p1) & ( p2 | p3)) +
    ((1-p3) & ( p4 | p5)) +
    ((1-p5) & ( p6 | p7)) +
    ((1-p7) & ( p0 | p1));
    if(C == 1 && c[]) {
                                        /// calculate N
      int N1 = (p0 | p1) +
      (p2 | p3) +
      (p4 | p5) +
      (p6 | p7);
      int N2 = (p1 | p2) +
      (p3 | p4) +
      (p5 | p6) +
      (p7 | p0);
      int N = min(N1,N2);
      if ((N == 2) || (N == 3)) {
                                                /// calculate criteria 3
        int c3 = ( p5 | p6 | (1-p0)) & p7;
        if(c3 == 0) {
          c[] = 0; 
          * modif = 1;
        }
      }
    }
  }
  boundary({c});
  restriction({c});
}




void thinning2D(struct thinningStruct p){
  // algo of Guo and Hall
  // c should have either a value of 0 or 1.
/**
Default max thickness is 30.
*/
  scalar c      = p.c;
  int thickness = p.thickness;

  if(thickness == 0) thickness = 30;

  bool modif = 1;
  bool parity = 0;
  int iteration = 0;
  do{
    modif = 0;
    if(parity)thinningIteration1(c, &modif);
    else thinningIteration2(c, &modif);
    parity = !parity;
    iteration++;
    if(iteration>thickness){ // thinning stops in zones where we have more than
     // a 100
     // cells.
      modif = 0;
    }
  }while(modif);
}

void filtered_thinning(scalar skeleton, scalar f, int layers, int thickness){
  foreach(){
    if(f[]> 1.e-3){
      skeleton[] = 1;
    }
    else{
     skeleton[] = 0;
   }
  }
  boundary({skeleton});
  restriction({skeleton});

  thinning2D(skeleton);
  /**
  Simple method to select only cells near the interface. Here the idea is to
  remove the skeleton in the bulk and only select the film. Once we have the
  skeleton we can remove cells which are further than `n` cells away from the
  interface using a `n` pass algorithm and tagging functions.
  */
  int n = layers;
  scalar tag[], tag2[];
  /**
  First we tag interfacial cells.
  */
  foreach(){
    if(f[]> 1.e-3 && f[] < 0.999){
      tag[] = 1;
    }
    else 
      tag[] = 0;
  }
  boundary({tag});

  /**
  We propagate the tag value.
  */
  for (int i = 0; i < n; i++){
    foreach(){
      tag2[] = tag[];
    }
    boundary({tag2});
    
    foreach(){
      if(f[] > 1.e-3 && tag[] == 0){
        double maxval = 0;
        foreach_neighbor(1)
          maxval = max(maxval,tag2[]);
        tag[] = max(maxval,tag[]); 
      }
    }
    boundary({tag});
  }

  /**
  After the propagation we remove the skeleton in the bulk.
  */
  foreach(){
    if(tag[] == 0) skeleton[] =0;
  }
  boundary({skeleton});
}


/**

##3D

Method from [Lee et al., 94](#lee1994building). I took the c++/java code of 
[Fiji](https://github.com/fiji/Skeletonize3D) and adapted it for basilisk.

The basic idea is to do 6 sub-iterations, (up down left right front back) and
select only cells that can be thinned. Such cells are categorized using Euler
invariancy of a 3x3 stencil and an associated lookup table (LUT) of the octree
that can be built from such a stencil.
One must also select cells that have at least to neighbours, so as to not "cut"
the skeleton that we want to build.

This algorithm has been tested only in sequential test cases, but it should work
in parallel (it was written for it). I may have done some bad translating from
c++ and will check it the near future.

*/

static const int LUT[256] = {
  0,  1,  0, -1,  0, -1,  0,  1, //  0..7
  0, -3,  0, -1,  0, -1,  0,  1, //  8..15
  0, -1,  0,  1,  0,  1,  0, -1, //  16..23
  0,  3,  0,  1,  0,  1,  0, -1, //  24..31
  0, -3,  0, -1,  0,  3,  0,  1, //  32..39
  0,  1,  0, -1,  0,  3,  0,  1, //  40..47
  0, -1,  0,  1,  0,  1,  0, -1, //  48..55
  0,  3,  0,  1,  0,  1,  0, -1, //  56..63
  0, -3,  0,  3,  0, -1,  0,  1, //  64..71
  0,  1,  0,  3,  0, -1,  0,  1, //  72..79
  0, -1,  0,  1,  0,  1,  0, -1, //  80..87
  0,  3,  0,  1,  0,  1,  0, -1, //  88..95
  0,  1,  0,  3,  0,  3,  0,  1, //  96..103
  0,  5,  0,  3,  0,  3,  0,  1, // 104..111
  0, -1,  0,  1,  0,  1,  0, -1, // 112..119
  0,  3,  0,  1,  0,  1,  0, -1, // 120..127
  0, -7,  0, -1,  0, -1,  0,  1, // 128..135
  0, -3,  0, -1,  0, -1,  0,  1, // 136..143
  0, -1,  0,  1,  0,  1,  0, -1, // 144..151
  0,  3,  0,  1,  0,  1,  0, -1, // 152..159
  0, -3,  0, -1,  0,  3,  0,  1, // 160..167
  0,  1,  0, -1,  0,  3,  0,  1, // 168..175
  0, -1,  0,  1,  0,  1,  0, -1, // 176..183
  0,  3,  0,  1,  0,  1,  0, -1, // 184..191
  0, -3,  0,  3,  0, -1,  0,  1, // 192..199
  0,  1,  0,  3,  0, -1,  0,  1, // 200..207
  0, -1,  0,  1,  0,  1,  0, -1, // 208..215
  0,  3,  0,  1,  0,  1,  0, -1, // 216..223
  0,  1,  0,  3,  0,  3,  0,  1, // 224..231
  0,  5,  0,  3,  0,  3,  0,  1, // 232..239
  0, -1,  0,  1,  0,  1,  0, -1, // 240..247
  0,  3,  0,  1,  0,  1,  0, -1  // 248..255
};



static inline void getNeighborhood(Point point, scalar c, int neighborhood[27]){
  neighborhood[ 0] = c[-1, -1, -1];
  neighborhood[ 1] = c[0 , -1, -1];
  neighborhood[ 2] = c[1 , -1, -1];

  neighborhood[ 3] = c[-1, 0, -1];
  neighborhood[ 4] = c[ 0, 0, -1];
  neighborhood[ 5] = c[ 1, 0, -1];

  neighborhood[ 6] = c[-1, 1, -1];
  neighborhood[ 7] = c[ 0, 1, -1];
  neighborhood[ 8] = c[ 1, 1, -1];

  neighborhood[ 9] = c[-1, -1, 0];
  neighborhood[10] = c[ 0, -1, 0];
  neighborhood[11] = c[ 1, -1, 0];

  neighborhood[12] = c[-1, 0, 0];
  neighborhood[13] = c[ ];
  neighborhood[14] = c[ 1, 0, 0];

  neighborhood[15] = c[-1, 1, 0];
  neighborhood[16] = c[ 0, 1, 0];
  neighborhood[17] = c[ 1, 1, 0];

  neighborhood[18] = c[-1, -1, 1];
  neighborhood[19] = c[ 0, -1, 1];
  neighborhood[20] = c[ 1, -1, 1];

  neighborhood[21] = c[ -1, 0, 1];
  neighborhood[22] = c[  0, 0, 1];
  neighborhood[23] = c[  1, 0, 1];

  neighborhood[24] = c[-1, 1, 1];
  neighborhood[25] = c[ 0, 1, 1];
  neighborhood[26] = c[ 1, 1, 1];

} 
/* end getNeighborhood */
  
static bool isEndPoint(scalar c, Point point)
{
  int numberOfNeighbors = -1;   // -1 and not 0 because the center pixel will be counted as well
  int neighbor[27];
  getNeighborhood(point, c, neighbor);
  for( int i = 0; i < 27; i++ ) // i =  0..26
  {                   
    if( neighbor[i] == 1 )
      numberOfNeighbors++;
  }
  return  numberOfNeighbors == 1;        
}

static bool isEulerInvariant(int neighbor[27]){ 
  // calculate Euler characteristic for each octant and sum up 
  int EulerChar = 0; 
  unsigned char n; 
  // Octant SWU 
  n = 1; 
  if( neighbor[24]==1 ) 
    n |= 128; 
  if( neighbor[25]==1 ) 
    n |=  64; 
  if( neighbor[15]==1 ) 
    n |=  32; 
  if( neighbor[16]==1 ) 
    n |=  16; 
  if( neighbor[21]==1 ) 
    n |=   8; 
  if( neighbor[22]==1 ) 
    n |=   4; 
  if( neighbor[12]==1 ) 
    n |=   2; 
  EulerChar += LUT[n]; 
  // Octant SEU 
  n = 1; 
  if( neighbor[26]==1 ) 
    n |= 128; 
  if( neighbor[23]==1 ) 
    n |=  64; 
  if( neighbor[17]==1 ) 
    n |=  32; 
  if( neighbor[14]==1 ) 
    n |=  16; 
  if( neighbor[25]==1 ) 
    n |=   8; 
  if( neighbor[22]==1 ) 
    n |=   4; 
  if( neighbor[16]==1 ) 
    n |=   2; 
  EulerChar += LUT[n]; 
  // Octant NWU 
  n = 1; 
  if( neighbor[18]==1 ) 
    n |= 128; 
  if( neighbor[21]==1 ) 
    n |=  64; 
  if( neighbor[9]==1 ) 
    n |=  32; 
  if( neighbor[12]==1 ) 
    n |=  16; 
  if( neighbor[19]==1 ) 
    n |=   8; 
  if( neighbor[22]==1 ) 
    n |=   4; 
  if( neighbor[10]==1 ) 
    n |=   2; 
  EulerChar += LUT[n]; 
  // Octant NEU 
  n = 1; 
  if( neighbor[20]==1 ) 
    n |= 128; 
  if( neighbor[23]==1 ) 
    n |=  64; 
  if( neighbor[19]==1 ) 
    n |=  32; 
  if( neighbor[22]==1 ) 
    n |=  16; 
  if( neighbor[11]==1 ) 
    n |=   8; 
  if( neighbor[14]==1 ) 
    n |=   4; 
  if( neighbor[10]==1 ) 
    n |=   2; 
  EulerChar += LUT[n]; 
  // Octant SWB 
  n = 1; 
  if( neighbor[6]==1 ) 
    n |= 128; 
  if( neighbor[15]==1 ) 
    n |=  64; 
  if( neighbor[7]==1 ) 
    n |=  32; 
  if( neighbor[16]==1 ) 
    n |=  16; 
  if( neighbor[3]==1 ) 
    n |=   8; 
  if( neighbor[12]==1 ) 
    n |=   4; 
  if( neighbor[4]==1 ) 
    n |=   2; 
  EulerChar += LUT[n]; 
  // Octant SEB 
  n = 1; 
  if( neighbor[8]==1 ) 
    n |= 128; 
  if( neighbor[7]==1 ) 
    n |=  64; 
  if( neighbor[17]==1 ) 
    n |=  32; 
  if( neighbor[16]==1 ) 
    n |=  16; 
  if( neighbor[5]==1 ) 
    n |=   8; 
  if( neighbor[4]==1 ) 
    n |=   4; 
  if( neighbor[14]==1 ) 
    n |=   2; 
  EulerChar += LUT[n]; 
  // Octant NWB 
  n = 1; 
  if( neighbor[0]==1 ) 
    n |= 128; 
  if( neighbor[9]==1 ) 
    n |=  64; 
  if( neighbor[3]==1 ) 
    n |=  32; 
  if( neighbor[12]==1 ) 
    n |=  16; 
  if( neighbor[1]==1 ) 
    n |=   8; 
  if( neighbor[10]==1 ) 
    n |=   4; 
  if( neighbor[4]==1 ) 
    n |=   2; 
  EulerChar += LUT[n]; 
  // Octant NEB 
  n = 1; 
  if( neighbor[2]==1 ) 
    n |= 128; 
  if( neighbor[1]==1 ) 
    n |=  64; 
  if( neighbor[11]==1 ) 
    n |=  32; 
  if( neighbor[10]==1 ) 
    n |=  16; 
  if( neighbor[5]==1 ) 
    n |=   8; 
  if( neighbor[4]==1 ) 
    n |=   4; 
  if( neighbor[14]==1 ) 
    n |=   2; 
  EulerChar += LUT[n]; 
  if( EulerChar == 0 ) 
    return true; 
  else 
    return false; 
} 

static void Octree_labeling(int octant, int label, int *cube) 
{ 
  // check if there are points in the octant with value 1 
  if( octant==1 ) 
  { 
    // set points in this octant to current label 
    // and recurseive labeling of adjacent octants 
    if( cube[0] == 1 ) 
      cube[0] = label; 
    if( cube[1] == 1 ) 
    { 
      cube[1] = label;         
      Octree_labeling( 2, label, cube); 
    } 
    if( cube[3] == 1 ) 
    { 
      cube[3] = label;         
      Octree_labeling( 3, label, cube); 
    } 
    if( cube[4] == 1 ) 
    { 
      cube[4] = label;         
      Octree_labeling( 2, label, cube); 
      Octree_labeling( 3, label, cube); 
      Octree_labeling( 4, label, cube); 
    } 
    if( cube[9] == 1 ) 
    { 
      cube[9] = label;         
      Octree_labeling( 5, label, cube); 
    } 
    if( cube[10] == 1 ) 
    { 
      cube[10] = label;         
      Octree_labeling( 2, label, cube); 
      Octree_labeling( 5, label, cube); 
      Octree_labeling( 6, label, cube); 
    } 
    if( cube[12] == 1 ) 
    { 
      cube[12] = label;         
      Octree_labeling( 3, label, cube); 
      Octree_labeling( 5, label, cube); 
      Octree_labeling( 7, label, cube); 
    } 
  } 
  if( octant==2 ) 
  { 
    if( cube[1] == 1 ) 
    { 
      cube[1] = label; 
      Octree_labeling( 1, label, cube); 
    } 
    if( cube[4] == 1 ) 
    { 
      cube[4] = label;         
      Octree_labeling( 1, label, cube); 
      Octree_labeling( 3, label, cube); 
      Octree_labeling( 4, label, cube); 
    } 
    if( cube[10] == 1 ) 
    { 
      cube[10] = label;         
      Octree_labeling( 1, label, cube); 
      Octree_labeling( 5, label, cube); 
      Octree_labeling( 6, label, cube); 
    } 
    if( cube[2] == 1 ) 
      cube[2] = label;         
    if( cube[5] == 1 ) 
    { 
      cube[5] = label;         
      Octree_labeling( 4, label, cube); 
    } 
    if( cube[11] == 1 ) 
    { 
      cube[11] = label;         
      Octree_labeling( 6, label, cube); 
    } 
    if( cube[13] == 1 ) 
    { 
      cube[13] = label;         
      Octree_labeling( 4, label, cube); 
      Octree_labeling( 6, label, cube); 
      Octree_labeling( 8, label, cube); 
    } 
  } 
  if( octant==3 ) 
  { 
    if( cube[3] == 1 ) 
    { 
      cube[3] = label;         
      Octree_labeling( 1, label, cube); 
    } 
    if( cube[4] == 1 ) 
    { 
      cube[4] = label;         
      Octree_labeling( 1, label, cube); 
      Octree_labeling( 2, label, cube); 
      Octree_labeling( 4, label, cube); 
    } 
    if( cube[12] == 1 ) 
    { 
      cube[12] = label;         
      Octree_labeling( 1, label, cube); 
      Octree_labeling( 5, label, cube); 
      Octree_labeling( 7, label, cube); 
    } 
    if( cube[6] == 1 ) 
      cube[6] = label;         
    if( cube[7] == 1 ) 
    { 
      cube[7] = label;         
      Octree_labeling( 4, label, cube); 
    } 
    if( cube[14] == 1 ) 
    { 
      cube[14] = label;         
      Octree_labeling( 7, label, cube); 
    } 
    if( cube[15] == 1 ) 
    { 
      cube[15] = label;         
      Octree_labeling( 4, label, cube); 
      Octree_labeling( 7, label, cube); 
      Octree_labeling( 8, label, cube); 
    } 
  } 
  if( octant==4 ) 
  { 
    if( cube[4] == 1 ) 
    { 
      cube[4] = label;         
      Octree_labeling( 1, label, cube); 
      Octree_labeling( 2, label, cube); 
      Octree_labeling( 3, label, cube); 
    } 
    if( cube[5] == 1 ) 
    { 
      cube[5] = label;         
      Octree_labeling( 2, label, cube); 
    } 
    if( cube[13] == 1 ) 
    { 
      cube[13] = label;         
      Octree_labeling( 2, label, cube); 
      Octree_labeling( 6, label, cube); 
      Octree_labeling( 8, label, cube); 
    } 
    if( cube[7] == 1 ) 
    { 
      cube[7] = label;         
      Octree_labeling( 3, label, cube); 
    } 
    if( cube[15] == 1 ) 
    { 
      cube[15] = label;         
      Octree_labeling( 3, label, cube); 
      Octree_labeling( 7, label, cube); 
      Octree_labeling( 8, label, cube); 
    } 
    if( cube[8] == 1 ) 
      cube[8] = label;         
    if( cube[16] == 1 ) 
    { 
      cube[16] = label;         
      Octree_labeling( 8, label, cube); 
    } 
  } 
  if( octant==5 ) 
  { 
    if( cube[9] == 1 ) 
    { 
      cube[9] = label;         
      Octree_labeling( 1, label, cube); 
    } 
    if( cube[10] == 1 ) 
    { 
      cube[10] = label;         
      Octree_labeling( 1, label, cube); 
      Octree_labeling( 2, label, cube); 
      Octree_labeling( 6, label, cube); 
    } 
    if( cube[12] == 1 ) 
    { 
      cube[12] = label;         
      Octree_labeling( 1, label, cube); 
      Octree_labeling( 3, label, cube); 
      Octree_labeling( 7, label, cube); 
    } 
    if( cube[17] == 1 ) 
      cube[17] = label;         
    if( cube[18] == 1 ) 
    { 
      cube[18] = label;         
      Octree_labeling( 6, label, cube); 
    } 
    if( cube[20] == 1 ) 
    { 
      cube[20] = label;         
      Octree_labeling( 7, label, cube); 
    } 
    if( cube[21] == 1 ) 
    { 
      cube[21] = label;         
      Octree_labeling( 6, label, cube); 
      Octree_labeling( 7, label, cube); 
      Octree_labeling( 8, label, cube); 
    } 
  } 
  if( octant==6 ) 
  { 
    if( cube[10] == 1 ) 
    { 
      cube[10] = label;         
      Octree_labeling( 1, label, cube); 
      Octree_labeling( 2, label, cube); 
      Octree_labeling( 5, label, cube); 
    } 
    if( cube[11] == 1 ) 
    { 
      cube[11] = label;         
      Octree_labeling( 2, label, cube); 
    } 
    if( cube[13] == 1 ) 
    { 
      cube[13] = label;         
      Octree_labeling( 2, label, cube); 
      Octree_labeling( 4, label, cube); 
      Octree_labeling( 8, label, cube); 
    } 
    if( cube[18] == 1 ) 
    { 
      cube[18] = label;         
      Octree_labeling( 5, label, cube); 
    } 
    if( cube[21] == 1 ) 
    { 
      cube[21] = label;         
      Octree_labeling( 5, label, cube); 
      Octree_labeling( 7, label, cube); 
      Octree_labeling( 8, label, cube); 
    } 
    if( cube[19] == 1 ) 
      cube[19] = label;         
    if( cube[22] == 1 ) 
    { 
      cube[22] = label;         
      Octree_labeling( 8, label, cube); 
    } 
  } 
  if( octant==7 ) 
  { 
    if( cube[12] == 1 ) 
    { 
      cube[12] = label;         
      Octree_labeling( 1, label, cube); 
      Octree_labeling( 3, label, cube); 
      Octree_labeling( 5, label, cube); 
    } 
    if( cube[14] == 1 ) 
    { 
      cube[14] = label;         
      Octree_labeling( 3, label, cube); 
    } 
    if( cube[15] == 1 ) 
    { 
      cube[15] = label;         
      Octree_labeling( 3, label, cube); 
      Octree_labeling( 4, label, cube); 
      Octree_labeling( 8, label, cube); 
    } 
    if( cube[20] == 1 ) 
    { 
      cube[20] = label;         
      Octree_labeling( 5, label, cube); 
    } 
    if( cube[21] == 1 ) 
    { 
      cube[21] = label;         
      Octree_labeling( 5, label, cube); 
      Octree_labeling( 6, label, cube); 
      Octree_labeling( 8, label, cube); 
    } 
    if( cube[23] == 1 ) 
      cube[23] = label;         
    if( cube[24] == 1 ) 
    { 
      cube[24] = label;         
      Octree_labeling( 8, label, cube); 
    } 
  } 
  if( octant==8 ) 
  { 
    if( cube[13] == 1 ) 
    { 
      cube[13] = label;         
      Octree_labeling( 2, label, cube); 
      Octree_labeling( 4, label, cube); 
      Octree_labeling( 6, label, cube); 
    } 
    if( cube[15] == 1 ) 
    { 
      cube[15] = label;         
      Octree_labeling( 3, label, cube); 
      Octree_labeling( 4, label, cube); 
      Octree_labeling( 7, label, cube); 
    } 
    if( cube[16] == 1 ) 
    { 
      cube[16] = label;         
      Octree_labeling( 4, label, cube); 
    } 
    if( cube[21] == 1 ) 
    { 
      cube[21] = label;         
      Octree_labeling( 5, label, cube); 
      Octree_labeling( 6, label, cube); 
      Octree_labeling( 7, label, cube); 
    } 
    if( cube[22] == 1 ) 
    { 
      cube[22] = label;         
      Octree_labeling( 6, label, cube); 
    } 
    if( cube[24] == 1 ) 
    { 
      cube[24] = label;         
      Octree_labeling( 7, label, cube); 
    } 
    if( cube[25] == 1 ) 
      cube[25] = label;         
  }  
}

static bool isSimplePoint(int neighbor[26]) 
{ 
  // copy neighbors for labeling 
  int cube[26]; 
  int i; 
  for( i = 0; i < 13; i++ )  // i =  0..12 -> cube[0..12] 
    cube[i] = neighbor[i]; 
  // i != 13 : ignore center pixel when counting (see [Lee94]) 
  for( i = 14; i < 27; i++ ) // i = 14..26 -> cube[13..25] 
    cube[i-1] = neighbor[i]; 
  // set initial label 
  int label = 2; 
  // for all points in the neighborhood 
  for( int i = 0; i < 26; i++ ){ 
    if( cube[i]==1 ){     // voxel has not been labelled yet  
      // start recursion with any octant that contains the point i 
      switch( i ) 
      { 
      case 0: 
      case 1: 
      case 3: 
      case 4: 
      case 9: 
      case 10: 
      case 12: 
        Octree_labeling(1, label, cube ); 
        break; 
      case 2: 
      case 5: 
      case 11: 
      case 13: 
        Octree_labeling(2, label, cube ); 
        break; 
      case 6: 
      case 7: 
      case 14: 
      case 15: 
        Octree_labeling(3, label, cube ); 
        break; 
      case 8: 
      case 16: 
        Octree_labeling(4, label, cube ); 
        break; 
      case 17: 
      case 18: 
      case 20: 
      case 21: 
        Octree_labeling(5, label, cube ); 
        break; 
      case 19: 
      case 22: 
        Octree_labeling(6, label, cube ); 
        break; 
      case 23: 
      case 24: 
        Octree_labeling(7, label, cube ); 
        break; 
      case 25: 
        Octree_labeling(8, label, cube ); 
        break; 
      } 
      label++; 
      if( label-2 >= 2 ) 
      { 
        return false; 
      } 
    } 
  } 
  //return label-2; in [Lee94] if the number of connected compontents would be needed 
  return true; 
}






void thinning3D(struct thinningStruct p){

  scalar c      = p.c;
  int thickness = p.thickness;

  if(thickness==0)thickness = 10;

  scalar c3[];
  int unchangedBorders = 0;
  int loopcount;
  for (loopcount = 0; loopcount < thickness; ++loopcount){
    unchangedBorders = 0;
    for( int currentBorder = 1; currentBorder <= 6; currentBorder++) 
    {
      scalar tag[];


      scalar cp[];
      c.restriction = cp.restriction = myrestrict;
      c.prolongation = cp.prolongation = myprolongation;

      restriction({c});
      foreach(){
        cp[] = c[];
        tag[] = 0.;
      }
      boundary({cp});
      restriction({cp});


      foreach(){

        if(cp[] != 1)continue; // action only if pixel is foreground
  
        bool isBorderPoint = false; 
        if( currentBorder == 1 && cp[0,-1,0]<=0 ) 
          isBorderPoint = true; 
        if( currentBorder == 2 && cp[0,1,0]<=0 ) 
          isBorderPoint = true; 
        if( currentBorder == 3 && cp[1,0,0]<=0 ) 
          isBorderPoint = true; 
        if( currentBorder == 4 && cp[-1,0,0]<=0 ) 
          isBorderPoint = true; 
        if( currentBorder == 5 && cp[0,0,1]<=0 ) 
          isBorderPoint = true; 
        if( currentBorder == 6 && cp[0,0,-1]<=0 ) 
          isBorderPoint = true; 
        if( !isBorderPoint ){ 
          continue;         // current point is not deletable 
        }         
        if( !isBorderPoint ){ 
          continue;         // current point is not deletable 
        }         
        // check if point is the end of an arc 
        if( isEndPoint( cp, point) ){
          continue;
        }

        int neighbor[27];
        getNeighborhood(point, cp, neighbor);
        // Check if point is simple (deletion does not change connectivity in the 3x3x3 neighborhood)
        // (conditions 2 and 3 in Lee[94])
        if( !isEulerInvariant( neighbor ) ){
          continue;         // current point is not deletable
        }

        if( !isSimplePoint( neighbor)){ 
          continue;         // current point is not deletable 
        }
        // add all simple border points to a cache for sequential re-checking 
        tag[] = 1.;
      }

      bool noChange = true;
      boundary({tag});
/**

*/
      foreach(){
        if(tag[]){
          c[] = 0.;
          noChange = false;
        }
      }  
      
      boundary({c});
      restriction({c});
      if(noChange)
        unchangedBorders++;
    } // end currentBorder for loop
  }
#ifdef thinDebug
  fprintf(stderr, "##Thinning3D, number of thinning operations %d\n", loopcount);
#endif
}

void filtered3Dthinning(struct thinningStruct p){
  scalar c      = p.c;
  int thickness = p.thickness;
  if(p.thickness<2)thickness =4;
  thinning3D(c,thickness-2);
  scalar copy[];
  foreach()
    copy[] = c[];
  boundary({copy});
  thinning3D(c,2);

  foreach()
    c[] = c[]*copy[];
  boundary({c});
}

/**
~~~bib 
#article{guo1989parallel,
  title={Parallel thinning with two-subiteration algorithms},
  author={Guo, Zicheng and Hall, Richard W},
  journal={Communications of the ACM},
  volume={32},
  number={3},
  pages={359--373},
  year={1989},
  publisher={ACM New York, NY, USA}
}

#article{lee1994building,
  title={Building skeleton models via 3-D medial surface axis thinning algorithms},
  author={Lee, Ta-Chih and Kashyap, Rangasami L and Chu, Chong-Nam},
  journal={CVGIP: Graphical Models and Image Processing},
  volume={56},
  number={6},
  pages={462--478},
  year={1994},
  publisher={Elsevier}
}


~~~

##Examples

[2D atomisation](1_test_cases/atomisation_skeleton.c)

[Scardovelli's reversed flow](1_test_cases/reversed_skeleton.c)

[Skeleton with a large range of characteristic sizes](1_test_cases/distance.c)

[Rayleigh-Plateau](1_test_cases/simili-RP-instab.c)


##To do list

- add a pruning method. [For instance erosion method](https://github.com/danielyan86129/ET) or [Here too](#https://github.com/danielyan86129/voxel_ma.git)



*/