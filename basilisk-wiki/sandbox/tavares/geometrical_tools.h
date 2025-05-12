/**
# Miscalleneous functions

## 2D Intersection routine -- this function tells if two segments intersect with their coordinates

*/

int IntersectB(coord A1, coord A2, coord B1, coord B2){
	coord t1, t2, t3; 
	foreach_dimension(){
		t1.x = B2.x-B1.x;
		t2.x = B1.x-A1.x;
		t3.x = B1.x-A2.x;
	}
  double det = ( t2.x*t1.y-t2.y*t1.x ) * ( t3.x*t1.y-t3.y*t1.x );
  int intersection = (det<=0.); //B1B2 cut the segment A1A2

  foreach_dimension(){
    t1.x=A2.x-A1.x;
    t2.x=A1.x-B1.x;
    t3.x=A1.x-B2.x;
  }	
  det=( t2.x*t1.y-t2.y*t1.x ) * ( t3.x*t1.y-t3.y*t1.x );
  return intersection = (intersection && (det<=0.) );// A1A2 cut the segment B1B2
}


static inline int Intersect(coord A, coord B, coord C, coord D, coord * I){
  coord t1, t2, t3; 
  foreach_dimension(){
    t1.x = D.x-C.x;
    t2.x = C.x-A.x;
    t3.x = C.x-B.x;
  }
  double det = ( t2.x*t1.y-t2.y*t1.x ) * ( t3.x*t1.y-t3.y*t1.x );
  int intersection = (det<=0.); //B1B2 cut the segment A1A2

  foreach_dimension(){
    t1.x=B.x-A.x;
    t2.x=A.x-C.x;
    t3.x=A.x-D.x;
  } 
  det=( t2.x*t1.y-t2.y*t1.x ) * ( t3.x*t1.y-t3.y*t1.x );
  intersection = (intersection && (det<=0.) );// A1A2 cut the segment B1B2
  if (intersection){ // We calculate the intersection between the two segments
  //fprintf(stdout, "in %g %g\n", I.x, I.y );
    foreach_dimension() 
      t2.x = D.x-C.x;
    det = t2.y*t1.x  - t2.x*t1.y; //We use the Cramer method
    double detAB = B.x*A.y - A.x*B.y;
    double detCD = D.x*C.y - C.x*D.y;
    //fprintf(stdout, "in %g %g %g %g %g\n", I.x, I.y, det, detAB, detCD );
    foreach_dimension() 
      I->x = detAB * t2.x - detCD * t1.x / det;
    //fprintf(stdout, "in %g %g\n", I->x, I->y );
    }
  return intersection;
}


/**
## 3D Intersection routine -- this function tells if a ray intersects a bounded box (see[http://www.basilisk.fr/sandbox/Antoonvh/bwatch.h#intersections]).
*/

// A ray class
typedef struct ray {
  coord O;
  coord dir;
  int   depth;   // for recasted rays
} ray;

// static inline int ray_box (ray r, coord bb[2]) { //bb is the extent of the box
//   double in, out;
//   int intersection;
//   in = -HUGE; out = HUGE;
//   foreach_dimension() {
//     if (r.dir.x > 0) {
//       in  = max(in,  (bb[0].x - r.O.x)/r.dir.x);
//       out = min(out, (bb[1].x - r.O.x)/r.dir.x);
//     } else if (r.dir.x < 0) {
//       in  = max(in,  (bb[1].x - r.O.x)/r.dir.x);
//       out = min(out, (bb[0].x - r.O.x)/r.dir.x);
//     }
//   }
//   if (in >= out || out < 0) intersection = 0;
//     //return HUGE;
//   //in = in < 0 ? 0 : in; //The origin is in the box.
//   //return in;
//   intersection = (in>0);
//   return intersection;
// }


//static inline double ray_box (ray r, coord bb[2]) { //bb is the extent of the box
static inline double ray_box (ray r, coord bb[2]) { //bb is the extent of the box
  
  double in, out;
  in = -HUGE; out = HUGE;
  foreach_dimension() {
    if (r.dir.x > 0) {
      in  = max(in,  (bb[0].x - r.O.x)/r.dir.x);
      out = min(out, (bb[1].x - r.O.x)/r.dir.x);
    } else if (r.dir.x < 0) {
      in  = max(in,  (bb[1].x - r.O.x)/r.dir.x);
      out = min(out, (bb[0].x - r.O.x)/r.dir.x);
    }
  }
  if (in >= out || out < 0);
    return HUGE;
  in = in < 0 ? 0 : in; //The origin is in the box.
  //fprintf(fout, "%g\n",in);
  return in;
}


/** This function returns *true* if there is an intersection between a segment and a facet  */
//static inline bool facet_embed_intersect (scalar c, face vector s, scalar f, coord I[2], coord n, coord m, Point point) { 
//static inline bool facet_embed_intersect (double c, double f, coord n, coord  m, coord I0, Point point) { 

//bool facet_embed_intersect (double c, double f, coord n, coord  m, coord I0, Point point) {   
  int face_embed_intersect (double c, double f, coord n, coord  m, coord It, Point point) {
  coord cc = {x, y, z}; 
  coord a[2];
  double alphan = plane_alpha (f, n);
  double alpham = plane_alpha (c, m);
  //fprintf(stdout, "in %g %g\n", It.x, It.y);
  if (facets (m, alpham, a) >=2) {
    double ALP = 0, ALP2 = 0;
    foreach_dimension() {
      ALP  += n.x*(a[0].x - cc.x)/Delta;
      ALP2 += n.x*(a[1].x - cc.x)/Delta;
    }
    fprintf(stdout, "in %g \n", (ALP2 - alphan)/(ALP - alphan));
    if ((ALP2 - alphan)/(ALP - alphan) > 0.05) // 3% gap filling
      return 0;
    double w = fabs((ALP2 - alphan)) / (fabs(ALP - alphan) + fabs(ALP2 - alphan));
    foreach_dimension() {
      It.x = w*a[0].x + (1 - w)*a[1].x;
    }
    fprintf(stdout, "out %g %g\n", It.x, It.y);
    return 1;
  }
  else 
    return 0;
}
