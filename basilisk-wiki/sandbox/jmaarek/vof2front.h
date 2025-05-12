#include "fractions.h"
#include "distance.h"

#define FRONTTHR 1e-8


//compute normal with height function if possible
coord interface_normal2(Point point, scalar c)
{
  coord n;
  if (!c.height.x.i || (n = height_normal (point, c, c.height)).x == nodata)
    n = mycs (point, c);
  else {
    double nn = 0.;
    foreach_dimension()
      nn += fabs(n.x);
    foreach_dimension()
      n.x /= nn;
  }
  return n;
}

void vof2dist(scalar c, vertex scalar phi){
  //if phi.height exists compute height function
  if (c.height.x.i)
    heights (c, c.height);
  //set magnitude of signed distance to be greater than maximum domain size
  foreach_vertex()
    phi[] = L0*sqrt(2)+1;

  boundary({phi});
  restriction({phi});

  foreach()
    if (fabs(c[] - 0.5) <= (0.5-FRONTTHR)){
      coord n = interface_normal2(point, c);
      //compute alpha and facets
      double alpha = plane_alpha (c[], n);
#if dimension == 2
      coord segment[2];
      coord temp1;
      coord temp2;
      if (facets (n, alpha, segment) == 2){
        //define segments in absolute location
        segment[0] = (coord){x + segment[0].x*Delta, y + segment[0].y*Delta, 0.0};
        segment[1] = (coord){x + segment[1].x*Delta, y + segment[1].y*Delta, 0.0};
        //define segment such that normal matches
        if ((sign(segment[1].y -segment[0].y) == sign(n.x)) && (sign(segment[0].x -segment[1].x) == sign(n.y))){
           temp1 = segment[1];
           temp2 = segment[0];
        }
        else{
          temp1 = segment[0];
          temp2 = segment[1];
        }
        //compute signed minimum distance on local stencil
        foreach_neighbor(2){
          coord r;
          coord c = (coord){x-0.5*Delta,y-0.5*Delta,0};
          double s, d2 = PointSegmentDistance (&c, &(temp1), &(temp2), &r, &s);
          if (sqrt(d2) < fabs(phi[])){
            phi[] = sqrt(d2)*((double)PointSegmentOrientation(&c, &(temp1), &(temp2)));}
          //if local cell is refine compute signed minimum distance for children
          if (is_refined (cell)){
            for (int ii = 0; ii<=2; ii++)
              for (int jj = 0; jj<=2; jj++){
                coord r;
                coord c = (coord){x+(-0.5+ii*0.5)*Delta,y+(-0.5+jj*0.5)*Delta,0};
                double s, d2 = PointSegmentDistance (&c, &(temp1), &(temp2), &r, &s);
                if (sqrt(d2) < fabs(fine(phi,ii,jj))){
                  fine(phi,ii,jj) = sqrt(d2)*((double)PointSegmentOrientation(&c, &(temp1), &(temp2)));}
              }

          }
          if (is_prolongation(cell)){
            //determine whether the vertex exists on the coarser grid in the case of prolongation cells
            if ((point.i + GHOSTS)%2 + (point.j + GHOSTS)%2 == 0){
              if (fabs(phi[]) < fabs(coarse(phi,0,0)))
                coarse(phi,0,0) = phi[];
              }
          }

        }
    }
  #else //dimension = 3
      coord v[12];
      int m = facets (n, alpha, v, 1.); //m is the number of coordinates
      //split facet into a set of triangles and follow same method
      for (int j = 0; j < m - 2; j++) {
        coord temp1 = {x + v[0].x*Delta  , y + v[0].y*Delta  , z + v[0].z*Delta};
        coord temp2 = {x + v[j+1].x*Delta, y + v[j+1].y*Delta, z + v[j+1].z*Delta};
        coord temp3 = {x + v[j+2].x*Delta, y + v[j+2].y*Delta, z + v[j+2].z*Delta};

        foreach_neighbor(2){
          coord r;
          coord c = (coord){x-0.5*Delta,y-0.5*Delta,z-0.5*Delta};
          double s, t, d2 = PointTriangleDistance (&c, &(temp1), &(temp2), &(temp3), &s, &t);

          if (sqrt(d2) < fabs(phi[])){
            phi[] = sqrt(d2)*((double)PointTriangleOrientation(&c, &(temp1), &(temp2), &(temp3)));}

          //if local cell is refine compute signed minimum distance for children
          if (is_refined (cell))
            for (int ii = 0; ii<=2; ii++)
              for (int jj = 0; jj<=2; jj++)
                for (int kk = 0; kk<=2; kk++){
                  coord r;
                  coord c = (coord){x+(-0.5+ii*0.5)*Delta,y+(-0.5+jj*0.5)*Delta,z+(-0.5+kk*0.5)*Delta};
                  double s, t, d2 = PointTriangleDistance (&c, &(temp1), &(temp2), &(temp3), &s, &t);
                  if (sqrt(d2) < fabs(fine(phi,ii,jj,kk))){
                    phi[] = sqrt(d2)*((double)PointTriangleOrientation(&c, &(temp1), &(temp2), &(temp3)));}
              }
          if (is_prolongation(cell)){
            //determine whether the vertex exists on the coarser grid
            if ((point.i + GHOSTS)%2 + (point.j + GHOSTS)%2 + (point.k + GHOSTS)%2 == 0){
              if (fabs(phi[]) < fabs(coarse(phi,0,0)))
                coarse(phi,0,0) = phi[];
              }
          }


        }
      }
  #endif


  }

  boundary({phi});
  restriction({phi});

  /*heuristic check to confirm sign of non interfacial cells
  foreach vertex check vof values of adjacent cells at most refined level.
  If any are non interfacial set the sign of the vertex such that it is consistent
  */

  #if dimension == 2
  foreach_vertex()
    for (int ii = -1; ii<=0; ii++)
      for (int jj = -1; jj<=0; jj++){
        if is_refined(neighbor(ii,jj)){
          if (fabs(fine(c,ii,jj) - 0.5) > (0.5-FRONTTHR)) //if !interfacial
            phi[] = sign(fine(c,ii,jj) - 0.5)*fabs(phi[]);}
        else{
          if (fabs(c[ii,jj] - 0.5) > (0.5-FRONTTHR)) //if !interfacial
            phi[] = sign(c[ii,jj] - 0.5)*fabs(phi[]);}
      }

  #else //dimension == 3
  foreach_vertex()
    for (int ii = -1; ii<=0; ii++)
      for (int jj = -1; jj<=0; jj++)
        for (int kk = -1; kk<=0; kk++){
          if is_refined(neighbor(ii,jj,kk)){
            if (fabs(fine(c,ii,jj,kk) - 0.5) > (0.5-FRONTTHR)) //if !interfacial
              phi[] = sign(fine(c,ii,jj,kk) - 0.5)*fabs(phi[]);}
          else{
            if (fabs(c[ii,jj,kk] - 0.5) > (0.5-FRONTTHR)) //if !interfacial
              phi[] = sign(c[ii,jj,kk] - 0.5)*fabs(phi[]);}
        }

  #endif
  boundary({phi});
  restriction({phi});
}
