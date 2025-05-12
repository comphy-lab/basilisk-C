/**
# Functions to extend or restrict a field around an interface

This file is poorly documented. Do not hesitate to ask me questions:
quentinmagdelaine@gmail.com */

#ifndef F_ERR
  #define F_ERR 1e-10
#endif

void compute_mycs_nodata (scalar f, vector normal_vector) {
  #if TREE
    foreach_dimension() {
      normal_vector.x.refine = normal_vector.x.prolongation = curvature_prolongation;
      normal_vector.x.restriction = curvature_restriction;
    }
  #endif
  foreach() {
    bool close_to_interface = false;
    foreach_neighbor(1)
      close_to_interface = (interfacial(point, f) ? true : close_to_interface);
    // coord m = (close_to_interface ? mycs (point, f) : {0., 0.}); // does not work
    coord m = mycs (point, f);
    foreach_dimension() 
      normal_vector.x[] = (close_to_interface ? m.x : nodata);
  }
  boundary((scalar*){normal_vector});
}

void compute_normal_nodata (scalar f, vector normal_vector) {
  #if TREE
    foreach_dimension() {
      normal_vector.x.prolongation = normal_vector.x.refine = curvature_prolongation;
      normal_vector.x.restriction = curvature_restriction;    
    }
  #endif
  foreach() {
    bool close_to_interface = false;
    foreach_neighbor(1)
      close_to_interface = (interfacial(point, f) ? true : close_to_interface);
    // coord m = (close_to_interface ? mycs (point, f) : {0., 0.}); // does not work
    coord m = normal (point, f);
    foreach_dimension() 
      normal_vector.x[] = (close_to_interface ? m.x : nodata);
  }
  boundary((scalar*){normal_vector});
}

/**
You should not do mycs() in a ghost cell and you can't do interfacial() in the 
second round of neighbors because it need neighbor cells. The solution
here is to store the value you need in scalars or vector and to apply
boundary to extend it the ghost cell. An other solution would be to test first
if a cell is not a ghost before calling mycs and interfacial. */

void compute_extended_mycs (scalar f, vector normal_vector) {
  scalar interfacial_cell[];
  foreach() {
    interfacial_cell[] = (interfacial(point, f) ? 1 : 0);
    coord m = mycs (point, f);
    foreach_dimension() 
      normal_vector.x[] = (interfacial_cell[] > 0 ? m.x : nodata);
  }
  #if TREE
    foreach_dimension() {
      normal_vector.x.prolongation = normal_vector.x.refine = curvature_prolongation;
      normal_vector.x.restriction = curvature_restriction;    
    }
  #endif
  boundary({interfacial_cell, normal_vector});
  foreach() {
    coord m = {0., 0.};
    if (interfacial_cell[] == 0) {
      foreach_neighbor() {
        if (interfacial_cell[] > 0) {
          foreach_dimension()
            m.x += normal_vector.x[];
        }
      }
      double nm = 0.;
      foreach_dimension()
        nm += fabs(m.x);
      foreach_dimension()
        normal_vector.x[] = (nm > 0.5 ? m.x/nm : nodata);
    }
  }
  boundary((scalar*){normal_vector});
}

/** normal norm 2 */

void compute_extended_normal (scalar f, vector normal_vector) {
  scalar interfacial_cell[];
  foreach() {
    interfacial_cell[] = (interfacial(point, f) ? 1 : 0);
    coord m = normal (point, f);
    foreach_dimension() 
      normal_vector.x[] = (interfacial_cell[] > 0 ? m.x : nodata);
  }
  #if TREE
    foreach_dimension() {
      normal_vector.x.prolongation = normal_vector.x.refine = curvature_prolongation;
      normal_vector.x.restriction = curvature_restriction;    
    }
  #endif
  boundary({interfacial_cell, normal_vector});
  foreach() {
    coord m = {0., 0.};
    if (interfacial_cell[] == 0) {
      foreach_neighbor() {
        if (interfacial_cell[] > 0) {
          foreach_dimension()
            m.x += normal_vector.x[];
        }
      }
      double nm = 0.;
      foreach_dimension()
        nm += sq(m.x);
      foreach_dimension()
        normal_vector.x[] = (nm > 0.5 ? m.x/sqrt(nm) : nodata);
    }
  }
  boundary((scalar*){normal_vector});
}

void compute_restricted_normal (scalar f, vector normal_vector) {
  foreach() {
    coord m = normal (point, f);
    foreach_dimension() 
      normal_vector.x[] = (interfacial(point, f) ? m.x : nodata);
  }
  #if TREE
    foreach_dimension() {
      normal_vector.x.prolongation = normal_vector.x.refine = curvature_prolongation;
      normal_vector.x.restriction = curvature_restriction;    
    }
  #endif
  boundary((scalar *){normal_vector});
}

void compute_restricted_mycs (scalar f, vector normal_vector) {
  foreach() {
    coord m = mycs (point, f);
    foreach_dimension() 
      normal_vector.x[] = (interfacial(point, f) ? m.x : nodata);
  }
  #if TREE
    foreach_dimension() {
      normal_vector.x.prolongation = normal_vector.x.refine = curvature_prolongation;
      normal_vector.x.restriction = curvature_restriction;    
    }
  #endif
  boundary((scalar *){normal_vector});
}

/**
We have to use the same version of the fraction field as the one used to define
the field to extend. If not, you may have some errors. */ 

void extend_face_vector (face vector v, face vector ev, scalar f) {
  face vector interfacial_face[];
  foreach_face() {
    interfacial_face.x[] = (interfacial(point, f) || interfacial(neighborp(-1), f) ? 1 : 0);
    ev.x[] = (interfacial_face.x[] > 0 ? v.x[] : nodata);
  }
  #if TREE
    foreach_dimension() {
      ev.x.prolongation = ev.x.refine = curvature_prolongation;
      ev.x.restriction = curvature_restriction;    
    }
  #endif
  boundary((scalar *){interfacial_face, ev});
  foreach_face() {
    double count = 0., v_sum = 0.;
    if (interfacial_face.x[] == 0) {
      foreach_neighbor(1) {
        if (interfacial_face.x[] > 0) {
            v_sum += ev.x[];
            count += 1.;
        }
      }
      ev.x[] = (count > 0. ? v_sum/count : nodata);
    }
  }
  boundary((scalar *){ev});
}

/**
Version of extend_face_vector, with 0 instead of nodata and no particular refine
functions. */

void extend_face_vector_w0 (face vector v, face vector ev, scalar f) {
  face vector interfacial_face[];
  foreach_face() {
    interfacial_face.x[] = (interfacial(point, f) || interfacial(neighborp(-1), f) ? 1 : 0);
    ev.x[] = (interfacial_face.x[] > 0 ? v.x[] : 0.);
  }
  boundary((scalar *){interfacial_face, ev});
  foreach_face() {
    double count = 0., v_sum = 0.;
    if (interfacial_face.x[] == 0) {
      foreach_neighbor(1) {
        if (interfacial_face.x[] > 0) {
            v_sum += ev.x[];
            count += 1.;
        }
      }
      ev.x[] = (count > 0. ? v_sum/count : 0.);
    }
  }
  boundary((scalar *){ev});
}

void restrict_face_vector (face vector v, face vector ev, scalar f) {
  #if TREE
    foreach_dimension() {
      ev.x.prolongation = ev.x.refine = curvature_prolongation;
      ev.x.restriction = curvature_restriction;    
    }
  #endif
  foreach_face()
    ev.x[] = (interfacial(point, f) || interfacial(neighborp(-1), f) ? v.x[] : nodata);
  boundary((scalar *){ev});
}
