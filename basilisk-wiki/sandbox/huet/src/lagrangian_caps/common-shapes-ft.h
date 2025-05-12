/**
# Common initial shapes for immersed membranes
*/

#ifndef NLP
  #define NLP 100
#endif
#ifndef RADIUS
  #define RADIUS 1.
#endif

/**
## 2D circular membrane
*/
struct _initialize_circular_capsule {
  lagMesh* mesh;
  double radius;
  int nln;
  int level;
  double inclination;
  coord shift;
  bool disregard_shift;
  bool verbose;
};

#if dimension < 3
void initialize_circular_capsule(struct _initialize_circular_capsule p) {
  double radius = (p.radius) ? p.radius : RADIUS;
  int nln = (p.nln) ? p.nln : NLP;
  coord shift;
  if (p.shift.x || p.shift.y || p.shift.z)
    {shift.x = p.shift.x; shift.y = p.shift.y; shift.z = p.shift.z;}
  else {shift.x = 0.; shift.y = 0.; shift.z = 0.;}
  p.mesh->nln = nln;
  p.mesh->nle = nln;
  p.mesh->nodes = malloc(nln*sizeof(lagNode));
  p.mesh->edges = malloc(nln*sizeof(Edge));

  double alpha = 2*pi/(nln);
  /** Fill the array of nodes */
  for(int i=0; i<nln; i++) {
    p.mesh->nodes[i].pos.x = radius*cos(alpha*i) + shift.x;
    p.mesh->nodes[i].pos.y = radius*sin(alpha*i) + shift.y;
    p.mesh->nodes[i].edge_ids[0] = -1;
    p.mesh->nodes[i].edge_ids[1] = -1;
    foreach_dimension() {
      p.mesh->nodes[i].lagForce.x = 0.;
      p.mesh->nodes[i].lagVel.x = 0.;
    }
  }
  /** Fill the array of edges.
  For the last edge, the next vertex id is 0 */
  for(int i=0; i<p.mesh->nle; i++) {
    p.mesh->edges[i].node_ids[0] = i;
    if (p.mesh->nodes[i].edge_ids[0] < 0)
      p.mesh->nodes[i].edge_ids[0] = i;
    else
      p.mesh->nodes[i].edge_ids[1] = i;
    int next_vertex_id = (i+1<nln) ? i+1 : 0;
    p.mesh->edges[i].node_ids[1] = next_vertex_id;
    if (p.mesh->nodes[next_vertex_id].edge_ids[0] < 0)
      p.mesh->nodes[next_vertex_id].edge_ids[0] = i;
    else
      p.mesh->nodes[next_vertex_id].edge_ids[1] = i;
  }
  /** The above procedure switches the two egde ids for the first node, which
  we correct below */
  p.mesh->nodes[0].edge_ids[0] = p.mesh->nle-1;
  p.mesh->nodes[0].edge_ids[1] = 0;
  correct_lag_pos(p.mesh);
  for(int i=0.; i<p.mesh->nle; i++) {
    p.mesh->edges[i].l0 = edge_length(p.mesh, i);
    p.mesh->edges[i].length = p.mesh->edges[i].l0;
  }

  #ifdef CAPS_VISCOSITY
    fraction(prevI, sq(radius) - sq(x - shift.x) - sq(y - shift.y));
  #endif

  comp_centroid(p.mesh);
  comp_volume(p.mesh);
  p.mesh->initial_volume = p.mesh->volume;
}

/**
## 2D biconcave membrane
*/
void initialize_biconcave_capsule(struct _initialize_circular_capsule p) {
  double radius = (p.radius) ? p.radius : RADIUS;
  int nln = (p.nln) ? p.nln : NLP;
  double inclination = (p.inclination) ? p.inclination : 0.;
  coord shift;
  if (p.shift.x || p.shift.y || p.shift.z)
    {shift.x = p.shift.x; shift.y = p.shift.y; shift.z = p.shift.z;}
  else {shift.x = 0.; shift.y = 0.; shift.z = 0.;}
  double c = 1.3858189;
  p.mesh->nln = nln;
  p.mesh->nle = nln;
  p.mesh->nodes = malloc(nln*sizeof(lagNode));
  p.mesh->edges = malloc(nln*sizeof(Edge));

  double alpha = 2*pi/(nln);
  /** Fill the array of nodes */
  for(int i=0; i<nln; i++) {
    double nrx = radius*c*cos(alpha*i);
    double nry = .5*radius*c*sin(alpha*i)*(0.207 +
      2.003*sq(cos(2*pi-alpha*i)) - 1.123*sq(sq(cos(2*pi-alpha*i))));
    p.mesh->nodes[i].pos.x = nrx*cos(inclination) - nry*sin(inclination) +
      shift.x;
    p.mesh->nodes[i].pos.y = nrx*sin(inclination) + nry*cos(inclination) +
      shift.y;
    p.mesh->nodes[i].edge_ids[0] = -1;
    p.mesh->nodes[i].edge_ids[1] = -1;
    foreach_dimension() {
      p.mesh->nodes[i].lagForce.x = 0.;
      p.mesh->nodes[i].lagVel.x = 0.;
    }
  }
  /** Fill the array of edges.
  For the last edge, the next vertex id is 0 */
  for(int i=0; i<p.mesh->nle; i++) {
    p.mesh->edges[i].node_ids[0] = i;
    if (p.mesh->nodes[i].edge_ids[0] < 0)
      p.mesh->nodes[i].edge_ids[0] = i;
    else
      p.mesh->nodes[i].edge_ids[1] = i;
    int next_vertex_id = (i+1<nln) ? i+1 : 0;
    p.mesh->edges[i].node_ids[1] = next_vertex_id;
    if (p.mesh->nodes[next_vertex_id].edge_ids[0] < 0)
      p.mesh->nodes[next_vertex_id].edge_ids[0] = i;
    else
      p.mesh->nodes[next_vertex_id].edge_ids[1] = i;
  }
  /** The above procedure switches the two egde ids for the first node, which
  we correct below */
  p.mesh->nodes[0].edge_ids[0] = p.mesh->nle-1;
  p.mesh->nodes[0].edge_ids[1] = 0;
  correct_lag_pos(p.mesh);
  for(int i=0.; i<p.mesh->nle; i++) {
    p.mesh->edges[i].l0 = edge_length(p.mesh, i);
    p.mesh->edges[i].length = p.mesh->edges[i].l0;
  }
  comp_centroid(p.mesh);
  comp_volume(p.mesh);
  p.mesh->initial_volume = p.mesh->volume;
}


/**
## 2D elliptic membrane
*/
struct _initialize_elliptic_capsule {
  lagMesh* mesh;
  double a;
  double b;
  double nln;
  double inclination;
};

void initialize_elliptic_capsule(struct _initialize_elliptic_capsule p) {
  double a = (p.a) ? p.a : RADIUS;
  double b = (p.b) ? p.b : RADIUS;
  int nln = (p.nln) ? p.nln : NLP;
  // double inclination = (p.inclination) ? p.inclination : 0.;
  p.mesh->nln = nln;
  p.mesh->nle = nln;
  p.mesh->nodes = malloc(nln*sizeof(lagNode));
  p.mesh->edges = malloc(nln*sizeof(Edge));

  double alpha = 2*pi/(nln);
  /** Fill the array of nodes */
  for(int i=0; i<nln; i++) {
    p.mesh->nodes[i].pos.x = a*cos(alpha*i);
    p.mesh->nodes[i].pos.y = b*sin(alpha*i);
    p.mesh->nodes[i].edge_ids[0] = -1;
    p.mesh->nodes[i].edge_ids[1] = -1;
    foreach_dimension() {
      p.mesh->nodes[i].lagForce.x = 0.;
      p.mesh->nodes[i].lagVel.x = 0.;
    }
  }
  /** Fill the array of edges.
  For the last edge, the next vertex id is 0 */
  for(int i=0; i<p.mesh->nle; i++) {
    p.mesh->edges[i].node_ids[0] = i;
    if (p.mesh->nodes[i].edge_ids[0] < 0)
      p.mesh->nodes[i].edge_ids[0] = i;
    else
      p.mesh->nodes[i].edge_ids[1] = i;
    int next_vertex_id = (i+1<nln) ? i+1 : 0;
    p.mesh->edges[i].node_ids[1] = next_vertex_id;
    if (p.mesh->nodes[next_vertex_id].edge_ids[0] < 0)
      p.mesh->nodes[next_vertex_id].edge_ids[0] = i;
    else
      p.mesh->nodes[next_vertex_id].edge_ids[1] = i;
  }
  /** The above procedure switches the two egde ids for the first node, which
  we correct below */
  p.mesh->nodes[0].edge_ids[0] = p.mesh->nle-1;
  p.mesh->nodes[0].edge_ids[1] = 0;
  correct_lag_pos(p.mesh);
  for(int i=0.; i<p.mesh->nle; i++) {
    p.mesh->edges[i].l0 = edge_length(p.mesh, i);
    p.mesh->edges[i].length = p.mesh->edges[i].l0;
  }

  #ifdef CAPS_VISCOSITY
    fraction(prevI, 1 - sq(x/a) - sq(y/b));
  #endif

  comp_centroid(p.mesh);
  comp_volume(p.mesh);
  p.mesh->initial_volume = p.mesh->volume;
}

#else // dimension > 2
/**
## 3D icosahedron
*/
void initialize_icosahedron(struct _initialize_circular_capsule p) {
  double radius = (p.radius) ? p.radius : RADIUS;
  p.mesh->nln = 12;
  p.mesh->nodes = malloc(p.mesh->nln*sizeof(lagNode));
  p.mesh->nle = 0;
  p.mesh->edges = malloc(30*sizeof(Edge));
  p.mesh->nlt = 0;
  p.mesh->triangles = malloc(20*sizeof(Triangle));

  /**
  * Create the Lagrangian nodes
  The nodes of a regular icosahedron are located at the vertices of [three
  mutually orthogonal golden rectangles](https://en.wikipedia.org/wiki/Regular_icosahedron#/media/File:Icosahedron-golden-rectangles.svg).
   */
  double sl, ll, r;
  r = sqrt(1. + sq(.5*(1. + sqrt(5))));
  sl = radius*1./r; // sl for "small length"
  ll = radius*.5*(1 + sqrt(5))/r; // ll for "large length"
  coord c[3] = {{0, sl, ll}, {sl, ll, 0}, {ll, 0, sl}};
  int s[9] = {0, 1, 1, -1, -1, -1, 1, -1, 1}; // s for "signs"
  for(int i=0; i<3; i++) {
    for(int j=0; j<4; j++) {
      p.mesh->nodes[i*4+j].nb_neighbors = 0;
      p.mesh->nodes[i*4+j].nb_triangles = i;
      foreach_dimension()
        p.mesh->nodes[i*4+j].pos.x = c[i].x*
          s[(1+j+4*((fabs(c[i].x - sl) < 1.e-8) ? 0 : 1))*
            ((fabs(c[i].x) < 1.e-8) ? 0 : 1)];
    }
  }

  /**
  * Create the edges
  */
  for(int i=0; i<p.mesh->nln; i++) {
    int my_gr = p.mesh->nodes[i].nb_triangles; // gr for "golden ratio"
    int my_ld_sign = GET_LD_SIGN(p.mesh->nodes[i]); // ld for "large distance"
    for(int j=0; j<p.mesh->nln; j++) {
      if (p.mesh->nodes[j].nb_triangles == my_gr && j != i) {
        if (my_ld_sign == GET_LD_SIGN(p.mesh->nodes[j])) {
          if (write_edge(p.mesh, p.mesh->nle, i, j, new_mesh = true))
            p.mesh->nle++;
        }
      }
      else if (GET_LD(p.mesh->nodes[j]) == GET_ZD(p.mesh->nodes[i])) {
        if (GET_SD_SIGN(p.mesh->nodes[j]) == GET_LD_SIGN(p.mesh->nodes[i])) {
          if (write_edge(p.mesh, p.mesh->nle, i, j, new_mesh = true))
            p.mesh->nle++;
        }
      }
    }
  }

  /**
  * Create the triangular faces
  */
  for(int i=0; i<p.mesh->nln; i++) {
    p.mesh->nodes[i].nb_triangles = 0;
  }
  for(int i=0; i<p.mesh->nle; i++) {
    int nid[2];
    for(int j=0; j<2; j++) nid[j] = p.mesh->edges[i].node_ids[j];
    for(int j=0; j<p.mesh->nodes[nid[0]].nb_neighbors; j++) {
      for(int k=0; k<p.mesh->nodes[nid[1]].nb_neighbors; k++) {
        if (p.mesh->nodes[nid[0]].neighbor_ids[j] ==
          p.mesh->nodes[nid[1]].neighbor_ids[k] && p.mesh->nodes[nid[0]].neighbor_ids[j] != -1)
          if (write_triangle(p.mesh, p.mesh->nlt, nid[0], nid[1],
            p.mesh->nodes[nid[0]].neighbor_ids[j]))
              p.mesh->nlt++;
      }
    }
  }

  comp_centroid(p.mesh);
  comp_volume(p.mesh);
  p.mesh->initial_volume = p.mesh->volume;
}

/**
## 3D spherical membrane
*/
/** The function below triangulates a sphere: it starts from an icosahedron,
subdivides each of its triangles into four smaller ones until the desired number
of Lagrangian nodes is reached or exceeded, and projects the resulting mesh
onto a sphere. */
void initialize_spherical_capsule(struct _initialize_circular_capsule p) {
  initialize_icosahedron(p);

  double radius = (p.radius) ? p.radius : RADIUS;
  int nln = (p.nln) ? p.nln : -1;
  int ns = (p.level) ? p.level : -1;
  coord shift;
  if (p.shift.x || p.shift.y || p.shift.z)
    {shift.x = p.shift.x; shift.y = p.shift.y; shift.z = p.shift.z;}
  else {shift.x = 0.; shift.y = 0.; shift.z = 0.;}
  bool verbose = p.verbose ? true : false;

  /** For each subdivision:
    * the number of additional nodes equals the number of edges
    * the number of edges is multiplied by 4: each edge is split in two, and
    each new nodes results in two new edges (one per neighboring triangles)
    * the number of triangles is multiplied by 4
  */
  int nn, ne, nt; // the numbers of nodes, edges and triangles.
  nn = 12; // at first, an icosahedron has 12 nodes
  if (nln > 0) {
    ns = 0;
    while (nn < nln) {
      ne = 30*pow(4,ns);
      nn += ne;
      ns++;
    }
  }
  if (ns > 0) {
    int cns = 0;
    while (cns < ns) {
      ne = 30*pow(4,cns);
      nn += ne;
      cns++;
    }
  }
  else {
    fprintf(stderr, "Error: number of Lagrangian nodes or Lagrangian mesh "
      "levels need to be specified.\n");
  }
  ne = 30*pow(4,ns);
  nt = 20*pow(4,ns);
  p.mesh->nodes = realloc(p.mesh->nodes, nn*sizeof(lagNode));
  p.mesh->edges = realloc(p.mesh->edges, ne*sizeof(Edge));
  p.mesh->triangles = realloc(p.mesh->triangles, nt*sizeof(Triangle));

  /** It's time to perform our $ns$ triangle subdivisions */
  for(int i=0; i<ns; i++) refine_mesh(p.mesh);

  /** At last, we project each node onto a sphere of desired radius */
  for(int i=0; i<p.mesh->nln; i++) {
    double cr = 0.;
    foreach_dimension() cr += sq(p.mesh->nodes[i].pos.x);
    cr = sqrt(cr);
    foreach_dimension() p.mesh->nodes[i].pos.x *= radius/cr;
  }
  
  if (verbose) {
    fprintf(stderr, "Number of triangle refinements: %d\n", ns);
    fprintf(stderr, "Number of Lagrangian nodes: %d\n", p.mesh->nln);
    fprintf(stderr, "Number of Lagrangian edges: %d\n", p.mesh->nle);
    fprintf(stderr, "Number of Lagrangian triangles: %d\n", p.mesh->nlt);
  }

  comp_initial_area_normals(p.mesh);
  if (!p.disregard_shift) {
    if (shift.x > 1.e-10 || shift.y > 1.e-10 || shift.z > 1.e-10) {
      for(int i=0; i<p.mesh->nln; i++)
        foreach_dimension()
          p.mesh->nodes[i].pos.x += shift.x;
    }

    #ifdef CAPS_VISCOSITY
      fraction(I, sq(RADIUS) - sq(x - shift.x) - sq(y - shift.y)
        - sq(z - shift.z));
      foreach() if (I[] > 1.e-6) prevI[] = I[];
    #endif
  }
  correct_lag_pos(p.mesh);
  comp_normals(p.mesh);
  comp_centroid(p.mesh);
  comp_volume(p.mesh);
  p.mesh->initial_volume = p.mesh->volume;
}

void initialize_rbc_capsule(struct _initialize_circular_capsule p) {
  initialize_spherical_capsule(mesh = p.mesh, radius = p.radius, nln = p.nln,
    level = p.level, disregard_shift = true);

  double c0, c1, c2;
  c0 = 0.2072; c1 = 2.0026; c2 = -1.1228;
  double radius = (p.radius) ? p.radius : RADIUS;
  for(int i=0; i<p.mesh->nln; i++) {
    double rho = sqrt(sq(p.mesh->nodes[i].pos.x) +
      sq(p.mesh->nodes[i].pos.z))/radius;
    rho = (rho > 1) ? 1 : rho;
    int sign = (p.mesh->nodes[i].pos.y > 0.) ? 1 : -1;
    p.mesh->nodes[i].pos.y = sign*.5*radius*sqrt(1 - sq(rho))*
      (c0 + c1*sq(rho) + c2*sq(sq(rho)));
  }
  correct_lag_pos(p.mesh);
  comp_normals(p.mesh);

  coord shift;
  if (p.shift.x || p.shift.y || p.shift.z)
    {shift.x = p.shift.x; shift.y = p.shift.y; shift.z = p.shift.z;}
  else {shift.x = 0.; shift.y = 0.; shift.z = 0.;}
  if (fabs(shift.x) > 1.e-10 || fabs(shift.y) > 1.e-10 ||
    fabs(shift.z) > 1.e-10) {
    for(int i=0; i<p.mesh->nln; i++)
      foreach_dimension()
        p.mesh->nodes[i].pos.x += shift.x;
  }

  #ifdef CAPS_VISCOSITY
    double a, c;
    c = 1.3858189;
    a = RADIUS/c;
    // We define below the local coordinates of the RBC and the parametric angle
    #define COSPHI2 ((sq(x - shift.x) + sq(z - shift.z))/sq(a*c))
    #define RHS (0.207 + 2.003*COSPHI2 - 1.123*sq(COSPHI2))
    fraction(I, 1. - sq((x - shift.x)/(a*c)) -
      sq(2*(y - shift.y)/(RHS*a*c)) - sq((z - shift.z)/(a*c)));
    foreach() if (I[] > 1.e-6) prevI[] = I[];
  #endif
  correct_lag_pos(p.mesh);
  comp_centroid(p.mesh);
  comp_volume(p.mesh);
  p.mesh->initial_volume = p.mesh->volume;
}

void activate_spherical_capsule(struct _initialize_circular_capsule p) {
  initialize_empty_capsule(p.mesh);
  p.mesh->isactive = true;
  initialize_spherical_capsule(p);
  initialize_capsule_stencils(p.mesh);
  generate_lag_stencils(no_warning = true);
  #if _ELASTICITY_FT
    store_initial_configuration(p.mesh);
  #endif
  #if _BENDING_FT
    initialize_refcurv_onecaps(p.mesh);
  #endif
}

void activate_biconcave_capsule(struct _initialize_circular_capsule p) {
  initialize_empty_capsule(p.mesh);
  p.mesh->isactive = true;
  initialize_rbc_capsule(p);
  initialize_capsule_stencils(p.mesh);
  generate_lag_stencils(no_warning = true);
  #if _ELASTICITY_FT
    store_initial_configuration(p.mesh);
  #endif
  #if _BENDING_FT
    initialize_refcurv_onecaps(p.mesh);
  #endif
}

#endif
