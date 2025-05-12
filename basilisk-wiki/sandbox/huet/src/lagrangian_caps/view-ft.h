/**
# View front-tracking

This file provides functions to display the Lagrangian mesh. */

#include "view.h"

macro draw_vertices (bview * view, float color[3], float ps)
{
  glMatrixMode (GL_PROJECTION);
  glPushMatrix();
  glTranslatef (0., 0., view->lc*view->fov/24.);
  glColor3f (color[0], color[1], color[2]);
  glEnable(GL_POINTS);
  glPointSize (view->samples*(ps > 0. ? ps : 8.));
  bool _reversed = view->reversed;
  view->reversed = false;
  {...}
  glMatrixMode (GL_PROJECTION);
  glPopMatrix();
  view->reversed = _reversed;
}

/**
The function `draw_lag` displays a triangulation in the current view. It does
not include shading effects, so displaying the edges of the triangles is
necessary to get a three-dimensional perspective.

The arguments are:

* `mesh`: the Lagrangian mesh, aka the triangulation to display
* `nodes`: a boolean which sets whether to display Lagrangian nodes or not. It
is false by default.
* `edges`: a boolean which sets whether to display Lagrangian edges or not. It
is false by default.
* `facets`: a boolean which sets whether to display triangles or not. It
is false by default.
* `nc`: the nodes color, an array of three doubles (RGB color convention).
* `ec`: the edges color, an array of three doubles (RGB color convention).
* `fc`: the facets color, an array of three doubles (RGB color convention).
* `ns`: a double setting the nodes size. If different than 0, it sets `nodes` to
true.
* `lw`: a double setting the line width of the edges. If different than 0, it 
sets `nodes` to true.
*/

struct _draw_lag {
  lagMesh * mesh; // Compulsory
  bool nodes;
  bool edges;
  bool facets;
  float fc[3], lc[3], nc[3], lw, ns;
};

void draw_lag (struct _draw_lag p)
{
  if (pid() == 0) {
    bool edges = p.edges;
    if (p.lw > 0 || (p.lc[0] > 0 || p.lc[1] > 0 || p.lc[2] > 0))
      edges = true;
    bool nodes = p.nodes;
    if (p.ns > 0)
      nodes = true;
    float my_color[3] = {0};
    bview * view = draw();
    if (edges) {
      p.lw = p.lw ? p.lw : 2.;
      if (p.lc[0])
	my_color[0] = p.lc[0], my_color[1] = p.lc[1], my_color[2] = p.lc[2];
      draw_lines (view, my_color, p.lw) {
        for (int i = 0; i < p.mesh->nle; i++) {
          int v1, v2;
          v1 = p.mesh->edges[i].node_ids[0];
          v2 = p.mesh->edges[i].node_ids[1];
	  coord node1 = correct_periodic_node_pos(p.mesh->nodes[v1].pos,
						  p.mesh->centroid);
	  coord node2 = correct_periodic_node_pos(p.mesh->nodes[v2].pos,
						  p.mesh->centroid);
	  glBegin(GL_LINES);
#if dimension < 3
	  glvertex2d(view, p.mesh->nodes[v1].pos.x,
		     p.mesh->nodes[v1].pos.y);
	  glvertex2d(view, p.mesh->nodes[v2].pos.x,
		     p.mesh->nodes[v2].pos.y);
#else
	  glvertex3d(view, node1.x, node1.y, node1.z);
	  glvertex3d(view, node2.x, node2.y, node2.z);
#endif
	  glEnd();
	  view->ni++;
        }
      }
    }
    if (nodes) {
      p.ns = p.ns ? p.ns : 8.;
      if (p.nc[0])
	my_color[0] = p.nc[0], my_color[1] = p.nc[1], my_color[2] = p.nc[2];
      draw_vertices(view, my_color, p.ns) {
        glBegin(GL_POINTS);
	for (int i = 0; i < p.mesh->nln; i++)
#if dimension < 3
	  glvertex2d(view, p.mesh->nodes[i].pos.x,
		     p.mesh->nodes[i].pos.y);
#else
	glvertex3d(view, p.mesh->nodes[i].pos.x,
		   p.mesh->nodes[i].pos.y, p.mesh->nodes[i].pos.z);
#endif
        glEnd();
        view->ni++;
      }
    }
#if dimension > 2
    bool facets = p.facets;
    if (facets) {
      if (p.fc[0])
	my_color[0] = p.fc[0], my_color[1] = p.fc[1], my_color[2] = p.fc[2];
      else
	my_color[0] = 1., my_color[1] = 1., my_color[2] = 1.;
      for (int i = 0; i < p.mesh->nlt; i++) {
	coord nodes[3];
	for(int j = 0; j < 3; j++)
	  nodes[j] = correct_periodic_node_pos(p.mesh->nodes[p.mesh->triangles[i].node_ids[j]].pos, p.mesh->centroid);
	glBegin (GL_TRIANGLE_FAN);
	for(int j=0; j<3; j++) {
	  glColor3f (my_color[0], my_color[1], my_color[2]);
	  glVertex3d(
		     nodes[j].x,
		     nodes[j].y,
		     nodes[j].z);
	}
	glEnd ();
	view->ni++;
      }
    }
#endif
  }
}

/**
The function below display all active capsules. Its optional arguments are the
same as the `draw_lag` function.
*/

void draw_lags (struct _draw_lag p)
{
  for (int k = 0; k < NCAPS; k++)
    if (CAPS(k).isactive)
      draw_lag (&CAPS(k), nodes = p.nodes, edges = p.edges,
		facets = p.facets,
		fc = {p.fc[0], p.fc[1], p.fc[2]},
		lc = {p.lc[0], p.lc[1], p.lc[2]},
		nc = {p.nc[0], p.nc[1], p.nc[2]},
		lw = p.lw, ns = p.ns);
}
