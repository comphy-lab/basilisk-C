/**
The interface reconstruction is modified here, in order to deal with the two-phase cells in the three-phase flow simulation. */

#include "geometry.h"
#if dimension == 1
coord mycs(Point point, scalar c)
{
    coord n = {1.};
    return n;
}
#elif dimension == 2
#include "myc2d.h"
#else // dimension == 3
#include "myc.h"
#endif

#if TREE

void fraction_refine(Point point, scalar c)
{
    double cc = c[];
    if (cc <= 0. || cc >= 1.)
        foreach_child()
            c[] = cc;
    else
    {
        coord n = mycs(point, c);
        double alpha = plane_alpha(cc, n);
        coord a, b;
        foreach_dimension()
        {
            a.x = 0.;
            b.x = 0.5;
        }
        foreach_child()
        {
            coord nc;
            foreach_dimension()
                nc.x = child.x * n.x;
            c[] = rectangle_fraction(nc, alpha, a, b);
        }
    }
}

attribute
{
    vector n;
}

static void alpha_refine(Point point, scalar alpha)
{
    vector n = alpha.n;
    double alphac = 2. * alpha[];
    coord m;
    foreach_dimension()
        m.x = n.x[];
    foreach_child()
    {
        alpha[] = alphac;
        foreach_dimension()
            alpha[] -= child.x * m.x / 2.;
    }
}

#endif // TREE

struct Fractions
{
    vertex scalar Phi; // compulsory
    scalar c;          // compulsory
    face vector s;     // optional
};

trace void fractions(struct Fractions a)
{
    vertex scalar Phi = a.Phi;
    scalar c = a.c;
    face vector s = automatic(a.s);

#if dimension == 3
    vector p[];
#else // dimension == 2
    vector p;
    p.x = s.y;
    p.y = s.x;
#endif

    foreach_edge()
    {
        if (Phi[] * Phi[1] < 0.)
        {
            p.x[] = Phi[] / (Phi[] - Phi[1]);
            if (Phi[] < 0.)
                p.x[] = 1. - p.x[];
        }
        else
            p.x[] = (Phi[] > 0. || Phi[1] > 0.);
    }

#if dimension == 3
    scalar s_x = s.x, s_y = s.y, s_z = s.z;
    foreach_face(z, x, y)
#else // dimension == 2
    boundary_flux({s});
    scalar s_z = c;
    foreach ()
#endif
    {
        coord n;
        double nn = 0.;
        foreach_dimension(2)
        {
            n.x = p.y[] - p.y[1];
            nn += fabs(n.x);
        }
        if (nn == 0.)
            s_z[] = p.x[];
        else
        {
            foreach_dimension(2)
                n.x /= nn;
            double alpha = 0., ni = 0.;
            for (int i = 0; i <= 1; i++)
            {
                foreach_dimension(2)
                if (p.x[0, i] > 0. && p.x[0, i] < 1.)
                {
                    double a = sign(Phi[0, i]) * (p.x[0, i] - 0.5);
                    alpha += n.x * a + n.y * (i - 0.5);
                    ni++;
                }
            }
            s_z[] = ni ? line_area(n.x, n.y, alpha / ni) : max(p.x[], p.y[]);
        }
    }

#if dimension == 3
    boundary_flux({s});
    foreach ()
    {
        coord n;
        double nn = 0.;
        foreach_dimension(3)
        {
            n.x = s.x[] - s.x[1];
            nn += fabs(n.x);
        }
        if (nn == 0.)
            c[] = s.x[];
        else
        {
            foreach_dimension(3)
                n.x /= nn;
            double alpha = 0., ni = 0.;
            for (int i = 0; i <= 1; i++)
                for (int j = 0; j <= 1; j++)
                    foreach_dimension(3)
                    if (p.x[0, i, j] > 0. && p.x[0, i, j] < 1.)
                    {
                        double a = sign(Phi[0, i, j]) * (p.x[0, i, j] - 0.5);
                        alpha += n.x * a + n.y * (i - 0.5) + n.z * (j - 0.5);
                        ni++;
                    }
            c[] = ni ? plane_volume(n, alpha / ni) : s.x[];
        }
    }
#endif
    boundary({c});
}

#define fraction(f, func)    \
    do                       \
    {                        \
        vertex scalar phi[]; \
        foreach_vertex()     \
            phi[] = func;    \
        fractions(phi, f);   \
    } while (0)

coord youngs_normal(Point point, scalar c)
{
    coord n;
    double nn = 0.;
    assert(dimension == 2);
    foreach_dimension()
    {
        n.x = (c[-1, 1] + 2. * c[-1, 0] + c[-1, -1] - c[+1, 1] - 2. * c[+1, 0] - c[+1, -1]);
        nn += fabs(n.x);
    }
    // normalize
    if (nn > 0.)
        foreach_dimension()
            n.x /= nn;
    else // this is a small fragment
        n.x = 1.;
    return n;
}

/**
### Normal approximation using MYC or face fractions
*/

coord facet_normal (Point point, scalar c, face vector s)
{
  coord n;
  if (s.x.i < 0)
    n = mycs (point, c);
  else { // compute normal from face fractions
    double nn = 0.;
    foreach_dimension() {
      n.x = s.x[] - s.x[1];
      nn += fabs(n.x);
    }
    assert (nn > 0.);
    foreach_dimension()
      n.x /= nn;
  }
  return n;
}

/**
Here is part that is modified for the three-phase flow simulation. 

The reconstruction function takes a volume fraction field `c` and
returns the corresponding normal vector field `n` and intercept field
$\alpha$. */

trace void reconstruction(const scalar c, vector n, scalar alpha)
{
    foreach ()
    {
        if (c[] <= 0. || c[] >= 1.)
        {
            alpha[] = 0.;
            foreach_dimension()
                n.x[] = 0.;
        }
        else
        {
/** Change to one interface for two-phase cells of the three-phase flow simulation, instead of having two not necessarily same interface for the two-phase cells.
It should be mentioned that only around the triple point we might have some two-phase flow cells with two different normals for the interface. */
            coord m1, m2, m1f, m2f, m, ma;
            double sum;
            char flag = 'n';
/** First the interface normal is determined as follow, */
            m = mycs(point, c);
/** Chesking which phases are in the cell and then calculate the interface normals,*/
            if (f1[] == 0. && f2[] > 0. && f2[] < 1. && f3[] > 0. && f3[] < 1.)
            {
                m1 = mycs(point, f2);
                m2 = mycs(point, f3);
/** Also, we are calculating the volumetric portion of each interface normal,*/
                foreach_dimension()
                {
                    m1f.x = m1.x * f2[];
                    m2f.x = m2.x * f3[];
                }
                flag = 'y';
            }
            else if (f2[] == 0. && f3[] > 0. && f3[] < 1. && f1[] > 0. && f1[] < 1.)
            {
                m1 = mycs(point, f3);
                m2 = mycs(point, f1);
                foreach_dimension()
                {
                    m1f.x = m1.x * f3[];
                    m2f.x = m2.x * f1[];
                }
                flag = 'y';
            }
            else if (f3[] == 0. && f1[] > 0. && f1[] < 1. && f2[] > 0. && f2[] < 1.)
            {
                m1 = mycs(point, f1);
                m2 = mycs(point, f2);
                foreach_dimension()
                {
                    m1f.x = m1.x * f1[];
                    m2f.x = m2.x * f2[];
                }
                flag = 'y';
            }
/** After all the combinations are checked, we will find which interface is now being calculated as *c* in the reconstruction function; and, we find the other interface and reverse it for the averaging process.  */
            if (flag == 'y')
            {
                if (fabs(m1.x - m.x) < 1.0e-6 && fabs(m1.y - m.y) < 1.0e-6)
                {
                    foreach_dimension()
                        m2f.x *= -1.0;
                }
                else if (fabs(m2.x - m.x) < 1.0e-6 && fabs(m2.y - m.y) < 1.0e-6)
                {
                    foreach_dimension()
                        m1f.x *= -1.0;
                }
                else
                {
                    printf ("error in finding the interface\r\n");
                }
/** Performing the volume average for the interface normals, */
                foreach_dimension()
                    ma.x = m1f.x + m2f.x;
                sum = 0.0;
                foreach_dimension()
                    sum += fabs(ma.x);
                foreach_dimension()
                    ma.x /= sum;
                 foreach_dimension()
                    m.x = ma.x;
            }
            foreach_dimension()
                n.x[] = m.x;
            alpha[] = plane_alpha(c[], m);
        }
    }

#if TREE
    foreach_dimension()
        n.x.refine = n.x.prolongation = refine_injection;
    alpha.n = n;
    alpha.refine = alpha.prolongation = alpha_refine;
#endif
    boundary({n, alpha});
}

struct OutputFacets
{
    scalar c;
    FILE *fp;      // optional: default is stdout
    face vector s; // optional: default is none
};

trace void output_facets(struct OutputFacets p)
{
    scalar c = p.c;
    face vector s = p.s;
    if (!p.fp)
        p.fp = stdout;

    foreach ()
        if (c[] > 1e-6 && c[] < 1. - 1e-6)
        {
            coord n;
            if (!s.x.i)
                // compute normal from volume fraction
                n = mycs(point, c);
            else
            {
                // compute normal from face fractions
                double nn = 0.;
                foreach_dimension()
                {
                    n.x = s.x[] - s.x[1];
                    nn += fabs(n.x);
                }
                assert(nn > 0.);
                foreach_dimension()
                    n.x /= nn;
            }
            double alpha = plane_alpha(c[], n);
#if dimension == 2
            coord segment[2];
            if (facets(n, alpha, segment) == 2)
                fprintf(p.fp, "%g %g\n%g %g\n\n",
                        x + segment[0].x * Delta, y + segment[0].y * Delta,
                        x + segment[1].x * Delta, y + segment[1].y * Delta);
#else // dimension == 3
            coord v[12];
            int m = facets(n, alpha, v, 1.);
            for (int i = 0; i < m; i++)
                fprintf(p.fp, "%g %g %g\n",
                        x + v[i].x * Delta, y + v[i].y * Delta, z + v[i].z * Delta);
            if (m > 0)
                fputc('\n', p.fp);
#endif
        }

    fflush(p.fp);
}

trace double interface_area(scalar c)
{
    double area = 0.;
    foreach ()
        if (c[] > 1e-6 && c[] < 1. - 1e-6)
        {
            coord n = mycs(point, c), p;
            double alpha = plane_alpha(c[], n);
            area += pow(Delta, dimension - 1) * plane_area_center(n, alpha, &p);
        }
    return area;
}
