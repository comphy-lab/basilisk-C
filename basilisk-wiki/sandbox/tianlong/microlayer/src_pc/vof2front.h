#include "fractions.h"
#include "distance.h"

#define FRONTTHR 1e-8


double getSignedDistanceCell(const coord * ifsg /*interface segments*/, const int num_seg,
                         const double color_center, const coord cc /*cell center*/)
{
    double distance = 1.0e10;
    for (int i = 0; i < num_seg; i++) {
        const coord * p0 = &ifsg[i*2];
        const coord * p1 = &ifsg[i*2 + 1];
        coord p_closest = {0.0};
        double sp = 0.0;

        double dnew = PointSegmentDistance(&cc, p0, p1, &p_closest, &sp);
        dnew = sqrt(dnew);

        if(dnew < distance)
        {
            distance = dnew;
        }
    }

    return color_center * distance;
}

void extendVariables(const double i_region, const scalar color_cc, const vector np, 
                     const scalar is_interf, const scalar is_cut, const bool fix_cut, scalar var_e)
{
    scalar dvar[];
    boundary({var_e});
    const int num_iter = 10;

    //we sue the neigbour values to set the initial values before extending
    foreach()
    {
#if USE_MY_SOLID
        if ((int)is_solid[] == 1)
            continue;
#endif
        bool is_extend = color_cc[] * i_region < 0.0 && is_interf[];
        int conf = (int)is_cut[];
        if (fix_cut)
        {
            is_extend = is_extend && conf == 0;
        }
        if (is_extend)
        {
            double extend_init = 0.0;
            int num_nonzero = 0;
            foreach_neighbor(1)
            {
#if USE_MY_SOLID
                if ((int)is_solid[] == 1)
                    continue;
#endif
                bool is_exist = color_cc[] * i_region > 0.0 && is_interf[];
                if (is_exist && var_e[] != 0.)
                {
                    extend_init += var_e[];
                    num_nonzero += 1;
                }
            }

            var_e[] = num_nonzero == 0 ? 0.0 : extend_init / (double) num_nonzero;            
        }
    }

    for (int it = 0; it < num_iter; it++)
    {
        // increment
        foreach ()
        {
            const double dtau = 0.5 * Delta; // psudo time step
            dvar[] = 0.0;
#if USE_MY_SOLID
            if((int)is_solid[] == 1)
                continue;
#endif
            const double cc = color_cc[];
            bool is_extend = cc * i_region < 0.0 && is_interf[];
            int conf = (int) is_cut[];
            if(fix_cut)
            {
                is_extend = is_extend && conf == 0;
            }
            if(is_extend)
            {
                foreach_dimension()
                {
                    double var = var_e[];
                    double var_m = var_e[-1];
                    double var_p = var_e[1];
                    double nd = - i_region * np.x[];
                    double dd = nd > 0.0 ? var_p - var : var - var_m;
                    dvar[] += dd * nd / Delta;
                }

                dvar[] *= dtau;
            }
        }

        // add
        foreach ()
        {
#if USE_MY_SOLID
            if((int)is_solid[] == 1)
                continue;
#endif
            if(is_interf[])
            {
                var_e[] += dvar[];
            }
        }

        boundary({var_e});
    }
}

void getSignedDistance(scalar color_cc, vector normal, scalar c, scalar phi_dis)
{
  foreach ()
  {
    phi_dis[] = color_cc[] * 10.0 * Delta;
#if USE_MY_SOLID
    if ((int)is_solid[] == 1 && (int)is_solid[1] != 0)
      continue;
#endif
    bool is_interf = false;

    coord *ifsg; // interface segments
    Array *arr;
    // make sure it's an interfacial cell
    foreach_neighbor(2)
    {
#if USE_MY_SOLID
      if ((int)is_solid[] == 1)
        continue;
#endif
      int conf = (c[] > 0.0 && c[] < 1.0);
      if (conf > 0)
      {
        // it's okay to set the is_interf multiple times
        is_interf = true;
      }
    }

    if (is_interf)
    {
      arr = array_new();
      foreach_neighbor(2)
      {
#if USE_MY_SOLID
        if ((int)is_solid[] == 1)
          continue;
#endif
        if (c[] > 0.0 && c[] < 1.0)
        {
          coord n = {normal.x[], normal.y[], 0};
          double alpha = plane_alpha(c[], n);
          coord segment[2];
          if (facets(n, alpha, segment) == 2)
          {
            segment[0] = (coord){x + segment[0].x * Delta, y + segment[0].y * Delta, 0.0};
            segment[1] = (coord){x + segment[1].x * Delta, y + segment[1].y * Delta, 0.0};
            coord temp1 = segment[0];
            coord temp2 = segment[1];
            array_append(arr, &temp1, sizeof(coord));
            array_append(arr, &temp2, sizeof(coord));
          }
        }
      }

      int num_points = arr->len / sizeof(coord);
      ifsg = (coord *)array_shrink(arr); // interface segments

      double color_center = color_cc[];
      coord pos_cc = {x, y, 0.0};

      phi_dis[] = getSignedDistanceCell(ifsg, num_points / 2, color_center, pos_cc);
      free(ifsg);
    }
  }
  boundary({phi_dis});
  return;
}