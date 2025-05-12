#ifndef _UTILS_PHASE_CHANGE_H
#define _UTILS_PHASE_CHANGE_H

/**
This file contains useful functions for simulation phase change flows with EBIT.
*/

#include "PointTriangle.h"
#include "intgrad.h"

static inline void restriction_EBIT_linear (Point point, scalar s)
{  

    extern scalar color_cc;
    double i_region = color_cc[];
    /**
     Otherwise, we use the average of the child cells which are defined
    (there is at least one). */
    
    double val = 0.0;
    int nv = 0;
    foreach_child()
    {
        if(color_cc[] * i_region > 0.0)
        {
            val += s[], nv++;
        }
    }
    if(nv == 0)
    {
        foreach_child()
        {
            val += s[], nv++;
        }
    }
    assert (nv > 0.);
    s[] = val/nv;
}

static inline void restriction_EBIT_face (Point point, scalar s)
{
    vector v = s.v;
    foreach_dimension()
    {
#if dimension == 1
        v.x[] = fine(v.x, 0);
        v.x[1] = fine(v.x, 2);
#elif dimension == 2
        if(fine(v.x,0,0) == 0 || fine(v.x, 0, 1) == 0)
        {
            v.x[] = 0.0;
        }
        else
        {
            v.x[] = (fine(v.x, 0, 0) + fine(v.x, 0, 1)) / 2.;
        }
        if(fine(v.x, 2, 0) == 0 || fine(v.x, 2, 1) == 0)
        {
            v.x[1] = 0.0;
        }
        else
        {
            v.x[1] = (fine(v.x, 2, 0) + fine(v.x, 2, 1)) / 2.;
        }
#else // dimension == 3
        v.x[] = (fine(v.x, 0, 0, 0) + fine(v.x, 0, 1, 0) +
                 fine(v.x, 0, 0, 1) + fine(v.x, 0, 1, 1)) /
                4.;
        v.x[1] = (fine(v.x, 2, 0, 0) + fine(v.x, 2, 1, 0) +
                  fine(v.x, 2, 0, 1) + fine(v.x, 2, 1, 1)) /
                 4.;
#endif
    }
}

static inline void restriction_EBIT_color (Point point, scalar color)
{  
  color[] = f[] > 0.5 ? 1.0 : -1.0;
}

void getInterfaceSegmentsinCell(Point point, Array * arr)
{
    //a strange bug. The funciton will be called at level 0 
    //by using foreach() iterator
    if(point.level == 0)
    {
      return;
    }
    int ii = 0, iip, ind_markers[4] = {-1, -1, -1, -1};
    double xx[4], yy[4], xy_edge[4][2];

    xx[0] = xx[1] = yy[0] = yy[1] = 0.;
    xx[2] = xx[3] = yy[2] = yy[3] = 0.;

    xy_edge[0][0] = x - Delta / 2.;
    xy_edge[0][1] = y - Delta / 2. + s.x[0, 0] * Delta;

    xy_edge[1][0] = x - Delta / 2. + s.y[0, 0] * Delta;
    xy_edge[1][1] = y - Delta / 2.;

    xy_edge[2][0] = x + Delta / 2.;
    xy_edge[2][1] = y - Delta / 2. + s.x[1, 0] * Delta;

    xy_edge[3][0] = x - Delta / 2. + s.y[0, 1] * Delta;
    xy_edge[3][1] = y + Delta / 2.;

    /** Start calculating the connections of markers basing on the
      color vertex (central point)*/
    int conf = (int)config_dict[];

    if (conf == 1)
    { // Connect the opposite edges, one interface
        ind_markers[0] = 1;
        ind_markers[1] = 3;
    }
    else if (conf == 2)
    {
        ind_markers[0] = 0;
        ind_markers[1] = 2;
    }
    else if (conf > 2 && conf < 7)
    { // Connect the consecutive edges, one interface
        ii = conf - 3;
        iip = (ii + 1) % 4;
        ind_markers[0] = ii;
        ind_markers[1] = iip;
    }
    else if (conf == 7)
    { // Connect the consecutive edges, two interfaces
        ind_markers[0] = 3;
        ind_markers[1] = 0;

        ind_markers[2] = 1;
        ind_markers[3] = 2;
    }
    else if (conf == 8)
    {
        ind_markers[0] = 0;
        ind_markers[1] = 1;

        ind_markers[2] = 2;
        ind_markers[3] = 3;
    }

    for (int iv = 0; iv < 4; iv += 2)
    {
        if (ind_markers[iv] != -1)
        {
            int ie1 = ind_markers[iv], ie2 = ind_markers[iv + 1];
            xx[0] = xy_edge[ie1][0];
            yy[0] = xy_edge[ie1][1];

            xx[1] = xy_edge[ie2][0];
            yy[1] = xy_edge[ie2][1];

            coord p = {xx[0], yy[0], 0.0};
            array_append(arr, &p, sizeof(coord));

            coord p2 = {xx[1], yy[1], 0.0};
            array_append(arr, &p2, sizeof(coord));
        }
    }
}

double weno3P(double *f)
{

	//assign value to v1, v2,...
	int k = 0;
	double v1 = *(f + k - 1);
	double v2 = *(f + k);
	double v3 = *(f + k + 1);

	//smoothness indicator
	double epsilon = 1.0e-16;
	double s1 = (v2 - v3)*(v2 - v3) + epsilon;
	double s2 = (v2 - v1)*(v2 - v1) + epsilon;

	//weights
	double a1 = 1.0e1/s1/s1;
	double a2 = 1.0e1/s2/s2;
	double tw1 = 1.0 / (a1 +a2);
	double w1 = tw1*a1;
	double w2 = tw1*a2;

	//return weighted average
	return  w1*(0.5*v2 + 0.5*v3)
		  + w2*(-0.5*v1 + 1.5*v2);
}

double weno3M(double *f)
{

	//assign value to v1, v2,...
	int k = 1;
	double v1 = *(f + k + 1);
	double v2 = *(f + k);
	double v3 = *(f + k - 1);

	//smoothness indicator
	double epsilon = 1.0e-16;
	double s1 = (v2 - v3)*(v2 - v3) + epsilon;
	double s2 = (v2 - v1)*(v2 - v1) + epsilon;

	//weights
	double a1 = 1.0e1/s1/s1;
	double a2 = 1.0e1/s2/s2;
	double tw1 = 1.0 / (a1 +a2);
	double w1 = tw1*a1;
	double w2 = tw1*a2;

	//return weighted average
	return  w1*(0.5*v2 + 0.5*v3)
		  + w2*(-0.5*v1 + 1.5*v2);
}

void setCenterVfromFaceV(vector u_c, face vector u_f)
{
    foreach ()
    {
        foreach_dimension()
        {
            u_c.x[] = (u_f.x[] + u_f.x[1]) / (fm.x[] + fm.x[1]);
        }
    }
}

void setVelPhaseChange(vector u_ns, vector u_pc)
{
    foreach ()
    {
        foreach_dimension()
        {
            u_pc.x[] = 0.0;
        }

        if (is_interfacial[])
        {
            double md = mdot[];
            double rho = color_cc[] > 0.0 ? rho1 : rho2;

            foreach_dimension()
            {
                u_pc.x[] = u_ns.x[];
                double vel_pc = md / rho;
#if !USE_DOUBLE_VEL
                double rhol = face_fraction.x[] > 0.5 ? rho1 : rho2;
                double rhor = face_fraction.x[1] > 0.5 ? rho1 : rho2;
                double fl = fm.x[];
                double fr = fm.x[1];
                vel_pc = (fl * md / rhol + fr * md / rhor) / (fl + fr + SEPS);
#endif             
                u_pc.x[] += vel_pc * n_pc.x[];
            }
        }
    }

    boundary({u_pc});
}

//this function is copied from Edorado's sandbox
coord normal (Point point, scalar c)
{
  coord n = mycs (point, c);
  double nn = 0.;
  foreach_dimension()
    nn += sq(n.x);
  nn = sqrt(nn);
  foreach_dimension()
    n.x /= nn;
  return n;
}

//set the source terms for the projection step
//this function is copied from Edorado's sandbox and is modified a little.
void setStefanFlow(scalar stefanflow)
{
    foreach ()
    {
        stefanflow[] = 0.;
        double alpha = 0.;
        double segment = 0.;
#if USE_MY_SOLID
        if ((int)is_solid[] == 1)
            continue;

#endif
        if ((int)config_dict[] > 0 && f[] > F_ERR && (f[] < 1. - F_ERR))
        {
            //coord m = normal(point, f);
            //this should recover exactly the interface segments obtained by semushin
            coord m = my_facet_normal(point, f, face_fraction);
            alpha = plane_alpha(f[], m);
            coord prel;
            segment = plane_area_center(m, alpha, &prel);
#ifdef AXI
            stefanflow[] = cm[] * mdot[] * segment * (y + prel.y * Delta) * (1. / rho2 - 1. / rho1) / (Delta * y);
#else
            stefanflow[] = mdot[] * segment * (1. / rho2 - 1. / rho1) / Delta * cm[];
#endif
        }
    }
}

void getColorExact(scalar color_cell)
{
    foreach ()
    {
        //-1 will be more convenient to use
        // Note that now the meaning of "-1" is differnt with before
        color_cell[] = f[] >= 0.5 ? 1.0 : -1.0;
    }

    boundary({color_cell});
}

//compute the signed distance for a narrow band near the interface
//as the interface segments are provided by EBIT, we can compute the distance directly.
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

void getSignedDistance(const scalar color_cc, scalar phi_dis)
{
    foreach()
    {
        phi_dis[] = color_cc[] * 10.0 * Delta;
#if USE_MY_SOLID
        if ((int)is_solid[] == 1)
            continue;
#endif
        bool is_interf = false;

        coord *ifsg; // interface segments
        Array *arr;
        //make sure it's an interfacial cell
        foreach_neighbor(2)
        {
            int conf = (int) config_dict[];
            if(conf > 0)
            {
                //it's okay to set the is_interf multiple times
                is_interf = true;
            }
        }

        if(is_interf)
        {
            arr = array_new();
            foreach_neighbor(2)
            {
                int conf = (int) config_dict[];
                if (conf > 0)
                {
                    getInterfaceSegmentsinCell(point, arr);
                }
            }

            int num_points = arr->len / sizeof(coord);
            ifsg = (coord *)array_shrink(arr); // interface segments

            double color_center = color_cc[];
            coord pos_cc = {x, y, 0.0};

            phi_dis[] = getSignedDistanceCell(ifsg,num_points/2, color_center, pos_cc);

            free(ifsg);
        }

    }
    boundary({phi_dis});
    return;
}
//compute the distance from the cell center to the interface segemet if the colors of two neighbours are different.
void getDistanceC2I(const double i_region, const scalar color_cell, const scalar phi, face vector dis_c2i)
{
    //as we will only use dict in user defined function, where the foreach_dimension doesn't work
    //we only neede it here.
    update_dict_x();
    int idim = 0;
    foreach_dimension()
    {
        idim++;
        foreach_face(x)
        {
            dis_c2i.x[] = 0.0;
#if USE_MY_SOLID
            if ((int) is_solid_face.x[] == 1)
                continue;          
#endif
            dis_c2i.x[] = 0.0;
            bool is_interface_exist = color_cell[] * color_cell[-1, 0] < 0.0;
            double color_0 = color_cell[];
            double color_1 = color_cell[-1, 0];

            if (is_interface_exist)
            {
                double u_frac = fabs(phi[] + 1.0e-16)/(fabs(phi[]) + fabs(phi[-1]) + 1.0e-16);
                dis_c2i.x[] = color_cell[] * i_region > 0.0 ? u_frac : 1.0 - u_frac;

                if(dis_c2i.x[] <= 0.0)
                {
                    printf("ERROR IN GETDISTANCEC2I");
                    exit(1);
                }
            }
        }
    }

    boundary((scalar *){dis_c2i});

    return;
}

void getDis2fromDis1(const face vector dis1, face vector dis2)
{
    foreach_face()
    {
        dis2.x[] = 0.0;
        if (dis1.x[] > 0.0)
        {
            dis2.x[] = 1.0 - dis1.x[];
        }
    }

    boundary((scalar *){dis2});
}

//the interfacial region where the ghost velocites are set to couple two phases
void tagInterfacialRegion(scalar is_interfacial)
{

    //find out interfacial cells
    foreach ()
    {
        is_interfacial[] = 0.0;
#if USE_MY_SOLID
        if ((int)is_solid[] == 1)
            continue;
#endif
        bool is_narrow = false;
        foreach_neighbor(1)
        {
            int conf = (int) config_dict[];
            if (conf > 0)
            {
                is_narrow = true;
            }
        }
        is_interfacial[] = is_narrow ? 1.0 : 0.0;
    }

    boundary({is_interfacial});
}

void getNormalbyDistance(const coord * ifsg /*interface segments*/, const int num_seg,
                         const double color_center, const coord cc /*cell center*/, coord * np)
{
    double distance = 1.0e10;
    coord p_final = {0.0};
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
            p_final = p_closest;
        }
    }

    //p_final -> cc
    *np = vecdiff(cc, p_final);

    foreach_dimension()
    {
        np->x = -np->x * color_center / distance;
    }
    
#if DEBUGTL
    double test = np->x * np->x + np->y * np->y;
    if(fabs(test - 1.0) > 1.0e-10)
    {
        printf("ERROR IN COMPUTING NORMAL!\n");
        exit(1);
    }
#endif
    return;
}

void getNormalinInterfacialRegion(const scalar color_cc, const scalar is_interfacial, vector np)
{
    //calculate vector in interfacial region
    foreach()
    {
        foreach_dimension()
            np.x[] = 0.0;

        coord n_temp = {0.0, 0.0, 0.0};
#if USE_MY_SOLID//It's not as good as embed boundary method,
        if ((int)is_solid[] == 1)
            continue;
#endif
        if(is_interfacial[])
        {
            //the signed distance is used to calculate the normal vector, which wll be smoother
            foreach_dimension()
            {
                n_temp.x = -(phi_dis[1] - phi_dis[-1] + 1.0e-16) / Delta;
            }

            normalize(&n_temp);
        }

        foreach_dimension()
        {
            np.x[] = n_temp.x;
        }
    }

    boundary((scalar *){np});
    return;
}

//Obtain the stencil for the finite difference scheme.
//For the cell which is located another phase, its ghost value is determined by extratpolation
//The temperature condition at the interface is taken into account.
foreach_dimension()
void getTempStencil_x(Point point, const double i_region, const scalar color_cc, const face vector dis_c2i, 
                        const scalar T, double Tst[5])
{
    const double epsilon = Delta;
    double Tc = T[];
    double Tm = T[-1];
    double Tmm = T[-2];
    double Tp = T[1];
    double Tpp = T[2];

    double cc = color_cc[] * i_region;
    double cm = color_cc[-1] * i_region;
    double cmm = color_cc[-2] * i_region;
    double cp = color_cc[1] * i_region;
    double cpp = color_cc[2] * i_region;

    if (cm < 0.0)
    {
        double dis = dis_c2i.x[];
#if DEBUGTL
        if (dis == 0.0)
        {
            printf("WRONG DISTANCE!\n");
            exit(1);
        }
#endif
        if (dis < epsilon)
        {
            Tm = cp > 0.0 ? 2.0 * Tc - Tp : Tsat;
        }
        else
        {
            Tm = Tc - (Tc - Tsat) / dis;
        }
    }

    if (cp < 0.0)
    {
        double dis = dis_c2i.x[1];
#if DEBUGTL
        if (dis == 0.0)
        {
            printf("WRONG DISTANCE!\n");
            exit(1);
        }
#endif
        if (dis < epsilon)
        {
            Tp = cm > 0.0 ? 2.0 * Tc - Tm : Tsat;
        }
        else
        {
            Tp = Tc - (Tc - Tsat) / dis;
        }
    }

    if (cmm < 0.0)
    {
        Tmm = 2.0 * Tm - Tc;
    }

    if (cpp < 0.0)
    {
        Tpp = 2.0 * Tp - Tc;
    }

    Tst[0] = Tmm;
    Tst[1] = Tm;
    Tst[2] = Tc;
    Tst[3] = Tp;
    Tst[4] = Tpp;
}

//Extrapolate T to the cell face
foreach_dimension()
void getTempStencilFace_x(Point point, const double i_region, const scalar color_cc, const face vector dis_c2i, 
                        const scalar T, double Tst[5])
{
    const double epsilon = Delta;
    double Tc = T[];
    double Tm = T[-1];
    double Tmm = T[-2];
    double Tp = T[1];

    double cc = color_cc[] * i_region;
    double cm = color_cc[-1] * i_region;
    double cmm = color_cc[-2] * i_region;
    double cp = color_cc[1] * i_region;

    if (cm < 0.0)
    {
        double dis = dis_c2i.x[];
#if DEBUGTL
        if (dis == 0.0)
        {
            printf("WRONG DISTANCE!\n");
            exit(1);
        }
#endif
        if (dis < epsilon)
        {
            Tm = cp > 0.0 ? 2.0 * Tc - Tp : Tsat;
        }
        else
        {
            Tm = Tc - (Tc - Tsat) / dis;
        }
    }

    if (cp < 0.0)
    {
        double dis = dis_c2i.x[1];
#if DEBUGTL
        if (dis == 0.0)
        {
            printf("WRONG DISTANCE!\n");
            exit(1);
        }
#endif
        if (dis < epsilon)
        {
            Tp = cm > 0.0 ? 2.0 * Tc - Tm : Tsat;
        }
        else
        {
            Tp = Tc - (Tc - Tsat) / dis;
        }
    }

    if (cmm < 0.0)
    {
        Tmm = 2.0 * Tm - Tc;
    }

    Tst[0] = Tmm;
    Tst[1] = Tm;
    Tst[2] = Tc;
    Tst[3] = Tp;
}

//Compute the normal derivative with WENO
void calNormalDerivative(const double i_region, const scalar color_cc, const vector np, 
                        const scalar is_interf, const face vector dis_c2i, const scalar T, 
                        scalar dTdn)
{
    //1.0 for fluid1 and -1.0 for fluid2
    //as the default normal is from fluid 1 (liquid) to fluid 2 (vapor)
    const double sign = i_region > 0.0 ? 1.0 : -1.0;
    foreach()
    {
        const double epsilon = Delta;
        dTdn[] = 0.0;
#if USE_MY_SOLID
        if((int)is_solid[] == 1)
            continue;
#endif
        bool is_in_region = color_cc[] * i_region > 0.0;
        int idim = 0; //debug reason
        if(is_interf[] && is_in_region)
        {
            idim++;
            double dT[4] = {0.0};
            double Tst[5] = {0.0};
            foreach_dimension()
            {
                double dTc = 0.0;
                double cc = color_cc[] * i_region;
                double cm = color_cc[-1] * i_region;
                double cp = color_cc[1] * i_region;

                if(cm < 0.0 && cp < 0.0)
                {
                    //do nothing
                }
                else if(cm > 0.0 && cp > 0.0)
                {
                    //standard central difference
                    dTc = (T[1] - T[-1])/2.0;
                }
                else if(cm < 0.0)
                {
                    getTempStencil_x(point, i_region, color_cc, dis_c2i, T, Tst);
                    for(int ist = 0; ist < 4; ist ++)
                    {
                        dT[ist] = Tst[ist + 1] - Tst[ist];
                    }
                    dTc = weno3M(&dT[1]);
                }
                else
                {
                    getTempStencil_x(point, i_region, color_cc, dis_c2i, T, Tst);
                    for(int ist = 0; ist < 4; ist ++)
                    {
                        dT[ist] = Tst[ist + 1] - Tst[ist];
                    }
                    dTc = weno3P(&dT[1]);
                }

                dTdn[] += sign * dTc / Delta * np.x[];
            }
        }
    }
    boundary({dTdn});
}

//extrapolation via solvig PDE
void extendVariables(const double i_region, const scalar color_cc, const vector np, 
                     const scalar is_interf, const bool fix_cut, scalar var_e)
{
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
        int conf = (int)config_dict[];
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
            int conf = (int) config_dict[];
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
//Extrapolation via solvig PDE
//This can make sure the 1st order derivative is continuious
void extendVariables1stDerivative(const double i_region, const scalar color_cc, const vector np, 
                                  const scalar is_interf, scalar var_e, scalar dvardn)
{
    boundary({var_e});
    const int num_iter = 10;

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

                dvar[] += dvardn[];
                dvar[] *= dtau;
            }
        }

        // add
        foreach ()
        {
            if(is_interf[])
            {
                var_e[] += dvar[];
            }
        }

        boundary({var_e});
    }
}

//Compute the mass fluxes and extending it into the whole narrow band
void calMassFluxes(const scalar color_cc, const vector np, const scalar is_interf, 
                   const scalar dTdnL, const scalar dTdnG, scalar mdot)
{
    foreach()
    {
        //cut cell
        //it's better to give an initial value
        mdot[] = 0.0;
#if USE_MY_SOLID
        if((int)is_solid[] == 1)
            continue;
#endif
        if(is_interfacial[])
        {
          mdot[] = (lambda1 * dTdnL[] + lambda2 * dTdnG[]) / dhev;
        }
    }
    boundary({mdot});
    //extending
    extendVariables(1.0, color_cc, np, is_interf, true, mdot);
    extendVariables(-1.0, color_cc, np, is_interf, true, mdot);
}

//The advection part of the energy equation
void getIncreTempAdv(const double i_region, const scalar color_cc, const face vector dis_c2i, 
                     const vector uadv, const scalar T, scalar dvar)
{
    foreach()
    {
        const double epsilon = Delta;
        dvar[] = 0.0;
#if USE_MY_SOLID
        if ((int)is_solid[] == 1)
            continue;
#endif
        double cc = color_cc[] * i_region;
        if(cc > 0.0)
        {
            double dT[4] = {0.0};
            foreach_dimension()
            {
                double vel = uadv.x[];
                double dTc = 0.0;

                double Tst[5] = {0.0};
                double dT[4] = {0.0};
                getTempStencil_x(point, i_region, color_cc, dis_c2i, T, Tst);
                for(int ist = 0; ist < 4; ist ++)
                {
                    dT[ist] = Tst[ist + 1] - Tst[ist];
                }
                // we can use different scheme, here we use 2nd upwind
                if (vel > 0.0)
                {
                    //dTc = vel / 2.0 / Delta * (3.0 * T[] - 4.0 * T[-1] + T[-2]);
                    dTc = weno3P(&dT[1]);
                    dTc = dTc / Delta * vel;
                }
                else
                {
                    //dTc = vel / 2.0 / Delta * (-3.0 * T[] + 4.0 * T[1] - T[2]);
                    dTc = weno3M(&dT[1]);
                    dTc = dTc / Delta * vel;
                }

                dvar[] += dTc;
            }
        }
    }
}

//solve the advection part.
void advectTemp(const vector uadv, const scalar color_cc, const double dt)
{
    //ADV_SCHEME 1: interpolate T at the cell face
    //ADV_SCHEME 2: interpolate \Delta T at the cell center
#if (ADV_SCHEME == 1)
    boundary({TL,TG});
    face vector TL_f[];
    foreach_face()
    {
        TL_f.x[] = 0.0;
#if USE_MY_SOLID
        if((int)is_solid_face.x[] == 1)
            continue;
#endif
        if(color_cc[] > 0.0 || color_cc[-1] > 0.0)
        {
            double Tst[5] = {0.0};
            getTempStencilFace_x(point, 1.0, color_cc, disL, TL, Tst);
            TL_f.x[] = uf.x[] > 0.0 ? weno3P(&Tst[1]) : weno3M(&Tst[1]);
        }
    }

    boundary({TL_f});

    foreach()
    {
        dvar[] = 0.0;
#if USE_MY_SOLID
        if ((int)is_solid[] == 1)
            continue;the 
#endif
        if(color_cc[] > 0.0)
        {
            double dT = 0.0;
            foreach_dimension()
                dT += uadv.x[] * (TL_f.x[]  - TL_f.x[1]) / (Delta);
            dvar[] = dT;
        }
    }

    foreach()
    {
#if USE_MY_SOLID
        if ((int)is_solid[] == 1)
            continue;
#endif
        if(color_cc[] > 0.0)
        {
            TL[] += dvar[] * dt;
        }
        else
        {
            //initial conditions for the exteding. To speedup the convergence
            TL[] = Tsat;
        }
    }
    extendVariables1stDerivative(1.0, color_cc, n_pc, is_interfacial, TL, dTdnL);

    face vector TG_f[];
    foreach_face()
    {
        TG_f.x[] = 0.0;
#if USE_MY_SOLID
        if((int)is_solid_face.x[] == 1)
            continue;
#endif
        if(color_cc[] < 0.0 || color_cc[-1] < 0.0)
        {
            double Tst[5] = {0.0};
            getTempStencilFace_x(point, -1.0, color_cc, disG, TG, Tst);
            TG_f.x[] = uf2.x[] > 0.0 ? weno3P(&Tst[1]) : weno3M(&Tst[1]);
        }
    }

    boundary({TG_f});

    foreach()
    {
        dvar[] = 0.0;
#if USE_MY_SOLID
        if ((int)is_solid[] == 1)
            continue;
#endif
        if(color_cc[] < 0.0)
        {
            double dT = 0.0;
            foreach_dimension()
                dT += uadv.x[] * (TG_f.x[]  - TG_f.x[1]) / (Delta);
            dvar[] = dT;
        }
    }

    foreach ()
    {
#if USE_MY_SOLID
        if ((int)is_solid[] == 1)
            continue;
#endif
        if (color_cc[] < 0.0)
        {
            TG[] += dvar[] * dt;
        }
        else
        {
            // initial conditions for the exteding. To speedup the convergence
            TG[] = Tsat;
        }
    }
    extendVariables1stDerivative(-1.0, color_cc, n_pc, is_interfacial, TG, dTdnG);
    boundary({TL,TG});
#else
    boundary({TL,TG});
    //get incremental
    getIncreTempAdv(1.0, color_cc, disL, uadv, TL, dvar);
    //update
    foreach()
    {
#if USE_MY_SOLID
        if ((int)is_solid[] == 1)
            continue;
#endif
        if(color_cc[] > 0.0)
        {
            TL[] -= dvar[] * dt;
        }
        else
        {
            //initial conditions for the exteding. To speedup the convergence
            TL[] = Tsat;
        }
    }
    extendVariables1stDerivative(1.0, color_cc, n_pc, is_interfacial, TL, dTdnL);

    //get incremental
    getIncreTempAdv(-1.0, color_cc, disG, uadv, TG, dvar);
    //update
    foreach()
    {
#if USE_MY_SOLID
        if ((int)is_solid[] == 1)
            continue;
#endif
        if(color_cc[] < 0.0)
        {
            TG[] -= dvar[] * dt;
        }
        else
        {
            //initial conditions for the exteding. To speedup the convergence
            TG[] = Tsat;
        }
    }
    extendVariables1stDerivative(-1.0, color_cc, n_pc, is_interfacial, TG, dTdnG);

    boundary({TL,TG});
#endif
}

//for cells whose centers are very close to the interface
//we set their tempeartures to saturation temperature during the solution of heat diffusion
void getInterfaceCell(const double i_region, const scalar color_cell, const face vector dis_c2i, 
                     scalar interface_cell)
{
    foreach()
    {
        const double EPS = 1.0e-6; //for real distance < dx*dx
        interface_cell[] = 0.0;
#if USE_MY_SOLID
        if ((int)is_solid[] == 1)
            continue;
#endif
        int conf = (int)config_dict[];
        double color_c = color_cell[];
        bool is_in_region = color_c * i_region > 0.0;

        if(is_in_region && conf > 0) //pay attention to the cells in the current region, 
        {
            //see if the segments is very close to the cell center
            //left
            bool is_left_close = dis_c2i.x[] > 0.0 && dis_c2i.x[] < EPS;
            
            //right
            bool is_right_close = dis_c2i.x[1,0] > 0.0 && dis_c2i.x[1,0] < EPS;

            //bottom
            bool is_bottom_close = dis_c2i.y[] > 0.0 && dis_c2i.y[] < EPS;

            //top
            bool is_top_close = dis_c2i.y[0,1] > 0.0 && dis_c2i.y[0,1] < EPS;

            bool is_close = is_left_close || is_right_close || is_bottom_close || is_top_close;
            
            if(is_close) 
            {
                interface_cell[] = 1.0;
            }
        }   
    }

    boundary({interface_cell});
}

void setPoissonTempDiffusionCoef(const double i_region, const scalar color_cell, const face vector dis_c2i, 
                                 const scalar is_interface_cell, const double dt, 
                                 face vector alpha, scalar lambda, scalar source)
{
    double coef_tm1 = dt * lambda1 / rho1 / cp1; 
    double coef_tm2 = dt * lambda2 / rho2 / cp2;
    const double coef_tm = i_region > 0.0 ? coef_tm1 : coef_tm2;
    //obtain alpha
    foreach_face()
    {
        bool is_neg_p = color_cell[] * i_region  < 0.0;
        bool is_neg_m = color_cell[-1] * i_region < 0.0;
        bool is_inter_p = is_interface_cell[] == 1.0;
        bool is_inter_m = is_interface_cell[-1] == 1.0;
        bool is_const = is_neg_p || is_neg_m || is_inter_p || is_inter_m;

        if(is_const)
        {
            alpha.x[] = 0.0;
        }
        else
        {
            alpha.x[] = coef_tm * fm.x[];
        }
    }

    //set lambda and source
    foreach()
    {
        if(color_cell[] * i_region < 0.0)
        {
            lambda[] = 1.0;
            source[] = Tsat;
        }
        else if(is_interface_cell[] == 1)
        {
            lambda[] = 1.0;
            source[] = Tsat;
        }
        else
        {
            lambda[] = -1.0 * cm[];
            source[] *= cm[];
            // left
            if (alpha.x[] == 0)
            {
                //when alpha.x[] = 0 and dis_c2i.x[] = 0, the neighbour cell is a interface cell
                double frac = dis_c2i.x[] > 0.0 ? dis_c2i.x[] : 1.0;
                double xm = x - frac * Delta;
                double ym = y;
                double T_sat = Tsat;

                double coef = coef_tm / frac / Delta / Delta * fm.x[];

                source[] -= T_sat * coef;
                lambda[] -= coef;
            }

            // right
            if (alpha.x[1, 0] == 0)
            {
                double frac = dis_c2i.x[1, 0] > 0.0 ? dis_c2i.x[1, 0] : 1.0;

                double xm = x + frac * Delta;
                double ym = y;
                double T_sat = Tsat;

                double coef = coef_tm / frac / Delta / Delta * fm.x[1,0];

                source[] -= T_sat * coef;
                lambda[] -= coef;
            }

            // bottom
            if (alpha.y[] == 0)
            {
                double frac = dis_c2i.y[] > 0.0 ? dis_c2i.y[] : 1.0;
                double xm = x;
                double ym = y - frac * Delta;
                double T_sat = Tsat;

                double coef = coef_tm / frac / Delta / Delta * fm.y[];

                source[] -= T_sat * coef;
                lambda[] -= coef;
            }

            // top
            if (alpha.y[0, 1] == 0)
            {
                double frac = dis_c2i.y[0, 1] > 0.0 ? dis_c2i.y[0, 1] : 1.0;
                double xm = x;
                double ym = y + frac * Delta;
                double T_sat = Tsat;

                double coef = coef_tm/frac/Delta/Delta * fm.y[0,1];

                source[] -= T_sat * coef;
                lambda[] -= coef;
            }
        }
    }

    boundary((scalar *){alpha});
    boundary({lambda, source});
}
// When a solid boundary is used, we always ensure that several layers of cells near it are located at the finest resolution.
// This function is used to create a field which can be accordingly used for mesh adaptation.
void setJumpVarAdaptation(scalar var, int num)
{
  if (num == 0)
    return;
  scalar new[];
  foreach ()
  {
    new[] = var[];
  }

  for (int ii = 0; ii < num * 2; ++ii)
  {
    foreach ()
    {
      if (var[] == 0. || var[] == 1.)
      {
        bool is_neigh = false;
        foreach_neighbor(1)
        {
          bool is_finest = point.level == grid->maxdepth;
          const double ferror = 1.0e-10;
          if (fabs(var[]) > ferror && fabs(var[]) < (1. - ferror))
          {
            is_neigh = true;
          }
          // if (fabs(var[]) > 0. && fabs(var[]) < 1. && is_finest)
          // {
          //   is_neigh = true;
          // }
        }
        if (is_neigh)
        {
          double rnumb = 0.7 * pow(-1, ii);
          new[] = rnumb;
        }
      }
    }
    foreach ()
    {
      var[] = new[];
    }
  }
}
//Refine grids near the solid boundries and the interface
void fillRefineVOFs(scalar f_refine, scalar solid_refine, scalar solid_refiney)
{
  extern const int num_refine;
#if USE_MY_SOLID
  extern scalar is_solid;
  extern scalar is_solid_y;
#endif
  foreach ()
  {
    f_refine[] = f[];
    solid_refine[] = 0.0;
    solid_refiney[] = 0.0;
#if USE_MY_SOLID
    solid_refine[] = is_solid[];
    solid_refiney[] = is_solid_y[];
#endif
  }

#if USE_MY_SOLID
  foreach ()
  {
    if (is_solid[] == 1.0 && is_solid[1] == 0.0)
    {
      solid_refine[] = 0.5;
    }
    if (is_solid[] == 0.0 && is_solid[-1] == 1.0)
    {
      solid_refine[] = 0.5;
    }
    if (is_solid_y[] == 0.0 && is_solid_y[0, -1] == 1.0)
    {
      solid_refiney[] = 0.5;
    }
    if (is_solid_y[] == 1.0 && is_solid_y[0, 1] == 0.0)
    {
      solid_refiney[] = 0.5;
    }
  }
#endif
  setJumpVarAdaptation(solid_refine, num_refine);
  setJumpVarAdaptation(solid_refiney, num_refine);
  setJumpVarAdaptation(f_refine, num_refine);
#if USE_MY_SOLID
  boundarySolidNeummanNoauto(f_refine);
#endif
}
#endif