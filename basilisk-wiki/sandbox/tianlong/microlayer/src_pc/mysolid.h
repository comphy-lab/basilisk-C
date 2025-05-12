#ifndef _MY_SOLID_H
#define _MY_SOLID_H

const bool is_slip_x = false; //slip wall or no-slip wall. 
const bool is_slip_y = true;

const double coef_slip_x = is_slip_x ? 1.0 : -1.0;
const double coef_slip_y = is_slip_y ? 1.0 : -1.0; 

double SOLID_LEN_x = 0.0;
double SOLID_LEN_y = 0.0;

bool IS_SOLID_x = false;
bool IS_SOLID_y = false;

#if USE_CONJUGATE_HEAT
double Rcc = 0.0;
double delta_q_zero(double x, double t) { return 0.0; }
double (*delta_q)(double x, double t) = delta_q_zero;
double (*delta_q_volume)(double x, double t) = delta_q_zero;
#endif //USE_CONJUGATE_HEAT

scalar is_solid[];
scalar is_solid_x[];
scalar is_solid_y[];
face vector is_solid_face[];
face vector is_solid_x_face[];
face vector is_solid_y_face[];

attribute {
  void (* boundarySolid_x) (scalar s, int dir);
  void (* boundarySolid_y) (scalar s, int dir);
}

//to use symmetry BC with solid
# define neumann_pressure_solid(i) (alpha.n[i] == 0.0 ? 0.0 : a.n[i]*fm.n[i]/alpha.n[i])

static inline void restriction_solid (Point point, scalar s)
{  
  s[] = (x <= 0.0 || y <= 0.0) ? 1.0 : 0.0;
}

static inline void restriction_xsolid (Point point, scalar s)
{  
  s[] = x <= 0.0 ? 1.0 : 0.0;
}

static inline void restriction_ysolid (Point point, scalar s)
{  
  s[] = y <= 0.0 ? 1.0 : 0.0;
}

static inline void refine_solid (Point point, scalar s)
{  
  foreach_child()
  {
    s[] = (x <= 0.0 || y <= 0.0) ? 1.0 : 0.0;
  }
}

static inline void refine_xsolid (Point point, scalar s)
{  
  foreach_child()
  {
    s[] = x <= 0.0 ? 1.0 : 0.0;
  }
}

static inline void refine_ysolid (Point point, scalar s)
{  
  foreach_child()
  {
    s[] = y <= 0.0 ? 1.0 : 0.0;
  }
}

static inline void restriction_solidface (Point point, scalar s)
{
    //fixme: 3D is not considered for now
    vector v = s.v;
    int idim = 0;
    foreach_dimension()
    {
      idim++;
#if dimension == 2
      double ic = idim == 1 ? x : y;
      double im = ic - Delta;
      double ip = ic + Delta;
      double cr = idim == 1 ? y : x;
      v.x[] = (im <= 0.0 || ic <= 0.0 || cr <= 0.0) ? 1.0 : 0.0;
      v.x[1] = (ic <= 0.0 || ip <= 0.0 || cr <= 0.0) ? 1.0 : 0.0;
#endif
    }
}

static inline void restriction_solidxface (Point point, scalar s)
{
    //fixme: 3D is not considered for now
    vector v = s.v;
    int idim = 0;
    foreach_dimension()
    {
      idim++;
#if dimension == 2
      double ic = x;
      double dd = idim == 1 ? Delta : 0.0;
      double im = ic - dd;
      double ip = ic + dd;
      v.x[] = (im <= 0.0 || ic <= 0.0) ? 1.0 : 0.0;
      v.x[1] = (ic <= 0.0 || ip <= 0.0) ? 1.0 : 0.0;
#endif
    }
}

static inline void restriction_solidyface (Point point, scalar s)
{
    //fixme: 3D is not considered for now
    vector v = s.v;
    int idim = 0;
    foreach_dimension()
    {
      idim++;
#if dimension == 2
      double ic = y;
      double dd = idim == 2 ? Delta : 0.0;
      double im = ic - dd;
      double ip = ic + dd;
      v.x[] = (im <= 0.0 || ic <= 0.0) ? 1.0 : 0.0;
      v.x[1] = (ic <= 0.0 || ip <= 0.0) ? 1.0 : 0.0;
#endif
    }
}

static inline void refine_solidface (Point point, scalar s)
{
    //fixme: 3D is not considered for now
    vector v = s.v;
    
    int idim = 0;
    foreach_dimension()
    {
      idim++;
      double ic = idim == 1 ? x : y;
      double dd = Delta / 2.0;
      ic = ic - dd / 2.0;
      double im = ic - dd;
      double ip = ic + dd;
      double ipp = ic + 2.0 * dd;

      double cr = idim == 1 ? y : x;
      cr = cr - dd / 2.0;
      double crp = cr + dd;

      if (!is_refined(neighbor(-1)) && (is_local(cell) || is_local(neighbor(-1))))
      {
        fine(v.x, 0, 0) = (im <= 0.0 || ic <= 0.0 || cr <= 0.0) ? 1.0 : 0.0;
        fine(v.x, 0, 1) = (im <= 0.0 || ic <= 0.0 || crp <= 0.0) ? 1.0 : 0.0;
      }

      if (is_local(cell))
      {
        fine(v.x, 1, 0) = (ic <= 0.0 || ip <= 0.0 || cr <= 0.0) ? 1.0 : 0.0;
        fine(v.x, 1, 1) = (ic <= 0.0 || ip <= 0.0 || crp <= 0.0) ? 1.0 : 0.0;
      }

      if (!is_refined(neighbor(1)) && neighbor(1).neighbors && (is_local(cell) || is_local(neighbor(1))))
      {
        fine(v.x, 2, 0) = (ip <= 0.0 || ipp <= 0.0 || cr <= 0.0) ? 1.0 : 0.0;
        fine(v.x, 2, 1) = (ip <= 0.0 || ipp <= 0.0 || crp <= 0.0) ? 1.0 : 0.0;
      }
    }
}

static inline void refine_solidxface (Point point, scalar s)
{
    //fixme: 3D is not considered for now
    vector v = s.v;
    
    int idim = 0;
    foreach_dimension()
    {
      idim++;

      //TODO: simplify
      double ic = x;
      ic = ic - Delta / 4.0;
      double dd = idim == 1 ? Delta / 2.0 : 0.0;
      double im = ic - dd;
      double ip = ic + dd;
      double ipp = ic + 2.0 * dd;

      double cr = ic;
      double dd_i = idim == 2 ? Delta / 2.0 : 0.0;
      double crp = cr + dd_i;
      
      if(idim == 1)
      {
        cr = crp = 1.0;
      }
      else if(idim == 2)
      {
        im = ic = ip = ipp = 1.0;
      }

      if (!is_refined(neighbor(-1)) && (is_local(cell) || is_local(neighbor(-1))))
      {
        fine(v.x, 0, 0) = (im <= 0.0 || ic <= 0.0 || cr <= 0.0) ? 1.0 : 0.0;
        fine(v.x, 0, 1) = (im <= 0.0 || ic <= 0.0 || crp <= 0.0) ? 1.0 : 0.0;
      }

      if (is_local(cell))
      {
        fine(v.x, 1, 0) = (ic <= 0.0 || ip <= 0.0 || cr <= 0.0) ? 1.0 : 0.0;
        fine(v.x, 1, 1) = (ic <= 0.0 || ip <= 0.0 || crp <= 0.0) ? 1.0 : 0.0;
      }

      if (!is_refined(neighbor(1)) && neighbor(1).neighbors && (is_local(cell) || is_local(neighbor(1))))
      {
        fine(v.x, 2, 0) = (ip <= 0.0 || ipp <= 0.0 || cr <= 0.0) ? 1.0 : 0.0;
        fine(v.x, 2, 1) = (ip <= 0.0 || ipp <= 0.0 || crp <= 0.0) ? 1.0 : 0.0;
      }

    }
}

static inline void refine_solidyface (Point point, scalar s)
{
    //fixme: 3D is not considered for now
    vector v = s.v;
    
    int idim = 0;
    foreach_dimension()
    {
      idim++;

      //TODO: simplify
      double ic = y;
      ic = ic - Delta / 4.0;
      double dd = idim == 2 ? Delta / 2.0 : 0.0;
      double im = ic - dd;
      double ip = ic + dd;
      double ipp = ic + 2.0 * dd;

      double cr = ic;
      double dd_i = idim == 1 ? Delta / 2.0 : 0.0;
      double crp = cr + dd_i;
      
      if(idim == 2)
      {
        cr = crp = 1.0;
      }
      else if(idim == 1)
      {
        im = ic = ip = ipp = 1.0;
      }

      if (!is_refined(neighbor(-1)) && (is_local(cell) || is_local(neighbor(-1))))
      {
        fine(v.x, 0, 0) = (im <= 0.0 || ic <= 0.0 || cr <= 0.0) ? 1.0 : 0.0;
        fine(v.x, 0, 1) = (im <= 0.0 || ic <= 0.0 || crp <= 0.0) ? 1.0 : 0.0;
      }

      if (is_local(cell))
      {
        fine(v.x, 1, 0) = (ic <= 0.0 || ip <= 0.0 || cr <= 0.0) ? 1.0 : 0.0;
        fine(v.x, 1, 1) = (ic <= 0.0 || ip <= 0.0 || crp <= 0.0) ? 1.0 : 0.0;
      }

      if (!is_refined(neighbor(1)) && neighbor(1).neighbors && (is_local(cell) || is_local(neighbor(1))))
      {
        fine(v.x, 2, 0) = (ip <= 0.0 || ipp <= 0.0 || cr <= 0.0) ? 1.0 : 0.0;
        fine(v.x, 2, 1) = (ip <= 0.0 || ipp <= 0.0 || crp <= 0.0) ? 1.0 : 0.0;
      }

    }
}

void setSolidFlag()
{

  foreach_dimension()
  {
    IS_SOLID_x = SOLID_LEN_x > 0.0;
  }

  int dep = depth();
  for (int il = 0; il <= dep; il++)
  {
    foreach_level(il)
    {
      is_solid[] = 0.0;
      is_solid_x[] = 0.0;
      is_solid_y[] = 0.0;
      
      //due to the AXI implementation in Basilisk
      //the solid will always be puted on the left side or the bottom side.
      if(IS_SOLID_x)
      {
        is_solid_x[] = x <= 0.0 ? 1.0 : 0.0;
      }
      if(IS_SOLID_y)
      {
        is_solid_y[] = y <= 0.0 ? 1.0 : 0.0;
      }

      is_solid[] = (is_solid_x[] + is_solid_y[]) == 0.0 ? 0.0 : 1.0;
    }
  }

  is_solid_x.refine = is_solid_x.prolongation = refine_xsolid;
  is_solid_x.restriction = restriction_xsolid;
  is_solid_y.refine = is_solid_y.prolongation = refine_ysolid;
  is_solid_y.restriction = restriction_ysolid;
  is_solid.refine = is_solid.prolongation = refine_solid;
  is_solid.restriction = restriction_solid;

  is_solid_x.dirty = true;
  is_solid_y.dirty = true;
  is_solid.dirty = true;
  boundary({is_solid, is_solid_x, is_solid_y});

  //as long as this face is next to a solid cell, it's a solid face
  foreach_face()
  {
    bool is_solid_p = (int) is_solid[] == 1;
    bool is_solid_m = (int) is_solid[-1] == 1;

    is_solid_face.x[] = (is_solid_p || is_solid_m) ? 1.0 : 0.0; 
  }

  //to avoid the automatic substitution of foreach_face()
  foreach_face(x)
  {
    bool is_solid_xp = (int) is_solid_x[] == 1;
    bool is_solid_xm = (int) is_solid_x[-1] == 1;

    is_solid_x_face.x[] = (is_solid_xp || is_solid_xm) ? 1.0 : 0.0; 

    bool is_solid_yp = (int) is_solid_y[] == 1;
    bool is_solid_ym = (int) is_solid_y[-1] == 1;

    is_solid_y_face.x[] = (is_solid_yp || is_solid_ym) ? 1.0 : 0.0; 
  }

  foreach_face(y)
  {
    bool is_solid_xp = (int) is_solid_x[] == 1;
    bool is_solid_xm = (int) is_solid_x[0,-1] == 1;

    is_solid_x_face.y[] = (is_solid_xp || is_solid_xm) ? 1.0 : 0.0; 

    bool is_solid_yp = (int) is_solid_y[] == 1;
    bool is_solid_ym = (int) is_solid_y[0,-1] == 1;

    is_solid_y_face.y[] = (is_solid_yp || is_solid_ym) ? 1.0 : 0.0; 
  }

  foreach_dimension()
  {
    is_solid_face.x.restriction = is_solid_face.x.refine = is_solid_face.x.prolongation = no_restriction;
    is_solid_x_face.x.restriction = is_solid_x_face.x.refine = is_solid_x_face.x.prolongation = no_restriction;
    is_solid_y_face.x.restriction = is_solid_y_face.x.refine = is_solid_y_face.x.prolongation = no_restriction;
  }
  is_solid_face.x.refine = is_solid_face.x.prolongation = refine_solidface;
  is_solid_face.x.restriction = restriction_solidface;

  is_solid_x_face.x.refine = is_solid_x_face.x.prolongation = refine_solidxface;
  is_solid_x_face.x.restriction = restriction_solidxface;

  is_solid_y_face.x.refine = is_solid_y_face.x.prolongation = refine_solidyface;
  is_solid_y_face.x.restriction = restriction_solidyface;
  boundary({is_solid_face, is_solid_x_face, is_solid_y_face});
}

foreach_dimension()
double getFluidScalar_x(Point point, scalar s, int ip)
{
  return s[ip];
}

foreach_dimension()
void boundarySolidNeumman_x(scalar s, int dir)
{
  //this function will be called in boundary();
  //use noauto to avoid the recursion
  foreach(noauto)
  {
    int sdc = (int)is_solid_x[];
    int sdp = (int)is_solid_x[1];
    int sdpp = (int)is_solid_x[2];

    if(sdc == 1)
    {
      s[] = 0.0;
      bool is_finest = point.level == grid->maxdepth;
      if (is_finest)
      {
        if (sdp == 0)
        {
          s[] = s[1];
        }
        else if (sdpp == 0)
        {
          Point target = point;
          if (dir == 1)
          {
            target.i += 3;
          }
          else if (dir == 2)
          {
            target.j += 3;
          }
          s[] = getFluidScalar_x(target, s, 0);
        }
      }
    }
  }
}

//this function shouldn't be called automatically in boundary();
//explicitly use it to impose the boundary conditions
void boundarySolidNeummanNoauto (scalar s)
{
  s.dirty = true;
  boundary({s});

  int idim = 0;
  foreach_dimension()
  {
    idim ++; 
    if(IS_SOLID_x)
    {
      boundarySolidNeumman_x(s, idim);
    }

    s.dirty = true;
    boundary({s});
  }
}

foreach_dimension()
void boundarySolidVelC_x (scalar s, int dir)
{
  vector v = s.v;
  foreach(noauto)
  {
    int sdc = (int)is_solid_x[];
    int sdp = (int)is_solid_x[1];
    int sdpp = (int)is_solid_x[2];

    if(sdc == 1)
    {
      v.x[] = 0.0;
      v.y[] = 0.0;
      bool is_finest = point.level == grid->maxdepth;
      if (is_finest)
      {
        if (sdp == 0)
        {
          v.x[] = -v.x[1];
          v.y[] = coef_slip_x * v.y[1];
        }
        else if (sdpp == 0)
        {
          Point target = point;
          if (dir == 1)
          {
            target.i += 3;
          }
          else if (dir == 2)
          {
            target.j += 3;
          }
          v.x[] = -getFluidScalar_x(target, v.x, 0);
          v.y[] = coef_slip_x * getFluidScalar_x(target, v.y, 0);
        }
      }
    }
  }
}

//this function shouldn't be called automatically in boundary();
//explicitly use it to impose the boundary conditions
void boundarySolidVelCNoauto (vector v)
{
  v.x.dirty = true;
  boundary((scalar *){v});

  int idim = 0;
  foreach_dimension()
  {
    idim++;
    if(IS_SOLID_x)
    {
      boundarySolidVelC_x(v.x,idim);
    }

    v.x.dirty = true;
    boundary((scalar *){v});
  }
}

foreach_dimension()
void boundarySolidVelF_x (scalar s, int dir)
{
  vector vf = s.v;
  foreach_face(x, noauto)
  {
    if((int)is_solid_face.x[] == 1)
    {
      vf.x[] = 0.0;
    }
  }

  foreach_face(y, noauto)
  {
    int sdc = (int)is_solid_x[];
    int sdp = (int)is_solid_x[1];
    int sdpp = (int)is_solid_x[2];
    if(sdc == 1)
    {
      vf.y[] = 0.0;
      bool is_finest = point.level == grid->maxdepth;
      if (is_finest)
      {
        if (sdp == 0)
        {
          vf.y[] = coef_slip_x * vf.y[1];
#if AXI
          if (dir == 2) // y = 0 will be the axis
          {
            vf.y[] = -vf.y[1];
          }
#endif
        }
        else if (sdpp == 0)
        {
          Point target = point;
          if (dir == 1)
          {
            target.i += 3;
          }
          else if (dir == 2)
          {
            target.j += 3;
          }
          vf.y[] = coef_slip_x * getFluidScalar_x(target, vf.y, 0);
#if AXI
          if (dir == 2) // y = 0 will be the axis
          {
            vf.y[] = -getFluidScalar_x(target, vf.y, 0);
          }
#endif
        }
      }
    }
  }

  //I think the boundary conditions for the other direction is not important
}

//this function shouldn't be called automatically in boundary();
//explicitly use it to impose the boundary conditions
void boundarySolidVelFNoauto (vector vf)
{
  vf.x.dirty = true;
  boundary((scalar *){vf});

  int idim = 0;
  foreach_dimension()
  {
    idim++;
    if(IS_SOLID_x)
    {
      boundarySolidVelF_x(vf.x, idim);
    }

    vf.x.dirty = true;
    boundary((scalar *){vf});
  }
}

foreach_dimension()
void boundarySolidHeight_x(scalar s, int dir)
{
  vector h = s.v;
  //this function will be called in boundary();
  //use noauto to avoid the recursion
  foreach(noauto)
  {
    int sdc = (int)is_solid_x[];
    int sdp = (int)is_solid_x[1];
    int sdpp = (int)is_solid_x[2];

    if(sdc == 1)
    {
      h.x[] = 0.0;
      h.y[] = 0.0;
      bool is_finest = point.level == grid->maxdepth;
      if (is_finest)
      {
        if (sdp == 0)
        {
          double hx = h.x[1];
          double hy = h.y[1];

          h.x[] = hx == nodata ? nodata : -hx; // wrong curvature will be obtianed without minus
          h.y[] = hy == nodata ? nodata : hy;  // the height_prolongation function
          // we only impose contact line boundary conditions for the x walls
#if USE_CONTACT_ANGLE
          if (dir == 1)
          {
            extern double theta0;
            double exact_val = (10e-6 - y) / Delta + 20.0 ;
            //h.y[] = hy == nodata ? nodata : -hy + 2. * exact_val;
            h.y[] = hy == nodata ? nodata : hy + (orientation(hy) ? -1. : 1.) / tan(theta0 * pi / 180.);
          }
#endif
        }
        else if (sdpp == 0)
        {
          Point target = point;
          if (dir == 1)
          {
            target.i += 2; // not the real symmetric condition, to be exactly the same as the bc before
          }
          else if (dir == 2)
          {
            target.j += 2;
          }
          double hx = getFluidScalar_x(target, h.x, 0);
          double hy = getFluidScalar_x(target, h.y, 0);
          h.x[] = hx == nodata ? nodata : -hx;
          h.y[] = hy == nodata ? nodata : hy;
#if USE_CONTACT_ANGLE
          if (dir == 1)
          {
            extern double theta0;
            double exact_val = (10e-6 - y) / Delta + 20.0;
            //h.y[] = hy == nodata ? nodata : -hy + 2. * exact_val;
            h.y[] = hy == nodata ? nodata : hy + (orientation(hy) ? -1. : 1.) / tan(theta0 * pi / 180.);
          }
#endif
        }
      }
    }
  }
}

foreach_dimension()
void boundarySolidVectorNeumann_x(scalar s, int dir)
{
  vector v = s.v;
  //this function will be called in boundary();
  //use noauto to avoid the recursion
  foreach(noauto)
  {
    int sdc = (int)is_solid_x[];
    int sdp = (int)is_solid_x[1];
    int sdpp = (int)is_solid_x[2];

    if(sdc == 1)
    {
      v.x[] = 0.0;
      v.y[] = 0.0;
      bool is_finest = point.level == grid->maxdepth;
      if (is_finest)
      {
        if (sdp == 0)
        {
          v.x[] = v.x[1];
          v.y[] = v.y[1];
        }
        else if (sdpp == 0)
        {
          Point target = point;
          if (dir == 1)
          {
            target.i += 3;
          }
          else if (dir == 2)
          {
            target.j += 3;
          }
          v.x[] = getFluidScalar_x(target, v.x, 0);
          v.y[] = getFluidScalar_x(target, v.y, 0);
        }
      }
    }
  }
}

#if USE_CONJUGATE_HEAT
void imposeTempBoundaryConjugateLevel(scalar TL, scalar TG, scalar TS, int l)
{
  extern scalar f;
  foreach_level(l)
  {
    if(is_solid[] == 1.0)
    {
      TL[] = TS[];
      TG[] = TS[];
    }
    else
    {
      double T_liq = (fL[] > F_ERR) ? TL[] * (fL[]) : 0.;
      double T_gas = (fG[] > F_ERR) ? TG[] * (1.0 - fL[]): 0.;
      TS[] = T_liq + T_gas;
    }
  }

  TL.dirty = TG.dirty = TS.dirty = true;
  boundary_level({TL, TG, TS}, l);
}

#endif //USE_CONJUGATE_HEAT

#if USE_CONTACT_ANGLE
foreach_dimension()
void boundarySolidVOFSmallAngles_x(scalar s, int dir)
{
  extern double theta0;
  foreach(noauto)
  {
    int sdc = (int)is_solid_x[];
    int sdp = (int)is_solid_x[1];
    int sdpp = (int)is_solid_x[2];

    if(sdc == 1)
    {
      s[] = 0.0;
      bool is_finest = point.level == grid->maxdepth;
      if (is_finest && ((sdp == 0) || sdpp == 0))
      {
        double f = 0.0;
        double fm = 0.0;
        double fp = 0.0;
        int ilayer = 0;
        if (sdp == 0)
        {
          //first layer
          ilayer = 1;
          f = s[1];
          fm = s[1,-1];
          fp = s[1,1];
        }
        else if (sdpp == 0)
        {
          //second layer 
          ilayer = 2;
          f = s[2];
          fm = s[2,-1];
          fp = s[2,1];
        }

        bool is_contact = false;
        double fv = 0.0;
        // 1. to decide wether it is the cell containing a triple point
        if (f > 0.0 && f < 1.0)
        {
          if (fm <= 0.0)
          {
            is_contact = true;
          }

          if (fp <= 0.0)
          {
            is_contact = true;
          }
        }

        coord n = {0.0,0.0};
        double val = theta0 / 180.0 * pi;

        if(dir == 1)
        {
          if(fp < fm)
          {
            n.x = cos(val); 
            n.y = sin(val);
          }
          else
          {
            n.x = cos(val);
            n.y = -sin(val);
          }
        }
        else
        {
          if(fp < fm)
          {
            n.x = sin(val);
            n.y = cos(val);
          }
          else
          {
            n.x = sin(val);
            n.y = -cos(val);
          }
        }
        double alpha = plane_alpha (f, n);
        double coefx = dir == 1 ? 1.0 : 0.0;
        double coefy = 1.0 - coefx;
        coefx *= abs(ilayer);
        coefy *= abs(ilayer);

        coord lf = {-1.0 / 2.0 - coefx, -1.0 / 2.0 - coefy, 0.0};
        coord hr = {1.0 / 2.0 - coefx, 1.0 / 2.0 - coefy, 0.0};

        if (is_contact)
        {
          fv = rectangle_fraction(n, alpha, lf, hr);
        }
        else
        {
          fv = f > 0.0;
        }

        s[] = fv;
      }
    }
  }
}

#endif


#endif //_MY_SOLID_H
