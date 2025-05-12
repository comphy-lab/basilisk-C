static inline void coarsen_quadratic (Point point, scalar s)
{
  double sum = 0.;
#if dimension == 1
  foreach_child()
    sum += 9*s[]-s[child.x,0,0];
  s[] = sum/(16);
  
#elif dimension==2
  foreach_child()
    sum += 10*s[]-s[child.x,0,0]-s[0,child.y,0];
  s[] = sum/(32);
  
#else //dimension ==3 
  foreach_child()
    sum += 11*s[]-s[child.x,0,0]-s[0,child.y,0]-s[0,0,child.z];
  s[] = sum/(64);
#endif
}

static inline void refine_quadratic (Point point, scalar s)
{
  #if dimension == 1
  foreach_child()
    s[]=(30*coarse(s,0)+5*coarse(s,child.x)-3*coarse(s,-child.x))/32;
  #elif dimension == 2
    foreach_child()
    s[]=(296.*coarse(s,0,0) +5.*coarse(s,child.x,child.y) +
	 146.*(coarse(s,child.x,0) + coarse(s,0,child.y)) +
	 98.*(coarse(s,-child.x, 0) + coarse(s,0,-child.y)) -
	 61.*(coarse(s,child.x,-child.y) + coarse(s,-child.x,child.y))-
	 91.*(coarse(s,-child.x,-child.y)))/576.;
  #else //dimension == 3
    foreach_child()
      s[]=(412*coarse(s,0,0,0)+
	   262.*(coarse(s,child.x,0,0)+coarse(s,0,child.y,0)+coarse(s,0,0,child.z))+
	   121.*(coarse(s,child.x,child.y,0)+coarse(s,child.x,0,child.z)+
		 coarse(s,0,child.y,child.z))+
	   214.*(coarse(s,-child.x,0,0)+coarse(s,0,-child.y,0)+coarse(s,0,0,-child.z))-
	   11.*(coarse(s,child.x,child.y,child.z))+
	   55.*(coarse(s,child.x,-child.y,0)+coarse(s,-child.x,child.y,0)+
		coarse(s,child.x,0,-child.z)+coarse(s,-child.x,0,child.z)+
		coarse(s,0,child.y,-child.z)+coarse(s,0,-child.y,child.z))-
	   95.*(coarse(s,child.x,child.y,-child.z)+coarse(s,child.x,-child.y,child.z)+
		coarse(s,-child.x,child.y,child.z))+
	   25.*(coarse(s,-child.x,0,-child.z)+coarse(s,0,-child.y,-child.z)+
		coarse(s,-child.x,-child.y,0))-
	   143.*(coarse(s,-child.x,-child.y,child.z)+coarse(s,-child.x,child.y,-child.z)+
		 coarse(s,child.x,-child.y,-child.z))-
	   155.*coarse(s,-child.x,-child.y,-child.z))/1728;
    #endif
}

static inline void refine_3x3_linear(Point point, scalar s){
  double a = 8.*(s[1,0]+s[1,-1]-s[-1,-1]-s[-1,0]-s[-1,1]+s[1,1]);
  double b = 8.*(-s[1,-1]-s[0,-1]-s[-1,-1]+s[-1,1]+s[0,1]+s[1,1]);
  double c = 3.*(-s[1,-1]+s[-1,-1]-s[-1,1]+s[1,1]);
  double d = s[];
  foreach_child()
    s[]=d+(a*child.x+b*child.y+c*(child.x*child.y))/192.;
}


# if dimension == 3
static inline void refine_3x3x3_linear(Point point, scalar s){
  double a = 64.*(s[1,0,0]+s[1,-1,0]-s[-1,-1,0]-s[-1,0,0]-s[-1,1,0]+s[1,1,0]
		  +s[1,0,-1]+s[1,-1,-1]-s[-1,-1,-1]-s[-1,0,-1]-s[-1,1,-1]+s[1,1,-1]
		  +s[1,0,1]+s[1,-1,1]-s[-1,-1,1]-s[-1,0,1]-s[-1,1,1]+s[1,1,1]);
  double b = 64.*(-s[1,-1,0]-s[0,-1,0]-s[-1,-1,0]+s[-1,1,0]+s[0,1,0]+s[1,1,0]
		  -s[1,-1,-1]-s[0,-1,-1]-s[-1,-1,-1]+s[-1,1,-1]+s[0,1,-1]+s[1,1,-1]
		  -s[1,-1,1]-s[0,-1,1]-s[-1,-1,1]+s[-1,1,1]+s[0,1,1]+s[1,1,1]);
  double c = 64.*(-s[0,0,-1]-s[1,1,-1]-s[-1,-1,-1]-s[1,0,-1]-s[-1,0,-1]-s[0,1,-1]
		  -s[0,-1,-1]-s[1,-1,-1]-s[-1,1,-1]+s[0,0,1]+s[1,1,1]+s[-1,-1,1]
		  +s[1,0,1]+s[-1,0,1]+s[0,1,1]+s[0,-1,1]+s[1,-1,1]+s[-1,1,1]);
  double d = 24.*(-s[1,-1,0]+s[-1,-1,0]-s[-1,1,0]+s[1,1,0]-s[1,-1,-1]+s[-1,-1,-1]
		  -s[-1,1,-1]+s[1,1,-1]-s[1,-1,1]+s[-1,-1,1]-s[-1,1,1]+s[1,1,1]) ; 
  double e = 24.*(s[1,-1,-1]+s[0,-1,-1]+s[-1,-1,-1]-s[-1,1,-1]-s[0,1,-1]-s[1,1,-1]
		  -s[1,-1,1]-s[0,-1,1]-s[-1,-1,1]+s[-1,1,1]+s[0,1,1]+s[1,1,1]) ; 
  double f = 24.*(-s[1,0,-1]-s[1,-1,-1]+s[-1,-1,-1]+s[-1,0,-1]+s[-1,1,-1]-s[1,1,-1]
		  +s[1,0,1]+s[1,-1,1]-s[-1,-1,1]-s[-1,0,1]-s[-1,1,1]+s[1,1,1]) ; 
  double g = 9.*(s[1,-1,-1]-s[-1,-1,-1]+s[-1,1,-1]-s[1,1,-1]-s[1,-1,1]+s[-1,-1,1]-s[-1,1,1]+s[1,1,1]); 
  double h = s[];
  foreach_child()
    s[]=h+(a*child.x+b*child.y+c*child.z+d*(child.x*child.y)+e*(child.y*child.z)+f*(child.x*child.z)+g*(child.x*child.y*child.z))/4608.;
}
static inline void refine_quad(Point point, scalar s){
  double a = 3.*448.*(s[1,0,0]+s[1,-1,0]-s[-1,-1,0]-s[-1,0,0]-s[-1,1,0]+s[1,1,0]
		  +s[1,0,-1]+s[1,-1,-1]-s[-1,-1,-1]-s[-1,0,-1]-s[-1,1,-1]+s[1,1,-1]
		  +s[1,0,1]+s[1,-1,1]-s[-1,-1,1]-s[-1,0,1]-s[-1,1,1]+s[1,1,1]);
  double b = 3.*448.*(-s[1,-1,0]-s[0,-1,0]-s[-1,-1,0]+s[-1,1,0]+s[0,1,0]+s[1,1,0]
		  -s[1,-1,-1]-s[0,-1,-1]-s[-1,-1,-1]+s[-1,1,-1]+s[0,1,-1]+s[1,1,-1]
		  -s[1,-1,1]-s[0,-1,1]-s[-1,-1,1]+s[-1,1,1]+s[0,1,1]+s[1,1,1]);
  double c = 3.*448.*(-s[0,0,-1]-s[1,1,-1]-s[-1,-1,-1]-s[1,0,-1]-s[-1,0,-1]-s[0,1,-1]
		  -s[0,-1,-1]-s[1,-1,-1]-s[-1,1,-1]+s[0,0,1]+s[1,1,1]+s[-1,-1,1]
		  +s[1,0,1]+s[-1,0,1]+s[0,1,1]+s[0,-1,1]+s[1,-1,1]+s[-1,1,1]);
  double d = 3.*168.*(-s[1,-1,0]+s[-1,-1,0]-s[-1,1,0]+s[1,1,0]-s[1,-1,-1]+s[-1,-1,-1]
		  -s[-1,1,-1]+s[1,1,-1]-s[1,-1,1]+s[-1,-1,1]-s[-1,1,1]+s[1,1,1]) ; 
  double e = 3.*168.*(s[1,-1,-1]+s[0,-1,-1]+s[-1,-1,-1]-s[-1,1,-1]-s[0,1,-1]-s[1,1,-1]
		  -s[1,-1,1]-s[0,-1,1]-s[-1,-1,1]+s[-1,1,1]+s[0,1,1]+s[1,1,1]) ; 
  double f = 3.*168.*(-s[1,0,-1]-s[1,-1,-1]+s[-1,-1,-1]+s[-1,0,-1]+s[-1,1,-1]-s[1,1,-1]
		   +s[1,0,1]+s[1,-1,1]-s[-1,-1,1]-s[-1,0,1]-s[-1,1,1]+s[1,1,1]) ; 
  double g = 3.*63.*(s[1,-1,-1]-s[-1,-1,-1]+s[-1,1,-1]-s[1,1,-1]-s[1,-1,1]+s[-1,-1,1]-s[-1,1,1]+s[1,1,1]);
  double h = 3.*(240.*(s[1,0,0]+s[-1,0,0])
	      +144.*(s[1,-1,0]+s[-1,-1,0]+s[1,1,0]+s[1,0,-1]+s[-1,0,-1]+s[1,0,1]+s[-1,0,1])
	      -96.*(s[0,-1,0]+s[0,1,0]+s[0,0,-1]+s[0,0,1])
	      +48.*(s[1,-1,-1]+s[-1,-1,-1]+s[-1,1,-1]+s[1,1,-1]+s[1,-1,1]+s[-1,-1,1]+s[-1,1,1]+s[1,1,1])
	      -192.*(s[0,-1,-1]+s[0,1,-1]+s[0,-1,1]+s[0,1,1]));
  double i = 0.;
  double j = 0. ;
  double k = 3584.*(7.*s[]
		    +4.*(s[1,0,0]+s[-1,0,0]+s[0,1,0]+s[0,-1,0]+s[0,0,1]+s[0,0,-1])
		    +(s[1,1,0]+s[-1,1,0]+s[1,-1,0]+s[-1,-1,0]
		      +s[1,0,1]+s[-1,0,1]+s[1,0,-1]+s[-1,0,-1]
		      +s[0,1,1]+s[0,-1,1]+s[0,1,-1]+s[0,-1,-1])
		    -2.*(s[1,1,1]+s[-1,1,1]+s[1,-1,1]+s[1,1,-1]
			 +s[-1,-1,1]+s[-1,1,-1]+s[1,-1,-1]+s[-1,-1,-1]));
						    
   
  foreach_child()
    s[]=(k+a*child.x+b*child.y+c*child.z+d*(child.x*child.y)+e*(child.y*child.z)
	 +f*(child.x*child.z)+g*(child.x*child.y*child.z)+h*(child.x*child.x)+
	 i*(child.y*child.y)+j*(child.z*child.z))/(32256.*3);
}

# endif
static inline void coarsen1D_twoleaf_quad (Point point, scalar s){
  double sum = 0.;
  foreach_child()
    sum += 16*s[];
  sum+=-(s[-1]+s[1]);
  s[] = sum/(30);
}

static inline void coarsen1D_leftleaf_quad (Point point, scalar s){
  double sum = 0.;
  foreach_child(){
    if (child.x==-1)
      sum += 591.2*s[];
    else
      sum+=(((490.5*s[])+(-32.3*s[child.x])));
  }
  sum+=-49.4*(s[-1]);
  s[] = sum/(1000);
}

static inline void coarsen1D_rightleaf_quad (Point point, scalar s){
  double sum = 0.;
  foreach_child(){
    if (child.x==1)
      sum += 591.2*s[];
    else
      sum+=(((490.5*s[])+(-32.3*s[child.x])));
  }
  sum+=-49.4*(s[1]);
  s[] = sum/(1000);
  }

static inline void coarsen2D_oneleaf_quad(Point point,scalar s,int n)
{
  double w[6]= {0.321645313553608, 0.282670262980445, -0.010788941335132,
		-0.052461227242077, -0.013486176668914,-0.055158462575860};
  if (n == 1){
    s[]=(w[0]*(fine(s,0,0)+fine(s,0,1))+
	 w[1]*(fine(s,1,0)+fine(s,1,1))+
	 w[2]*(fine(s,2,1)+fine(s,2,0))+
	 w[3]*(fine(s,1,2)+fine(s,1,-1))+
	 w[4]*(fine(s,0,2)+fine(s,0,-1))+
	 w[5]*s[-1]);
  }
  else if (n == 2){
    s[]=(w[0]*(fine(s,0,0)+fine(s,1,0))+
	 w[1]*(fine(s,0,1)+fine(s,1,1))+
	 w[2]*(fine(s,1,2)+fine(s,0,2))+
	 w[3]*(fine(s,2,1)+fine(s,-1,1))+
	 w[4]*(fine(s,2,0)+fine(s,-1,0))+
	 w[5]*s[0,-1]);
  }
  else if (n == 4)
    s[]=(w[0]*(fine(s,1,1)+fine(s,1,0))+
	 w[1]*(fine(s,0,0)+fine(s,0,1))+
	 w[2]*(fine(s,-1,0)+fine(s,-1,1))+
	 w[3]*(fine(s,0,-1)+fine(s,0,2))+
	 w[4]*(fine(s,1,-1)+fine(s,1,2))+
	 w[5]*s[1]);
  else // (n == 8)
    s[]=(w[0]*(fine(s,1,1)+fine(s,0,1))+
	 w[1]*(fine(s,0,0)+fine(s,1,0))+
	 w[2]*(fine(s,0,-1)+fine(s,1,-1))+
	 w[3]*(fine(s,-1,0)+fine(s,2,0))+
	 w[4]*(fine(s,-1,1)+fine(s,2,1))+
	 w[5]*s[0,1]);
}

static inline void coarsen2D_twoleaf_quad(Point point,scalar s,int n)
{
  if (n==5 || n==10)
    {
      if (n==5)
	s[]=(9*(fine(s,0,0)+fine(s,0,1)+
		fine(s,1,1)+fine(s,1,0))-
	     (s[-1]+s[1]+
	      fine(s,1,2)+fine(s,0,2)+
	      fine(s,1,-1)+fine(s,0,-1)))/30;
      
      else if (n == 10)
	s[]=(9*(fine(s,0,0)+fine(s,0,1)+
		fine(s,1,1)+fine(s,1,0))-
	     (s[0,-1]+s[0,1]+
	      fine(s,2,1)+fine(s,2,0)+
	      fine(s,-1,1)+fine(s,-1,0)))/30;
    }
  else
    {
      double w[6]={0.343774069319640, 0.287034659820282, 0.256482670089859,
		   -0.016174582798460, -0.020539152759949, -0.050449293966624};
      if (n == 3)
	s[]=(w[0]*fine(s,0,0)+
	     w[1]*(fine(s,1,0)+fine(s,0,1))+
	     w[2]*fine(s,1,1)+
	     w[3]*(fine(s,0,2)+fine(s,2,0))+
	     w[4]*(fine(s,1,2)+fine(s,2,1))+
	     w[5]*(s[-1]+s[0,-1]));
      else if (n == 6)
	s[]=(w[0]*fine(s,1,0)+
	     w[1]*(fine(s,1,1)+fine(s,0,0))+
	     w[2]*fine(s,0,1)+
	     w[3]*(fine(s,1,2)+fine(s,-1,0))+
	     w[4]*(fine(s,0,2)+fine(s,-1,1))+
	     w[5]*(s[1]+s[0,-1]));
      else if (n== 9)
	s[]=(w[0]*fine(s,0,1)+
	     w[1]*(fine(s,1,1)+fine(s,0,0))+
	     w[2]*fine(s,1,0)+
	     w[3]*(fine(s,0,-1)+fine(s,2,1))+
	     w[4]*(fine(s,1,-1)+fine(s,2,0))+
	     w[5]*(s[-1]+s[0,1]));
      else // (n==12)
	s[]=(w[0]*fine(s,1,1)+
	     w[1]*(fine(s,0,1)+fine(s,1,0))+
	     w[2]*fine(s,0,0)+
	     w[3]*(fine(s,1,-1)+fine(s,-1,1))+
	     w[4]*(fine(s,0,-1)+fine(s,-1,0))+
	     w[5]*(s[1]+s[0,1]));
	}
}

static inline void coarsen2D_threeleaf_quad(Point point,scalar s,int n)
{
  double w[5]= {0.312119914346895, 0.2674946466809427,
		-0.050963597430407, -0.035032119914347, -0.01910064239828};
  if (n == 7)
    s[]=(w[0]*(fine(s,0,0)+fine(s,1,0))+
	 w[1]*(fine(s,0,1)+fine(s,1,1))+
	 w[3]*(s[-1]+s[1])+w[2]*s[0,-1]+
	 w[4]*(fine(s,0,2)+fine(s,1,2)));
  else if (n == 11)
    s[]=(w[0]*(fine(s,0,0)+fine(s,0,1))+
	 w[1]*(fine(s,1,1)+fine(s,1,0))+
	 w[3]*(s[0,-1]+s[0,1])+w[2]*s[-1]+
	 w[4]*(fine(s,2,0)+fine(s,2,1)));
  else if (n == 13)
    s[]=(w[0]*(fine(s,0,1)+fine(s,1,1))+
	 w[1]*(fine(s,1,0)+fine(s,0,0))+
	 w[3]*(s[-1]+s[1])+w[2]*s[0,1]+
	 w[4]*(fine(s,1,-1)+fine(s,0,-1)));
  else // (n == 14)
    s[]=(w[0]*(fine(s,1,1)+fine(s,1,0))+
	 w[1]*(fine(s,0,0)+fine(s,0,1))+
	 w[3]*(s[0,-1]+s[0,1])+w[2]*s[1]+
	 w[4]*(fine(s,-1,1)+fine(s,-1,0)));
}


static inline void coarsen2D_fourleaf_quad(Point point, scalar s){
  double sum = 0;
  foreach_child()
    sum+=8*s[];
  sum-=(s[-1]+s[1]+s[0,1]+s[0,-1]);
  s[] = sum/28; 
}

static inline void coarsen3D_oneleaf_quad(Point point, scalar s, int n){
}
static inline void coarsen3D_twoleaf_quad(Point point, scalar s, int n){
}
static inline void coarsen3D_threeleaf_quad(Point point, scalar s, int n){
}
static inline void coarsen3D_fourleaf_quad(Point point, scalar s, int n){
}

static inline void coarsen3D_fiveleaf_quad(Point point, scalar s,int n){
  double w[5]={0.181093394077449, 0.119589977220957, 0.018223234624146,
	       -0.039863325740319, -0.116173120728929};
  fprintf(ferr,"Lets use this variable %g\n",w[2]);
}

static inline void coarsen3D_sixleaf_quad(Point point, scalar s){
  double sum = 0;
  foreach_child()
    sum+=4*s[];
  sum-=(s[-1]+s[1]+s[0,1]+s[0,-1]+s[0,0,1]+s[0,0,-1]);
  s[] = sum/28; 
}

static inline int caser(Point point)
{
  int n=0;
  int j=0;
  for (int i = -1;i<=1;i+=2){
    if (is_leaf(neighbor(i,0,0)))
      n+=pow(2,j);
    j++;
#if dimension > 1 
    if (is_leaf(neighbor(0,i,0)))
      n+=pow(2,j);
    j++;
#endif
#if dimension > 2 
    if (is_leaf(neighbor(0,0,i)))
      n+=pow(2,j);
    j++;
#endif
  }
  return n;
}


static inline void quadratic_coarsening(Point point, scalar s)
{
  int n = caser(point);
  if (n==0)
    coarsen_quadratic(point, s); 
#if dimension == 1
  else if (n==1) 
    coarsen1D_leftleaf_quad(point,s);
  else if (n==2)
    coarsen1D_rightleaf_quad(point,s);
  else if (n==3)
    coarsen1D_twoleaf_quad(point,s);
  
#elif dimension == 2
  else if (n==1||n==2||n==4||n==8)
    coarsen2D_oneleaf_quad(point,s,n);
  else if (n==3||n==5||n==6||n==9||n==10||n==12)  
    coarsen2D_twoleaf_quad(point,s,n); 
  else if(n==7||n==11||n==13||n==14) 
    coarsen2D_threeleaf_quad(point,s,n);
  else if (n==15)
    coarsen2D_fourleaf_quad(point,s);
  
#else  // dimension == 2
  else if (n==1||n==2||n==4||n==8||n==16||n==32)
    coarsen3D_oneleaf_quad(point,s,n);
  else if (n==3||n==5||n==9||n==17||n==33||n==6||n==10||n==18||
	   n==34||n==12||n==20||n==36||n==24||n==40||n==48)
    coarsen3D_twoleaf_quad(point,s,n);
  else if(n==7||n==11||n==19||n==35||n==13||n==21||n==37||
	  n==25||n==41||n==49||n==14||n==22||n==38||n==26||
	  n==42||n==50||n==28||n==44||n==56)
    coarsen3D_threeleaf_quad(point,s,n);
  else if(n==15||n==23||n==27||n==29||n==30||n==39||n==43||
	  n==45||n==46||n==51||n==53||n==54||n==57||n==58||n==60) 
    coarsen3D_fourleaf_quad(point,s,n);
  else if (n==31||n==47||n==55||n==59||n==61||n==62)
    coarsen3D_fiveleaf_quad(point,s,n);
  else if (n==63)
    coarsen3D_sixleaf_quad(point,s);
#endif
}

static inline void quadratic_boundary_restriction(Point point, scalar s)
{
  int n = caser(point);
  if (n==0)
    restriction_average(point, s); 
#if dimension == 1
  else if (n==1) 
    coarsen1D_leftleaf_quad(point,s);
  else if (n==2)
    coarsen1D_rightleaf_quad(point,s);
  else if (n==3)
    coarsen1D_twoleaf_quad(point,s);
  
#elif dimension == 2
  else if (n==1||n==2||n==4||n==8)
    coarsen2D_oneleaf_quad(point,s,n);
  else if (n==3||n==5||n==6||n==9||n==10||n==12)  
    coarsen2D_twoleaf_quad(point,s,n); 
  else if(n==7||n==11||n==13||n==14) 
    coarsen2D_threeleaf_quad(point,s,n);
  else if (n==15)
    coarsen2D_fourleaf_quad(point,s);
#else 
  else if (n==1||n==2||n==4||n==8||n==16||n==32)
    coarsen3D_oneleaf_quad(point,s,n);
  else if (n==3||n==5||n==9||n==17||n==33||n==6||n==10||n==18||
	   n==34||n==12||n==20||n==36||n==24||n==40||n==48)
    coarsen3D_twoleaf_quad(point,s,n);
  else if(n==7||n==11||n==19||n==35||n==13||n==21||n==37||
	  n==25||n==41||n==49||n==14||n==22||n==38||n==26||
	  n==42||n==50||n==28||n==44||n==56)
    coarsen3D_threeleaf_quad(point,s,n);
  else if(n==15||n==23||n==27||n==29||n==30||n==39||n==43||n==45||
	  n==46||n==51||n==53||n==54||n==57||n==58||n==60) 
    coarsen3D_fourleaf_quad(point,s,n);
  else if (n==31||n==47||n==55||n==59||n==61||n==62)
    coarsen3D_fiveleaf_quad(point,s,n);
  else if (n==63)
    coarsen3D_sixleaf_quad(point,s);
#endif
}
