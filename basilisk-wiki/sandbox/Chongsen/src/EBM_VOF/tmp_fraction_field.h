

#include "TPR2D.h"


static inline bool contact_embed_cell (Point point, scalar c, scalar contact_angle, scalar cs, scalar mark)
{

if (mark[]==10){
  if (c[] > VFTL && c[] < cs[]-VFTL){
    if (contact_angle[] <= 90){
      for (int i = -1; i <= 1; i += 2){
        foreach_dimension(){
          if (mark[i]==10 && c[i] <= VFTL)
            return true;
        }
        if ((mark[i,i]==10 && c[i,i] <= VFTL) || (mark[i,-i]==10 && c[i,-i] <= VFTL))
          return true;
      }
    }
    else if (contact_angle[] > 90){
      for (int i = -1; i <= 1; i += 2){
        foreach_dimension(){
          if (mark[i]==10 && c[i] >= cs[i]-VFTL)
            return true;
        }
        if ((mark[i,i]==10 && c[i,i] >= cs[i,i]-VFTL) || (mark[i,-i]==10 && c[i,-i] >= cs[i,-i]-VFTL))
          return true;
      }
    }
  }
  else if (c[] >= cs[]-VFTL){
    for (int i = -1; i <= 1; i += 2){
      foreach_dimension(){
        if (mark[i]==10 && c[i] <= VFTL)
          return true;
      }

    }
  }
  else if (c[] <= VFTL){
    for (int i = -1; i <= 1; i += 2){
      foreach_dimension(){
        if (mark[i]==10 && c[i] >= cs[i]-VFTL)
          return true;
      }

    }
  }
    return false;
}
else{
      return false;
}


}

static inline bool i_w (Point point, scalar cs)
{
  if (cs[] >= 1.) {
    for (int i = -1; i <= 1; i += 2)
      foreach_dimension()
	if (cs[i] <= 0.)
	  return true;
  }
  else if (cs[] > 0. && cs[] < 1) 
	  return true;
  return false;
}

static inline bool interfacial_w (Point point, scalar cs)
{
  if (cs[] >= 1.) {
    for (int i = -1; i <= 1; i += 2)
      foreach_dimension()
	      if (cs[i] <= 0.)
	        return true;
  }
  else if (cs[] <= 0.) {
    for (int i = -1; i <= 1; i += 2)
      foreach_dimension()
	      if (cs[i] >= 1.)
	        return true;
  }
  else // cs[] > 0. && cs[] < 1.
    return true;
  return false;
}

static inline bool interfacial_i (Point point, scalar c)   //  judge the interface of 'c' 
{
  if (c[] >= 1.) {
    for (int i = -1; i <= 1; i += 2)
      foreach_dimension()
	      if (c[i] <= 0.)
	        return true;
  }
  else if (c[] <= 0.) {
    for (int i = -1; i <= 1; i += 2)
      foreach_dimension()
	      if (c[i] >= 1.)
	        return true;
  }
  else // c[] > 0. && c[] < 1.
    return true;
  return false;
}

coord facet_normal_fs (Point point, face vector s)
{
  if (s.x.i >= 0) {
    coord n;
    double nn = 0.;
    foreach_dimension() {
      n.x = s.x[] - s.x[1];
      nn += fabs(n.x);
    }
    if (nn > 0.)
      foreach_dimension()
	n.x /= nn;
    else
      foreach_dimension()
	n.x = 1./dimension;
    return n;
  }
  coord n = {0,0};
  return n;
}

void reconstruction_cs (const scalar cs, face vector fs, vector ms, scalar alphacs)
{
  foreach() {
    if (cs[] <= 0. || cs[] >= 1.){
      alphacs[] = 0.;
      foreach_dimension()
	      ms.x[] = 0.;
    }
    if ((cs[] <= 0. || cs[] >= 1.) && interfacial_w(point, cs)){
      alphacs[] = 10.;
      coord m = facet_normal_fs (point, fs);
      foreach_dimension()
	      ms.x[] = m.x;
    }
    else {
      coord m = facet_normal (point, cs, fs);
      foreach_dimension()
	      ms.x[] = m.x;
      alphacs[] = plane_alpha (cs[], m);
    }
  }
}

void reconstruction_mc (const scalar c, vector mc)
{
  foreach() {
    if (c[] <= 0. || c[] >= 1.) {
      foreach_dimension()
	  mc.x[] = 0.;
    }
    else {
      coord m = interface_normal (point, c);
      foreach_dimension()
      mc.x[] = m.x;
    }
  }
}

void reconstruction_mc_myc (const scalar c, vector mc)
{
  foreach() {
    if (c[] <= 0. || c[] >= 1.) {
      foreach_dimension()
	  mc.x[] = 0.;
    }
    else {
      coord m = mycs (point, c);
      foreach_dimension()
      mc.x[] = m.x;
    }
  }
}


trace
void reconstruction_test (const scalar c, scalar tmp_c, vector n, scalar alpha)
{
  foreach() {

    if (c[] <= 0. || c[] >= 1.) {
      alpha[] = 0.;
      foreach_dimension()
	n.x[] = 0.;
    }
    else {
      coord m = interface_normal (point, tmp_c);
      foreach_dimension()
	      n.x[] = m.x;
      alpha[] = plane_alpha (tmp_c[], m);
    }
  }
}


/**
Construct the extended volume field, similar to Popinet/contact sandbox.*/

trace
void reconstruction_tmp_embed_fraction_field (const scalar c, scalar cs, vector vmc, vector vms, scalar alphacs, scalar contact_angle, scalar tmp_c, scalar mark)
{

    vector nc[];
    scalar alphac[];

    foreach(){
      coord mc={10,10};
      coord ms={10,10};
      coord p_mof[2]={{10,10},{10,10}};
      coord pp[5]={{10,10},{10,10},{10,10},{10,10},{10,10}};
        mc.x = vmc.x[];
        mc.y = vmc.y[];
        ms.x = vms.x[];
        ms.y = vms.y[];
        if (mark[] == 2){
          tmp_c[] = 0;
        }
        else if (mark[] == 3){
          tmp_c[] = 1;
        }
        else if ( mark[]==4 && (mc.x!=0||mc.y!=0)){
            if (cs[]>=1 || cs[]<=0){
                tmp_c[] = c[];
            }
            else {
                coord nnc = normal_contact(ms, mc, contact_angle[]*pi/180);

                double aalphac = polygon_alpha(c[], nnc, ms, alphacs[], p_mof, pp);
                
                tmp_c[] = line_area(nnc.x, nnc.y, aalphac);
            }
        }
        else if (mark[] == 6 || mark[] == 7 || mark[] == 8){
          tmp_c[] = c[];
        }
        else if (mark[] == 5){
            if (cs[]>=1 || cs[]<=0){
                tmp_c[] = c[];
            }
            else {
              
              coord nnc = normal_contact(ms, mc, contact_angle[]*pi/180);

              
              coord mst = {ms.y, -ms.x};
              int a = 0;
              int b = 0;
              if (mst.x*nnc.x + mst.y*nnc.y >=0){
                a = mst.x == 0 ? 0 :  sign(mst.x);
                b = mst.y == 0 ? 0 :  sign(mst.y);
              }
              else {
                a = mst.x == 0 ? 0 : -sign(mst.x);
                b = mst.y == 0 ? 0 : -sign(mst.y);
              }
              
              if (mark[a,0]!=4 && mark[0,b]!=4 && mark[a,b]!=4 && mark[2*a,0]!=4 && mark[0,2*b]!=4 && mark[2*a,2*b]!=4 && mark[a,2*b]!=4 && mark[2*a,b]!=4)
                tmp_c[] = c[];
              else 
                tmp_c[] = c[] + 1 - cs[];

              
            }
        }
    }

    foreach(){
      if (mark[]==5 && tmp_c[]==c[]){
        foreach_dimension(){
          if (cs[-1]>0 && cs[-1]<1 && tmp_c[-1]==c[-1]+1-cs[-1])
            tmp_c[]=c[]+1-cs[];
          if (cs[ 1]>0 && cs[ 1]<1 && tmp_c[ 1]==c[ 1]+1-cs[ 1])
            tmp_c[]=c[]+1-cs[];
        }
      }
    }

    foreach(){  // cs[] = 0
        if (mark[] == 1){

            double fc1 = 0., fc2 = 0., sfc1 = 0., sfc2 = 0., coef = 0.;
            coord o = {x, y};


            foreach_neighbor(2){
                

                if (mark[]==3){

              
                    coord mc={10,10};
                    mc.x = - vms.x[];
                    mc.y = - vms.y[];
                    double alpha = - alphacs[];
                    coef = 1;
                    sfc2 += coef;
                    coord a = {x, y}, b;
                    foreach_dimension()
	                    a.x = (o.x - a.x)/Delta - 0.5, b.x = a.x + 1.;

                    fc2 = max(fc2, coef*rectangle_fraction (mc, alpha, a, b)); 
                    if ((mc.x != 0 || mc.y != 0) && alphacs[]==10)

                      fc2 = 1;


                }
                 
            }
            if (sfc2 > 0.){
                  tmp_c[] = fc2;

                    }
                    else {
                    tmp_c[] = 0;

                }
            foreach_neighbor(2){

             
                if (mark[]==4)  {

                    
                    coord mc={10,10};
                    coord ms={10,10};
                    coord p_mof[2]={{10,10},{10,10}};
                    coord pp[5]={{10,10},{10,10},{10,10},{10,10},{10,10}};
                    mc.x = vmc.x[];
                    mc.y = vmc.y[];
                    ms.x = vms.x[];
                    ms.y = vms.y[];
                    coord nnc = normal_contact(ms, mc, contact_angle[]*pi/180);
                    double aalphac = polygon_alpha(c[], nnc, ms, alphacs[], p_mof, pp);
                    

                      coef = 1;
	                    sfc1 += coef;

	                    coord a = {x, y}, b;
	                    foreach_dimension()
	                      a.x = (o.x - a.x)/Delta - 0.5, b.x = a.x + 1.;

                      fc1 += coef*rectangle_fraction (nnc, aalphac, a, b);  

	              }
 
            }

                if (sfc1 > 0.){
	                tmp_c[] = fc1/sfc1;
                }


        }
        
        
        else if ((mark[] == 2 || mark[] == 3) && cs[]>0 && cs[]<1){

            double fc1 = 0., sfc1 = 0., coef = 0.;
            coord o = {x, y};
            
            coord ms={10,10};
            coord nc={10,10};
            foreach_neighbor(2){
              if (mark[]==4)  {
                coord mc={10,10};
                mc.x = vmc.x[];
                mc.y = vmc.y[];
                ms.x = vms.x[];
                ms.y = vms.y[];
                nc = normal_contact(ms, mc, contact_angle[]*pi/180);
              }
            }
            
            coord mst = {nc.y, -nc.x};
            int aa = 0;
            int bb = 0;
            if (mst.x*ms.x + mst.y*ms.y >=0){
              aa = mst.x == 0 ? 0 :  sign(mst.x);
              bb = mst.y == 0 ? 0 :  sign(mst.y);
            }
            else {
              aa = mst.x == 0 ? 0 : -sign(mst.x);
              bb = mst.y == 0 ? 0 : -sign(mst.y);
            }

            foreach_neighbor(2) {
                if (mark[]==4) {

                  coord a = {x, y};
                  if ((o.x - a.x)/Delta*aa>=0 && (o.y - a.y)/Delta*bb>=0){

                    coord mc={10,10};
                    coord ms={10,10};
                    coord p_mof[2]={{10,10},{10,10}};
                    coord pp[5]={{10,10},{10,10},{10,10},{10,10},{10,10}};
                    mc.x = vmc.x[];
                    mc.y = vmc.y[];
                    ms.x = vms.x[];
                    ms.y = vms.y[];
                    coord nnc = normal_contact(ms, mc, contact_angle[]*pi/180);
                    double aalphac = polygon_alpha(c[], nnc, ms, alphacs[], p_mof, pp);

                      coef = 1;
	                    sfc1 += coef;

	                    coord a = {x, y}, b;
	                    foreach_dimension()
	                      a.x = (o.x - a.x)/Delta - 0.5, b.x = a.x + 1.;
	                    fc1 += coef*rectangle_fraction (nnc, aalphac, a, b);
                    
                  }

              }
            }
                if (sfc1 > 0.){
	                tmp_c[] = fc1/sfc1;
                }
            
        }
        
        
        
    boundary ({tmp_c});
}}


/**
4 is the contact angle cell;*/

trace
void sort_cell(scalar c, scalar cs, scalar contact_angle, scalar mark)
{


  scalar mark2[];
  foreach(){
    if (interfacial_w(point,cs) && cs[] > 0)                          
      mark2[]=10;
    
    else
      mark2[]=0;
  }
  foreach(){
    mark[] = 0;                                                       
    if (cs[] <= 0)
      mark[] = 1;                                                     

    else if (interfacial_w(point, cs)){                               

      if (c[] <= VFTL)  
        mark[] = 2;                                                   
      else if (c[] >= (cs[]-VFTL))  
        mark[] = 3;                                                   
      else {
  
        mark[] = 5;                                                   
        
        coord ns = facet_normal(point, cs, fs);
        coord nf = mycs(point, c);
        coord nc = normal_contact(ns, nf, contact_angle[]*pi/180);
        coord mnc = {0,0};
        if (ns.x*nc.y - ns.y*nc.x > 0){
          mnc.x =  nc.y;
          mnc.y = -nc.x;
        }
        else {
          mnc.x = -nc.y;
          mnc.y =  nc.x;
        }
        
        int a = sign2(mnc.x), b = sign2(mnc.y), markk4 = 1;
        
        if (cs[a,0]>0 && cs[a,0]<1 && c[a,0]>=VFTL && c[a,0]<=cs[a,0]-VFTL)
          markk4 = 0;
        if (cs[0,b]>0 && cs[0,b]<1 && c[0,b]>=VFTL && c[0,b]<=cs[0,b]-VFTL)
          markk4 = 0;
        if (cs[a,b]>0 && cs[a,b]<1 && c[a,b]>=VFTL && c[a,b]<=cs[a,b]-VFTL)
          markk4 = 0;
        if ((cs[a,0]<=0 || cs[a,0]>=1) && (cs[0,b]<=0 || cs[0,b]>=1) && (cs[a,b]<=0 || cs[a,b]>=1))
          markk4 = 1;
        
        if (markk4)
          mark[] = 4;
        
      
      }
        
    }

    else if (cs[] >= 1.){
      if (interfacial_i(point, c))
        mark[] = 6;                                                   
      else if (c[] <= 0.)
        mark[] = 7;                                                   
      else if (c[] >= 1.)
        mark[] = 8;                                                   
    }

  }


}



/**
Calculate polygon area*/


#define nann 10
trace
double area_polygon (coord a[6]){
    coord o_xy[1] = {{0, 0}};
    int num = 0;
    for (int i = 0; i <= 5; i++){
        if (a[i].x != nann && a[i].y != nann){
            foreach_dimension()
            o_xy[0].x += a[i].x;
            num++;
        }
    }
    if (num > 0){
        foreach_dimension()
        o_xy[0].x /= num;
        double angle[num];
        for (int j = 0; j <= num-1; j++)
            angle[j] = atan2(a[j].y - o_xy[0].y, a[j].x - o_xy[0].x);
    
        for (int ii = 0; ii <= num-2; ii++) {
            for (int jj = 0; jj <= num-2-ii; jj++) {
                if(angle[ii] > angle[ii+jj+1]){
                    swap(double, angle[ii], angle[ii+jj+1]);
                    swap(coord, a[ii], a[ii+jj+1]);
                }
        }
    }

        double area = 0;


        for (int i = 0 ; i < num ; i++){
          int b = i - 1 <    0    ? i - 1 + num : i - 1;
          int f = i + 1 > num - 1 ? i + 1 - num : i + 1;
          area += a[i].x*(a[f].y-a[b].y)/2;
        }

        return area;
    }
    else{
        fprintf (stderr,"something wrong in area_polygon");
        exit (1);
    }
}



/**
Calculate the flux in mixed cells.*/

trace
double polygon_fraction(coord pp[5], double un, double s, coord p_mof[2], coord mc, double alpha)
{

    int ppf_number = 0;  
    int np_number = 0;
    coord ppf[5]={{nann, nann},{nann, nann},{nann, nann},{nann, nann},{nann, nann}}; 
    coord pp_mof[2]={{nann, nann},{nann, nann}};
    coord mmc={nann, nann};
    int np[5]={nann,nann,nann,nann,nann};    
    for (int pri = 0; pri <= 4; pri++){
        if (fabs(pp[pri].x) < 1e-30)
        pp[pri].x = 0;
        if (fabs(pp[pri].y) < 1e-30)
        pp[pri].y = 0;
    }

    if (s > 0){
        for (int i = 0; i <= 4; i++)

            if (pp[i].x != nann || pp[i].y != nann){
                ppf[ppf_number].x = - pp[i].x;
                ppf[ppf_number].y =   pp[i].y;
                ppf_number++;
            }
        pp_mof[0].x = - p_mof[0].x;
        pp_mof[0].y =   p_mof[0].y;
        pp_mof[1].x = - p_mof[1].x;
        pp_mof[1].y =   p_mof[1].y;
        mmc.x = - mc.x;
        mmc.y =   mc.y;
    }
    else{
        for (int ii = 0; ii <= 4; ii++)
            if (pp[ii].x != nann || pp[ii].y != nann){
                ppf[ppf_number].x =   pp[ii].x; 
                ppf[ppf_number].y =   pp[ii].y;
                ppf_number++;
            }
        pp_mof[0].x =   p_mof[0].x;
        pp_mof[0].y =   p_mof[0].y;
        pp_mof[1].x =   p_mof[1].x;
        pp_mof[1].y =   p_mof[1].y;
        mmc.x =   mc.x;
        mmc.y =   mc.y;
    }

    for (int iii = 0; iii <= ppf_number - 1; iii++){
        if (ppf[iii].x <= (s*un - 0.5)){
                np[np_number] = iii;
                np_number++;
        }
    }

    if (np_number != 0){
        coord fp[5]={{nann, nann},{nann, nann},{nann, nann},{nann, nann},{nann, nann}};




    if (np_number == 1){
        int tmp_front = np[0] + 1 > (ppf_number - 1) ? np[0] - ppf_number + 1 : np[0] + 1;
        int tmp_back  = np[0] - 1 <        0         ? np[0] + ppf_number - 1 : np[0] - 1;
        fp[0].x = ppf[np[0]].x;
        fp[0].y = ppf[np[0]].y;
        fp[1].x = s*un - 0.5;

        fp[1].y = fabs(ppf[tmp_front].y - ppf[np[0]].y) < 1e-10 ? ppf[np[0]].y : ppf[np[0]].y + ((s*un - 0.5) - ppf[np[0]].x)*(ppf[tmp_front].y - ppf[np[0]].y)/(ppf[tmp_front].x - ppf[np[0]].x);
        fp[2].x = s*un - 0.5;

        fp[2].y = fabs(ppf[tmp_back].y - ppf[np[0]].y) < 1e-10 ? ppf[np[0]].y : ppf[np[0]].y + ((s*un - 0.5) - ppf[np[0]].x)*(ppf[tmp_back].y - ppf[np[0]].y)/(ppf[tmp_back].x - ppf[np[0]].x);
    }
    else if (np_number == 2){
        int tmp_front = np[0] + 1 > (ppf_number - 1) ? np[0] - ppf_number + 1 : np[0] + 1;
        int tmp_back  = np[0] - 1 <        0         ? np[0] + ppf_number - 1 : np[0] - 1;


        fp[0].x = ppf[np[0]].x;
        fp[0].y = ppf[np[0]].y;
        if (ppf[tmp_back].x <= (s*un - 0.5)){
            fp[1].x = s*un - 0.5;

            fp[1].y = fabs(ppf[tmp_front].y - ppf[np[0]].y) < 1e-10 ? ppf[np[0]].y : ppf[np[0]].y + ((s*un - 0.5) - ppf[np[0]].x)*(ppf[tmp_front].y - ppf[np[0]].y)/(ppf[tmp_front].x - ppf[np[0]].x);
            fp[3].x = ppf[tmp_back].x;
            fp[3].y = ppf[tmp_back].y;
            int tmp_bback  = tmp_back - 1 < 0 ? tmp_back + ppf_number - 1 : tmp_back - 1;
            fp[2].x = s*un - 0.5;

            fp[2].y = fabs(ppf[tmp_bback].y - ppf[tmp_back].y) < 1e-10 ? ppf[tmp_back].y : ppf[tmp_back].y + ((s*un - 0.5) - ppf[tmp_back].x)*(ppf[tmp_bback].y - ppf[tmp_back].y)/(ppf[tmp_bback].x - ppf[tmp_back].x);
        }
        else{
            fp[1].x = s*un - 0.5;

            fp[1].y = fabs(ppf[tmp_back].y - ppf[np[0]].y) < 1e-10 ? ppf[np[0]].y : ppf[np[0]].y + ((s*un - 0.5) - ppf[np[0]].x)*(ppf[tmp_back].y - ppf[np[0]].y)/(ppf[tmp_back].x - ppf[np[0]].x);
            fp[3].x = ppf[tmp_front].x;
            fp[3].y = ppf[tmp_front].y;
            int tmp_ffront  = tmp_front + 1 > (ppf_number - 1) ? tmp_front - ppf_number + 1 : tmp_front + 1;

            fp[2].x = s*un - 0.5;

            fp[2].y = fabs(ppf[tmp_ffront].y - ppf[tmp_front].y) < 1e-10 ? ppf[tmp_front].y : ppf[tmp_front].y + ((s*un - 0.5) - ppf[tmp_front].x)*(ppf[tmp_ffront].y - ppf[tmp_front].y)/(ppf[tmp_ffront].x - ppf[tmp_front].x);
        }
    }
    else if (np_number == 3){
        int tmp_f = np[0] + 1 > (ppf_number - 1) ? np[0] - ppf_number + 1 : np[0] + 1;
        int tmp_b = np[0] - 1 <        0         ? np[0] + ppf_number - 1 : np[0] - 1;
        fp[0].x = ppf[np[0]].x;
        fp[0].y = ppf[np[0]].y;
        if (ppf[tmp_f].x <= (s*un - 0.5)){
            int tmp_ff  = tmp_f + 1 > (ppf_number - 1) ? tmp_f - ppf_number + 1 : tmp_f + 1;
            fp[1].x = ppf[tmp_f].x;
            fp[1].y = ppf[tmp_f].y;
            if (ppf[tmp_ff].x <= (s*un - 0.5)){
                fp[2].x = ppf[tmp_ff].x;
                fp[2].y = ppf[tmp_ff].y;
                int tmp_fff  = tmp_ff + 1 > (ppf_number - 1) ? tmp_ff - ppf_number + 1 : tmp_ff + 1;
                fp[3].x = s*un - 0.5;

                fp[3].y = fabs(ppf[tmp_fff].y - ppf[tmp_ff].y) < 1e-10 ? ppf[tmp_ff].y : ppf[tmp_ff].y + ((s*un - 0.5) - ppf[tmp_ff].x)*(ppf[tmp_fff].y - ppf[tmp_ff].y)/(ppf[tmp_fff].x - ppf[tmp_ff].x);
                fp[4].x = s*un - 0.5;

                fp[4].y = fabs(ppf[tmp_b].y - ppf[np[0]].y)  < 1e-10 ? ppf[np[0]].y : ppf[np[0]].y + ((s*un - 0.5) - ppf[np[0]].x)*(ppf[tmp_b].y - ppf[np[0]].y)/(ppf[tmp_b].x - ppf[np[0]].x);
            }
            else{
                fp[2].x = s*un - 0.5;

                fp[2].y = fabs(ppf[tmp_ff].y - ppf[tmp_f].y)  < 1e-10 ? ppf[tmp_f].y : ppf[tmp_f].y + ((s*un - 0.5) - ppf[tmp_f].x)*(ppf[tmp_ff].y - ppf[tmp_f].y)/(ppf[tmp_ff].x - ppf[tmp_f].x);
                int tmp_bb  = tmp_b - 1 < 0 ? tmp_b + ppf_number - 1 : tmp_b - 1;
                fp[4].x = ppf[tmp_b].x;
                fp[4].y = ppf[tmp_b].y;
                fp[3].x = s*un - 0.5;

                fp[3].y = fabs(ppf[tmp_bb].y - ppf[tmp_b].y) < 1e-10 ? ppf[tmp_b].y : ppf[tmp_b].y + ((s*un - 0.5) - ppf[tmp_b].x)*(ppf[tmp_bb].y - ppf[tmp_b].y)/(ppf[tmp_bb].x - ppf[tmp_b].x);
            }
        }
        else {
            fp[1].x = s*un - 0.5;

            fp[1].y = fabs(ppf[tmp_f].y - ppf[np[0]].y)  < 1e-10 ? ppf[np[0]].y : ppf[np[0]].y + ((s*un - 0.5) - ppf[np[0]].x)*(ppf[tmp_f].y - ppf[np[0]].y)/(ppf[tmp_f].x - ppf[np[0]].x);
            fp[4].x = ppf[tmp_b].x;
            fp[4].y = ppf[tmp_b].y;
            int tmp_bb  = tmp_b - 1 < 0 ? tmp_b + ppf_number - 1 : tmp_b - 1;
            fp[3].x = ppf[tmp_bb].x;
            fp[3].y = ppf[tmp_bb].y;
            int tmp_bbb  = tmp_bb - 1 < 0 ? tmp_bb + ppf_number - 1 : tmp_bb - 1;
            fp[2].x = s*un - 0.5;

            fp[2].y = fabs(ppf[tmp_bbb].y - ppf[tmp_bb].y) < 1e-10 ? ppf[tmp_bb].y : ppf[tmp_bb].y + ((s*un - 0.5) - ppf[tmp_bb].x)*(ppf[tmp_bbb].y - ppf[tmp_bb].y)/(ppf[tmp_bbb].x - ppf[tmp_bb].x);
        }
    }

    int jk = 0;
    coord ap[6]={{nann, nann},{nann, nann},{nann, nann},{nann, nann},{nann, nann},{nann, nann}};
    coord lp[1]={{nann,nann}};
    for (int iiii = 0; iiii <= 1; iiii++){

        if(fabs(pp_mof[iiii].x) != nann && pp_mof[iiii].x <= s*un-0.5){
            ap[jk].x = pp_mof[iiii].x;
            ap[jk].y = pp_mof[iiii].y;
            jk++;
        }
        else {
            lp[0].x = pp_mof[iiii].x;
            lp[0].y = pp_mof[iiii].y;
        }
    }

    int oriant = 0;
    for (int j = 0; j <= np_number + 1; j++){
        if (fp[j].x * mmc.x + fp[j].y * mmc.y - alpha < 0){
            ap[jk + oriant].x = fp[j].x;
            ap[jk + oriant].y = fp[j].y;
            oriant++;
        }
    }

    if (jk == 0)
        if (oriant == 0)
            return 0;
        else
            return area_polygon(ap);
    else if(jk == 2)
        return area_polygon(ap);
    else if(jk == 1){
        ap[jk + oriant].x = s*un - 0.5;
        ap[jk + oriant].y = ap[0].y + ((s*un - 0.5) - ap[0].x)*(lp[0].y - ap[0].y)/(lp[0].x - ap[0].x);


        return area_polygon(ap);
    }
    else{
        fprintf (stderr,"something wrong in ploy_fraction");
        exit (1);
    }

    }
    else {

        return 0;
    }

}

trace
void reconstruction_tmp_c (scalar c, scalar contact_angle, scalar cs, face vector fs, scalar tmp_c, scalar mark)
{
    vector ms[], tmp_mc[];
    scalar alphacs[];

    sort_cell(c, cs, contact_angle, mark);
    
    reconstruction_cs(cs, fs, ms, alphacs);

    reconstruction_mc_myc(c, tmp_mc);

    reconstruction_tmp_embed_fraction_field (c, cs, tmp_mc, ms, alphacs, contact_angle, tmp_c, mark);


}


/**
Clean too small mixed cells.*/


foreach_dimension()
  static double interface_fractions_x (coord m, double alpha, int droite){
  double vback;
  double ss = droite ? -0.5 : 0.5;

#if dimension <= 2  
  if (fabs(m.y) > 1.e-4){
    double aa = (alpha - ss * m.x)/m.y;
    if (aa >= - 0.5 && aa <= 0.5)
      vback = m.y < 0 ? 0.5 - aa : 0.5 + aa;
    else if (aa < -0.5)
      vback = m.y < 0 ? 1. : 0.;
    else
      vback = m.y < 0 ? 0. : 1.;
  }
  else if (m.x > 0.)
    vback = droite ? 1. : 0.;
  else
    vback = droite ? 0. : 1.;
#else
  if (fabs (m.y) > 1.e-4 || fabs (m.z) > 1.e-4) {
    alpha -= ss * m.x;
    double nsum = fabs(m.y) + fabs(m.z);
    vback = line_area (m.y/nsum, m.z/nsum, alpha/nsum);
  }
  else if (m.x > 0.)
    vback = droite ? 1. : 0.;
  else
    vback = droite ? 0. : 1.;
#endif  
  return vback;
}

void myface_fraction_C (scalar cs, face vector fs, double VTOL){

  boundary({cs});
  vector normal_vector[];
  foreach() {
    coord m = mycs (point, cs);
    foreach_dimension() 
      normal_vector.x[] = m.x;
  }
  boundary((scalar*){normal_vector});
  foreach_face() {
    if (cs[-1] <= VTOL || cs[] <= VTOL)
      fs.x[] = 0.;

    else if (cs[-1] >= 1. || cs[] >= 1.)
      fs.x[] = 1.;
    else {
      double vleft = 1., vright = 1.;
      coord n1;
      foreach_dimension()
        n1.x = normal_vector.x[];
      double alpha = plane_alpha (cs[], n1);
      vleft = interface_fractions_x(n1, alpha, 1);
      
      coord n2;

      n2.x = normal_vector.x[-1];
      n2.y = normal_vector.y[-1];
#if dimension > 2
      n2.z = normal_vector.z[-1];
#endif
      double alpha2 = plane_alpha (cs[-1], n2);
      vright = interface_fractions_x(n2, alpha2, 0);      
      fs.x[] = (vleft + vright)/2;
      if(fs.x[]==0 && (cs[] > 1e-4 && cs[-1] > 1e-4) && cs[] != cs[-1])
        fs.x[] = 1e-10;
    }    
  }
}

void cleansmallcell (scalar cs, face vector fs, double csTL){


    foreach(){
      if (cs[] < csTL)
        cs[] = 0;
    }


    foreach(){
      if (cs[] >= 1. && (cs[1,1] < csTL || cs[1,-1] < csTL || cs[-1,1] < csTL || cs[-1,-1] < csTL)){
        cs[] = 1. - 1.1*csTL;
      }
    }

    myface_fraction_C (cs, fs, csTL);

    foreach(){

      double markface = 0;
      double special_treat = 0;
      foreach_dimension(){
        for (int i=0;i<=1;i++)
          if (fs.x[i] > 0 && fs.x[i] < 1){
            markface++;
          }
      }

      double x1 = 0, y1 = 0;
      if ((fs.x[]*fs.y[]==0 && fs.x[]+fs.y[]>=1)){
        markface++;
        special_treat++;
        x1 = -0.5;
        y1 = -0.5;
      }
      if ((fs.x[]*fs.y[0,1]==0 && fs.x[]+fs.y[0,1]>=1)){
        markface++;
        special_treat++;
        x1 = -0.5;
        y1 =  0.5;
      }
      if((fs.x[1]*fs.y[]==0 && fs.x[1]+fs.y[]>=1)){
        markface++;
        special_treat++;
        x1 =  0.5;
        y1 = -0.5;
      }
      if ((fs.x[1]*fs.y[0,1]==0 && fs.x[1]+fs.y[0,1]>=1)){
        markface++;
        special_treat++;
        x1 =  0.5;
        y1 =  0.5;
      }

      if(markface != 2){
        if ((fs.x[] == 1 && fs.x[1,0] == 1 && fs.y[] == 1 && fs.y[0,1] == 1) || (fs.x[]*fs.x[1,0]*fs.y[]*fs.y[0,1] == 0 && fs.x[]+fs.x[1,0]+fs.y[]+fs.y[0,1] >= 3 - csTL))
            cs[] = 1;
        else if ((fs.x[] == 0 && fs.x[1,0] == 0 && fs.y[] == 0 && fs.y[0,1] == 0))
            cs[] = 0;
      }
      else {

        coord n;
        double nn = 0.;
        foreach_dimension() {
          n.x = fs.x[] - fs.x[1];
          nn += fabs(n.x);
        }
        if (nn > 0.)
          foreach_dimension()
	          n.x /= nn;
        else
          foreach_dimension()
          	n.x = 1./dimension;

        if (special_treat){
          if ((x1!=0 && y1!=0)){
          
            double alphacs = n.x*x1 + n.y*y1;
            cs[] = plane_volume(n, alphacs);

          }
        }
        else{
          int a = 1;
          if (fs.x[]>0 && fs.x[]<1){
            x1 = -0.5;
            y1 = n.y < 0 ? 0.5 - fs.x[] : fs.x[] - 0.5;
            a = 0;
          }

          if (fs.x[1,0]>0 && fs.x[1,0]<1 && a){
            x1 = 0.5;
            y1 = n.y < 0 ? 0.5 - fs.x[1,0] : fs.x[1,0] - 0.5;
            a = 0;
          }

          if (fs.y[]>0 && fs.y[]<1 && a){
            y1 = -0.5;
            x1 = n.x < 0 ? 0.5 - fs.y[] : fs.y[] - 0.5;
            a = 0;
          }

          if (fs.y[0,1]>0 && fs.y[0,1]<1 && a){
            y1 = 0.5;
            x1 = n.x < 0 ? 0.5 - fs.y[0,1] : fs.y[0,1] - 0.5;
          }

          if ((x1!=0 && y1!=0)){
          
            double alphacs = n.x*x1 + n.y*y1;
            cs[] = plane_volume(n, alphacs);

          }
        }

        
            
      }
        
    }
    
}

