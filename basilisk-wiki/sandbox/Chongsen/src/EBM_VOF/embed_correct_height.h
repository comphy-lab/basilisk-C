#include "tmp_fraction_field.h"

/**
Calculate the height function using parabolic fitting, please refer to the article ~ for algorithm details.*/

trace
void embed_h(scalar c, scalar cs, vector mc, vector ms, scalar alphacs, scalar contact_angle, scalar mark, vector tmp_h, vector o_interface, int add){
  vector s[],scs[];
  vector s1[],s2[],cs1[],cs2[];

  foreach_dimension()
    for (int i = 0; i < nboundary; i++){
       s1.x.boundary[i] =  c.boundary[i];
       s2.x.boundary[i] =  c.boundary[i];
      cs1.x.boundary[i] = cs.boundary[i];
      cs2.x.boundary[i] = cs.boundary[i];
    }
  foreach()
    foreach_dimension(){
        s.x[] =  0;
      scs.x[] =  0;
       s1.x[] =  c[2];
       s2.x[] =  c[-2];
      cs1.x[] = cs[2];
      cs2.x[] = cs[-2];
    }

  foreach(){
    foreach_dimension(){
      for (int i = -2 ; i < 0 ; i++){
        if (mark[i]==4){
            s.x[] =  s1.x[];
          scs.x[] = cs1.x[];
        }
      }
      for (int i = 2 ; i > 0 ; i--){
        if (mark[i]==4){
            s.x[] =  s2.x[];
          scs.x[] = cs2.x[];
        }
      }
    }
  }


  vector secp[];
  vector vhcs_n[];
  foreach(){
    secp.x[]=nodata;
    secp.y[]=nodata;
    vhcs_n.x[]=nodata;
    vhcs_n.y[]=nodata;
    if (mark[]==4){

      foreach_dimension(){
        

        int hna = tmp_h.x[]==nodata ? 0 : -1;
        foreach_neighbor(1){
          if (tmp_h.x[]!=nodata)
            hna++;
        }

        if (hna>0){

        int a = sign(mc.x[])==0 ? 0 : sign(mc.x[]);

        int hcs_n = 10, h_pos = 10, intf=0, entf=0;
        double h = 0, hcs[9]={0,0,0,0,0,0,0,0,0}, fc=0;


        if (a!=0){

          int begin1 = 0;
          for (int i=-4*a; a*i<=4;i+=a){
            int j = sign(i)*(fabs(i)-2);
            double  ci = fabs(i) > 2 ?   s.x[j] :  c[i];
            double csi = fabs(i) > 2 ? scs.x[j] : cs[i];

            ci = ci > csi ? csi : ci;

            if (ci==csi)
              begin1 = 1;
            

            if (begin1){
              if (i/a+4==0)
                hcs[0] += csi;
              else
                hcs[i/a+4] = csi + hcs[i/a+4-1];

              if (a*i < 0 && ci==csi)
                intf++;
              if (a*i >= 0 && intf == 0){
                h_pos = 10;
                break;
              }

              if (a*i < 0 && ci==0 && csi != 0){
                entf++;
              }
                
              if (a*i >= 0 && entf != 0){
                h_pos = 100;
                break;
              }


              if (a*i >= 0 && ci>0 && ci==csi){
                h_pos = 10;
                break;
              }

              if ((a*i > 0) && (csi > 0) && (ci <= 1e-10)){
                h_pos = i - a;
                break;
              }
              else if ((a*i > 0) && (csi == 0)){
                h_pos = i - a;
                break;
              }
              else if ((csi > 0) && (ci > 0) && (ci <= csi))
                h+=ci;
            }

            
          }
          
          if (h_pos==10)
            hcs_n=10;
          else if (h_pos==100){
            fc = h;
            hcs_n = 0;
          }
          else{
            for (int i = -4*a; i*a <= 4;i+=a){
              if ((hcs[i/a+4] < h) && (hcs[i/a+4+1] > h)){
                fc = h-hcs[i/a+4];
                hcs_n = i + a;
                break;
              }
            }
          }

          double h0 = nodata;
          double ox = 0;
          if (fabs(hcs_n) > 2 || h_pos == 10){
            h0 = nodata;
            secp.x[] = 2;
          }


          else {

            coord mms={10,10};
            mms.x = ms.x[hcs_n];
            mms.y = ms.y[hcs_n];
            coord p_mof[2]={{10,10},{10,10}};
            coord pp[5]={{10,10},{10,10},{10,10},{10,10},{10,10}};

            double alphac = 10;


              if (cs[hcs_n]>0 && cs[hcs_n]<1 ){

                coord nc = {0,0};
                nc.x = a;
                alphac = polygon_alpha (fc, nc, mms, alphacs[hcs_n], p_mof, pp);
                ox = (p_mof[0].y + p_mof[1].y)/2;
                o_interface.x[] = ox;
                
              }
              else{
                ox = 0;
                o_interface.x[] = ox;
                coord nc = {0,0};
                nc.x = a;
                alphac = plane_alpha (fc, nc);
              }
              if (a>0 && alphac != 10)
                h0 = alphac + hcs_n;
              else if (a<0  && alphac != 10)
                h0 = 20 - alphac + hcs_n;
              
              int ori = orientation(h0);
              vhcs_n.x[] = hcs_n;
              secp.x[]=2;
            

            if (cs[hcs_n]>=1 && alphacs[hcs_n]!=10){

              secp.x[] = 2;
              ox = 0;
              o_interface.x[]=0;

            }
            else{
              if (fabs(ms.y[hcs_n])<1e-10)
                secp.x[] = 2;
              else if (tmp_h.x[0,1]==nodata && tmp_h.x[0,-1]==nodata)
                secp.x[] = 2;
              else{
                
                secp.x[] = alphacs[hcs_n]!=10 ? (alphacs[hcs_n]-ms.x[hcs_n]*(alphac)/sign(mc.x[]))/ms.y[hcs_n] : (0.5-ms.x[hcs_n]*(alphac)/sign(mc.x[]))/ms.y[hcs_n];
              }
                
              if (fabs(secp.x[]) > 0.5+1e-10){
              
                secp.x[] = 2;

              }
              else{
                for (int i = -1 ;i <= 1;i+=2){
                  if (tmp_h.x[0,i]!=nodata && orientation(tmp_h.x[0,i]) == ori){

                    if ((secp.x[] < i && secp.x[] > ox) || (secp.x[] > i && secp.x[] < ox))
                      secp.x[] = 2;

                  }
                }

              }

            }
            
            if (add){

              if (secp.x[]!=2 && secp.x[]!=nodata){
                  
              double k_thetac = -mc.y[hcs_n]/mc.x[hcs_n];

                for (int i = -1 ;i <= 1;i+=2){
                  if (tmp_h.x[0,i]!=nodata && ((i - ox)*(i + ox - 2*secp.x[]))!=0 && orientation(tmp_h.x[0,i]) == ori){
              tmp_h.x[] = (-k_thetac*ox + h0 + k_thetac*i*pow(ox,2) - 2*h0*secp.x[]*i - tmp_h.x[0,i]*pow(ox,2) + 2*tmp_h.x[0,i]*secp.x[]*ox)/((i - ox)*(i + ox - 2*secp.x[]));
              }
              }
              }
              else if (ox != 0 && (((c[-1,1] + c[0,1] + c[1,1]) > 0 && ox>0) || ((c[-1,-1] + c[0,-1] + c[1,-1]) > 0 && ox<0)))
                tmp_h.x[] += ox*mc.y[hcs_n]/mc.x[hcs_n];


          }
          else {
            tmp_h.x[] = h0;
          }

        }
        

        }
        
        }
        }
      }

    }

  foreach (serial)
    for (int i = -2; i <= 2; i++)
      foreach_dimension(){
        if (mark[i]==4){
          if (fabs(height(tmp_h.x[i])) <= 3.5 && tmp_h.x[i] != nodata){
            if (tmp_h.x[]==nodata || ((tmp_h.x[]!=nodata) && (orientation(tmp_h.x[])==orientation(tmp_h.x[i]))))
	            tmp_h.x[] = tmp_h.x[i] + i;
          }
          else if (tmp_h.x[i]==nodata){
            tmp_h.x[]=50;
          }
        }
        
        
      }
  
  foreach(){
    for (int i = -2; i <= 2; i++)
      foreach_dimension(){
        if (tmp_h.x[i]==50)
          tmp_h.x[] = 50;
      }
  }

  foreach(){
    foreach_dimension()
      if (tmp_h.x[]==50)
        tmp_h.x[]=nodata;
  }
  


if (add){
  foreach(){
    foreach_dimension(){
      if (tmp_h.x[]==nodata && (mark[]==1 || mark[]==2 || mark[]==3)){
        for (int i=-1;i<=1;i+=2){
          if (mark[0,i]==4 && tmp_h.x[0,i]!=nodata){

            if (secp.x[0,i]!=nodata && secp.x[0,i]<1 && tmp_h.x[0,2*i]!=nodata && ((i*(i - 2*secp.x[0,i]))!=0) && (orientation(tmp_h.x[0,2*i]) == orientation(tmp_h.x[0,i]))){

int hcsn = vhcs_n.x[0,i]==nodata? 0 : vhcs_n.x[0,i];

double k_thetac = -mc.y[hcsn,i]/mc.x[hcsn,i];




tmp_h.x[] = -(tmp_h.x[0,2*i] - 4*tmp_h.x[0,i] + 2*i*k_thetac + 2*tmp_h.x[0,i]*i*(2*i + 2*secp.x[0,i]) - tmp_h.x[0,2*i]*i*(2*i + 2*secp.x[0,i]))/(i*(i - 2*secp.x[0,i]));


if (mark[]<100)
mark[]*=100;

            }

            else if (secp.x[0,i]==2 || secp.x[0,i]==nodata){
int hcsn = vhcs_n.x[0,i]==nodata? 0 : vhcs_n.x[0,i];
              if (mark[]==1)          
tmp_h.x[]=(mc.x[hcsn,i] == 0 ? tmp_h.x[0,i] : tmp_h.x[0,i] + i*mc.y[hcsn,i]/mc.x[hcsn,i]);

              else {
                int sy = (ms.x[])==0 ? 0 : sign(ms.x[]);
                int aa = (mc.x[])==0 ? 0 : sign(mc.x[]);
                int na = (1+i)/2;
                 double alpha = (-0.5+fs.y[0,na])*sy*aa ;
                double hy1 = aa<0 ? 20 - alpha : alpha;
                double dc = i;  dc/=2;
                
                hy1 += dc*mc.y[]/mc.x[];                
                double hy2 = (mc.x[hcsn,i] == 0 ? tmp_h.x[0,i] : tmp_h.x[0,i] + i*mc.y[hcsn,i]/mc.x[hcsn,i]);
                if (aa<0)
                  tmp_h.x[] = mark[]==2 ? max(hy1,hy2) : min(hy1,hy2);
                else
                  tmp_h.x[] = mark[]==2 ? min(hy1,hy2) : max(hy1,hy2);
              }

if (mark[]<100)
mark[]*=100;
              
          }



          }
        }
      }
    }
  }

  foreach (serial)
    if (mark[] != 4){
      for (int i = -2; i <= 2; i++)
        foreach_dimension()
	    if (mark[i]>=100 && tmp_h.x[i] != nodata){
        if (tmp_h.x[]==nodata || ((tmp_h.x[]!=nodata) && (orientation(tmp_h.x[])==orientation(tmp_h.x[i]))))
	        tmp_h.x[] = tmp_h.x[i] + i;
      }
  }
  
  foreach(){
    if (mark[] >= 100){
      int a=mark[]/100;
      mark[]=a;
    }

  }
}
  

}




trace
void recompute_h (scalar f, scalar tmp_c, scalar cs, scalar mark, scalar contact_angle, vector mc, vector tmp_h, vector oxyi, int add){

  vector ms[];
  scalar alphacs[];
  reconstruction_cs(cs, fs, ms, alphacs);
  reconstruction_mc_myc(tmp_c, mc);

    foreach(){

      if (mark[]==4||mark[]==5){
        coord mms;
        coord mmc;
        mms.x = ms.x[];
        mms.y = ms.y[];
        mmc.x = mc.x[];
        mmc.y = mc.y[];
        coord nnc = normal_contact(mms, mmc, contact_angle[]*pi/180);
        foreach_dimension()
          mc.x[] = nnc.x;
          
      }
      else if (mark[]==2 || mark[]==3){
        double me = 0;
        coord nnc;
        coord mms;
        mms.x = ms.x[];
        mms.y = ms.y[];
        double anglemy = contact_angle[];
        foreach_neighbor(1){
          if (mark[]==4||mark[]==5){
            me++;
            coord mmc;
            mmc.x = mc.x[];
            mmc.y = mc.y[];
            nnc = normal_contact(mms, mmc, anglemy*pi/180);
          }
        }
        if (me){
          foreach_dimension()
            mc.x[] = nnc.x;
        }
      }
      else{
        double me = 0;
        coord nnc;
        double anglemy = contact_angle[];
        foreach_neighbor(1){
          if (mark[]==4||mark[]==5){
            me++;
            coord mmc;
            coord mms;
            mms.x = ms.x[];
            mms.y = ms.y[];
            mmc.x = mc.x[];
            mmc.y = mc.y[];
            nnc = normal_contact(mms, mmc, anglemy*pi/180);
          }
        }
        if (me){
          foreach_dimension()
            mc.x[] = nnc.x;
        }
      }
    }

embed_h(f, cs, mc, ms, alphacs, contact_angle, mark, tmp_h, oxyi, add);
    
    
}




    
attribute {
  vector oxyi;
  vector nc;
  vector hnew1;
}