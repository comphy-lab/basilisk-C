
#include "fractions.h"

/**
Interface reconstruction within unstructured meshes (triangles, quadrilaterals, pentagons).*/


int myfacets (coord n, double alpha, coord pp[2])
{
  int i = 0;
  int ii = 0;
  coord p[4];

  for (double s = -0.5; s <= 0.5; s += 1.){
    foreach_dimension(){

	  if (fabs (n.y) > 1e-4) {
	    double a = (alpha - s*n.x)/n.y;
	    if (a >= -0.5 && a <= 0.5) {
          
	      p[i].x   = s;
	      p[i++].y = a;

    	}
      }

    }
  }

  if (i==2){
    pp[0].x = p[0].x;
    pp[0].y = p[0].y;
    pp[1].x = p[1].x;
    pp[1].y = p[1].y;
    ii = 2;
  }
  else if (i==3){
    pp[0].x = p[0].x;
    pp[0].y = p[0].y;
    ii = 2;
    double d1 = sqrt(pow(p[0].x-p[1].x,2) + pow(p[0].y-p[1].y,2));
    double d2 = sqrt(pow(p[0].x-p[2].x,2) + pow(p[0].y-p[2].y,2));
    if (d1 > d2){
      pp[1].x = p[1].x;
      pp[1].y = p[1].y;
    }
    else {
      pp[1].x = p[2].x;
      pp[1].y = p[2].y;
    }
  }
  else if (i==4){
    pp[0].x = p[0].x;
    pp[0].y = p[0].y;
    ii = 2;
    double d1 = sqrt(pow(p[0].x-p[1].x,2) + pow(p[0].y-p[1].y,2));
    double d2 = sqrt(pow(p[0].x-p[2].x,2) + pow(p[0].y-p[2].y,2));
    double d3 = sqrt(pow(p[0].x-p[3].x,2) + pow(p[0].y-p[3].y,2));
    if (d1 < d2){
        swap(coord,p[1],p[2]);
    }
    if (d1 < d3){
        swap(coord,p[1],p[3]);
    }
    if (d2 < d3){
        swap(coord,p[2],p[3]);
    }
    pp[1].x = p[1].x;
    pp[1].y = p[1].y;
  }

  return ii;
}


trace
double icx (coord a, coord b, double y){
    return a.x + (y - a.y)*(b.x - a.x)/(b.y - a.y);
}

trace
coord theta_n (coord nc, coord e){
    coord a;
    if(nc.x*e.y - nc.y*e.x > 0){
        a.x = (nc.x*e.x+nc.y*e.y)/(sqrt(sq(nc.x)+sq(nc.y))*sqrt(sq(e.x)+sq(e.y)));
        a.y = sqrt(1-sq(a.x));
    }

    else{
        a.x = (nc.x*e.x+nc.y*e.y)/(sqrt(sq(nc.x)+sq(nc.y))*sqrt(sq(e.x)+sq(e.y)));
        a.y = - sqrt(1-sq(a.x));
    }

    return a;
}

double yk (coord a, coord b){
    if (fabs(a.x - b.x) > 1e-10)
        return (a.x-b.x)/(a.y-b.y);
    else
        return 0;
        
}

double ysolve (double a, double b, double c, double min, double max){
    double y1 = (-b + sqrt(sq(b) - 4*a*c))/(2*a) ;
    double y2 = (-b - sqrt(sq(b) - 4*a*c))/(2*a) ;
    if (y1 >= 0 && y1 <= (max - min))
        return y1;
    else
        return y2;
}

trace
double mof_3points(double c, coord pp[3], coord nc, coord p_mof[2])
{
    coord e = {0,1};
    coord rot = theta_n(nc, e);

    double cosa = rot.x;
    double sina = rot.y;
    double cosb = cosa;
    double sinb = -sina;

    coord np[3];
    for (int i = 0; i <= 2; i++) {
        np[i].x = pp[i].x * cosa - pp[i].y * sina;
        np[i].y = pp[i].x * sina + pp[i].y * cosa;
    }

    coord sp[3];
    for(int k=0; k<=2; k++){
        foreach_dimension()
        sp[k].x = np[k].x;
    }
    if(sp[0].y > sp[1].y)
    swap(coord, sp[0], sp[1]);
    if(sp[0].y > sp[2].y)
    swap(coord, sp[0], sp[2]);
    if(sp[1].y > sp[2].y)
    swap(coord, sp[1], sp[2]);

    double shift = sp[0].y;
    for(int k=0; k<=2; k++){
        sp[k].y = sp[k].y - shift;
    }

    double d = fabs(icx(sp[0], sp[2], sp[1].y) - sp[1].x);
    double c1 = fabs(d*(sp[1].y - sp[0].y)/2);
    double c2 = c1 + fabs(d*(sp[2].y - sp[1].y)/2);

    coord inp[2] = {{0,0},{0,0}};
    if(c >= 0 && c < c1 - 1e-10){
        inp[0].y = sqrt(sq(sp[0].y) + c*(sq(sp[1].y) - sq(sp[0].y))/c1);
        inp[1].y = inp[0].y;
        inp[0].x = icx(sp[0], sp[1], inp[0].y);
        inp[1].x = icx(sp[0], sp[2], inp[1].y);
    }
    else if(c > c1 + 1e-10 && c < c2 - 1e-10){
        double tmp_c = -(c-c1);
        double k12 = yk(sp[1],sp[2]);
        double k02 = yk(sp[0],sp[2]);
        if (sp[1].x > icx(sp[0], sp[2], sp[1].y))
            inp[0].y = ysolve((k12-k02)/2, d, tmp_c, sp[1].y, sp[2].y) + sp[1].y;
        else
            inp[0].y = ysolve((k02-k12)/2, d, tmp_c, sp[1].y, sp[2].y) + sp[1].y;
        inp[1].y = inp[0].y;
        inp[0].x = icx(sp[1], sp[2], inp[0].y);
        inp[1].x = icx(sp[0], sp[2], inp[1].y);
    }
    else if(fabs(c - c1) <= 1e-10){
        inp[0].y = sp[1].y;
        inp[1].y = inp[0].y;
        inp[0].x = sp[1].x;
        inp[1].x = icx(sp[0], sp[2], inp[1].y);
    }
    else if(fabs(c - c2) <= 1e-10){
        inp[0].y = sp[2].y;
        inp[1].y = inp[0].y;
        inp[0].x = sp[2].x;
        inp[1].x = inp[0].x;
    }

    inp[0].y = inp[0].y + shift;
    inp[1].y = inp[1].y + shift;

    for(int j=0; j<=1; j++){
            p_mof[j].x = inp[j].x*cosb - inp[j].y*sinb;
            p_mof[j].y = inp[j].x*sinb + inp[j].y*cosb;
    }
    double alphac = nc.x*p_mof[0].x + nc.y*p_mof[0].y;
    return alphac;
}
trace
double mof_4points(double c, coord pp[4], coord nc, coord p_mof[2])
{

    coord e = {0,1};
    coord rot = theta_n(nc, e);

    double cosa = rot.x;
    double sina = rot.y;
    double cosb = cosa;
    double sinb = -sina;

    coord np[4];
    for (int i = 0; i <= 3; i++) {
        np[i].x = pp[i].x * cosa - pp[i].y * sina;
        np[i].y = pp[i].x * sina + pp[i].y * cosa;
    }
    coord sp[4];
    int tmp_p[4];
    int number_p[4];
    int invnumber_p[4];

    for (int k = 0; k <= 3; k++){
        foreach_dimension()
        sp[k].x = np[k].x;
    }

    for (int i = 0; i <= 3; i++) {
        tmp_p[i] = i;
    }
    for (int i = 0; i <= 2; i++) {
        for (int j = 0; j <= 2-i; j++) {
            if(sp[i].y > sp[i+j+1].y){
            swap(coord, sp[i], sp[i+j+1]);
            swap(int, tmp_p[i], tmp_p[i+j+1]);
            }

        }
    }


    for (int k = 0; k <= 3; k++){
        number_p[k] = tmp_p[k];
    }



    for (int i = 0; i <= 3; i++) {
        invnumber_p[i] = i;
    }
    for (int i = 0; i <= 2; i++) {
        for (int j = 0; j <= 2-i; j++) {
            if(tmp_p[i] > tmp_p[i+j+1]){
            swap(int, tmp_p[i], tmp_p[i+j+1]);
            swap(int, invnumber_p[i], invnumber_p[i+j+1]);
            }
        }
    }
    

    double shift = sp[0].y;
    for(int k = 0; k <= 3; k++){
        sp[k].y = sp[k].y - shift;
    }

    int sort1  = number_p[2] + 1 > 3 ? number_p[2] - 3 : number_p[2] + 1;
    int sort_1 = number_p[2] - 1 < 0 ? number_p[2] + 3 : number_p[2] - 1;

    double d1, d2, c1, c2, c3; 
    coord inp[2] = {{0,0},{0,0}};

    if (invnumber_p[sort1] == 0 || invnumber_p[sort_1] == 0) {

        d1 = fabs(icx(sp[0], sp[2], sp[1].y) - sp[1].x);
        d2 = fabs(icx(sp[1], sp[3], sp[2].y) - sp[2].x);  
        c1 = fabs(d1 * (sp[1].y - sp[0].y)/2);   
        c2 = c1 + fabs((d1 + d2) * (sp[2].y - sp[1].y)/2);
        c3 = c2 + fabs(d2 * (sp[3].y - sp[2].y)/2);

        if(c >= 0 && c < c1 - 1e-10){
            inp[0].y = sqrt(sq(sp[0].y) + c*(sq(sp[1].y) - sq(sp[0].y))/c1);
            inp[1].y = inp[0].y;
            inp[0].x = icx(sp[0], sp[1], inp[0].y);
            inp[1].x = icx(sp[0], sp[2], inp[1].y);
        }

        else if(c > c1 + 1e-10 && c < c2 - 1e-10){
            if (fabs(d1 - d2) < 1e-10){
                inp[0].y = sp[1].y + (c - c1)*(sp[2].y - sp[1].y)/(c2 -c1);
            }
            else {
                double tmp_c = -(c-c1);
                double k13 = yk(sp[1],sp[3]);
                double k02 = yk(sp[0],sp[2]);
                if (sp[1].x > icx(sp[0], sp[2], sp[1].y))
                    inp[0].y = ysolve((k13-k02)/2, d1, tmp_c, sp[1].y, sp[2].y) + sp[1].y;
                else
                    inp[0].y = ysolve((k02-k13)/2, d1, tmp_c, sp[1].y, sp[2].y) + sp[1].y;
            }
            inp[1].y = inp[0].y;
            inp[0].x = icx(sp[1], sp[3], inp[0].y);
            inp[1].x = icx(sp[0], sp[2], inp[1].y);
        }
        else if(c >= c2 + 1e-10 && c < c3 - 1e-10){
            double tmp_c = -(c-c2);
            double k13 = yk(sp[1],sp[3]);
            double k23 = yk(sp[2],sp[3]);
            if (icx(sp[1], sp[3], sp[2].y) > sp[2].x)
                inp[0].y = ysolve((k13-k23)/2, d2, tmp_c, sp[2].y, sp[3].y) + sp[2].y;
            else
                inp[0].y = ysolve((k23-k13)/2, d2, tmp_c, sp[2].y, sp[3].y) + sp[2].y;
            inp[1].y = inp[0].y;
            inp[0].x = icx(sp[1], sp[3], inp[0].y);
            inp[1].x = icx(sp[2], sp[3], inp[1].y);
        }
        else if(fabs(c - c1) <= 1e-10){
            inp[0].y = sp[1].y;
            inp[1].y = inp[0].y;
            inp[0].x = sp[1].x;
            inp[1].x = icx(sp[0], sp[2], inp[1].y);
        }
        else if(fabs(c - c2) <= 1e-10){
            inp[0].y = sp[2].y;
            inp[1].y = inp[0].y;
            inp[0].x = sp[2].x;
            inp[1].x = icx(sp[1], sp[3], inp[1].y);
        }
        else if(fabs(c - c3) <= 1e-10){
            inp[0].y = sp[3].y;
            inp[1].y = inp[0].y;
            inp[0].x = sp[3].x;
            inp[1].x = inp[0].x;
        }
    }     


    else {

        d1 = fabs(icx(sp[0], sp[3], sp[1].y) - sp[1].x);
        d2 = fabs(icx(sp[0], sp[3], sp[2].y) - sp[2].x);  
        c1 = fabs(d1 * (sp[1].y - sp[0].y)/2);   
        c2 = c1 + fabs((d1 + d2) * (sp[2].y - sp[1].y)/2);
        c3 = c2 + fabs(d2 * (sp[3].y - sp[2].y)/2);  


        if(c >= 0 && c < c1 - 1e-10){
            inp[0].y = sqrt(sq(sp[0].y) + c*(sq(sp[1].y) - sq(sp[0].y))/c1);
            inp[1].y = inp[0].y;
            inp[0].x = icx(sp[0], sp[1], inp[0].y);
            inp[1].x = icx(sp[0], sp[3], inp[1].y);

        }
        else if(c > c1 + 1e-10 && c < c2 - 1e-10){
            if (fabs(d1 - d2) < 1e-10){
                inp[0].y = sp[1].y + (c - c1)*(sp[2].y - sp[1].y)/(c2 -c1);
            }
            else {
                double tmp_c = -(c-c1);
                double k03 = yk(sp[0],sp[3]);
                double k12 = yk(sp[1],sp[2]);
                if (icx(sp[0], sp[3], sp[1].y) > sp[1].x)
                    inp[0].y = ysolve((k03-k12)/2, d1, tmp_c, sp[1].y, sp[2].y) + sp[1].y;
                else
                    inp[0].y = ysolve((k12-k03)/2, d1, tmp_c, sp[1].y, sp[2].y) + sp[1].y;
                
            }
            inp[1].y = inp[0].y;
            inp[0].x = icx(sp[1], sp[2], inp[0].y);
            inp[1].x = icx(sp[0], sp[3], inp[1].y);

        }
        else if(c >= c2 + 1e-10 && c < c3 - 1e-10){
            double tmp_c = -(c-c2);
            double k03 = yk(sp[0],sp[3]);
            double k23 = yk(sp[2],sp[3]);
            if (icx(sp[0], sp[3], sp[2].y) > sp[2].x)
                inp[0].y = ysolve((k03-k23)/2, d2, tmp_c, sp[2].y, sp[3].y) + sp[2].y;
            else
                inp[0].y = ysolve((k23-k03)/2, d2, tmp_c, sp[2].y, sp[3].y) + sp[2].y;
            inp[1].y = inp[0].y;
            inp[0].x = icx(sp[2], sp[3], inp[0].y);
            inp[1].x = icx(sp[0], sp[3], inp[1].y);
        }
        else if(fabs(c - c1) <= 1e-10){
            inp[0].y = sp[1].y;
            inp[1].y = inp[0].y;
            inp[0].x = sp[1].x;
            inp[1].x = icx(sp[0], sp[3], inp[1].y);
        }
        else if(fabs(c - c2) <= 1e-10){
            inp[0].y = sp[2].y;
            inp[1].y = inp[0].y;
            inp[0].x = sp[2].x;
            inp[1].x = icx(sp[0], sp[3], inp[1].y);
        }
        else if(fabs(c - c3) <= 1e-10){
            inp[0].y = sp[3].y;
            inp[1].y = inp[0].y;
            inp[0].x = sp[3].x;
            inp[1].x = inp[0].x;
        }
    }


    inp[0].y = inp[0].y + shift;
    inp[1].y = inp[1].y + shift;

    for(int j=0; j<=1; j++){
            p_mof[j].x = inp[j].x*cosb - inp[j].y*sinb;
            p_mof[j].y = inp[j].x*sinb + inp[j].y*cosb;
    }
    double alphac = nc.x*p_mof[0].x + nc.y*p_mof[0].y;
    return alphac;

}
trace
double mof_5points(double c, coord pp[5], coord nc, coord p_mof[2])
{
    coord e = {0,1};
    coord rot = theta_n(nc, e);

    double cosa = rot.x;
    double sina = rot.y;
    double cosb = cosa;
    double sinb = -sina;

    coord np[5];
    for (int i = 0; i <= 4; i++) {
        np[i].x = pp[i].x * cosa - pp[i].y * sina;
        np[i].y = pp[i].x * sina + pp[i].y * cosa;
    }
    coord sp[5];

    for (int k = 0; k <= 4; k++){
        foreach_dimension()
        sp[k].x = np[k].x;
    }

    for (int i = 0; i <= 3; i++) {
        for (int j = 0; j <= 3-i; j++) {
            if(sp[i].y > sp[i+j+1].y)
            swap(coord, sp[i], sp[i+j+1]);

        }
    }
    
   

    double shift = sp[0].y;
    for(int k = 0; k <= 4; k++){
        sp[k].y = sp[k].y - shift;
    }

    coord o_xy[1] = {{0, 0}};
    for (int k = 0; k <= 4; k++){
        foreach_dimension()
        o_xy[0].x += sp[k].x;
    }
    foreach_dimension()
    o_xy[0].x /= 5;

    double l_angle[5];
    int ta[5];
    for (int k = 0; k <= 4; k++){
        l_angle[k] = atan2(sp[k].y - o_xy[0].y, sp[k].x - o_xy[0].x);
        if (l_angle[k] < l_angle[0])
            l_angle[k] += 2*pi;
    }


    for (int i = 0; i <= 4; i++) {
        ta[i] = i;
    }
    for (int i = 0; i <= 3; i++) {
        for (int j = 0; j <= 3-i; j++) {
            if(l_angle[i] > l_angle[i+j+1]){
            swap(double, l_angle[i], l_angle[i+j+1]);
            swap(int, ta[i], ta[i+j+1]);
            }
        }
    }
    
    


    int case_n = 0;
    double d1, d2, d3, c1, c2, c3, c4; 
    coord inp[2] = {{0,0},{0,0}};

    
    if ((ta[1] == 1 && ta[2] == 2 && ta[3] == 4) || (ta[4] == 1 && ta[3] == 2 && ta[2] == 4)){
        case_n = 1;
    }
    else if ((ta[1] == 1 && ta[2] == 3) || (ta[4] == 1 && ta[3] == 3)){
        case_n = 2;
    }
    else if ((ta[1] == 1 && ta[2] == 4) || (ta[4] == 1 && ta[3] == 4)){
        case_n = 3;
    }
    else if ((ta[1] == 1 && ta[4] == 4) || (ta[1] == 4 && ta[4] == 1)){
        case_n = 4;
    }


    if (case_n == 1){
        d1 = fabs(icx(sp[0], sp[3], sp[1].y) - sp[1].x);
        d2 = fabs(icx(sp[0], sp[3], sp[2].y) - sp[2].x);
        d3 = fabs(icx(sp[2], sp[4], sp[3].y) - sp[3].x);
        c1 = fabs(d1 * (sp[1].y - sp[0].y)/2);   
        c2 = c1 + fabs((d1 + d2) * (sp[2].y - sp[1].y)/2);
        c3 = c2 + fabs((d2 + d3) * (sp[3].y - sp[2].y)/2);
        c4 = c3 + fabs(d3 * (sp[4].y - sp[3].y)/2);

        if(c >= 0 && c < c1 - 1e-10){
            inp[0].y = sqrt(sq(sp[0].y) + c*(sq(sp[1].y) - sq(sp[0].y))/c1);
            inp[1].y = inp[0].y;
            inp[0].x = icx(sp[0], sp[1], inp[0].y);
            inp[1].x = icx(sp[0], sp[3], inp[1].y);
        }

        else if(c >= c1 + 1e-10 && c < c2 - 1e-10){
            if (fabs(d1 - d2) < 1e-10){
                inp[0].y = sp[1].y + (c - c1)*(sp[2].y - sp[1].y)/(c2 - c1);
            }
            else {
                double tmp_c = -(c-c1);
                double k03 = yk(sp[0],sp[3]);
                double k12 = yk(sp[1],sp[2]);
                if (icx(sp[0], sp[3], sp[1].y) > sp[1].x)
                    inp[0].y = ysolve((k03-k12)/2, d1, tmp_c, sp[1].y, sp[2].y) + sp[1].y;
                else
                    inp[0].y = ysolve((k12-k03)/2, d1, tmp_c, sp[1].y, sp[2].y) + sp[1].y;
            }
            inp[1].y = inp[0].y;
            inp[0].x = icx(sp[1], sp[2], inp[0].y);
            inp[1].x = icx(sp[0], sp[3], inp[1].y);
        }

        else if(c >= c2 + 1e-10 && c < c3 - 1e-10){
            if (fabs(d2 - d3) < 1e-10){
                inp[0].y = sp[2].y + (c - c2)*(sp[3].y - sp[2].y)/(c3 - c2);
            }
            else {
                double tmp_c = -(c-c2);
                double k03 = yk(sp[0],sp[3]);
                double k24 = yk(sp[2],sp[4]);
                if (icx(sp[0], sp[3], sp[2].y) > sp[2].x)
                    inp[0].y = ysolve((k03-k24)/2, d2, tmp_c, sp[2].y, sp[3].y) + sp[2].y;
                else
                    inp[0].y = ysolve((k24-k03)/2, d2, tmp_c, sp[2].y, sp[3].y) + sp[2].y;
            }
            inp[1].y = inp[0].y;
            inp[0].x = icx(sp[2], sp[4], inp[0].y);
            inp[1].x = icx(sp[0], sp[3], inp[1].y);
        }

        else if(c >= c3 + 1e-10 && c < c4 - 1e-10){
                double tmp_c = -(c-c3);
                double k34 = yk(sp[3],sp[4]);
                double k24 = yk(sp[2],sp[4]);
                if (sp[3].x > icx(sp[2], sp[4], sp[3].y))
                    inp[0].y = ysolve((k34-k24)/2, d3, tmp_c, sp[3].y, sp[4].y) + sp[3].y;
                else
                    inp[0].y = ysolve((k24-k34)/2, d3, tmp_c, sp[3].y, sp[4].y) + sp[3].y;
            inp[1].y = inp[0].y;
            inp[0].x = icx(sp[2], sp[4], inp[0].y);
            inp[1].x = icx(sp[3], sp[4], inp[1].y);
        }

        else if(fabs(c - c1) <= 1e-10){
            inp[0].y = sp[1].y;
            inp[1].y = inp[0].y;
            inp[0].x = sp[1].x;
            inp[1].x = icx(sp[0], sp[3], inp[1].y);
        }

        else if(fabs(c - c2) <= 1e-10){
            inp[0].y = sp[2].y;
            inp[1].y = inp[0].y;
            inp[0].x = sp[2].x;
            inp[1].x = icx(sp[0], sp[3], inp[1].y);
        }
        
        else if(fabs(c - c3) <= 1e-10){
            inp[0].y = sp[3].y;
            inp[1].y = inp[0].y;
            inp[0].x = sp[3].x;
            inp[1].x = icx(sp[2], sp[4], inp[1].y);
        }
        
        else if(fabs(c - c4) <= 1e-10){
            inp[0].y = sp[4].y;
            inp[1].y = inp[0].y;
            inp[0].x = sp[4].x;
            inp[1].x = inp[0].x;
        }

        }

    else if (case_n == 2){
        d1 = fabs(icx(sp[0], sp[2], sp[1].y) - sp[1].x);
        d2 = fabs(icx(sp[1], sp[3], sp[2].y) - sp[2].x);
        d3 = fabs(icx(sp[2], sp[4], sp[3].y) - sp[3].x);
        c1 = fabs(d1 * (sp[1].y - sp[0].y)/2);   
        c2 = c1 + fabs((d1 + d2) * (sp[2].y - sp[1].y)/2);
        c3 = c2 + fabs((d2 + d3) * (sp[3].y - sp[2].y)/2);
        c4 = c3 + fabs(d3 * (sp[4].y - sp[3].y)/2);

        if(c >= 0 && c < c1 - 1e-10){
            inp[0].y = sqrt(sq(sp[0].y) + c*(sq(sp[1].y) - sq(sp[0].y))/c1);
            inp[1].y = inp[0].y;
            inp[0].x = icx(sp[0], sp[1], inp[0].y);
            inp[1].x = icx(sp[0], sp[2], inp[1].y);
        }

        else if(c >= c1 + 1e-10 && c < c2 - 1e-10){
            if (fabs(d1 - d2) < 1e-10){
                inp[0].y = sp[1].y + (c - c1)*(sp[2].y - sp[1].y)/(c2 - c1);
            }
            else {
                double tmp_c = -(c-c1);
                double k02 = yk(sp[0],sp[2]);
                double k13 = yk(sp[1],sp[3]);
                if (icx(sp[0], sp[2], sp[1].y) > sp[1].x)
                    inp[0].y = ysolve((k02-k13)/2, d1, tmp_c, sp[1].y, sp[2].y) + sp[1].y;
                else
                    inp[0].y = ysolve((k13-k02)/2, d1, tmp_c, sp[1].y, sp[2].y) + sp[1].y;
            }
            inp[1].y = inp[0].y;
            inp[0].x = icx(sp[0], sp[2], inp[0].y);
            inp[1].x = icx(sp[1], sp[3], inp[1].y);
        }

        else if(c >= c2 + 1e-10 && c < c3 - 1e-10){
            if (fabs(d2 - d3) < 1e-10){
                inp[0].y = sp[2].y + (c - c2)*(sp[3].y - sp[2].y)/(c3 - c2);
            }
            else {
                double tmp_c = -(c-c2);
                double k24 = yk(sp[2],sp[4]);
                double k13 = yk(sp[1],sp[3]);
                if (sp[2].x > icx(sp[1], sp[3], sp[2].y))
                    inp[0].y = ysolve((k24-k13)/2, d2, tmp_c, sp[2].y, sp[3].y) + sp[2].y;
                else
                    inp[0].y = ysolve((k13-k24)/2, d2, tmp_c, sp[2].y, sp[3].y) + sp[2].y;
            }
            inp[1].y = inp[0].y;
            inp[0].x = icx(sp[1], sp[3], inp[0].y);
            inp[1].x = icx(sp[2], sp[4], inp[1].y);
        }

        else if(c >= c3 + 1e-10 && c < c4 - 1e-10){
                double tmp_c = -(c-c3);
                double k24 = yk(sp[2],sp[4]);
                double k34 = yk(sp[3],sp[4]);
                if (icx(sp[2], sp[4], sp[3].y) > sp[3].x)
                    inp[0].y = ysolve((k24-k34)/2, d3, tmp_c, sp[3].y, sp[4].y) + sp[3].y;
                else
                    inp[0].y = ysolve((k34-k24)/2, d3, tmp_c, sp[3].y, sp[4].y) + sp[3].y;
            inp[1].y = inp[0].y;
            inp[0].x = icx(sp[3], sp[4], inp[0].y);
            inp[1].x = icx(sp[2], sp[4], inp[1].y);
        }

        else if(fabs(c - c1) <= 1e-10){
            inp[0].y = sp[1].y;
            inp[1].y = inp[0].y;
            inp[0].x = sp[1].x;
            inp[1].x = icx(sp[0], sp[2], inp[1].y);
        }

        else if(fabs(c - c2) <= 1e-10){
            inp[0].y = sp[2].y;
            inp[1].y = inp[0].y;
            inp[0].x = sp[2].x;
            inp[1].x = icx(sp[1], sp[3], inp[1].y);
        }
        
        else if(fabs(c - c3) <= 1e-10){
            inp[0].y = sp[3].y;
            inp[1].y = inp[0].y;
            inp[0].x = sp[3].x;
            inp[1].x = icx(sp[2], sp[4], inp[1].y);
        }
        
        else if(fabs(c - c4) <= 1e-10){
            inp[0].y = sp[4].y;
            inp[1].y = inp[0].y;
            inp[0].x = sp[4].x;
            inp[1].x = inp[0].x;
        }
        }
    
    else if (case_n == 3){

        d1 = fabs(icx(sp[0], sp[2], sp[1].y) - sp[1].x);
        d2 = fabs(icx(sp[1], sp[4], sp[2].y) - sp[2].x);
        d3 = fabs(icx(sp[1], sp[4], sp[3].y) - sp[3].x);
        c1 = fabs(d1 * (sp[1].y - sp[0].y)/2);   
        c2 = c1 + fabs((d1 + d2) * (sp[2].y - sp[1].y)/2);
        c3 = c2 + fabs((d2 + d3) * (sp[3].y - sp[2].y)/2);
        c4 = c3 + fabs(d3 * (sp[4].y - sp[3].y)/2);
  
        if(c >= 0 && c < c1 - 1e-10){

            inp[0].y = sqrt(sq(sp[0].y) + c*(sq(sp[1].y) - sq(sp[0].y))/c1);
            inp[1].y = inp[0].y;
            inp[0].x = icx(sp[0], sp[1], inp[0].y);
            inp[1].x = icx(sp[0], sp[2], inp[1].y);
        }

        else if(c >= c1 + 1e-10 && c < c2 - 1e-10){

            if (fabs(d1 - d2) < 1e-10){
                inp[0].y = sp[1].y + (c - c1)*(sp[2].y - sp[1].y)/(c2 - c1);
            }
            else {
                double tmp_c = -(c-c1);
                double k02 = yk(sp[0],sp[2]);
                double k14 = yk(sp[1],sp[4]);
                if (icx(sp[0], sp[2], sp[1].y) > sp[1].x)
                    inp[0].y = ysolve((k02-k14)/2, d1, tmp_c, sp[1].y, sp[2].y) + sp[1].y;
                else
                    inp[0].y = ysolve((k14-k02)/2, d1, tmp_c, sp[1].y, sp[2].y) + sp[1].y;
            }
            inp[1].y = inp[0].y;
            inp[0].x = icx(sp[1], sp[4], inp[0].y);
            inp[1].x = icx(sp[0], sp[2], inp[1].y);
        }

        else if(c >= c2 + 1e-10 && c < c3 - 1e-10){

            if (fabs(d2 - d3) < 1e-10){
                inp[0].y = sp[2].y + (c - c2)*(sp[3].y - sp[2].y)/(c3 - c2);
            }
            else {
                double tmp_c = -(c-c2);
                double k23 = yk(sp[2],sp[3]);
                double k14 = yk(sp[1],sp[4]);
                if (sp[2].x > icx(sp[1], sp[4], sp[2].y))
                    inp[0].y = ysolve((k23-k14)/2, d2, tmp_c, sp[2].y, sp[3].y) + sp[2].y;
                else
                    inp[0].y = ysolve((k14-k23)/2, d2, tmp_c, sp[2].y, sp[3].y) + sp[2].y;
            }
            inp[1].y = inp[0].y;
            inp[0].x = icx(sp[2], sp[3], inp[0].y);
            inp[1].x = icx(sp[1], sp[4], inp[1].y);
        }

        else if(c >= c3 + 1e-10 && c < c4 - 1e-10){

                double tmp_c = -(c-c3);
                double k34 = yk(sp[3],sp[4]);
                double k14 = yk(sp[1],sp[4]);
                if (sp[3].x > icx(sp[1], sp[4], sp[3].y))
                    inp[0].y = ysolve((k34-k14)/2, d3, tmp_c, sp[3].y, sp[4].y) + sp[3].y;
                else
                    inp[0].y = ysolve((k14-k34)/2, d3, tmp_c, sp[3].y, sp[4].y) + sp[3].y;
            inp[1].y = inp[0].y;
            inp[0].x = icx(sp[3], sp[4], inp[0].y);
            inp[1].x = icx(sp[1], sp[4], inp[1].y);
        }

        else if(fabs(c - c1) <= 1e-10){

            inp[0].y = sp[1].y;
            inp[1].y = inp[0].y;
            inp[0].x = sp[1].x;
            inp[1].x = icx(sp[0], sp[2], inp[1].y);
        }

        else if(fabs(c - c2) <= 1e-10){

            inp[0].y = sp[2].y;
            inp[1].y = inp[0].y;
            inp[0].x = sp[2].x;
            inp[1].x = icx(sp[1], sp[4], inp[1].y);
        }
        
        else if(fabs(c - c3) <= 1e-10){

            inp[0].y = sp[3].y;
            inp[1].y = inp[0].y;
            inp[0].x = sp[3].x;
            inp[1].x = icx(sp[1], sp[4], inp[1].y);
        }
        
        else if(fabs(c - c4) <= 1e-10){

            inp[0].y = sp[4].y;
            inp[1].y = inp[0].y;
            inp[0].x = sp[4].x;
            inp[1].x = inp[0].x;
        }
    }

    else if (case_n == 4){
        d1 = fabs(icx(sp[0], sp[4], sp[1].y) - sp[1].x);
        d2 = fabs(icx(sp[0], sp[4], sp[2].y) - sp[2].x);
        d3 = fabs(icx(sp[0], sp[4], sp[3].y) - sp[3].x);
        c1 = fabs(d1 * (sp[1].y - sp[0].y)/2);   
        c2 = c1 + fabs((d1 + d2) * (sp[2].y - sp[1].y)/2);
        c3 = c2 + fabs((d2 + d3) * (sp[3].y - sp[2].y)/2);
        c4 = c3 + fabs(d3 * (sp[4].y - sp[3].y)/2);

        if(c >= 0 && c < c1 - 1e-10){
            inp[0].y = sqrt(sq(sp[0].y) + c*(sq(sp[1].y) - sq(sp[0].y))/c1);
            inp[1].y = inp[0].y;
            inp[0].x = icx(sp[0], sp[1], inp[0].y);
            inp[1].x = icx(sp[0], sp[4], inp[1].y);
        }

        else if(c >= c1 + 1e-10 && c < c2 - 1e-10){
            if (fabs(d1 - d2) < 1e-10){
                inp[0].y = sp[1].y + (c - c1)*(sp[2].y - sp[1].y)/(c2 - c1);
            }
            else {
                double tmp_c = -(c-c1);
                double k12 = yk(sp[1],sp[2]);
                double k04 = yk(sp[0],sp[4]);
                if (sp[1].x > icx(sp[0], sp[4], sp[1].y))
                    inp[0].y = ysolve((k12-k04)/2, d1, tmp_c, sp[1].y, sp[2].y) + sp[1].y;
                else
                    inp[0].y = ysolve((k04-k12)/2, d1, tmp_c, sp[1].y, sp[2].y) + sp[1].y;
            }
            inp[1].y = inp[0].y;
            inp[0].x = icx(sp[0], sp[4], inp[0].y);
            inp[1].x = icx(sp[1], sp[2], inp[1].y);
        }

        else if(c >= c2 + 1e-10 && c < c3 - 1e-10){
            if (fabs(d2 - d3) < 1e-10){
                inp[0].y = sp[2].y + (c - c2)*(sp[3].y - sp[2].y)/(c3 - c2);
            }
            else {
                double tmp_c = -(c-c2);
                double k23 = yk(sp[2],sp[3]);
                double k04 = yk(sp[0],sp[4]);
                if (sp[2].x > icx(sp[0], sp[4], sp[2].y))
                    inp[0].y = ysolve((k23-k04)/2, d2, tmp_c, sp[2].y, sp[3].y) + sp[2].y;
                else
                    inp[0].y = ysolve((k04-k23)/2, d2, tmp_c, sp[2].y, sp[3].y) + sp[2].y;
            }
            inp[1].y = inp[0].y;
            inp[0].x = icx(sp[2], sp[3], inp[0].y);
            inp[1].x = icx(sp[0], sp[4], inp[1].y);
        }

        else if(c >= c3 + 1e-10 && c < c4 - 1e-10){
                double tmp_c = -(c-c3);
                double k34 = yk(sp[3],sp[4]);
                double k04 = yk(sp[0],sp[4]);
                if (sp[3].x > icx(sp[0], sp[4], sp[3].y))
                    inp[0].y = ysolve((k34-k04)/2, d3, tmp_c, sp[3].y, sp[4].y) + sp[3].y;
                else
                    inp[0].y = ysolve((k04-k34)/2, d3, tmp_c, sp[3].y, sp[4].y) + sp[3].y;
            inp[1].y = inp[0].y;
            inp[0].x = icx(sp[3], sp[4], inp[0].y);
            inp[1].x = icx(sp[0], sp[4], inp[1].y);
        }

        else if(fabs(c - c1) <= 1e-10){
            inp[0].y = sp[1].y;
            inp[1].y = inp[0].y;
            inp[0].x = sp[1].x;
            inp[1].x = icx(sp[0], sp[4], inp[1].y);
        }

        else if(fabs(c - c2) <= 1e-10){
            inp[0].y = sp[2].y;
            inp[1].y = inp[0].y;
            inp[0].x = sp[2].x;
            inp[1].x = icx(sp[0], sp[4], inp[1].y);
        }
        
        else if(fabs(c - c3) <= 1e-10){
            inp[0].y = sp[3].y;
            inp[1].y = inp[0].y;
            inp[0].x = sp[3].x;
            inp[1].x = icx(sp[0], sp[4], inp[1].y);
        }
        
        else if(fabs(c - c4) <= 1e-10){
            inp[0].y = sp[4].y;
            inp[1].y = inp[0].y;
            inp[0].x = sp[4].x;
            inp[1].x = inp[0].x;
        }
        }

        inp[0].y = inp[0].y + shift;
        inp[1].y = inp[1].y + shift;
        for(int j=0; j<=1; j++){
            p_mof[j].x = inp[j].x*cosb - inp[j].y*sinb;
            p_mof[j].y = inp[j].x*sinb + inp[j].y*cosb;
    }
    double alphac = nc.x*p_mof[0].x + nc.y*p_mof[0].y;

    return alphac;
}

/**
Copy from Popinet/contact sandbox.*/

trace
static inline coord normal_contact (coord ns, coord nf, double angle)
{
  assert (dimension == 2);
  coord n;
  if (- ns.x*nf.y + ns.y*nf.x > 0) {
    n.x = - ns.x*cos(angle) + ns.y*sin(angle);
    n.y = - ns.x*sin(angle) - ns.y*cos(angle);
  }
  else {
    n.x = - ns.x*cos(angle) - ns.y*sin(angle);
    n.y =   ns.x*sin(angle) - ns.y*cos(angle);
  }

  return n;
}



trace
double polygon_alpha (double c, coord mc, coord ms, double alphacs, 
                    coord p_mof[2], coord pp[5])
{

  if (mc.x == 0 && mc.y == 0){
    p_mof[0].x = 10, p_mof[0].y = 10;
    p_mof[1].x = 10, p_mof[1].y = 10;
    pp[0].x = 10, pp[0].y = 10;
    pp[1].x = 10, pp[1].y = 10;
    pp[2].x = 10, pp[2].y = 10;
    pp[3].x = 10, pp[3].y = 10;
    pp[4].x = 10, pp[4].y = 10;   
    return 0;
  }

  else{

    if ((ms.x != 0 || ms.y != 0) && alphacs != 10){

        coord pcs[2];
        int nfacets;
        nfacets = myfacets (ms, alphacs, pcs);
if(nfacets == 2){

            int aorient[4];
            int a = 0, an = 0;
            coord vex[4];
            vex[0].x = 0.5, vex[0].y = 0.5, vex[1].x = -0.5, vex[1].y = 0.5;
            vex[2].x = -0.5, vex[2].y = -0.5, vex[3].x = 0.5, vex[3].y = -0.5;  
            for(int i = 0; i <= 3; i++){

                if(ms.x * vex[i].x + ms.y * vex[i].y - alphacs <= 0){
                  if ((fabs(vex[i].x-pcs[0].x)<1e-12 && fabs(vex[i].y-pcs[0].y)<1e-12) || ((fabs(vex[i].x-pcs[1].x)<1e-12 && fabs(vex[i].y-pcs[1].y)<1e-12))){
                    aorient[a++] = 0; 

                  }
                  else{
                    aorient[a++] = 1; 
                    an++;
                  }
                }
                else aorient[a++] = 0; 
            }
            

            if(nfacets + an == 3){

                foreach_dimension(){
                    pp[0].x = pcs[0].x;
                    pp[1].x = pcs[1].x;
                }
                for(int j = 0; j <= 3; j++){
                    if(aorient[j] == 1){
                        foreach_dimension(){
                            pp[2].x = vex[j].x;
                        }
                    } 
                }


                double alphac = mof_3points (c, pp, mc, p_mof);
                return alphac;


            }

            else if(nfacets + an == 4){

                foreach_dimension(){
                    pp[0].x = pcs[0].x;
                    pp[1].x = pcs[1].x;
                }

                if (fabs(fabs(pp[1].x)-0.5)<1e-12 && fabs(fabs(pp[1].y)-0.5)<1e-12){
                    coord change = {0,0};
                    change.x = pp[0].x;
                    change.y = pp[0].y;
                    pp[0].x  = pp[1].x;
                    pp[0].y  = pp[1].y;
                    pp[1].x  = change.x;
                    pp[1].y  = change.y;
                }

                if (fabs(pp[1].y - 0.5) < 1e-6)
                    foreach_dimension(){
                        pp[2].x = aorient[0] == 1 ? vex[0].x : vex[1].x;
                        pp[3].x = aorient[0] == 1 ? vex[3].x : vex[2].x;
                    }

                else if (fabs(pp[1].y + 0.5) < 1e-6)
                    foreach_dimension(){
                        pp[2].x = aorient[3] == 1 ? vex[3].x : vex[2].x;
                        pp[3].x = aorient[3] == 1 ? vex[0].x : vex[1].x;
                    }
                        
                else if (fabs(pp[1].x - 0.5) < 1e-6)
                    foreach_dimension(){
                        pp[2].x = aorient[0] == 1 ? vex[0].x : vex[3].x;
                        pp[3].x = aorient[0] == 1 ? vex[1].x : vex[2].x;
                }

                else if (fabs(pp[1].x + 0.5) < 1e-6)
                    foreach_dimension(){
                        pp[2].x = aorient[1] == 1 ? vex[1].x : vex[2].x;
                        pp[3].x = aorient[1] == 1 ? vex[0].x : vex[3].x;
                }


                double alphac = mof_4points (c, pp, mc, p_mof);
                return alphac;

            }

            else if(nfacets + an == 5){

                foreach_dimension(){
                    pp[0].x = pcs[0].x;
                    pp[1].x = pcs[1].x;
                }

                if (fabs(fabs(pp[1].x)-0.5)<1e-12 && fabs(fabs(pp[1].y)-0.5)<1e-12){
                    coord change = {0,0};
                    change.x = pp[0].x;
                    change.y = pp[0].y;
                    pp[0].x  = pp[1].x;
                    pp[0].y  = pp[1].y;
                    pp[1].x  = change.x;
                    pp[1].y  = change.y;
                }

                int tmp_ori = -1;

                for(int j = 0; j <= 3; j++){
                    if(aorient[j] == 0){
                        foreach_dimension(){
                            pp[3].x = - vex[j].x;
                        }
                        tmp_ori = j;
                    } 
                }

                if (fabs(pp[1].y - 0.5) < 1e-6 && tmp_ori >= 0){
                    int tmp_f = tmp_ori + 1 > 3 ? tmp_ori - 3 : tmp_ori + 1;
                    int tmp_b = tmp_ori - 1 < 0 ? tmp_ori + 3 : tmp_ori - 1;
                    if (vex[tmp_f].y == 0.5){
                        foreach_dimension(){
                            pp[2].x = vex[tmp_f].x;
                            pp[4].x = - vex[tmp_f].x;
                        }
                    }
                    else if (vex[tmp_b].y == 0.5){
                            foreach_dimension(){
                                pp[2].x = vex[tmp_b].x;
                                pp[4].x = - vex[tmp_b].x;
                            }
                    }
                    else
                        fprintf(stderr,"mof_5points grid error\n");
                }

                if (fabs(pp[1].y + 0.5) < 1e-6 && tmp_ori >= 0){
                    int tmp_f = tmp_ori + 1 > 3 ? tmp_ori - 3 : tmp_ori + 1;
                    int tmp_b = tmp_ori - 1 < 0 ? tmp_ori + 3 : tmp_ori - 1;
                    if (vex[tmp_f].y == - 0.5){
                        foreach_dimension(){
                            pp[2].x = vex[tmp_f].x;
                            pp[4].x = - vex[tmp_f].x;
                        }
                    }
                    else if (vex[tmp_b].y == - 0.5){
                            foreach_dimension(){
                                pp[2].x = vex[tmp_b].x;
                                pp[4].x = - vex[tmp_b].x;
                            }
                    }
                    else
                        fprintf(stderr,"mof_5points grid error\n");
                }

                if (fabs(pp[1].x - 0.5) < 1e-6 && tmp_ori >= 0){
                    int tmp_f = tmp_ori + 1 > 3 ? tmp_ori - 3 : tmp_ori + 1;
                    int tmp_b = tmp_ori - 1 < 0 ? tmp_ori + 3 : tmp_ori - 1;
                    if (vex[tmp_f].x == 0.5){
                        foreach_dimension(){
                            pp[2].x = vex[tmp_f].x;
                            pp[4].x = - vex[tmp_f].x;
                        }
                    }
                    else if (vex[tmp_b].x == 0.5){
                            foreach_dimension(){
                                pp[2].x = vex[tmp_b].x;
                                pp[4].x = - vex[tmp_b].x;
                            }
                    }
                    else
                        fprintf(stderr,"mof_5points grid error\n");
                }
                        
                if (fabs(pp[1].x + 0.5) < 1e-6 && tmp_ori >= 0){
                    int tmp_f = tmp_ori + 1 > 3 ? tmp_ori - 3 : tmp_ori + 1;
                    int tmp_b = tmp_ori - 1 < 0 ? tmp_ori + 3 : tmp_ori - 1;
                    if (vex[tmp_f].x == - 0.5){
                        foreach_dimension(){
                            pp[2].x = vex[tmp_f].x;
                            pp[4].x = - vex[tmp_f].x;
                        }
                    }
                    else if (vex[tmp_b].x == - 0.5){
                            foreach_dimension(){
                                pp[2].x = vex[tmp_b].x;
                                pp[4].x = - vex[tmp_b].x;
                            }
                    }
                    else
                        fprintf(stderr,"mof_5points grid error\n");
                }

                double alphac = mof_5points (c, pp, mc, p_mof);
                
                return alphac;

            }

      }
      else if(nfacets == 0 || nfacets == 1){

        double alphac = plane_alpha(c, mc);
        if (myfacets(mc, alphac, p_mof)==2){
            pp[0].x =  0.5, pp[0].y =  0.5;
            pp[1].x =  0.5, pp[1].y = -0.5;
            pp[2].x = -0.5, pp[2].y = -0.5;
            pp[3].x = -0.5, pp[3].y =  0.5;
            pp[4].x =   10, pp[4].y =   10;   
        }

        else{
            p_mof[0].x = 10, p_mof[0].y = 10;
            p_mof[1].x = 10, p_mof[1].y = 10;

            pp[0].x =  0.5, pp[0].y =  0.5;
            pp[1].x =  0.5, pp[1].y = -0.5;
            pp[2].x = -0.5, pp[2].y = -0.5;
            pp[3].x = -0.5, pp[3].y =  0.5;
            pp[4].x =   10, pp[4].y =   10;   

        }
        return alphac;
      }

    }
    
  }
    p_mof[0].x = 10, p_mof[0].y = 10;
    p_mof[1].x = 10, p_mof[1].y = 10;
    return 0;
}