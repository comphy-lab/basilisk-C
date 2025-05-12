#include "higher-order.h"

face vector uff;
double dtf2;
double p_left (Point point, Point neighbor, scalar _s, bool * data) { 
  return (layer_nr_x == 1 ? val(_s,0,0,0) - (val(uff.x,0,0,0)/dtf2)*Delta :
	  val(_s,1,0,0) - 3.*(val(uff.x,0,0,0)/dtf2)*Delta);
}

double p_left_h (Point point, Point neighbor, scalar _s, bool * data) { 
  return layer_nr_x == 1 ? val(_s,0,0,0) : val(_s,1,0,0);
}

double p_right (Point point, Point neighbor, scalar _s, bool * data) { 
  return (layer_nr_x == 1 ? val(_s,0,0,0) - (-val(uff.x,1,0,0)/dtf2)*Delta :
	  val(_s,-1,0,0) - 3.*(-val(uff.x,1,0,0)/dtf2)*Delta);
}

double p_right_h (Point point, Point neighbor, scalar _s, bool * data) { 
  return layer_nr_x == 1 ? val(_s,0,0,0) : val(_s,-1,0,0);
}

double p_bottom (Point point, Point neighbor, scalar _s, bool * data) {
  return  (layer_nr_y == 1 ? val(_s,0,0,0) - (val(uff.y,0,0,0)/dtf2)*Delta :
	   val(_s,0,1,0) - 3.*(val(uff.y,0,0,0)/dtf2)*Delta);
}

static double p_bottom_h (Point point, Point neighbor, scalar _s, bool * data) { 
  return  layer_nr_y == 1 ? val(_s,0,0,0) : val(_s,0,1,0);
}

double p_top (Point point, Point neighbor, scalar _s, bool * data) { 
  return  (layer_nr_y == 1 ? val(_s,0,0,0) + (val(uff.y,0,1,0)/dtf2)*Delta :
	   val(_s,0,-1,0) + 3.*(val(uff.y,0,1,0)/dtf2)*Delta);
}

static double p_top_h (Point point, Point neighbor, scalar _s, bool * data) { 
  return  layer_nr_y == 1 ? val(_s,0,0,0) : val(_s,0,-1,0);
}
