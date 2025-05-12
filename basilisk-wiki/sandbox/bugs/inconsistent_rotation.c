/**
# Inconsistent rotation of doubly indexed variables

The code... 
 */
tensor a[];

@define str(x) #x  //string replacement before expansion of foreach_dimension()
  
int main() {
  foreach_dimension() {
    puts (a.x.x.name);
    puts (str(b_x_x));
    puts (str(c.x.x));
  }
}

/**
... [prints](inconsistent_rotation/out): 

~~~literatec
a.x.x  
b_x_x  
c.x.x   
a.y.y  
b_x_y  <-- _x_y not _y_y 
c.y.y  
~~~

Is seams _ and . suffixes are not treated the same.
*/


 

