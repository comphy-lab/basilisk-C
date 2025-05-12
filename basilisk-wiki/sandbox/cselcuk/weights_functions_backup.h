/** 
# Set of functions that computes the weights associated to a Lagrange
  multiplier
*/

/** 
## Some definitions
*/
# define sup2deltax (relnl.x < -1.9*delta) && (relnl.x > -2.1*delta)
# define sup1deltax (relnl.x < -0.9*delta) && (relnl.x > -1.1*delta)
# define sup0deltax fabs(relnl.x) < 0.1*delta

# define sup2deltay (relnl.y < -1.9*delta) && (relnl.y > -2.1*delta)
# define sup1deltay (relnl.y < -0.9*delta) && (relnl.y > -1.1*delta)
# define sup0deltay fabs(relnl.y) < 0.1*delta

# define sup2deltaz (relnl.z < -1.9*delta) && (relnl.z > -2.1*delta)
# define sup1deltaz (relnl.z < -0.9*delta) && (relnl.z > -1.1*delta)
# define sup0deltaz fabs(relnl.z) < 0.1*delta

# define inf2deltax (relnl.x < 2.1*delta) && (relnl.x > 1.9*delta)
# define inf1deltax (relnl.x < 1.1*delta) && (relnl.x > 0.9*delta)
  
# define inf2deltay (relnl.y < 2.1*delta) && (relnl.y > 1.9*delta)
# define inf1deltay (relnl.y < 1.1*delta) && (relnl.y > 0.9*delta)

# define inf2deltaz (relnl.z < 2.1*delta) && (relnl.z > 1.9*delta)
# define inf1deltaz (relnl.z < 1.1*delta) && (relnl.z > 0.9*delta)

void compute_relative_vector (coord vec1, coord vec2, coord * rel) {
  foreach_dimension() {
    rel->x = vec2.x - vec1.x;
  }
}

void assign_dial (coord rel, int * CX) {
#if dimension == 2
  if (rel.x > 0 ) {
    if (rel.y > 0 ) {
      /* We are on the lower left cell */
      *CX = 1;
    }
    else if(rel.y < 0 ) {
      /* We are on the upper left cell */
      *CX = 4;
    }
  }
  
  if (rel.x < 0 ) {
    if (rel.y > 0 ) {
      /* we are on the lower right cell */
      *CX = 2;
    }
    else if(rel.y < 0 ) {
      /* We are on the upper right cell */
      *CX = 3;
    }
  }

#elif dimension == 3 
  if (rel.z > 0){
    /* one has to pick the forward cells + the middle ones according to their 2D counterpart */
    if (rel.x > 0 ){
      if (rel.y > 0 ){
	/* We are on the lower left cell */
	*CX = 10;
      }
      else if(rel.y < 0 ) {
	/* We are on the upper left cell */
	*CX = 40;
      }
    }
    if (rel.x < 0 ){
      if (rel.y > 0 ){
	/* we are on the lower right cell */
	*CX = 20;
      }
      else if(rel.y < 0 ) {
	/* We are on the upper right cell */
	*CX = 30;
      }
    }
  }
  
  if (rel.z < 0){
    /* one has to pick the backward cells + the middle ones according to their 2D counterpart */

    if (rel.x > 0 ){
      if (rel.y > 0 ){
	/* We are on the lower left cell */
	*CX = 11;
      }
      else if(rel.y < 0 ) {
	/* We are on the upper left cell */
	*CX = 41;
      }
    }
    if (rel.x < 0 ){
      if (rel.y > 0 ){
	/* We are on the lower right cell */
	*CX = 21;
      }
      else if(rel.y < 0 ) {
	/* We are on the upper right cell */
	*CX = 31;
      }
    }
  }
#endif 
}

void NCX1_CX1 (const coord relnl, const double delta, const int weight_numbers[9], int * weight_id, size_t * goflag) {

  // weight_numbers[2] -> right column top
  // weight_numbers[5] -> right column middle
  // weight_numbers[1] -> right column bottom

  // weight_numbers[6] -> middle column top
  // weight_numbers[8] -> middle colum middle
  // weight_numbers[4] -> middle column bottom

  // weight_numbers[3] -> left column top
  // weight_numbers[7] -> left column middle
  // weight_numbers[0] -> left column bottom

  /* Stencil config, O is the cell containning the Lagrange multiplier */
  /* |_|_|x|x|x| --> |_|_|3|6|2| */
  /* |_|_|x|x|x| --> |_|_|7|8|5| */
  /* |_|_|O|x|x| --> |_|_|0|4|1| */ 
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  // Right column // sup2deltax
  if ((relnl.x < -1.9*delta) && (relnl.x > -2.1*delta)) {
     
    if ((relnl.y < -1.9*delta) && (relnl.y > -2.1*delta)) {
      *weight_id = weight_numbers[2]; *goflag = 1; // right top cell // sup2deltay
    }
    if ((relnl.y < -0.9*delta) && (relnl.y > -1.1*delta)) {
      *weight_id = weight_numbers[5]; *goflag = 1; // right middle cell // sup1deltay
    }
    if (fabs(relnl.y) < 0.1*delta) {
      *weight_id = weight_numbers[1]; *goflag = 1; // right bottom cell // sup0deltay
    }
  }

  
  // middle column // sup1deltax
  if ((relnl.x < -0.9*delta) && (relnl.x > -1.1*delta)) {

    if ((relnl.y < -1.9*delta) && (relnl.y > -2.1*delta)) {
      *weight_id = weight_numbers[6]; *goflag = 1; // middle top cell // sup2deltay
    }
    if ((relnl.y < -0.9*delta) && (relnl.y > -1.1*delta)) {
      *weight_id = weight_numbers[8]; *goflag = 1; // middle cell // sup1deltay
    }
    if (fabs(relnl.y) < 0.1*delta) {
      *weight_id = weight_numbers[4]; *goflag = 1; // middle bottom cell // sup0deltay
    }
  }
  
 
  // left column // sup0deltax
  if (fabs(relnl.x) < 0.1*delta) {
    
    if((relnl.y < -1.9*delta) && (relnl.y > -2.1*delta)) {
      *weight_id = weight_numbers[3]; *goflag = 1; // left top cell // sup2deltay
    }
    if((relnl.y < -0.9*delta) && (relnl.y > -1.1*delta)) {
      *weight_id = weight_numbers[7]; *goflag = 1; // left middle cell // sup1deltay
    }
    if(fabs(relnl.y) < 0.1*delta) {
      *weight_id = weight_numbers[0]; *goflag = 1; // left bottom cell // sup0deltay
    }
  }

  
}
//recoded on 28-02-2019

void NCX1_CX2 (const coord relnl, const double delta, const int weight_numbers[9], int * weight_id, size_t * goflag) {

  // weight_numbers[2] -> right column top
  // weight_numbers[5] -> right column middle
  // weight_numbers[1] -> right column bottom

  // weight_numbers[6] -> middle column top
  // weight_numbers[8] -> middle colum middle
  // weight_numbers[4] -> middle column bottom

  // weight_numbers[3] -> left column top
  // weight_numbers[7] -> left column middle
  // weight_numbers[0] -> left column bottom

  /* Stencil config, O is the cell containning the Lagrange multiplier */
  /* |_|x|x|x|_| -->  |_|3|6|2|_| */
  /* |_|x|x|x|_| -->  |_|7|8|5|_| */
  /* |_|x|O|x|_| -->  |_|0|4|1|_| */ 
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */

  // Right column sup1deltax
  if ((relnl.x < -0.9*delta) && (relnl.x > -1.1*delta)) {
    
      if ((relnl.y < -1.9*delta) && (relnl.y > -2.1*delta)) { // top sup2deltay
      *weight_id = weight_numbers[2]; *goflag = 1;
    }
      if ((relnl.y < -0.9*delta) && (relnl.y > -1.1*delta)) { // middle sup1deltay
      *weight_id = weight_numbers[5]; *goflag = 1;
    }
      if (fabs(relnl.y) < 0.1*delta) { // bottom sup0deltay
      *weight_id = weight_numbers[1]; *goflag = 1;
    }
  }

  // middle column
  if (fabs(relnl.x) < 0.1*delta) {
    if ((relnl.y < -1.9*delta) && (relnl.y > -2.1*delta)) {
      *weight_id = weight_numbers[6]; *goflag = 1;
    }
    if ((relnl.y < -0.9*delta) && (relnl.y > -1.1*delta)) {
      *weight_id = weight_numbers[8]; *goflag = 1;
    }
    if (fabs(relnl.y) < 0.1*delta) {
      *weight_id = weight_numbers[4]; *goflag = 1;
    }
  }

  // left column
  if ((relnl.x < 1.1*delta) && (relnl.x > 0.9*delta)) {
    
    if ((relnl.y < -1.9*delta) && (relnl.y > -2.1*delta)) {
      *weight_id = weight_numbers[3]; *goflag = 1;
    }
    if ((relnl.y < -0.9*delta) && (relnl.y > -1.1*delta)) {
      *weight_id = weight_numbers[7]; *goflag = 1;
    }
    if (fabs(relnl.y) < 0.1*delta) {
      *weight_id = weight_numbers[0]; *goflag = 1;
    }
  }
}
//recoded on 28-02-2019

void NCX_centred (const coord relnl, const double delta, const int weight_numbers[9], int * weight_id, size_t * goflag) {
  /* Stencil config, O is the cell containning the Lagrange multiplier */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |_|x|x|x|_| -->  |_|3|6|2|_| */
  /* |_|x|O|x|_| -->  |_|7|8|5|_| */ 
  /* |_|x|x|x|_| -->  |_|0|4|1|_| */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */

  
  // right column / sup1deltax
  if ((relnl.x < -0.9*delta) && (relnl.x > -1.1*delta)) {
    
    if ((relnl.y < -0.9*delta) && (relnl.y > -1.1*delta)) { // top row : sup1deltay
      *weight_id = weight_numbers[2]; *goflag = 1;
    }
    if (fabs(relnl.y) <  0.1*delta) { // middle row : sup0deltay
      *weight_id = weight_numbers[5]; *goflag = 1;
    }
    if ((relnl.y < 1.1*delta) && (relnl.y > 0.9*delta)) { //bottom row : inf1deltay
      *weight_id = weight_numbers[1]; *goflag = 1;
    }
  }
  // middle column : sup0deltax
  if (fabs(relnl.x) < 0.1*delta) {
    
    if ((relnl.y < -0.9*delta) && (relnl.y > -1.1*delta)) { // top row : sup1deltay
      *weight_id = weight_numbers[6]; *goflag = 1;
    }
    if (fabs(relnl.y) <  0.1*delta) { // middle row : sup0deltay
      *weight_id = weight_numbers[8]; *goflag = 1;
    }
    if((relnl.y < 1.1*delta) && (relnl.y > 0.9*delta)) { //bottom row : inf1deltay
      *weight_id = weight_numbers[4]; *goflag = 1;
    }
  }
  // left column : inf1deltax
  if ((relnl.x < 1.1*delta) && (relnl.x > 0.9*delta)) {

    if ((relnl.y < -0.9*delta) && (relnl.y > -1.1*delta)) { // top row : sup1deltay
      *weight_id = weight_numbers[3]; *goflag = 1;
    }
    if (fabs(relnl.y) <  0.1*delta) { // middle row : sup0deltay
      *weight_id = weight_numbers[7]; *goflag = 1; 
    }
    if ((relnl.y < 1.1*delta) && (relnl.y > 0.9*delta)) { //bottom row : inf1deltay
      *weight_id = weight_numbers[0]; *goflag = 1;
    }
  }

}
//recoded on 28-02-2019

void NCX1_CX4 (const coord relnl, const double delta, const int weight_numbers[9], int * weight_id, size_t * goflag) {
  // weight_numbers[2] -> right column top
  // weight_numbers[5] -> right column middle
  // weight_numbers[1] -> right column bottom

  // weight_numbers[6] -> middle column top
  // weight_numbers[8] -> middle colum middle
  // weight_numbers[4] -> middle column bottom

  // weight_numbers[3] -> left column top
  // weight_numbers[7] -> left column middle
  // weight_numbers[0] -> left column bottom


  /* Stencil config, O is the cell containning the Lagrange multiplier */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |_|_|x|x|x| -->  |_|_|3|6|2| */
  /* |_|_|O|x|x| -->  |_|_|7|8|5| */ 
  /* |_|_|x|x|x| -->  |_|_|0|4|1| */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */

  
  // right column : sup2deltax
  if ((relnl.x < -1.9*delta) && (relnl.x > -2.1*delta)) {

    if ((relnl.y < -0.9*delta) && (relnl.y > -1.1*delta)) { // top row : sup1deltay
      *weight_id = weight_numbers[2]; *goflag = 1;
    }
    if (fabs(relnl.y) <  0.1*delta) { // middle row : sup0deltay
      *weight_id = weight_numbers[5]; *goflag = 1;
    }
    if ((relnl.y < 1.1*delta) && (relnl.y > 0.9*delta)) { // bottom row : inf1deltay
      *weight_id = weight_numbers[1]; *goflag = 1;
    }
  }
  
  // middle column : sup1deltax
  if ((relnl.x < -0.9*delta) && (relnl.x > -1.1*delta)) {
    
    if ((relnl.y < -0.9*delta) && (relnl.y > -1.1*delta)) { // top row : sup1deltay
      *weight_id = weight_numbers[6]; *goflag = 1;
    }
    if (fabs(relnl.y) <  0.1*delta) { // middle row : sup0deltay
      *weight_id = weight_numbers[8]; *goflag = 1;
    }
    if ((relnl.y < 1.1*delta) && (relnl.y > 0.9*delta)) { // bottom row : inf1deltay
      *weight_id = weight_numbers[4]; *goflag = 1;
    }
  }
  // left column : sup0deltax
  if (fabs(relnl.x) < 0.1*delta) {

    if ((relnl.y < -0.9*delta) && (relnl.y > -1.1*delta)) { // top row : sup1deltay
      *weight_id = weight_numbers[3]; *goflag = 1;
    }
    if(fabs(relnl.y) <  0.1*delta) { // middle row : sup0deltay
      *weight_id = weight_numbers[7]; *goflag = 1;
    }
    if((relnl.y < 1.1*delta) && (relnl.y > 0.9*delta)) { // bottom row : inf1deltay
      *weight_id = weight_numbers[0]; *goflag = 1;
    }
  }
}
//recoded on 28-02-2019

void NCX2_CX2 (const coord relnl, const double delta, const int weight_numbers[9], int * weight_id, size_t * goflag) {
  // weight_numbers[2] -> right column top
  // weight_numbers[5] -> right column middle
  // weight_numbers[1] -> right column bottom

  // weight_numbers[6] -> middle column top
  // weight_numbers[8] -> middle colum middle
  // weight_numbers[4] -> middle column bottom

  // weight_numbers[3] -> left column top
  // weight_numbers[7] -> left column middle
  // weight_numbers[0] -> left column bottom


  /* Stencil config, O is the cell containning the Lagrange multiplier */
  /* |x|x|x|_|_| -->  |3|6|2|_|_| */
  /* |x|x|x|_|_| -->  |7|8|5|_|_| */
  /* |x|x|O|_|_| -->  |0|4|1|_|_| */ 
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */

  
  // right column : sup0deltax
  if (fabs(relnl.x) < 0.1*delta) {

    if ((relnl.y < -1.9*delta) && (relnl.y > -2.1*delta)) { // top row : sup2deltay
      *weight_id = weight_numbers[2]; *goflag = 1;
    }
    if ((relnl.y < -0.9*delta) && (relnl.y > -1.1*delta)) { // middle row : sup1deltay
      *weight_id = weight_numbers[5]; *goflag = 1;
    }
    if (fabs(relnl.y) < 0.1*delta) {
      *weight_id = weight_numbers[1]; *goflag = 1; // bottom row : sup0deltay
    }
  }
  
  // middle column : inf1deltax 
  if ((relnl.x < 1.1*delta) && (relnl.x > 0.9*delta)) {

    if ((relnl.y < -1.9*delta) && (relnl.y > -2.1*delta)) { // top row : sup2deltay
      *weight_id = weight_numbers[6]; *goflag = 1;
    }
    if ((relnl.y < -0.9*delta) && (relnl.y > -1.1*delta)) { // middle row : sup1deltay
      *weight_id = weight_numbers[8]; *goflag = 1;
    }
    if (fabs(relnl.y) < 0.1*delta) { // bottom row : sup0deltay
      *weight_id = weight_numbers[4]; *goflag = 1;
    }
  }
  
  // left column inf2deltax
  if ((relnl.x < 2.1*delta) && (relnl.x > 1.9*delta)) {

    if ((relnl.y < -1.9*delta) && (relnl.y > -2.1*delta)) { // top row : sup2deltay
      *weight_id = weight_numbers[3]; *goflag = 1;
    }
    if ((relnl.y < -0.9*delta) && (relnl.y > -1.1*delta)) { // middle row : sup1deltay
      *weight_id = weight_numbers[7]; *goflag = 1;
    }
    if (fabs(relnl.y) < 0.1*delta) { // bottom row : sup0deltay
      *weight_id = weight_numbers[0]; *goflag = 1;
    }
  }
}
//recoded on 28-02-2019

void NCX2_CX3 (const coord relnl, const double delta, const int weight_numbers[9], int * weight_id, size_t * goflag) {
  // weight_numbers[2] -> right column top
  // weight_numbers[5] -> right column middle
  // weight_numbers[1] -> right column bottom

  // weight_numbers[6] -> middle column top
  // weight_numbers[8] -> middle colum middle
  // weight_numbers[4] -> middle column bottom

  // weight_numbers[3] -> left column top
  // weight_numbers[7] -> left column middle
  // weight_numbers[0] -> left column bottom

  /* Stencil config, O is the cell containning the Lagrange multiplier */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |x|x|x|_|_| -->  |3|6|2|_|_| */
  /* |x|x|O|_|_| -->  |7|8|5|_|_| */ 
  /* |x|x|x|_|_| -->  |0|4|1|_|_| */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */

  
  // right column : sup0deltax
  if (fabs(relnl.x) < 0.1*delta) {

    if ((relnl.y < -0.9*delta) && (relnl.y > -1.1*delta)) { // top row : sup1deltay
      *weight_id = weight_numbers[2]; *goflag = 1;
    }
    if (fabs(relnl.y) <  0.1*delta) { // middle row : sup0deltay
      *weight_id = weight_numbers[5]; *goflag = 1;
    }
    if ((relnl.y < 1.1*delta) && (relnl.y > 0.9*delta)) { // bottom row : inf1deltay
      *weight_id = weight_numbers[1]; *goflag = 1;
    }
  }
  
  // middle column : inf1deltax
  if ((relnl.x < 1.1*delta) && (relnl.x > 0.9*delta)) {

    if ((relnl.y < -0.9*delta) && (relnl.y > -1.1*delta)) { // top row : sup1deltay
      *weight_id = weight_numbers[6]; *goflag = 1;
    }
    if (fabs(relnl.y) <  0.1*delta) { // middle row : sup0deltay
      *weight_id = weight_numbers[8]; *goflag = 1;
    }
    if ((relnl.y < 1.1*delta) && (relnl.y > 0.9*delta)) { // bottom row : inf1deltay
      *weight_id = weight_numbers[4]; *goflag = 1;
    }
  }
  
  // left column : inf2deltax
  if ((relnl.x < 2.1*delta) && (relnl.x > 1.9*delta)) {
    if ((relnl.y < -0.9*delta) && (relnl.y > -1.1*delta)) { // top row : sup1deltay
      *weight_id = weight_numbers[3]; *goflag = 1;
    }
    if (fabs(relnl.y) <  0.1*delta) { // middle row : sup0deltay
      *weight_id = weight_numbers[7]; *goflag = 1;
    }
    if ((relnl.y < 1.1*delta) && (relnl.y > 0.9*delta)) { // bottom row : inf1deltay
      *weight_id = weight_numbers[0]; *goflag = 1;
    }	
  }
}
//recoded on 28-02-2019

void NCX3_CX3 (const coord relnl, const double delta, const int weight_numbers[9], int * weight_id, size_t * goflag) {
  // weight_numbers[2] -> right column top
  // weight_numbers[5] -> right column middle
  // weight_numbers[1] -> right column bottom

  // weight_numbers[6] -> middle column top
  // weight_numbers[8] -> middle colum middle
  // weight_numbers[4] -> middle column bottom

  // weight_numbers[3] -> left column top
  // weight_numbers[7] -> left column middle
  // weight_numbers[0] -> left column bottom

  /* Stencil config, O is the cell containning the Lagrange multiplier */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |x|x|O|_|_| -->  |3|6|2|_|_| */ 
  /* |x|x|x|_|_| -->  |7|8|5|_|_| */
  /* |x|x|x|_|_| -->  |0|4|1|_|_| */
  
  // Right column : sup0deltax
  if (fabs(relnl.x) < 0.1*delta) {
      
    if (fabs(relnl.y) < 0.1*delta) { // top row : sup0deltay
      *weight_id = weight_numbers[2]; *goflag = 1;
    }
    if ((relnl.y < 1.1*delta) && (relnl.y > 0.9*delta)) { // middle row : inf1deltay
      *weight_id = weight_numbers[5]; *goflag = 1; 
    }
    if ((relnl.y < 2.1*delta) && (relnl.y > 1.9*delta)) { // bottom row : inf2deltay
      *weight_id = weight_numbers[1]; *goflag = 1;
    }
  }
  
  // middle column : inf1deltax
  if ((relnl.x < 1.1*delta) && (relnl.x > 0.9*delta)) {

    if (fabs(relnl.y) < 0.1*delta) { // top row : sup0deltay
      *weight_id = weight_numbers[6]; *goflag = 1;
    }
    if ((relnl.y < 1.1*delta) && (relnl.y > 0.9*delta)) { // middle row : inf1deltay
      *weight_id = weight_numbers[8]; *goflag = 1; 
    }
    if ((relnl.y < 2.1*delta) && (relnl.y > 1.9*delta)) { // bottom row : inf2deltay
      *weight_id = weight_numbers[4]; *goflag = 1;
    }
  }
  
  // left column : inf2deltax
  if ((relnl.x < 2.1*delta) && (relnl.x > 1.9*delta)) {
	
    if (fabs(relnl.y) < 0.1*delta) { // top row : sup0deltay
      *weight_id = weight_numbers[3]; *goflag = 1;
    }
    if ((relnl.y < 1.1*delta) && (relnl.y > 0.9*delta)) { // middle row : inf1deltay
      *weight_id = weight_numbers[7]; *goflag = 1; 
    }
    if ((relnl.y < 2.1*delta) && (relnl.y > 1.9*delta)) { // bottom row : inf2deltay
      *weight_id = weight_numbers[0]; *goflag = 1;
    }
  }

}
//recoded on 28-02-2019

void NCX3_CX4 (const coord relnl, const double delta, const int weight_numbers[9], int * weight_id, size_t * goflag) {
  // weight_numbers[2] -> right column top
  // weight_numbers[5] -> right column middle
  // weight_numbers[1] -> right column bottom

  // weight_numbers[6] -> middle column top
  // weight_numbers[8] -> middle colum middle
  // weight_numbers[4] -> middle column bottom

  // weight_numbers[3] -> left column top
  // weight_numbers[7] -> left column middle
  // weight_numbers[0] -> left column bottom

  /* Stencil config, O is the cell containning the Lagrange multiplier */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |_|x|O|x|_| -->  |_|3|6|2|_| */ 
  /* |_|x|x|x|_| -->  |_|7|8|5|_| */
  /* |_|x|x|x|_| -->  |_|0|4|1|_| */

  
  // right column : sup1deltax
  if ((relnl.x < -0.9*delta) && (relnl.x > -1.1*delta)) {

    if (fabs(relnl.y) < 0.1*delta) { // top row : sup0deltay
      *weight_id = weight_numbers[2]; *goflag = 1;
    }
    if ((relnl.y < 1.1*delta) && (relnl.y > 0.9*delta)) { // middle row : inf1deltay
      *weight_id = weight_numbers[5]; *goflag = 1; 
    }
    if ((relnl.y < 2.1*delta) && (relnl.y > 1.9*delta)) { // bottom row : inf2deltay
      *weight_id = weight_numbers[1]; *goflag = 1;
    }
  }

  // middle column : sup0deltax
  if (fabs(relnl.x) < 0.1*delta) {
    
    if (fabs(relnl.y) < 0.1*delta) { // top row : sup0deltay
      *weight_id = weight_numbers[6]; *goflag = 1;
    }
    if ((relnl.y < 1.1*delta) && (relnl.y > 0.9*delta)) { // middle row : inf1deltay
      *weight_id = weight_numbers[8]; *goflag = 1; 
    }
    if ((relnl.y < 2.1*delta) && (relnl.y > 1.9*delta)) { // bottom row : inf2deltay
      *weight_id = weight_numbers[4]; *goflag = 1;
    }
  }
  
  // left column : inf1deltax
  if ((relnl.x < 1.1*delta) && (relnl.x > 0.9*delta)) {

    if (fabs(relnl.y) < 0.1*delta) { // top row : sup0deltay
      *weight_id = weight_numbers[3]; *goflag = 1;
    }
    if ((relnl.y < 1.1*delta) && (relnl.y > 0.9*delta)) { // middle row : inf1deltay
      *weight_id = weight_numbers[7]; *goflag = 1; 
    }
    if ((relnl.y < 2.1*delta) && (relnl.y > 1.9*delta)) { // bottom row : inf2deltay
      *weight_id = weight_numbers[0]; *goflag = 1;
    }
  }
}
//recoded on 28-02-2019

void NCX4_CX4 (const coord relnl, const double delta, const int weight_numbers[9], int * weight_id, size_t * goflag) {
  /* Stencil config, O is the cell containning the Lagrange multiplier */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |_|_|O|x|x| -->  |_|_|3|6|2| */ 
  /* |_|_|x|x|x| -->  |_|_|7|8|5| */
  /* |_|_|x|x|x| -->  |_|_|0|4|1| */
  
  // right column : sup2deltax
  if ((relnl.x < -1.9*delta) && (relnl.x > -2.1*delta)) {

    if (fabs(relnl.y) < 0.1*delta) { // top row : sup0deltay
      *weight_id = weight_numbers[2]; *goflag = 1;
    }
    if ((relnl.y < 1.1*delta) && (relnl.y > 0.9*delta)) { // middle row : inf1deltay
      *weight_id = weight_numbers[5]; *goflag = 1; 
    }
    if ((relnl.y < 2.1*delta) && (relnl.y > 1.9*delta)) { // bottom row : inf2deltay
      *weight_id = weight_numbers[1]; *goflag = 1;
    }
  }
  
  // middle column : sup1deltax
  if ((relnl.x < -0.9*delta) && (relnl.x > -1.1*delta)) {

    if (fabs(relnl.y) < 0.1*delta) { // top row : sup0deltay
      *weight_id = weight_numbers[6]; *goflag = 1;
    }
    if ((relnl.y < 1.1*delta) && (relnl.y > 0.9*delta)) { // middle row : inf1deltay
      *weight_id = weight_numbers[8]; *goflag = 1; 
    }
    if ((relnl.y < 2.1*delta) && (relnl.y > 1.9*delta)) { // bottom row : inf2deltay
      *weight_id = weight_numbers[4]; *goflag = 1;
    }
  }

  // left column : sup0deltax
  if (fabs(relnl.x) < 0.1*delta) {
    
    if (fabs(relnl.y) < 0.1*delta) { // top row : sup0deltay
      *weight_id = weight_numbers[3]; *goflag = 1;
    }
    if ((relnl.y < 1.1*delta) && (relnl.y > 0.9*delta)) { // middle row : inf1deltay
      *weight_id = weight_numbers[7]; *goflag = 1; 
    }
    if ((relnl.y < 2.1*delta) && (relnl.y > 1.9*delta)) { // bottom row : inf2deltay
      *weight_id = weight_numbers[0]; *goflag = 1;
    }
  }

}
//recoded on 28-02-2019

void fill_weight_numbers (const int wisz, int * weight_numbers) {

  // In the plane:
  // weight_numbers[2] -> right column top
  // weight_numbers[5] -> right column middle
  // weight_numbers[1] -> right column bottom

  // weight_numbers[6] -> middle column top
  // weight_numbers[8] -> middle colum middle
  // weight_numbers[4] -> middle column bottom

  // weight_numbers[3] -> left column top
  // weight_numbers[7] -> left column middle
  // weight_numbers[0] -> left column bottom

  for (int i = 0; i < 9; i ++) {
    weight_numbers[i] = -1;
  }

  
  // wisz == -1 backward plane 
  if(wisz == -1) {
    weight_numbers[3] = 6; weight_numbers[6] = 7; weight_numbers[2] = 8;
    weight_numbers[7] = 5; weight_numbers[8] = 4; weight_numbers[5] = 3;
    weight_numbers[0] = 0; weight_numbers[4] = 1; weight_numbers[1] = 2;
  }
  
  // wisz == 0 middle plane 
  if(wisz == 0) {
    weight_numbers[3] = 15; weight_numbers[6] = 16; weight_numbers[2] = 17;
    weight_numbers[7] = 14; weight_numbers[8] = 13; weight_numbers[5] = 12;
    weight_numbers[0] = 9;  weight_numbers[4] = 10; weight_numbers[1] = 11;
  }

  // wisz == 1 forward plane 
  if(wisz == 1) {
    weight_numbers[3] = 24; weight_numbers[6] = 25; weight_numbers[2] = 26;
    weight_numbers[7] = 23; weight_numbers[8] = 22; weight_numbers[5] = 21;
    weight_numbers[0] = 18; weight_numbers[4] = 19; weight_numbers[1] = 20;
  }
}
//recoded on 28-02-2019

#if dimension == 3 
//NCX10_CX Bloc, 8 functions
void NCX10_CX10 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containning the Lagrange multiplier */
  
  /* z = localplane (plane containning the multiplier) */
  /* |_|_|x|x|x| --> |_|_|6|7|8| */
  /* |_|_|x|x|x| --> |_|_|5|4|3| */
  /* |_|_|O|x|x| --> |_|_|0|1|2| */ 
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  // backward z plane: z = zlocal : sup0deltaz
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (-1,  &weight_numbers[0]);
    NCX1_CX1 (relnl, delta, weight_numbers, weight_id, goflag);
  }

  /* z = localplane + delta */
  /* |_|_|x|x|x| --> |_|_|15|16|18| */
  /* |_|_|x|x|x| --> |_|_|14|13|12| */
  /* |_|_|x|x|x| --> |_|_|_9|10|11| */ 
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  
   // middle z plane: z = zlocal + Delta : sup1deltaz
  if ((relnl.z < -0.9*delta) && (relnl.z > -1.1*delta)) {
    fill_weight_numbers (0,  &weight_numbers[0]);
    NCX1_CX1 (relnl, delta, weight_numbers, weight_id, goflag);
  }

  /* z = localplane + 2delta */
  /* |_|_|x|x|x| --> |_|_|24|25|26| */
  /* |_|_|x|x|x| --> |_|_|23|22|21| */
  /* |_|_|x|x|x| --> |_|_|18|19|20| */ 
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  
  // forward z plane: z = zlocal + 2*Delta : sup2deltaz
  if ((relnl.z < -1.9*delta) && (relnl.z > -2.1*delta)) {
    fill_weight_numbers (1,  &weight_numbers[0]);
    NCX1_CX1 (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

void NCX10_CX20 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] =  {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
 /* Stencil config, O is the cell containning the Lagrange multiplier */
  
  /* backward z plane: z = localplane : sup0deltaz (plane containning the multiplier) */
  /* |_|x|x|x|_| --> |_|6|7|8|_| */
  /* |_|x|x|x|_| --> |_|5|4|3|_| */
  /* |_|x|O|x|_| --> |_|0|1|2|_| */ 
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (-1,  &weight_numbers[0]);
    NCX1_CX2 (relnl, delta, weight_numbers, weight_id, goflag);
  }

  /* middle z plane: z = zlocal + Delta : sup1deltaz */
  /* |_|x|x|x|_| --> |_|15|16|17|_| */
  /* |_|x|x|x|_| --> |_|14|13|12|_| */
  /* |_|x|x|x|_| --> |_|_9|10|11|_| */ 
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  if ((relnl.z < -0.9*delta) && (relnl.z > -1.1*delta)) {
    fill_weight_numbers (0,  &weight_numbers[0]);
    NCX1_CX2 (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* forward z plane: z = zlocal + 2Delta : sup2deltaz */
  /* |_|x|x|x|_| --> |_|24|25|26|_| */
  /* |_|x|x|x|_| --> |_|23|22|21|_| */
  /* |_|x|x|x|_| --> |_|18|19|20|_| */ 
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  if ((relnl.z < -1.9*delta) && (relnl.z > -2.1*delta)) {
    fill_weight_numbers(1,  &weight_numbers[0]);
    NCX1_CX2 (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

void NCX10_CX30 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};

  /* Stencil config, O is the cell containning the Lagrange multiplier */
  
  /* backward z plane: z = localplane : sup0deltaz (plane containning the multiplier) */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|x|x|x|_| --> |_|6|7|8|_| */
  /* |_|x|O|x|_| --> |_|5|4|3|_| */ 
  /* |_|x|x|x|_| --> |_|0|1|2|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */

  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX_centred (relnl, delta, weight_numbers, weight_id, goflag);
  }

  /* middle z plane: z = zlocal + Delta : sup1deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|15|16|17|_| */
  /* |_|x|x|x|_| --> |_|14|13|12|_| */ 
  /* |_|x|x|x|_| --> |_|_9|10|11|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  if ((relnl.z < -0.9*delta) && (relnl.z > -1.1*delta)) {
    fill_weight_numbers (0, &weight_numbers[0]);
    NCX_centred (relnl, delta, weight_numbers, weight_id, goflag);
  }

  /* middle z plane: z = zlocal + 2Delta : sup2deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|24|25|26|_| */
  /* |_|x|x|x|_| --> |_|23|22|21|_| */ 
  /* |_|x|x|x|_| --> |_|18|19|20|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  if ((relnl.z < -1.9*delta) && (relnl.z > -2.1*delta)) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX_centred (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

void NCX10_CX40 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};

  /* Stencil config, O is the cell containning the Lagrange multiplier */
  /* backward z plane: z = localplane : sup0deltaz (plane containning the multiplier) */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|O|x|x| --> |_|6|7|8|_| */ 
  /* |_|_|x|x|x| --> |_|5|4|3|_| */
  /* |_|_|x|x|x| --> |_|0|1|2|_| */

  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX1_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }

  /* middle z plane: z = zlocal + Delta : sup1deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|x|x|x| --> |_|15|16|17|_| */ 
  /* |_|_|x|x|x| --> |_|14|13|12|_| */
  /* |_|_|x|x|x| --> |_|_9|10|11|_| */

  if ((relnl.z < -0.9*delta) && (relnl.z > -1.1*delta)) {
    fill_weight_numbers (0, &weight_numbers[0]);
    NCX1_CX4(relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* forward z plane: z = zlocal + 2*Delta : sup2deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|x|x|x| --> |_|24|25|26|_| */ 
  /* |_|_|x|x|x| --> |_|23|22|21|_| */
  /* |_|_|x|x|x| --> |_|18|19|20|_| */
  
  if ((relnl.z < -1.9*delta) && (relnl.z > -2.1*delta)) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX1_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

void NCX10_CX11 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};

  /* Stencil config, O is the cell containning the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - Delta : inf1deltaz*/
  /* |_|_|x|x|x| --> |_|_|6|7|8| */
  /* |_|_|x|x|x| --> |_|_|5|4|3| */
  /* |_|_|x|x|x| --> |_|_|0|1|2| */ 
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if ((relnl.z < 1.1*delta) && (relnl.z > 0.9*delta)) {
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX1_CX1 (relnl, delta, weight_numbers, weight_id, goflag);
  }

  /* middle z plane: z = zlocal : sup0deltaz */
  /* |_|_|x|x|x| --> |_|_|15|16|17| */
  /* |_|_|x|x|x| --> |_|_|14|13|12| */
  /* |_|_|O|x|x| --> |_|_|_9|10|11| */ 
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
    
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (0, &weight_numbers[0]);
    NCX1_CX1 (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* forward z plane: z = zlocal + Delta : sup1deltaz */
  /* |_|_|x|x|x| --> |_|_|24|25|26| */
  /* |_|_|x|x|x| --> |_|_|23|22|21| */
  /* |_|_|x|x|x| --> |_|_|18|19|20| */ 
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */

  if ((relnl.z < -0.9*delta) && (relnl.z > -1.1*delta)) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX1_CX1 (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

void NCX10_CX21 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};

  /* Stencil config, O is the cell containning the Lagrange multiplier */
  /* backward z plane: z = zlocal - Delta : inf1deltaz*/
  /* |_|x|x|x|_| -->  |_|6|7|8|_| */
  /* |_|x|x|x|_| -->  |_|5|4|3|_| */
  /* |_|x|x|x|_| -->  |_|0|1|2|_| */ 
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */

  if ((relnl.z < 1.1*delta) && (relnl.z > 0.9*delta)) {
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX1_CX2 (relnl, delta, weight_numbers, weight_id, goflag);
  }
  /* middle z plane: z = zlocal : sup0deltaz*/
  /* |_|x|x|x|_| -->  |_|15|16|17|_| */
  /* |_|x|x|x|_| -->  |_|14|13|12|_| */
  /* |_|x|O|x|_| -->  |_|_9|10|11|_| */ 
  /* |_|_|_|_|_| -->  |_|__|__|__|_| */
  /* |_|_|_|_|_| -->  |_|__|__|__|_| */
  
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (0, &weight_numbers[0]);
    NCX1_CX2 (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* forward z plane: z = zlocal + Delta : sup1deltaz */
  /* |_|x|x|x|_| -->  |_|24|25|26|_| */
  /* |_|x|x|x|_| -->  |_|23|22|21|_| */
  /* |_|x|x|x|_| -->  |_|18|19|20|_| */ 
  /* |_|_|_|_|_| -->  |_|__|__|__|_| */
  /* |_|_|_|_|_| -->  |_|__|__|__|_| */
  
  if ((relnl.z < -0.9*delta) && (relnl.z > -1.1*delta)) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX1_CX2 (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

void NCX10_CX31 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};

  /* Stencil config, O is the cell containning the Lagrange multiplier */
  /* backward z plane: z = zlocal - Delta : inf1deltaz*/
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |_|x|x|x|_| -->  |_|6|7|8|_| */
  /* |_|x|x|x|_| -->  |_|5|4|3|_| */ 
  /* |_|x|x|x|_| -->  |_|0|1|2|_| */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
 
  if ((relnl.z < 1.1*delta) && (relnl.z > 0.9*delta)) {
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX_centred (relnl, delta, weight_numbers, weight_id, goflag);
  }

  /* middle z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| -->  |_|__|__|__|_| */
  /* |_|x|x|x|_| -->  |_|15|16|17|_| */
  /* |_|x|O|x|_| -->  |_|14|13|12|_| */ 
  /* |_|x|x|x|_| -->  |_|_9|10|11|_| */
  /* |_|_|_|_|_| -->  |_|__|__|__|_| */
  
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (0,  &weight_numbers[0]);
    NCX_centred (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* forward z plane: z = zlocal + Delta : sup1deltaz */
  /* |_|_|_|_|_| -->  |_|__|__|__|_| */
  /* |_|x|x|x|_| -->  |_|24|25|26|_| */
  /* |_|x|x|x|_| -->  |_|23|22|21|_| */ 
  /* |_|x|x|x|_| -->  |_|18|19|20|_| */
  /* |_|_|_|_|_| -->  |_|__|__|__|_| */
  
  if ((relnl.z < -0.9*delta) && (relnl.z > -1.1*delta)) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX_centred (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

void NCX10_CX41 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containning the Lagrange multiplier */
  /* backward z plane: z = zlocal - Delta : inf1deltaz */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |_|_|x|x|x| -->  |_|_|6|7|8| */
  /* |_|_|x|x|x| -->  |_|_|5|4|3| */ 
  /* |_|_|x|x|x| -->  |_|_|0|1|2| */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */

  if ((relnl.z < 1.1*delta) && (relnl.z > 0.9*delta)) {
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX1_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }

  /* middle z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| -->  |_|_|__|__|__| */
  /* |_|_|x|x|x| -->  |_|_|15|16|17| */
  /* |_|_|O|x|x| -->  |_|_|14|13|12| */ 
  /* |_|_|x|x|x| -->  |_|_|_9|10|11| */
  /* |_|_|_|_|_| -->  |_|_|__|__|__| */

  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (0, &weight_numbers[0]);
    NCX1_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  // forward z plane: z = zlocal + Delta : sup1deltaz
  /* |_|_|_|_|_| -->  |_|_|__|__|__| */
  /* |_|_|x|x|x| -->  |_|_|24|25|26| */
  /* |_|_|x|x|x| -->  |_|_|23|22|21| */ 
  /* |_|_|x|x|x| -->  |_|_|18|19|20| */
  /* |_|_|_|_|_| -->  |_|_|__|__|__| */
  
  if ((relnl.z < -0.9*delta) && (relnl.z > -1.1*delta)) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX1_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

// NCX11_CX Bloc, 4 functions

/* NCX11_CX20 = NCX10_CX21 */

/* NCX11_CX30 = NCX10_CX31 */

/* NCX11_CX40 = NCX10_CX41 */

void NCX11_CX11 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containning the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - 2Delta : inf2deltaz */
  /* |_|_|x|x|x| --> |_|_|6|7|8| */
  /* |_|_|x|x|x| --> |_|_|5|4|3| */
  /* |_|_|x|x|x| --> |_|_|0|1|2| */ 
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if ((relnl.z < 2.1*delta) && (relnl.z > 1.9*delta)) {
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX1_CX1 (relnl, delta, weight_numbers, weight_id, goflag);
  }

  /* middle z plane: z = zlocal - Delta : inf1deltaz */
  /* |_|_|x|x|x| --> |_|_|15|16|17| */
  /* |_|_|x|x|x| --> |_|_|14|13|12| */
  /* |_|_|x|x|x| --> |_|_|_9|10|11| */ 
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  
  if ((relnl.z < 1.1*delta) && (relnl.z > 0.9*delta)) {
    fill_weight_numbers (0, &weight_numbers[0]);
    NCX1_CX1 (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* forward z plane: z = zlocal : sup0deltaz */
  /* |_|_|x|x|x| --> |_|_|24|25|26| */
  /* |_|_|x|x|x| --> |_|_|23|22|21| */
  /* |_|_|O|x|x| --> |_|_|18|19|20| */ 
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (1,  &weight_numbers[0]);
    NCX1_CX1 (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

void NCX11_CX21 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  /* Stencil config, O is the cell containning the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - 2Delta : inf2deltaz */
  /* |_|x|x|x|_| --> |_|6|7|8|_| */
  /* |_|x|x|x|_| --> |_|5|4|3|_| */
  /* |_|x|x|x|_| --> |_|0|1|2|_| */ 
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if ((relnl.z < 2.1*delta) && (relnl.z > 1.9*delta)) {
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX1_CX2 (relnl, delta, weight_numbers, weight_id, goflag);
  }

  /* middle z plane: z = zlocal - Delta : inf1deltaz */
  /* |_|x|x|x|_| --> |_|15|16|17|_| */
  /* |_|x|x|x|_| --> |_|14|13|12|_| */
  /* |_|x|x|x|_| --> |_|_9|10|11|_| */ 
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  if ((relnl.z < 1.1*delta) && (relnl.z > 0.9*delta)) {
    fill_weight_numbers (0, &weight_numbers[0]);
    NCX1_CX2 (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* forward z plane: z = zlocal : sup0deltaz */
  /* |_|x|x|x|_| --> |_|24|25|26|_| */
  /* |_|x|x|x|_| --> |_|23|22|21|_| */
  /* |_|x|O|x|_| --> |_|18|19|20|_| */ 
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX1_CX2 (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

void NCX11_CX31 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  /* Stencil config, O is the cell containning the Lagrange multiplier */

  /* backward z plane: z = zlocal - 2Delta : inf2deltaz */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|x|x|x|_| --> |_|6|7|8|_| */
  /* |_|x|x|x|_| --> |_|5|4|3|_| */ 
  /* |_|x|x|x|_| --> |_|0|1|2|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if ((relnl.z < 2.1*delta) && (relnl.z > 1.9*delta)) { 
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX_centred (relnl, delta, weight_numbers, weight_id, goflag);
  }

  /* middle z plane: z = zlocal - Delta : inf1deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|15|16|17|_| */
  /* |_|x|x|x|_| --> |_|14|13|12|_| */ 
  /* |_|x|x|x|_| --> |_|_9|10|11|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  if ((relnl.z < 1.1*delta) && (relnl.z > 0.9*delta)) {
    fill_weight_numbers (0, &weight_numbers[0]);
    NCX_centred (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* forward z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|24|25|26|_| */
  /* |_|x|x|x|_| --> |_|23|22|21|_| */ 
  /* |_|x|x|x|_| --> |_|18|19|20|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX_centred (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

void NCX11_CX41 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  /* Stencil config, O is the cell containning the Lagrange multiplier */

  /* backward z plane: z = zlocal - 2Delta : inf2deltaz */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|x|x|x| --> |_|_|6|7|8| */
  /* |_|_|x|x|x| --> |_|_|5|4|3| */ 
  /* |_|_|x|x|x| --> |_|_|0|1|2| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if ((relnl.z < 2.1*delta) && (relnl.z > 1.9*delta)) {
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX1_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }

  /* middle z plane: z = zlocal - Delta : inf1deltaz */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|x|x|x| --> |_|_|15|16|17| */
  /* |_|_|x|x|x| --> |_|_|14|13|12| */ 
  /* |_|_|x|x|x| --> |_|_|_9|10|11| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  
  if ((relnl.z < 1.1*delta) && (relnl.z > 0.9*delta)) {
    fill_weight_numbers (0, &weight_numbers[0]);
    NCX1_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* forward z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|x|x|x| --> |_|_|24|25|26| */
  /* |_|_|x|x|x| --> |_|_|23|22|21| */ 
  /* |_|_|x|x|x| --> |_|_|18|19|20| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX1_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

//NCX20_CX Bloc, 6 functions
void NCX20_CX20 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containning the Lagrange multiplier */
  
  /* z = localplane (plane containning the multiplier) */
  /* backward z plane: z = zlocal : sup0deltaz */
  /* |x|x|x|_|_| --> |6|7|8|_|_| */
  /* |x|x|x|_|_| --> |5|4|3|_|_| */
  /* |x|x|O|_|_| --> |0|1|2|_|_| */ 
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX2_CX2 (relnl, delta, weight_numbers, weight_id, goflag);
  }

  /* middle z plane: z = zlocal + Delta : sup1delta*/
  /* |x|x|x|_|_| --> |15|16|17|_|_| */
  /* |x|x|x|_|_| --> |14|13|12|_|_| */
  /* |x|x|x|_|_| --> |_9|10|11|_|_| */ 
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */

  if ((relnl.z < -0.9*delta) && (relnl.z > -1.1*delta)) {
    fill_weight_numbers (0, &weight_numbers[0]);
    NCX2_CX2 (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* forward z plane: z = zlocal + 2*Delta : sup2delta*/
  /* |x|x|x|_|_| --> |24|25|26|_|_| */
  /* |x|x|x|_|_| --> |23|22|21|_|_| */
  /* |x|x|x|_|_| --> |18|19|20|_|_| */ 
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  
  if ((relnl.z < -1.9*delta) && (relnl.z > -2.1*delta)) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX2_CX2 (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

void NCX20_CX30 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containning the Lagrange multiplier */
  
  /* z = localplane (plane containning the multiplier) */
  /* backward z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |x|x|x|_|_| --> |6|7|8|_|_| */
  /* |x|x|O|_|_| --> |5|4|3|_|_| */ 
  /* |x|x|x|_|_| --> |0|1|2|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX2_CX3 (relnl, delta, weight_numbers, weight_id, goflag);
  }

  /* middle z plane: z = zlocal + Delta : sup1deltaz */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |x|x|x|_|_| --> |15|16|17|_|_| */
  /* |x|x|x|_|_| --> |14|13|12|_|_| */ 
  /* |x|x|x|_|_| --> |_9|10|11|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  
  if ((relnl.z < -0.9*delta) && (relnl.z > -1.1*delta)) {
    fill_weight_numbers (0, &weight_numbers[0]);
    NCX2_CX3 (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* forward z plane: z = zlocal + 2*Delta : sup2deltaz */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |x|x|x|_|_| --> |24|25|26|_|_| */
  /* |x|x|x|_|_| --> |23|22|21|_|_| */ 
  /* |x|x|x|_|_| --> |18|19|20|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  
  if ((relnl.z < -1.9*delta) && (relnl.z > -2.1*delta)) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX2_CX3 (relnl, delta, weight_numbers, weight_id, goflag);
  }  
}
//recoded on 28-02-2019

void NCX20_CX40 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containning the Lagrange multiplier */
  
  /* z = localplane (plane containning the multiplier) */
  /* backward z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|x|x|x|_| --> |_|6|7|8|_| */
  /* |_|x|O|x|_| --> |_|5|4|3|_| */ 
  /* |_|x|x|x|_| --> |_|0|1|2|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX_centred (relnl, delta, weight_numbers, weight_id, goflag);
  }

  // middle z plane: z = zlocal + Delta : sup1deltaz
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|15|16|17|_| */
  /* |_|x|x|x|_| --> |_|14|13|12|_| */ 
  /* |_|x|x|x|_| --> |_|_9|10|11|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  if ((relnl.z < -0.9*delta) && (relnl.z > -1.1*delta)) {
    fill_weight_numbers (0, &weight_numbers[0]);
    NCX_centred (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* forward z plane: z = zlocal + 2*Delta : sup2deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|24|25|26|_| */
  /* |_|x|x|x|_| --> |_|23|22|21|_| */ 
  /* |_|x|x|x|_| --> |_|18|19|20|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */

  if ((relnl.z < -1.9*delta) && (relnl.z > -2.1*delta)) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX_centred (relnl, delta, weight_numbers, weight_id, goflag);
  }  
}
//recoded on 28-02-2019

void NCX20_CX21 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};

  /* Stencil config, O is the cell containning the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - Delta : inf1deltaz*/
  /* |x|x|x|_|_| --> |6|7|8|_|_| */
  /* |x|x|x|_|_| --> |5|4|3|_|_| */
  /* |x|x|x|_|_| --> |0|1|2|_|_| */ 
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if ((relnl.z < 1.1*delta) && (relnl.z > 0.9*delta)) {
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX2_CX2 (relnl, delta, weight_numbers, weight_id, goflag);
  }

  /* middle z plane: z = zlocal : sup0deltaz */
  /* |x|x|x|_|_| --> |15|16|17|_|_| */
  /* |x|x|x|_|_| --> |14|13|12|_|_| */
  /* |x|x|O|_|_| --> |_9|10|11|_|_| */ 
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (0, &weight_numbers[0]);
    NCX2_CX2 (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* forward z plane: z = zlocal + Delta : sup1deltaz */
  /* |x|x|x|_|_| --> |24|25|26|_|_| */
  /* |x|x|x|_|_| --> |23|22|21|_|_| */
  /* |x|x|x|_|_| --> |18|19|20|_|_| */ 
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  
  if ((relnl.z < -0.9*delta) && (relnl.z > -1.1*delta)) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX2_CX2 (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

void NCX20_CX31 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};

  /* Stencil config, O is the cell containning the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - Delta : inf1deltaz*/
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |x|x|x|_|_| --> |6|7|8|_|_| */
  /* |x|x|x|_|_| --> |5|4|3|_|_| */ 
  /* |x|x|x|_|_| --> |0|1|2|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if ((relnl.z < 1.1*delta) && (relnl.z > 0.9*delta)) {
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX2_CX3 (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* middle z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |x|x|x|_|_| --> |15|16|17|_|_| */
  /* |x|x|x|_|_| --> |14|13|12|_|_| */ 
  /* |x|x|x|_|_| --> |_9|10|11|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (0, &weight_numbers[0]);
    NCX2_CX3 (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* forward z plane: z = zlocal + Delta : sup1deltaz */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |x|x|x|_|_| --> |24|25|26|_|_| */
  /* |x|x|x|_|_| --> |23|22|21|_|_| */ 
  /* |x|x|x|_|_| --> |18|19|20|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  
  if ((relnl.z < -0.9*delta) && (relnl.z > -1.1*delta)) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX2_CX3 (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

void NCX20_CX41 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};

  /* Stencil config, O is the cell containning the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - Delta : inf1deltaz*/
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|x|x|x|_| --> |_|6|7|8|_| */
  /* |_|x|x|x|_| --> |_|5|4|3|_| */ 
  /* |_|x|x|x|_| --> |_|0|1|2|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if ((relnl.z < 1.1*delta) && (relnl.z > 0.9*delta)) {
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX_centred (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* middle z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|15|16|17|_| */
  /* |_|x|x|x|_| --> |_|14|13|12|_| */ 
  /* |_|x|x|x|_| --> |_|_9|10|11|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (0, &weight_numbers[0]);
    NCX_centred (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* forward z plane: z = zlocal + Delta : sup1deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|24|25|26|_| */
  /* |_|x|x|x|_| --> |_|23|22|21|_| */ 
  /* |_|x|x|x|_| --> |_|18|19|20|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  if ((relnl.z < -0.9*delta) && (relnl.z > -1.1*delta)) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX_centred (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

// NCX21_CX Bloc 3 functions

/* NCX21_CX30 = NCX20_CX31 */

/* NCX21_CX40 = NCX20_CX41 */

void NCX21_CX21 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containning the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - 2Delta : inf2deltaz */
  /* |x|x|x|_|_| --> |6|7|8|_|_| */
  /* |x|x|x|_|_| --> |5|4|3|_|_| */
  /* |x|x|x|_|_| --> |0|1|2|_|_| */ 
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */

  if ((relnl.z < 2.1*delta) && (relnl.z > 1.9*delta)) {
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX2_CX2 (relnl, delta, weight_numbers, weight_id, goflag);
  }

  /* middle z plane: z = zlocal - Delta : inf1deltaz */
  /* |x|x|x|_|_| --> |15|16|17|_|_| */
  /* |x|x|x|_|_| --> |14|13|12|_|_| */
  /* |x|x|x|_|_| --> |_9|10|11|_|_| */ 
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  
  if ((relnl.z < 1.1*delta) && (relnl.z > 0.9*delta)) {
    fill_weight_numbers (0,  &weight_numbers[0]);
    NCX2_CX2 (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* forward z plane: z = zlocal : sup0deltaz */
  /* |x|x|x|_|_| --> |24|25|26|_|_| */
  /* |x|x|x|_|_| --> |23|22|21|_|_| */
  /* |x|x|O|_|_| --> |18|19|20|_|_| */ 
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX2_CX2 (relnl, delta, weight_numbers, weight_id, goflag);
  }
}

void NCX21_CX31 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containning the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - 2Delta : inf2deltaz */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |x|x|x|_|_| --> |6|7|8|_|_| */
  /* |x|x|x|_|_| --> |5|4|3|_|_| */ 
  /* |x|x|x|_|_| --> |0|1|2|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if ((relnl.z < 2.1*delta) && (relnl.z > 1.9*delta)) {
    fill_weight_numbers (-1,  &weight_numbers[0]);
    NCX2_CX3(relnl, delta, weight_numbers, weight_id, goflag);
  }

  // middle z plane: z = zlocal - Delta : inf1deltaz
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |x|x|x|_|_| --> |15|16|17|_|_| */
  /* |x|x|x|_|_| --> |14|13|12|_|_| */ 
  /* |x|x|x|_|_| --> |_9|10|11|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  
  if ((relnl.z < 1.1*delta) && (relnl.z > 0.9*delta)) {
    fill_weight_numbers (0,  &weight_numbers[0]);
    NCX2_CX3 (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* forward z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |x|x|x|_|_| --> |24|25|26|_|_| */
  /* |x|x|O|_|_| --> |23|22|21|_|_| */ 
  /* |x|x|x|_|_| --> |18|19|20|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX2_CX3 (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

void NCX21_CX41 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containning the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - 2Delta : inf2deltaz */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|x|x|x|_| --> |_|6|7|8|_| */
  /* |_|x|x|x|_| --> |_|5|4|3|_| */ 
  /* |_|x|x|x|_| --> |_|0|1|2|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if ((relnl.z < 2.1*delta) && (relnl.z > 1.9*delta)) {
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX_centred (relnl, delta, weight_numbers, weight_id, goflag);
  }

  // middle z plane: z = zlocal - Delta : inf1deltaz
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|15|16|17|_| */
  /* |_|x|x|x|_| --> |_|14|13|12|_| */ 
  /* |_|x|x|x|_| --> |_|_9|10|11|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  if ((relnl.z < 1.1*delta) && (relnl.z > 0.9*delta)) {
    fill_weight_numbers (0, &weight_numbers[0]);
    NCX_centred (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  // forward z plane: z = zlocal : sup0deltaz
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|24|25|26|_| */
  /* |_|x|O|x|_| --> |_|23|22|21|_| */ 
  /* |_|x|x|x|_| --> |_|18|19|20|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX_centred (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

//NCX30_CX Bloc, 4 functions
void NCX30_CX30 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containning the Lagrange multiplier */
  
  /* z = localplane (plane containning the multiplier) */
  /* backward z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |x|x|O|_|_| --> |6|7|8|_|_| */ 
  /* |x|x|x|_|_| --> |5|4|3|_|_| */
  /* |x|x|x|_|_| --> |0|1|2|_|_| */
  
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX3_CX3 (relnl, delta, weight_numbers, weight_id, goflag);
  }

  /* middle z plane: z = zlocal + Delta : sup1delta */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |x|x|x|_|_| --> |15|16|17|_|_| */ 
  /* |x|x|x|_|_| --> |14|13|12|_|_| */
  /* |x|x|x|_|_| --> |_9|10|11|_|_| */
  
  if ((relnl.z < -0.9*delta) && (relnl.z > -1.1*delta)) {
    fill_weight_numbers (0, &weight_numbers[0]);
    NCX3_CX3 (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* forward z plane: z = zlocal + 2*Delta : sup2delta*/
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |x|x|x|_|_| --> |24|25|26|_|_| */ 
  /* |x|x|x|_|_| --> |23|22|21|_|_| */
  /* |x|x|x|_|_| --> |18|19|20|_|_| */
  
  if ((relnl.z < -1.9*delta) && (relnl.z > -2.1*delta)) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX3_CX3 (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

void NCX30_CX40 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containning the Lagrange multiplier */
  
  /* z = localplane (plane containning the multiplier) */
  /* backward z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|x|O|x|_| --> |_|6|7|8|_| */ 
  /* |_|x|x|x|_| --> |_|5|4|3|_| */
  /* |_|x|x|x|_| --> |_|0|1|2|_| */

  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX3_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }
 
  /* middle z plane: z = zlocal + Delta : sup1deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|15|16|17|_| */ 
  /* |_|x|x|x|_| --> |_|14|13|12|_| */
  /* |_|x|x|x|_| --> |_|_9|10|11|_| */
  
  if ((relnl.z < -0.9*delta) && (relnl.z > -1.1*delta)) {
    fill_weight_numbers (0, &weight_numbers[0]);
    NCX3_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* forward z plane: z = zlocal + 2*Delta : sup2deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|24|25|26|_| */ 
  /* |_|x|x|x|_| --> |_|23|22|21|_| */
  /* |_|x|x|x|_| --> |_|18|19|20|_| */
  
  if ((relnl.z < -1.9*delta) && (relnl.z > -2.1*delta)) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX3_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

void NCX30_CX31 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};

  /* Stencil config, O is the cell containning the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - Delta : inf1deltaz*/
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |x|x|x|_|_| --> |6|7|8|_|_| */ 
  /* |x|x|x|_|_| --> |5|4|3|_|_| */
  /* |x|x|x|_|_| --> |0|1|2|_|_| */
  
  if ((relnl.z < 1.1*delta) && (relnl.z > 0.9*delta)) {
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX3_CX3 (relnl, delta, weight_numbers, weight_id, goflag);
  }

  /* middle z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |x|x|x|_|_| --> |15|16|17|_|_| */ 
  /* |x|x|x|_|_| --> |14|13|12|_|_| */
  /* |x|x|x|_|_| --> |_9|10|11|_|_| */
  
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (0, &weight_numbers[0]);
    NCX3_CX3 (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* forward z plane: z = zlocal + Delta : sup1deltaz */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |x|x|x|_|_| --> |24|25|26|_|_| */ 
  /* |x|x|x|_|_| --> |23|22|21|_|_| */
  /* |x|x|x|_|_| --> |18|19|20|_|_| */
  
  if ((relnl.z < -0.9*delta) && (relnl.z > -1.1*delta)) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX3_CX3 (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

void NCX30_CX41 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};

  /* Stencil config, O is the cell containning the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - Delta : inf1deltaz*/
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|x|x|x|_| --> |_|6|7|8|_| */ 
  /* |_|x|x|x|_| --> |_|5|4|3|_| */
  /* |_|x|x|x|_| --> |_|0|1|2|_| */
  
  if ((relnl.z < 1.1*delta) && (relnl.z > 0.9*delta)) {
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX3_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }

  /* middle z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|15|16|17|_| */ 
  /* |_|x|O|x|_| --> |_|14|13|12|_| */
  /* |_|x|x|x|_| --> |_|_9|10|11|_| */
  
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (0,  &weight_numbers[0]);
    NCX3_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* forward z plane: z = zlocal + Delta : sup1deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|24|25|26|_| */ 
  /* |_|x|x|x|_| --> |_|23|22|21|_| */
  /* |_|x|x|x|_| --> |_|18|19|20|_| */
  
  if ((relnl.z < -0.9*delta) && (relnl.z > -1.1*delta)) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX3_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

//NCX31_CX Bloc 3 functions

/* NCX31_CX40 = NCX30_CX41 */

void NCX31_CX31 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containning the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - 2Delta : inf2deltaz */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |x|x|x|_|_| --> |6|7|8|_|_| */ 
  /* |x|x|x|_|_| --> |5|4|3|_|_| */
  /* |x|x|x|_|_| --> |0|1|2|_|_| */
  
  if ((relnl.z < 2.1*delta) && (relnl.z > 1.9*delta)) {
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX3_CX3 (relnl, delta, weight_numbers, weight_id, goflag);
  }

  /* middle z plane: z = zlocal - Delta : inf1deltaz */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |x|x|x|_|_| --> |15|16|17|_|_| */ 
  /* |x|x|x|_|_| --> |14|13|12|_|_| */
  /* |x|x|x|_|_| --> |_9|10|11|_|_| */
  
  if ((relnl.z < 1.1*delta) && (relnl.z > 0.9*delta)) {
    fill_weight_numbers (0, &weight_numbers[0]);
    NCX3_CX3 (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  // forward z plane: z = zlocal 
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |x|x|O|_|_| --> |24|25|26|_|_| */ 
  /* |x|x|x|_|_| --> |23|22|21|_|_| */
  /* |x|x|x|_|_| --> |18|19|20|_|_| */
  
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX3_CX3 (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

void NCX31_CX41 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containning the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - 2Delta : inf2deltaz */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|x|x|x|_| --> |_|6|7|8|_| */ 
  /* |_|x|x|x|_| --> |_|5|4|3|_| */
  /* |_|x|x|x|_| --> |_|0|1|2|_| */
  
  if ((relnl.z < 2.1*delta) && (relnl.z > 1.9*delta)) {
    fill_weight_numbers(-1, &weight_numbers[0]);
    NCX3_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }

  /* middle z plane: z = zlocal - Delta : inf1deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|15|16|17|_| */ 
  /* |_|x|x|x|_| --> |_|14|13|12|_| */
  /* |_|x|x|x|_| --> |_|_9|10|11|_| */
  
  if ((relnl.z < 1.1*delta) && (relnl.z > 0.9*delta)) {
    fill_weight_numbers (0, &weight_numbers[0]);
    NCX3_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* forward z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|24|25|26|_| */ 
  /* |_|x|x|x|_| --> |_|23|22|21|_| */
  /* |_|x|x|x|_| --> |_|18|19|20|_| */
  
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX3_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

//NCX40_CX Bloc 2 functions
void NCX40_CX40 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containning the Lagrange multiplier */
  
  /* z = localplane (plane containning the multiplier) */
  /* backward z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|O|x|x| --> |_|_|6|7|8| */ 
  /* |_|_|x|x|x| --> |_|_|5|4|3| */
  /* |_|_|x|x|x| --> |_|_|0|1|2| */
  
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX4_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }

  // middle z plane: z = zlocal + Delta : sup1deltaz 
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|O|x|x| --> |_|_|15|16|17| */ 
  /* |_|_|x|x|x| --> |_|_|14|13|12| */
  /* |_|_|x|x|x| --> |_|_|_9|10|11| */
  
  if ((relnl.z < -0.9*delta) && (relnl.z > -1.1*delta)) {
    fill_weight_numbers (0, &weight_numbers[0]);
    NCX4_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* forward z plane: z = zlocal + 2*Delta : sup2deltaz */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|O|x|x| --> |_|_|24|25|26| */ 
  /* |_|_|x|x|x| --> |_|_|23|22|21| */
  /* |_|_|x|x|x| --> |_|_|18|19|20| */
  
  if ((relnl.z < -1.9*delta) && (relnl.z > -2.1*delta)) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX4_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

void NCX40_CX41 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};

  /* Stencil config, O is the cell containning the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - Delta : inf1deltaz*/
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|x|x|x| --> |_|_|6|7|8| */ 
  /* |_|_|x|x|x| --> |_|_|5|4|3| */
  /* |_|_|x|x|x| --> |_|_|0|1|2| */
 
  if ((relnl.z < 1.1*delta) && (relnl.z > 0.9*delta)) {
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX4_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }

  /* middle z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|O|x|x| --> |_|_|15|16|17| */ 
  /* |_|_|x|x|x| --> |_|_|14|13|12| */
  /* |_|_|x|x|x| --> |_|_|_9|10|11| */
  
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (0, &weight_numbers[0]);
    NCX4_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  // forward z plane: z = zlocal + Delta : sup1deltaz
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|x|x|x| --> |_|_|24|25|26| */ 
  /* |_|_|x|x|x| --> |_|_|23|22|21| */
  /* |_|_|x|x|x| --> |_|_|18|19|20| */
  
  if ((relnl.z < -0.9*delta) && (relnl.z > -1.1*delta)) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX4_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019

//NCX41_CX Bloc 1 function
void NCX41_CX41 (const coord relnl, const double delta, int * weight_id, size_t * goflag) {
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containning the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - 2Delta : inf2deltaz */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|x|x|x| --> |_|_|6|7|8| */ 
  /* |_|_|x|x|x| --> |_|_|5|4|3| */
  /* |_|_|x|x|x| --> |_|_|0|1|2| */
  
  if ((relnl.z < 2.1*delta) && (relnl.z > 1.9*delta)) {
    fill_weight_numbers (-1, &weight_numbers[0]);
    NCX4_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }

  /* middle z plane: z = zlocal - Delta : inf1deltaz */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|x|x|x| --> |_|_|15|16|17| */ 
  /* |_|_|x|x|x| --> |_|_|14|13|12| */
  /* |_|_|x|x|x| --> |_|_|_9|10|11| */
  
  if ((relnl.z < 1.1*delta) && (relnl.z > 0.9*delta)) {
    fill_weight_numbers (0, &weight_numbers[0]);
    NCX4_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  /* forward z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|O|x|x| --> |_|_|24|25|26| */ 
  /* |_|_|x|x|x| --> |_|_|23|22|21| */
  /* |_|_|x|x|x| --> |_|_|18|19|20| */
  
  if (fabs(relnl.z) < 0.1*delta) {
    fill_weight_numbers (1, &weight_numbers[0]);
    NCX4_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }
}
//recoded on 28-02-2019
#endif

void assign_weight_id_quad_outward (const int NCX, const int CX, const coord relnl, const double Delta, int * weight_id, size_t * goflag) {
  double delta = Delta;
  
#if dimension == 2
  // Table taken from Peligriff
  /* (*Q2numb)(0,0) = 0 ;  */
  /* (*Q2numb)(1,0) = 4 ;  */
  /* (*Q2numb)(2,0) = 1 ;  */
  
  /* (*Q2numb)(0,1) = 7 ;  */
  /* (*Q2numb)(1,1) = 8 ;  */
  /* (*Q2numb)(2,1) = 5 ;  */
  
  /* (*Q2numb)(0,2) = 3 ;  */
  /* (*Q2numb)(1,2) = 6 ;  */
  /* (*Q2numb)(2,2) = 2 ;  */

  // weight_numbers[2] -> right column top
  // weight_numbers[5] -> right column middle
  // weight_numbers[1] -> right column bottom

  // weight_numbers[6] -> middle column top
  // weight_numbers[8] -> middle colum middle
  // weight_numbers[4] -> middle column bottom

  // weight_numbers[3] -> left column top
  // weight_numbers[7] -> left column middle
  // weight_numbers[0] -> left column bottom

  int weight_numbers[9];
  weight_numbers[3] = 3; weight_numbers[6] = 6; weight_numbers[2] = 2;
  weight_numbers[7] = 7; weight_numbers[8] = 8; weight_numbers[5] = 5;
  weight_numbers[0] = 0; weight_numbers[4] = 4; weight_numbers[1] = 1;

  if (NCX == 1) {
    // fictitous-domain-s boundary is oriented (x+,y+)
    if (CX == 1)
      // rel vector is oriented (x+,y+)
      NCX1_CX1 (relnl, delta, weight_numbers, weight_id, goflag);      

    if(CX == 2)
      // rel vector is oriented (x-,y+)
      NCX1_CX2 (relnl, delta, weight_numbers, weight_id, goflag);

    if(CX == 3)
      // rel vector is oriented (x-,y-)
      NCX_centred (relnl, delta, weight_numbers, weight_id, goflag);
    
    if(CX == 4)
      // rel vector is oriented (x-,y-)
      NCX1_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  } 

  if (NCX == 2) {
    // fictitous-domain-s boundary is oriented (x-,y+)
    if (CX == 1)
      // rel vector is oriented (x+,y+)
      // NCX2_CX1 == NCX1_CX2
      NCX1_CX2 (relnl, delta, weight_numbers, weight_id, goflag);      
    if (CX == 2)
      // rel vector is oriented (x-,y+)
      NCX2_CX2 (relnl, delta, weight_numbers, weight_id, goflag);
    if (CX == 3)
      // rel vector is oriented (x-,y-)
      NCX2_CX3 (relnl, delta, weight_numbers, weight_id, goflag);
    if (CX == 4)
      // rel vector is oriented (x-,y-)
      NCX_centred (relnl, delta, weight_numbers, weight_id, goflag);
  } 

  if (NCX == 3) {
    // fictitous-domain-s boundary is oriented (x-,y-)
    if(CX == 1)
      // rel vector is oriented (x+,y+)
      // NCX3_CX1 == NCX1_CX3 == NCX_centred
      NCX_centred (relnl, delta, weight_numbers, weight_id, goflag);

    if(CX == 2)
      // rel vector is oriented (x-,y+)
      // NCX3_CX2 == NCX2_CX3
      NCX2_CX3(relnl, delta, weight_numbers, weight_id, goflag);
    
    if(CX == 3)
      // rel vector is oriented (x-,y-)
      NCX3_CX3 (relnl, delta, weight_numbers, weight_id, goflag);
    
    if(CX == 4)
      // rel vector is oriented (x+,y-)
      NCX3_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }
  
  if (NCX == 4) {
    // fictitous-domain-s boundary is oriented (x+,y-)
    if (CX == 1)
      // rel vector is oriented (x+,y+)
      // NCX4_CX1 = NCX1_CX4
      NCX1_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  
    if (CX == 2)
      // rel vector is oriented (x-,y+)
      // NCX4_CX2 = NCX2_CX4 = NCX_centred
      NCX_centred (relnl, delta, weight_numbers, weight_id, goflag);
    
    if (CX == 3)
      // rel vector is oriented (x-,y-)
      // NCX4_CX3 == NCX3_CX4
      NCX3_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
    
    if (CX == 4)
      // rel vector is oriented (x+,y-)
      NCX4_CX4 (relnl, delta, weight_numbers, weight_id, goflag);
  }
#endif
  
#if dimension == 3

  if (NCX == 10) {
    // fictitous-domain-s boundary is oriented (x+,y+,z+)

    if (CX == 10) {
      // rel vector is oriented (x+,y+,z+)
      NCX10_CX10 (relnl, delta, weight_id, goflag);
    }
    if (CX == 20) {
      // rel vector is oriented (x-,y+,z+)
      NCX10_CX20 (relnl, delta, weight_id, goflag);
    }  
    if (CX == 30) {
      // rel vector is oriented (x-,y-,z+)
      NCX10_CX30 (relnl, delta, weight_id, goflag);
    }
    if (CX == 40) {
      // rel vector is oriented (x+,y-,z+)
      NCX10_CX40 (relnl, delta, weight_id, goflag);
    }

    if (CX == 11) {
      // rel vector is oriented (x+,y+,z-)
      NCX10_CX11 (relnl, delta, weight_id, goflag);
    }
    if (CX == 21) {
      // rel vector is oriented (x-,y+,z-)
      NCX10_CX21 (relnl, delta, weight_id, goflag);
    }  
    if (CX == 31) {
      // rel vector is oriented (x-,y-,z-)
      NCX10_CX31 (relnl, delta, weight_id, goflag);
    }
    if (CX == 41) {
      // rel vector is oriented (x+,y-,z-)
      NCX10_CX41 (relnl, delta, weight_id, goflag);
    }
  }
  // NCX = 10, ok no bug
  // ====================================================//
  
  if (NCX == 11) {
    // fictitous-domain-s boundary is oriented (x+,y+,z-)

    if (CX == 10) {
      // rel vector is oriented (x+,y+,z+)
      //NCX11_CX10 = NCX10_CX11
      NCX10_CX11 (relnl, delta, weight_id, goflag);
    }
    if (CX == 20) {
      // rel vector is oriented (x-,y+,z+)
      /* NCX11_CX20 = NCX10_CX21 */
      NCX10_CX21 (relnl, delta, weight_id, goflag);
    }
    if (CX == 30) {
      // rel vector is oriented (x-,y-,z+)
      /* NCX11_CX30 = NCX10_CX31 */
      NCX10_CX31 (relnl, delta, weight_id, goflag);
    }
    if (CX == 40) {
      // rel vector is oriented (x+,y-,z+)
      /* NCX11_CX40 = NCX10_CX41 */
      NCX10_CX41 (relnl, delta, weight_id, goflag);
    }

    if (CX == 11) {
      // rel vector is oriented (x+,y+,z-)
      NCX11_CX11 (relnl, delta, weight_id, goflag);
    }
    if (CX == 21) {
      // rel vector is oriented (x-,y+,z-)
      NCX11_CX21 (relnl, delta, weight_id, goflag);
    }
    if (CX == 31) {
      // rel vector is oriented (x-,y-,z-)
      NCX11_CX31 (relnl, delta, weight_id, goflag);
    }
    if (CX == 41) {
      // rel vector is oriented (x+,y-,z-)
      NCX11_CX41 (relnl, delta, weight_id, goflag);
    }
  }
  // NCX = 11, ok no bug
  // ====================================================//
  
  if (NCX == 20) {
    // fictitous-domain-s boundary is oriented (x-,y+,z+)

    if (CX == 10) {
      // rel vector is oriented (x+,y+,z+)
      // NCX20_CX10 = NCX10_CX20
      NCX10_CX20 (relnl, delta, weight_id, goflag);
    }
    if (CX == 20) {
      // rel vector is oriented (x-,y+,z+)
      NCX20_CX20 (relnl, delta, weight_id, goflag);
    }  
    if (CX == 30) {
      // rel vector is oriented (x-,y-,z+)
      NCX20_CX30 (relnl, delta, weight_id, goflag);
    }
    if (CX == 40) {
      // rel vector is oriented (x+,y-,z+)
      NCX20_CX40 (relnl, delta, weight_id, goflag);
    }

    if (CX == 11) {
      // rel vector is oriented (x+,y+,z-)
      //NCX20_CX11 = NCX11_CX20 = NCX10_CX21
       NCX10_CX21 (relnl, delta, weight_id, goflag);
    }
    if (CX == 21) {
      // rel vector is oriented (x-,y+,z-)
      NCX20_CX21 (relnl, delta, weight_id, goflag);
    }
    if (CX == 31) {
      // rel vector is oriented (x-,y-,z-)
      NCX20_CX31 (relnl, delta, weight_id, goflag);
    }
    if (CX == 41) {
      // rel vector is oriented (x+,y-,z-)
      NCX20_CX41 (relnl, delta, weight_id, goflag);
    }
  }
  // NCX = 20, ok no bug
  // ====================================================//
  
  if (NCX == 21) {
    // fictitous-domain-s boundary is oriented (x-,y+,z-)
    if (CX == 10) {
      // rel vector is oriented (x+,y+,z+)
      //NCX21_CX10 = NCX10_CX21
      NCX10_CX21 (relnl, delta, weight_id, goflag);
    }
    if (CX == 20) {
      // rel vector is oriented (x-,y+,z+)
      //NCX21_CX20 = NCX20_CX21
      NCX20_CX21 (relnl, delta, weight_id, goflag);
    }
    if (CX == 30) {
      // rel vector is oriented (x-,y-,z+)
      /* NCX21_CX30 = NCX20_CX31 */
      NCX20_CX31 (relnl, delta, weight_id, goflag);
    }
    if (CX == 40) {
      // rel vector is oriented (x+,y-,z+)
      /* NCX21_CX40 = NCX20_CX41 */
      NCX20_CX41 (relnl, delta, weight_id, goflag);
    }
    if (CX == 11) {
      // rel vector is oriented (x+,y+,z-)
      // NCX21_CX11 = NCX11_CX21
      NCX11_CX21 (relnl, delta, weight_id, goflag);
    }
    if (CX == 21) {
      // rel vector is oriented (x-,y+,z-)
      NCX21_CX21 (relnl, delta, weight_id, goflag);
    }
    if (CX == 31) {
      // rel vector is oriented (x-,y-,z-)
      NCX21_CX31 (relnl, delta, weight_id, goflag);
    }
    if (CX == 41) {
      // rel vector is oriented (x+,y-,z-)
      NCX21_CX41 (relnl, delta, weight_id, goflag);
    }
  }
  // NCX = 21, ok no bug
  // ====================================================//
  
  if (NCX == 30) {
    // fictitous-domain-s boundary is oriented (x-,y-,z+)

    if (CX == 10) {
      // rel vector is oriented (x+,y+,z+)
      // NCX30_CX10 = NCX10_CX30
      NCX10_CX30 (relnl, delta, weight_id, goflag);
    }
    if (CX == 20) {
      // rel vector is oriented (x-,y+,z+)
      // NCX30_CX20 = NCX20_CX30
      NCX20_CX30 (relnl, delta, weight_id, goflag);
    }
    if (CX == 30) {
      // rel vector is oriented (x-,y-,z+)
      NCX30_CX30 (relnl, delta, weight_id, goflag);
    }
    if (CX == 40) {
      // rel vector is oriented (x+,y-,z+)
      NCX30_CX40 (relnl, delta, weight_id, goflag);
    }
    
    if (CX == 11) {
      // rel vector is oriented (x+,y+,z-)
      // NCX30_CX11 = NCX11_CX30 = NCX10_CX31
      NCX10_CX31 (relnl, delta, weight_id, goflag);
    }
    if (CX == 21) {
      // rel vector is oriented (x-,y+,z-)
      //NCX30_CX21 = NCX21_CX30 = NCX20_CX31
      NCX20_CX31 (relnl, delta, weight_id, goflag);
    }
    if (CX == 31) {
      // rel vector is oriented (x-,y-,z-)
      NCX30_CX31 (relnl, delta,  weight_id, goflag);
    }
    if (CX == 41) {
      // rel vector is oriented (x+,y-,z-)
      NCX30_CX41 (relnl, delta, weight_id, goflag);
    }
  }
  // NCX = 30, ok no bug
  // ====================================================//
  
  if (NCX == 31) {
    // fictitous-domain-s boundary is oriented (x-,y-,z-)

    if (CX == 10) {
      // rel vector is oriented (x+,y+,z+)
      //NCX31_CX10 = NCX10_CX31
      NCX10_CX31 (relnl, delta, weight_id, goflag);
    }
    if (CX == 20) {
      // rel vector is oriented (x-,y+,z+)
      //NCX31_CX20  = NCX20_CX31
      NCX20_CX31 (relnl, delta, weight_id, goflag);
    }
    if (CX == 30) {
      // rel vector is oriented (x-,y-,z+)
      //NCX31_CX30 = NCX30_CX31
      NCX30_CX31 (relnl, delta, weight_id, goflag);
    }
    if (CX == 40) {
      // rel vector is oriented (x+,y-,z+)
      /* NCX31_CX40 = NCX30_CX41 */
      NCX30_CX41 (relnl, delta, weight_id, goflag);
    }

    if (CX == 11) {
      // rel vector is oriented (x+,y+,z-)
      // NCX31_CX11 = NCX11_CX31
      NCX11_CX31 (relnl, delta, weight_id, goflag);
    }
    if (CX == 21) {
      // rel vector is oriented (x-,y+,z-)
      // NCX31_CX21 = NCX21_CX31
      NCX21_CX31 (relnl, delta, weight_id, goflag);
    }
    if (CX == 31) {
      // rel vector is oriented (x-,y-,z-)
      NCX31_CX31 (relnl, delta, weight_id, goflag);
    }
    if (CX == 41) {
      // rel vector is oriented (x+,y-,z-)
      NCX31_CX41 (relnl, delta, weight_id, goflag);
    }
  }
  // NCX = 31, ok no bug
  // ====================================================//
  
  if (NCX == 40) {
    // fictitous-domain-s boundary is oriented (x+,y-,z+)

    if (CX == 10) {
      // rel vector is oriented (x+,y+,z+)
      // NCX40_CX10 = NCX10_CX40
      NCX10_CX40 (relnl, delta, weight_id, goflag);
    }
    if (CX == 20) {
      // rel vector is oriented (x-,y+,z+)
      // NCX40_CX20 = NCX20_CX40
      NCX20_CX40 (relnl, delta, weight_id, goflag);
    }
    if (CX == 30) {
      // rel vector is oriented (x-,y-,z+)
      // NCX40_CX30 = NCX30_CX40
      NCX30_CX40 (relnl, delta, weight_id, goflag);
    }
    if (CX == 40) {
      // rel vector is oriented (x+,y-,z+)
      NCX40_CX40 (relnl, delta, weight_id, goflag);
    }

    if (CX == 11) {
      // rel vector is oriented (x+,y+,z-)
      // NCX40_CX11 = NCX11_CX40 = NCX10_CX41
      NCX10_CX41 (relnl, delta, weight_id, goflag);
    }
    if (CX == 21) {
      // rel vector is oriented (x-,y+,z-)
      // NCX40_CX21 = NCX21_CX40 = NCX20_CX41
      NCX20_CX41 (relnl, delta, weight_id, goflag);
    }
    if (CX == 31) {
      // rel vector is oriented (x-,y-,z-)
      // NCX40_CX31 = NCX31_CX40 = NCX30_CX41
      NCX30_CX41 (relnl, delta, weight_id, goflag);
    }
    if (CX == 41) {
      // rel vector is oriented (x+,y-,z-)
      NCX40_CX41 (relnl, delta, weight_id, goflag);
    }
  }
  // NCX = 40, ok no bug
  // ====================================================//
  
  if (NCX == 41) {
    // fictitous-domain-s boundary is oriented (x+,y-,z-)

    if (CX == 10) {
      // rel vector is oriented (x+,y+,z+)
      // NCX41_CX10 = NCX10_CX41
       NCX10_CX41 (relnl, delta, weight_id, goflag);
    }
    if (CX == 20) {
      // rel vector is oriented (x-,y+,z+)
      // NCX41_CX20 = NCX20_CX41
      NCX20_CX41 (relnl, delta, weight_id, goflag);
    }
    if (CX == 30) {
      // rel vector is oriented (x-,y-,z+)
      // NCX41_CX30 = NCX30_CX41
      NCX30_CX41 (relnl, delta, weight_id, goflag);
    }
    if (CX == 40) {
      // rel vector is oriented (x+,y-,z+)
      // NCX41_CX40 = NCX40_CX41
      NCX40_CX41 (relnl, delta, weight_id, goflag);
    }

    if (CX == 11) {
      // rel vector is oriented (x+,y+,z-)
      // NCX41_CX11 = NCX11_CX41
      NCX11_CX41 (relnl, delta, weight_id, goflag);
    }
    if (CX == 21) {
      // rel vector is oriented (x-,y+,z-)
      // NCX41_CX21 = NCX21_CX41
      NCX21_CX41 (relnl, delta, weight_id, goflag);
    }
    if (CX == 31) {
      // rel vector is oriented (x-,y-,z-)
      // NCX41_CX31 = NCX31_CX41
      NCX31_CX41(relnl, delta, weight_id, goflag);
    }
    if (CX == 41) {
      // rel vector is oriented (x+,y-,z-)
      NCX41_CX41 (relnl, delta, weight_id, goflag);
    }
  }
#endif
}

double Q2weighting (size_t  i, double  x, double y, double z ) { 
  double result = 0. ;

#if dimension == 2
  
  switch( i )
  {
      case 0 :
         result = (1.0-x) * (1.0-2.0*x) * (1.0-y) * (1.0-2.0*y) ;
         break ;
      case 1 :
         result = x * (2.0*x-1.0) * (1.0-y) * (1.0-2.0*y) ;
         break ;
      case 2 :
         result = x * (2.0*x-1.0) * y * (2.0*y-1.0) ;
         break ;
      case 3 :
         result = (1.0-x) * (1.0-2.0*x) * y * (2.0*y-1.0) ;
         break ;
      case 4 :
         result = 4.0 * x * (1.0-x) * (1.0-y) * (1.0-2.0*y) ;
         break ;
      case 5 :
         result = x * (2.0*x-1.0) * 4.0 * y * (1.0-y) ;
         break ;
      case 6 :
         result = 4.0 * x * (1.0-x) * y * (2.0*y-1.0) ;
         break ;
      case 7 :
         result = (1.0-x) * (1.0-2.0*x) * 4.0 * y * (1.0-y) ;
         break ;
      case 8 :
         result = 4.0 * x * (1.0-x) * 4.0 * y * (1.0-y) ;
         break ;
  }
#endif

#if dimension == 3
  
  switch( i )
  {
      case 0 :
         result = 8.0 * (1.0-x)*(0.5-x) * (1.0-y)*(0.5-y) * (1.0-z)*(0.5-z);
         break ;
      case 1 :
         result = 16.0 * x*(1.0-x) * (1.0-y)*(0.5-y) * (1.0-z)*(0.5-z);
         break ;
      case 2 :
         result = -8.0 * x*(0.5-x) * (1.0-y)*(0.5-y) * (1.0-z)*(0.5-z) ;
         break ;
      case 3 :
         result = -16.0 * x*(0.5-x) * y*(1.0-y) * (1.0-z)*(0.5-z);
         break ;
      case 4 :
         result = 32.0 * x*(1.0-x) * y*(1.0-y) * (1.0-z)*(0.5-z) ;
         break ;
      case 5 :
         result = 16.0 * (1.0-x)*(0.5-x) * y*(1.0-y) * (1.0-z)*(0.5-z) ;
         break ;
      case 6 :
         result =-8.0 * (1.0-x)*(0.5-x) * y*(0.5-y) * (1.0-z)*(0.5-z) ;
         break ;
      case 7 :
         result = -16.0 * x*(1.0-x) * y*(0.5-y) * (1.0-z)*(0.5-z);
         break ;
      case 8 :
         result = 8.0 * x*(0.5-x) * y*(0.5-y) * (1.0-z)*(0.5-z);
         break ;
      case 9 :
         result = 16.0 * (1.0-x)*(0.5-x) * (1.0-y)*(0.5-y) * z*(1.0-z) ;
         break ;
      case 10 :
         result = 32.0 * x*(1.0-x) * (1.0-y)*(0.5-y) *  z*(1.0-z);
         break ;
      case 11 :
         result = -16.0 * x*(0.5-x) * (1.0-y)*(0.5-y) *  z*(1.0-z) ;
         break ;
      case 12 :
         result = -32.0 * x*(0.5-x) * y*(1.0-y) * z*(1.0-z) ;
         break ;
      case 13 :
         result = 64.0 * x*(1.0-x) * y*(1.0-y) *  z*(1.0-z) ;
         break ;
      case 14 :
         result = 32.0 * (1.0-x)*(0.5-x) * y*(1.0-y) *  z*(1.0-z) ;
         break ;
      case 15 :
         result =-16.0 * (1.0-x)*(0.5-x) * y*(0.5-y) *  z*(1.0-z) ;
         break ;
      case 16 :
         result = -32.0 * x*(1.0-x) * y*(0.5-y) *  z*(1.0-z) ;
         break ;
      case 17 :
         result = 16.0 * x*(0.5-x) * y*(0.5-y) * z*(1.0-z);
         break ;
      case 18 :
         result = -8.0 * (1.0-x)*(0.5-x) * (1.0-y)*(0.5-y) * z*(0.5-z) ;
         break ;
      case 19 :
         result = -16.0 * x*(1.0-x) * (1.0-y)*(0.5-y) *  z*(0.5-z);
         break ;
      case 20 :
         result = 8.0 * x*(0.5-x) * (1.0-y)*(0.5-y) *  z*(0.5-z) ;
         break ;
      case 21 :
         result = 16.0 * x*(0.5-x) * y*(1.0-y) * z*(0.5-z) ;
         break ;
      case 22 :
         result = -32.0 * x*(1.0-x) * y*(1.0-y) *  z*(0.5-z) ;
         break ;
      case 23 :
         result = -16.0 * (1.0-x)*(0.5-x) * y*(1.0-y) *  z*(0.5-z) ;
         break ;
      case 24 :
         result = 8.0 * (1.0-x)*(0.5-x) * y*(0.5-y) *  z*(0.5-z) ;
         break ;
      case 25 :
         result = 16.0 * x*(1.0-x) * y*(0.5-y) *  z*(0.5-z) ;
         break ;
      case 26 :
         result = -8.0 * x*(0.5-x) * y*(0.5-y) * z*(0.5-z);
         break ; 
  }
   
#endif
  return(result);  
}

void weight_NCX1_CX1 (const coord poslocal, const double Delta, double * x1, double * y1) {
  // local cell is the bottom left cell of the stencil
  *x1 = poslocal.x;
  *y1 = poslocal.y;
}

void weight_NCX1_CX2 (const coord poslocal, const double Delta, double * x1, double * y1) {
  // local cell is the bottom center cell of the stencil
  *x1 = poslocal.x - Delta;
  *y1 = poslocal.y;
}

void weight_centred (const coord poslocal, const double Delta, double * x1, double * y1) {
  // local cell is the center cell of the stencil
  *x1 = poslocal.x - Delta;
  *y1 = poslocal.y - Delta;
}

void weight_NCX1_CX4 (const coord poslocal, const double Delta, double * x1, double * y1) {
  // local cell is the center left cell of the stencil
  *x1 = poslocal.x;
  *y1 = poslocal.y - Delta;
}

void weight_NCX2_CX2 (const coord poslocal, const double Delta, double * x1, double * y1) {
  // local cell is the bottom right cell of the stencil
  *x1 = poslocal.x - 2.*Delta;
  *y1 = poslocal.y;
}

void weight_NCX2_CX3 (const coord poslocal, const double Delta, double * x1, double * y1) {
  // local cell is the center right cell of the stencil
  *x1 = poslocal.x - 2.*Delta;
  *y1 = poslocal.y - Delta;
}

void weight_NCX3_CX3 (const coord poslocal, const double Delta, double * x1, double * y1) {
  // local cell is the top right cell of the stencil
  *x1 = poslocal.x - 2.*Delta; 
  *y1 = poslocal.y - 2.*Delta; 
}

void weight_NCX3_CX4 (const coord poslocal, const double Delta, double * x1, double * y1) {
  // local cell is the top center cell of the stencil
  *x1 = poslocal.x - Delta; 
  *y1 = poslocal.y - 2.*Delta; 
}

void weight_NCX4_CX4 (const coord poslocal, const double Delta, double * x1, double * y1) {
  // local cell is at the top left cell of the stencil 
  *x1 = poslocal.x; 
  *y1 = poslocal.y - 2.*Delta;
}

double compute_weight_Quad (const int weight_id, const coord posb ,const coord poslocal, const int NCX, const int CX, const double Delta) {
  double x1 = 0., y1 = 0., weight = 0.;
  coord posref;
  
#if dimension == 2
  if (NCX == 1) {

    if (CX == 1)
      weight_NCX1_CX1 (poslocal, Delta, &x1, &y1);
    
    if (CX == 2)
      weight_NCX1_CX2 (poslocal, Delta, &x1, &y1);
      
    if (CX == 3)
      weight_centred (poslocal, Delta, &x1, &y1);
    
    if (CX == 4)
      weight_NCX1_CX4 (poslocal, Delta, &x1, &y1);
  }
  
  if (NCX == 2) {

    if (CX == 1)
      // weight_NCX2_CX1 = weight_NCX1_CX2
      weight_NCX1_CX2 (poslocal, Delta, &x1, &y1);
    
    if (CX == 2)
      weight_NCX2_CX2 (poslocal, Delta, &x1, &y1);
      

    if (CX == 3)
      weight_NCX2_CX3 (poslocal, Delta, &x1, &y1);
    

    if (CX == 4)
      weight_centred (poslocal, Delta, &x1, &y1);
  }

  if (NCX == 3) {

    if (CX == 1)
      // weight_NCX3_CX1 = weight_NCX1_CX3 = weight_centred
      weight_centred (poslocal, Delta, &x1, &y1);
    
    if (CX == 2)
      // weight_NCX3_CX2 = weight_NCX2_CX3
      weight_NCX2_CX3 (poslocal, Delta, &x1, &y1);
      
    if (CX == 3)
      weight_NCX3_CX3 (poslocal, Delta, &x1, &y1);
    
    if (CX == 4)
      weight_NCX3_CX4 (poslocal, Delta, &x1, &y1);  
  }

  if (NCX == 4) {
    if (CX == 1)
      // weight_NCX4_CX1 = weight_NCX1_CX4
      weight_NCX1_CX4 (poslocal, Delta, &x1, &y1);
      
    if (CX == 2)
      // weight_NCX4_CX2 = weight_NCX2_CX4 = weight_centred
      weight_centred (poslocal, Delta, &x1, &y1);
    
    if (CX == 3)
      // weight_NCX4_CX3 = weight_NCX3_CX4
      weight_NCX3_CX4 (poslocal, Delta, &x1, &y1);

    if (CX == 4)
      weight_NCX4_CX4 (poslocal, Delta, &x1, &y1);
  }
  
  posref.x = (posb.x - x1)/(2.*Delta);
  posref.y = (posb.y - y1)/(2.*Delta);
    
  return weight = Q2weighting(weight_id, posref.x, posref.y, 0. );
  
#elif dimension == 3
  double z1 = 0.;
  
  if (NCX == 10) {

    if (CX == 10) {
      //NCX10_CX10
      z1 = poslocal.z;
      weight_NCX1_CX1 (poslocal, Delta, &x1, &y1);}

    if (CX == 20) {
      //NCX10_CX20
      z1 = poslocal.z;
      weight_NCX1_CX2 (poslocal, Delta, &x1, &y1);}

    if (CX == 30) {
      //NCX10_CX30
      z1 = poslocal.z;
      weight_centred (poslocal, Delta, &x1, &y1);}

    if (CX == 40) {
      //NCX10_CX40
      z1 = poslocal.z;
      weight_NCX1_CX4 (poslocal, Delta, &x1, &y1);}
    
    if (CX == 11) {
      //NCX10_CX11
      z1 = poslocal.z - Delta;
      weight_NCX1_CX1 (poslocal, Delta, &x1, &y1);}

    if (CX == 21) {
      //NCX10_CX21
      z1 = poslocal.z - Delta;
      weight_NCX1_CX2 (poslocal, Delta, &x1, &y1);}

    if (CX == 31) {
      //NCX10_CX31
      z1 = poslocal.z - Delta;
      weight_centred (poslocal, Delta, &x1, &y1);}

    if(CX == 41) {
      //NCX10_CX41
      z1 = poslocal.z - Delta;
      weight_NCX1_CX4 (poslocal, Delta, &x1, &y1);}
  }
  // ok for NCX = 10 
  
  if (NCX == 11) {

    if (CX == 10) {
      //NCX11_CX10 = NCX10_CX11 
      z1 = poslocal.z - Delta;
      weight_NCX1_CX1 (poslocal, Delta, &x1, &y1);}

    if(CX == 20) {
      //NCX11_CX20
      z1 = poslocal.z - Delta;
      weight_NCX1_CX2 (poslocal, Delta, &x1, &y1);}

    if (CX == 30) {
      //NCX11_CX30
      z1 = poslocal.z - Delta;
      weight_centred (poslocal, Delta, &x1, &y1);}

    if (CX == 40) {
      //NCX11_CX40
      z1 = poslocal.z - Delta;
      weight_NCX1_CX4 (poslocal, Delta, &x1, &y1);}

    if (CX == 11) {
      //NCX11_CX11
      z1 = poslocal.z - 2*Delta; 
      weight_NCX1_CX1 (poslocal, Delta, &x1, &y1);}
    
    if (CX == 21) {
      //NCX11_CX21
      z1 = poslocal.z - 2*Delta;
      weight_NCX1_CX2 (poslocal, Delta, &x1, &y1);}
    
    if (CX == 31) {
      //NCX11_CX31
      z1 = poslocal.z - 2*Delta;
      weight_centred (poslocal, Delta, &x1, &y1);}
    
    if (CX == 41) {
      //NCX11_CX41
      z1 = poslocal.z - 2*Delta;
      weight_NCX1_CX4 (poslocal, Delta, &x1, &y1);}
  }
  // ok z1 for NCX = 11

  if (NCX == 20) {

    if (CX == 10) {
      //NCX20_CX10 = NCX10_CX20
      z1 = poslocal.z;
      weight_NCX1_CX2 (poslocal, Delta, &x1, &y1);}

    if (CX == 20) {
      //NCX20_CX20
      z1 = poslocal.z;
      weight_NCX2_CX2 (poslocal, Delta, &x1, &y1);}

    if (CX == 30) {
      //NCX20_CX30
      z1 = poslocal.z;
      weight_NCX2_CX3 (poslocal, Delta, &x1, &y1);}

    if (CX == 40) {
      //NCX20_CX40
      z1 = poslocal.z;
      weight_centred (poslocal, Delta, &x1, &y1);}

    if (CX == 11) {
      //NCX20_CX11 = NCX11_CX20
      z1 = poslocal.z - Delta;
      weight_NCX1_CX2 (poslocal, Delta, &x1, &y1);}

    if (CX == 21) {
      //NCX20_CX21
      z1 = poslocal.z - Delta;
      weight_NCX2_CX2 (poslocal, Delta, &x1, &y1);}

    if (CX == 31) {
      //NCX20_CX31
      z1 = poslocal.z - Delta;
      weight_NCX2_CX3 (poslocal, Delta, &x1, &y1);}
    
    if (CX == 41) {
      //NCX20_CX41
      z1 = poslocal.z - Delta;
      weight_centred (poslocal, Delta, &x1, &y1);}
  }
  // ok for NCX = 20
  
  if (NCX == 21) {

    if (CX == 10) {
      //NCX21_CX10 = NCX10_CX21
      z1 = poslocal.z - Delta;
      weight_NCX1_CX2 (poslocal, Delta, &x1, &y1);}

    if (CX == 20) {
      //NCX21_CX20 = NCX20_CX21
      z1 = poslocal.z - Delta;
      weight_NCX2_CX2 (poslocal, Delta, &x1, &y1);}

    if (CX == 30) {
      //NCX21_CX30
      z1 = poslocal.z - Delta;
      weight_NCX2_CX3 (poslocal, Delta, &x1, &y1);}

    if (CX == 40) {
      //NCX21_CX40
      z1 = poslocal.z - Delta;
      weight_centred (poslocal, Delta, &x1, &y1);}
    
    if (CX == 11) {
      //NCX21_CX11
      z1 = poslocal.z - 2*Delta;
      weight_NCX1_CX2 (poslocal, Delta, &x1, &y1);}
    
    if (CX == 21) {
      //NCX21_CX21
      z1 = poslocal.z - 2*Delta;
      weight_NCX2_CX2 (poslocal, Delta, &x1, &y1);}
    
    if (CX == 31) {
      //NCX21_CX31
      z1 = poslocal.z - 2*Delta;
      weight_NCX2_CX3 (poslocal, Delta, &x1, &y1);}
    
    if (CX == 41) {
      //NCX21_CX41
      z1 = poslocal.z - 2*Delta;
      weight_centred (poslocal, Delta, &x1, &y1);}
  }
  // ok z1 for NCX = 21
  
  if (NCX == 30) {

    if (CX == 10) {
      //NCX30_CX10
      z1 = poslocal.z;
      weight_centred (poslocal, Delta, &x1, &y1);}

    if (CX == 20) {
      //NCX30_CX20
      z1 = poslocal.z;
      weight_NCX2_CX3 (poslocal, Delta, &x1, &y1);}

    if (CX == 30) {
      //NCX30_CX30
      z1 = poslocal.z;
      weight_NCX3_CX3 (poslocal, Delta, &x1, &y1);}

    if (CX == 40) {
      //NCX30_CX40
      z1 = poslocal.z;
      weight_NCX3_CX4 (poslocal, Delta, &x1, &y1);}
    
    if (CX == 11) {
      //NCX30_CX11
      z1 = poslocal.z - Delta;
      weight_centred (poslocal, Delta, &x1, &y1);}
    
    if (CX == 21) {
      //NCX30_CX21
      z1 = poslocal.z - Delta;
      weight_NCX2_CX3 (poslocal, Delta, &x1, &y1);}
    
    if (CX == 31) {
      //NCX30_CX31
      z1 = poslocal.z - Delta;
      weight_NCX3_CX3 (poslocal, Delta, &x1, &y1);}
    
    if (CX == 41) {
      //NCX30_CX41
      z1 = poslocal.z - Delta;
      weight_NCX3_CX4 (poslocal, Delta, &x1, &y1);}
  }
  // ok for NCX = 30

  if (NCX == 31) {

    if (CX == 10) {
      //NCX31_CX10
      z1 = poslocal.z - Delta;
      weight_centred (poslocal, Delta, &x1, &y1);}
    
    if (CX == 20) {
      //NCX31_CX20
      z1 = poslocal.z - Delta;
      weight_NCX2_CX3 (poslocal, Delta, &x1, &y1);}
    
    if (CX == 30) {
      //NCX31_CX30
      z1 = poslocal.z - Delta;
      weight_NCX3_CX3 (poslocal, Delta, &x1, &y1);}
    
    if (CX == 40) {
      //NCX31_CX40
      z1 = poslocal.z - Delta;
      weight_NCX3_CX4 (poslocal, Delta, &x1, &y1);}
    

    if (CX == 11) {
      //NCX31_CX11
      z1 = poslocal.z - 2*Delta;
      weight_centred (poslocal, Delta, &x1, &y1);}
    
    if (CX == 21) {
      //NCX31_CX21
      z1 = poslocal.z - 2*Delta;
      weight_NCX2_CX3 (poslocal, Delta, &x1, &y1);}
    
    if (CX == 31) {
      //NCX31_CX31
      z1 = poslocal.z - 2*Delta;
      weight_NCX3_CX3 (poslocal, Delta, &x1, &y1);}
    
    if (CX == 41) {
      //NCX31_CX41
      z1 = poslocal.z - 2*Delta;
      weight_NCX3_CX4 (poslocal, Delta, &x1, &y1);}
    
  }
  // ok z1 for NCX = 31
  
  if (NCX == 40) {

    if (CX == 10) {
      //NCX40_CX10
      z1 = poslocal.z;
      weight_NCX1_CX4 (poslocal, Delta, &x1, &y1);}
    
    if (CX == 20) {
      //NCX40_CX20
      z1 = poslocal.z;
      weight_centred (poslocal, Delta, &x1, &y1);}
    
    if (CX == 30) {
      //NCX40_CX30
      z1 = poslocal.z;
      weight_NCX3_CX4 (poslocal, Delta, &x1, &y1);}
    
    if (CX == 40) {
      //NCX40_CX40
      z1 = poslocal.z;
      weight_NCX4_CX4 (poslocal, Delta, &x1, &y1);}
    
    if (CX == 11) {
      //NCX40_CX11
      z1 = poslocal.z - Delta;
      weight_NCX1_CX4 (poslocal, Delta, &x1, &y1);}
    
    if (CX == 21) {
      //NCX40_CX21
      z1 = poslocal.z - Delta;
      weight_centred (poslocal, Delta, &x1, &y1);}
    
    if (CX == 31) {
      //NCX40_CX31
      z1 = poslocal.z - Delta;
      weight_NCX3_CX4 (poslocal, Delta, &x1, &y1);}
    
    if (CX == 41) {
      //NCX40_CX41
      z1 = poslocal.z - Delta;
      weight_NCX4_CX4 (poslocal, Delta, &x1, &y1);}
  }
  // ok z1 for NCX = 40

  if (NCX == 41) {

    if (CX == 10) {
      //NCX41_CX10
      z1 = poslocal.z - Delta;
      weight_NCX1_CX4 (poslocal, Delta, &x1, &y1);}
    
    if (CX == 20) {
      //NCX41_CX20
      z1 = poslocal.z - Delta;
     weight_centred (poslocal, Delta, &x1, &y1);}
    
    if (CX == 30) {
      //NCX41_CX30
      z1 = poslocal.z - Delta;
      weight_NCX3_CX4 (poslocal, Delta, &x1, &y1);}

    if (CX == 40) {
      //NCX41_CX40
      z1 = poslocal.z - Delta;
      weight_NCX4_CX4 (poslocal, Delta, &x1, &y1);}
    
    if (CX == 11) {
      //NCX41_CX11
      z1 = poslocal.z - 2*Delta;
      weight_NCX1_CX4 (poslocal, Delta, &x1, &y1);}

    if (CX == 21) {
      //NCX41_CX21
      z1 = poslocal.z - 2*Delta;
      weight_centred (poslocal, Delta, &x1, &y1);}

    if (CX == 31) {
      //NCX41_CX31
      z1 = poslocal.z - 2*Delta;
      weight_NCX3_CX4 (poslocal, Delta, &x1, &y1);}

    if (CX == 41) {
      //NCX41_CX41
      z1 = poslocal.z - 2*Delta;
      weight_NCX4_CX4 (poslocal, Delta, &x1, &y1);}
  }
  // ok z1 for NCX = 41
  
  posref.x = (posb.x - x1)/(2*Delta); 
  posref.y = (posb.y - y1)/(2*Delta);
  posref.z = (posb.z - z1)/(2*Delta);
  
  return weight = Q2weighting(weight_id, posref.x, posref.y, posref.z); 
#endif
}

size_t is_it_in_sphere (const double x, const double y, const double z, const GeomParameter gp) {
  size_t isin = 0;

#if dimension == 2
  if (sqrt (sq (x - gp.center.x) + sq (y - gp.center.y)) < gp.radius)
    isin = 1;
#elif dimension == 3
  if (sqrt (sq (x - gp.center.x) + sq (y - gp.center.y) + sq (z - gp.center.z)) < gp.radius)
    isin = 1;
#endif
 return isin;
}

size_t is_it_in_cube (const double x, const double y, const double z, const GeomParameter gp) {
  size_t isin = 0;

  double xmax = gp.center.x + 0.5*gp.radius;
  double xmin = gp.center.x - 0.5*gp.radius;
  double ymax = gp.center.y + 0.5*gp.radius;
  double ymin = gp.center.y - 0.5*gp.radius;

#if dimension > 2  
  double zmax = gp.center.z + 0.5*gp.radius;
  double zmin = gp.center.z - 0.5*gp.radius;
#endif

  if ((x < xmax) && (x > xmin) && (y < ymax) && (y > ymin)
#if dimension > 2
      && (z< zmax) && (z > zmin)
#endif
      ) {
    isin = 1;
  }
      
  return isin;
}

void assign_dial_fd_boundary (particle * p, const coord posb, const GeomParameter gp, const double Delta,  int * NCX) {
// find the normal to the fictitious domain and assign the dial accordingly
  
  size_t  isin_xp = 0, isin_yp = 0;
  double RDelta = Delta / 100.;
  
  if (!(p->iscube) && !(p->iswall)) {
    isin_xp = is_it_in_sphere (posb.x + RDelta, posb.y, posb.z, gp);
    isin_yp = is_it_in_sphere (posb.x, posb.y + RDelta, posb.z, gp);
  }

  if (p->iscube) {
    coord checkpt = {posb.x + RDelta, posb.y, posb.z};
    isin_xp = is_it_in_cube_v2 (&(p->g.u1), &(p->g.v1), &(p->g.w1), &(p->g.mins), &(p->g.maxs), &checkpt);

    checkpt.x = posb.x;
    checkpt.y = posb.y + RDelta;
    checkpt.z = posb.z;
    isin_yp = is_it_in_cube_v2 (&(p->g.u1), &(p->g.v1), &(p->g.w1), &(p->g.mins), &(p->g.maxs), &checkpt);
  }
  
#if dimension == 2
  if (isin_xp) {
    if(isin_yp){
      // fictitious-boundary's normal is oriented x- and y-. 
      *NCX = 3;
    }
    else{
	// fictitious-boundary's normal is oriented x- and y+.
	*NCX = 2;
      }
  }
  else {
    if (isin_yp) {
      // fictitious-boundary's normal is oriented x+ and y-. 
      *NCX = 4;}
    else {
      // fictitious-boundary's normal is oriented x+ and y+.
      *NCX = 1;
    }
  }
      
#elif dimension == 3
  size_t isin_zp = 0;

  if (!(p->iscube) && !(p->iswall)) {
    isin_zp = is_it_in_sphere (posb.x, posb.y, posb.z + RDelta, gp);
  }
  coord checkpt;
  checkpt.x = posb.x;
  checkpt.y = posb.y;
  checkpt.z  = posb.z + RDelta;
  
  if (p->iscube)
    isin_zp = is_it_in_cube_v2 (&(p->g.u1), &(p->g.v1), &(p->g.w1), &(p->g.mins), &(p->g.maxs), &checkpt);
  
  if (isin_zp) {
    // fictitious-boundary's normal is oriented z-.  
    if (isin_xp) {
      if (isin_yp) {
	// fictitious-boundary's normal is oriented x-, y-, z- . 
	*NCX = 31;
      }
      else{
	// fictitious-boundary's normal is oriented x-, y+, z-.
	*NCX = 21;
      }
    }
    else{
      if(isin_yp){
	// fictitious-boundary's normal is oriented x+, y-, z-. 
	*NCX = 41;}
      else{
	// fictitious-boundary's normal is oriented x+, y+, z-.
	*NCX = 11;
      }
    }
  }
  else {
    // fictitious-boundary's normal is oriented z+.  
    if (isin_xp) {
      if (isin_yp) {
	// fictitious-boundary's normal is oriented x-, y-, z+ . 
	*NCX = 30;
      }
      else {
	// fictitious-boundary's normal is oriented x-, y+, z+.
	*NCX = 20;
      }
    }
    else {
      if (isin_yp){ 
	// fictitious-boundary's normal is oriented x+, y-, z+. 
	*NCX = 40;}
      else {
	// fictitious-boundary's normal is oriented x+, y+, z+.
	*NCX = 10;
      }
    }
  }
  if (*NCX == 0) {
    fprintf(stderr, "NCX = 0, isin_xp = %zu, isin_yp = %zu, isin_zp = %zu \n",  isin_xp, isin_yp, isin_zp);
  }
#endif
}
