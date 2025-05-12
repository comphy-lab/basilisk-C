#include "PointTriangle.h"
void draw_space_curve(int n_seg, coord *p){
  // // Set the clear color to a shade of blue (RGBA)
  // glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

  // // Clear the color buffer
  // glClear(GL_COLOR_BUFFER_BIT);

  // Set the line width to a thicker value
  glLineWidth(3.0f);

  // Change the color
  glColor3f(1.0f, 0.0f, 0.0f);
  for (int i = 0; i < n_seg-1; i++){
    glBegin(GL_LINES);
      glVertex3f(p[i  ].x, p[i  ].y, p[i  ].z);
      glVertex3f(p[i+1].x, p[i+1].y, p[i+1].z);
    glEnd();
  }
}

void draw_space_curve_with_vectors(int n_seg, coord *p, coord *t, coord *n, coord *b, float scale=1.0) {
  // // Set the clear color to a shade of blue (RGBA)
  // glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

  // // Clear the color buffer
  // glClear(GL_COLOR_BUFFER_BIT);

  glLineWidth(5.0f);
  for (int i = 0; i < n_seg-1; i++) {
    // Start point of the vector is the curve point p[i]
    // End point of the vector is p[i] + t[i]
    glColor3f(1.0f, 0.0f, 0.0f);
    glBegin(GL_LINES);
      glVertex3f(p[i].x,                p[i].y,                p[i].z);
      glVertex3f(p[i].x + scale*t[i].x, p[i].y + scale*t[i].y, p[i].z + scale*t[i].z);
    glEnd();
    
    // End point of the vector is p[i] + n[i]
    glColor3f(0.0f, 1.0f, 0.0f);
    glBegin(GL_LINES);
      glVertex3f(p[i].x,                p[i].y,                p[i].z);
      glVertex3f(p[i].x + scale*n[i].x, p[i].y + scale*n[i].y, p[i].z + scale*n[i].z);
    glEnd();

    // End point of the vector is p[i] + b[i]
    glColor3f(0.0f, 0.0f, 1.0f);
    glBegin(GL_LINES);
      glVertex3f(p[i].x,                p[i].y,                p[i].z);
      glVertex3f(p[i].x + scale*b[i].x, p[i].y + scale*b[i].y, p[i].z + scale*b[i].z);
    glEnd();
  }
}

void draw_tube_along_curve(int n, coord *p, double *a, int segments=16) {
  // // Set the clear color to a shade of blue (RGBA)
  // glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

  // // Clear the color buffer
  // glClear(GL_COLOR_BUFFER_BIT);

  if (n < 2) return; // Need at least two points to form a segment

  for (int i = 0; i < n - 1; ++i) {
    // Define the two endpoints of the current segment
    coord p1 = p[i];
    float a1 = a[i];
    coord p2 = p[i + 1];    
    float a2 = a[i + 1];

    // Calculate the tangent vector of the curve segment
    coord tangent = vecdiff(p2, p1);
    
    // Normalize the tangent vector
    float tangentLength = sqrtf(vecdot(tangent,tangent));
    if (tangentLength < 1e-6f) continue; // Avoid division by zero for degenerate segments
    foreach_dimension()
      tangent.x /= tangentLength;

    // Find a vector orthogonal to the tangent    
    coord up = {0.0f, 1.0f, 0.0f}; // A default up vector
    // If the tangent is close to the up vector, use a different reference
    if (fabsf(tangent.y) > 0.99f) {
      up = (coord){1.0f, 0.0f, 0.0f};
    }

    coord normal;
    foreach_dimension()
      normal.x = tangent.y * up.z - tangent.z * up.y;

    // Normalize the normal vector
    float normalLength = sqrtf(vecdot(normal,normal));
    foreach_dimension()
      normal.x /= normalLength;
    
    // Calculate the binormal vector (orthogonal to both tangent and normal)
    coord binormal;
    foreach_dimension()
      binormal.x = tangent.y * normal.z - tangent.z * normal.y;
    
    // Generate vertices for the tube segment using GL_QUADS
    for (int j = 0; j < segments; ++j) {
      float angle1 = 2.0f * pi * j / segments;
      float angle2 = 2.0f * pi * (j + 1) / segments;

      float cosAngle1 = cosf(angle1);
      float sinAngle1 = sinf(angle1);
      float cosAngle2 = cosf(angle2);
      float sinAngle2 = sinf(angle2);

      // Calculate the coordinates of the vertices on the circles at p1 and p2
      coord v1p1 = {p1.x + a1 * (cosAngle1 * normal.x + sinAngle1 * binormal.x),
                    p1.y + a1 * (cosAngle1 * normal.y + sinAngle1 * binormal.y),
                    p1.z + a1 * (cosAngle1 * normal.z + sinAngle1 * binormal.z)};

      coord v2p1 = {p1.x + a1 * (cosAngle2 * normal.x + sinAngle2 * binormal.x),
                    p1.y + a1 * (cosAngle2 * normal.y + sinAngle2 * binormal.y),
                    p1.z + a1 * (cosAngle2 * normal.z + sinAngle2 * binormal.z)};

      coord v1p2 = {p2.x + a2 * (cosAngle1 * normal.x + sinAngle1 * binormal.x),
                    p2.y + a2 * (cosAngle1 * normal.y + sinAngle1 * binormal.y),
                    p2.z + a2 * (cosAngle1 * normal.z + sinAngle1 * binormal.z)};

      coord v2p2 = {p2.x + a2 * (cosAngle2 * normal.x + sinAngle2 * binormal.x),
                    p2.y + a2 * (cosAngle2 * normal.y + sinAngle2 * binormal.y),
                    p2.z + a2 * (cosAngle2 * normal.z + sinAngle2 * binormal.z)};

      // Draw the white panels (GL_QUADS)
      glColor3f(1.0f, 1.0f, 1.0f);
      glBegin(GL_QUADS);
        glVertex3f(v1p1.x, v1p1.y, v1p1.z);
        glVertex3f(v2p1.x, v2p1.y, v2p1.z);
        glVertex3f(v2p2.x, v2p2.y, v2p2.z);
        glVertex3f(v1p2.x, v1p2.y, v1p2.z);
      glEnd();

      // Draw the black edges (GL_LINES)
      glColor3f(0.0f, 0.0f, 0.0f);
      glBegin(GL_LINES);
        // Around the first circle
        glVertex3f(v1p1.x, v1p1.y, v1p1.z);
        glVertex3f(v2p1.x, v2p1.y, v2p1.z);

        // Around the second circle
        glVertex3f(v1p2.x, v1p2.y, v1p2.z);
        glVertex3f(v2p2.x, v2p2.y, v2p2.z);

        // Connecting the two circles
        glVertex3f(v1p1.x, v1p1.y, v1p1.z);
        glVertex3f(v1p2.x, v1p2.y, v1p2.z);

        glVertex3f(v2p1.x, v2p1.y, v2p1.z);
        glVertex3f(v2p2.x, v2p2.y, v2p2.z);
      glEnd();
    }
  }
}