#include <math.h>
#include <assert.h>

typedef double real;

typedef struct { real x, y; } vec2;

static inline vec2 vec2_sub (vec2 v1, vec2 v2)
{
  return (vec2){v1.x - v2.x, v1.y - v2.y};
}

static inline real vec2_cross (vec2 v1, vec2 v2)
{
  return v1.x*v2.y - v1.y*v2.x;
}

static inline real vec2_norm2 (vec2 v)
{
  return v.x*v.x + v.y*v.y;
}

static inline real vec2_norm  (vec2 v)
{
  return sqrt (vec2_norm2 (v));
}

static inline vec2 vec2_normalized (vec2 v)
{
  real norm = vec2_norm (v);
  if (norm)
    return (vec2){v.x/norm, v.y/norm};
  return v;
}

typedef struct { real x, y, z; } vec3;

static inline real vec3_norm2 (vec3 v)
{
  return v.x*v.x + v.y*v.y + v.z*v.z;
}

static inline real vec3_norm  (vec3 v)
{
  return sqrt (vec3_norm2 (v));
}

static inline vec3 vec3_normalized (vec3 v)
{
  real norm = vec3_norm (v);
  if (norm)
    return (vec3){v.x/norm, v.y/norm, v.z/norm};
  return v;
}

static inline vec3 vec3_embed (vec2 v)
{
  return (vec3) { v.x, v.y, 1. };
}

static inline vec3 vec3_sub (vec3 v1, vec3 v2)
{
  return (vec3){v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
}

static inline vec3 vec3_div (vec3 v, real a)
{
  assert (a);
  return (vec3){v.x/a, v.y/a, v.z/a};
}

static inline vec3 vec3_mul (vec3 v, real a)
{
  return (vec3){v.x*a, v.y*a, v.z*a};
}

static inline real vec3_scalar (vec3 v1, vec3 v2)
{
  return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

typedef struct { real x, y, z, t; } vec4;

static inline vec4 vec4_div (vec4 v, real a)
{
  assert (a);
  return (vec4){v.x/a, v.y/a, v.z/a, v.t/a};  
}

static inline vec2 vec4_proj2 (vec4 v)
{
  return (vec2) { v.x, v.y };
}

static inline vec3 vec4_proj3 (vec4 v)
{
  return (vec3) { v.x, v.y, v.z };
}

static inline vec4 vec4_embed (vec3 v, real fill)
{
  return (vec4) { v.x, v.y, v.z, fill };
}

typedef struct { vec3 x, y, z; } mat3;

static inline mat3 mat3_transpose (mat3 m)
{
  return (mat3) {
    { m.x.x, m.y.x, m.z.x },
    { m.x.y, m.y.y, m.z.y },
    { m.x.z, m.y.z, m.z.z }
  };
}

static inline real mat3_det (mat3 m)
{
  return (m.x.x*(m.y.y*m.z.z - m.z.y*m.y.z) - 
	  m.x.y*(m.y.x*m.z.z - m.z.x*m.y.z) + 
	  m.x.z*(m.y.x*m.z.y - m.z.x*m.y.y));
}

static inline mat3 mat3_invert (mat3 m)
{
  real det = mat3_det (m);
  assert (det);

  mat3 mi;
  mi.x.x = (m.y.y*m.z.z - m.y.z*m.z.y)/det; 
  mi.x.y = (m.z.y*m.x.z - m.x.y*m.z.z)/det;
  mi.x.z = (m.x.y*m.y.z - m.y.y*m.x.z)/det; 
  mi.y.x = (m.y.z*m.z.x - m.y.x*m.z.z)/det; 
  mi.y.y = (m.x.x*m.z.z - m.z.x*m.x.z)/det; 
  mi.y.z = (m.y.x*m.x.z - m.x.x*m.y.z)/det; 
  mi.z.x = (m.y.x*m.z.y - m.z.x*m.y.y)/det; 
  mi.z.y = (m.z.x*m.x.y - m.x.x*m.z.y)/det; 
  mi.z.z = (m.x.x*m.y.y - m.x.y*m.y.x)/det; 
  return mi;
}

static inline mat3 mat3_invert_transpose (const mat3 m)
{
  return mat3_transpose (mat3_invert (m));
}

static inline vec3 mat3_mul (mat3 m, vec3 v)
{
  return (vec3) {
    m.x.x*v.x + m.x.y*v.y + m.x.z*v.z,
    m.y.x*v.x + m.y.y*v.y + m.y.z*v.z,
    m.z.x*v.x + m.z.y*v.y + m.z.z*v.z
  };
}

static inline vec3 mat3_col (mat3 m, int idx)
{
  assert (idx >= 0 && idx < 3);
  return (vec3) { *(&m.x.x + idx), *(&m.y.x + idx), *(&m.z.x + idx) };
}

static inline int set_col3 (mat3 * m, int idx, vec3 v)
{
  assert (idx >= 0 && idx < 3);
  *(&m->x.x + idx) = v.x;
  *(&m->y.x + idx) = v.y;
  *(&m->z.x + idx) = v.z;
  return 0; // fixme
}

typedef struct { vec4 x, y, z, t; } mat4;

static inline vec4 mat4_mul (mat4 m, vec4 v)
{
  return (vec4) {
    m.x.x*v.x + m.x.y*v.y + m.x.z*v.z + m.x.t*v.t,
    m.y.x*v.x + m.y.y*v.y + m.y.z*v.z + m.y.t*v.t,
    m.z.x*v.x + m.z.y*v.y + m.z.z*v.z + m.z.t*v.t,
    m.t.x*v.x + m.t.y*v.y + m.t.z*v.z + m.t.t*v.t
  };
}

static inline mat4 mat4_mul4 (mat4 m1, mat4 m2)
{
  mat4 ret;
  const real * p1 = &m1.x.x, * p2 = &m2.x.x;
  real * p = &ret.x.x;
  for (int i = 0; i < 4; i++, p1 += 4, p += 4)
    for (int j = 0; j < 4; j++) {
      real a = 0.;
      for (int k = 0; k < 4; k++)
	a += p1[k]*(p2 + 4*k)[j];
      p[j] = a;
    }
  return ret;
}

static inline mat4 mat4_invert (mat4 M)
{
  mat4 I;
  real * m = &M.x.x, * i = &I.x.x;
  
  i[0] = m[5]*m[10]*m[15] - m[5]*m[11]*m[14] - m[9]*m[6]*m[15] + m[9]*m[7]*m[14] + m[13]*m[6]*m[11] - m[13]*m[7]*m[10];
  i[4] = -m[4]*m[10]*m[15] + m[4]*m[11]*m[14] + m[8]*m[6]*m[15] - m[8]*m[7]*m[14] - m[12]*m[6]*m[11] + m[12]*m[7]*m[10];
  i[8] = m[4]*m[9]*m[15] - m[4]*m[11]*m[13] - m[8]*m[5]*m[15] + m[8]*m[7]*m[13] + m[12]*m[5]*m[11] - m[12]*m[7]*m[9];
  i[12] = -m[4]*m[9]*m[14] + m[4]*m[10]*m[13] +m[8]*m[5]*m[14] - m[8]*m[6]*m[13] - m[12]*m[5]*m[10] + m[12]*m[6]*m[9];
  i[1] = -m[1]*m[10]*m[15] + m[1]*m[11]*m[14] + m[9]*m[2]*m[15] - m[9]*m[3]*m[14] - m[13]*m[2]*m[11] + m[13]*m[3]*m[10];
  i[5] = m[0]*m[10]*m[15] - m[0]*m[11]*m[14] - m[8]*m[2]*m[15] + m[8]*m[3]*m[14] + m[12]*m[2]*m[11] - m[12]*m[3]*m[10];
  i[9] = -m[0]*m[9]*m[15] + m[0]*m[11]*m[13] + m[8]*m[1]*m[15] - m[8]*m[3]*m[13] - m[12]*m[1]*m[11] + m[12]*m[3]*m[9];
  i[13] = m[0]*m[9]*m[14] - m[0]*m[10]*m[13] - m[8]*m[1]*m[14] + m[8]*m[2]*m[13] + m[12]*m[1]*m[10] - m[12]*m[2]*m[9];
  i[2] = m[1]*m[6]*m[15] - m[1]*m[7]*m[14] - m[5]*m[2]*m[15] + m[5]*m[3]*m[14] + m[13]*m[2]*m[7] - m[13]*m[3]*m[6];
  i[6] = -m[0]*m[6]*m[15] + m[0]*m[7]*m[14] + m[4]*m[2]*m[15] - m[4]*m[3]*m[14] - m[12]*m[2]*m[7] + m[12]*m[3]*m[6];
  i[10] = m[0]*m[5]*m[15] - m[0]*m[7]*m[13] - m[4]*m[1]*m[15] + m[4]*m[3]*m[13] + m[12]*m[1]*m[7] - m[12]*m[3]*m[5];
  i[14] = -m[0]*m[5]*m[14] + m[0]*m[6]*m[13] + m[4]*m[1]*m[14] - m[4]*m[2]*m[13] - m[12]*m[1]*m[6] + m[12]*m[2]*m[5];
  i[3] = -m[1]*m[6]*m[11] + m[1]*m[7]*m[10] + m[5]*m[2]*m[11] - m[5]*m[3]*m[10] - m[9]*m[2]*m[7] + m[9]*m[3]*m[6];
  i[7] = m[0]*m[6]*m[11] - m[0]*m[7]*m[10] - m[4]*m[2]*m[11] + m[4]*m[3]*m[10] + m[8]*m[2]*m[7] - m[8]*m[3]*m[6];
  i[11] = -m[0]*m[5]*m[11] + m[0]*m[7]*m[9] + m[4]*m[1]*m[11] - m[4]*m[3]*m[9] - m[8]*m[1]*m[7] + m[8]*m[3]*m[5];
  i[15] = m[0]*m[5]*m[10] - m[0]*m[6]*m[9] - m[4]*m[1]*m[10] + m[4]*m[2]*m[9] + m[8]*m[1]*m[6] - m[8]*m[2]*m[5];

  real det = m[0]*i[0] + m[1]*i[4] + m[2]*i[8] + m[3]*i[12];
  if (det == 0)
    return M;

  for (int j = 0; j < 16; j++)
    i[j] /= det;

  return I;
}

static inline mat4 mat4_transpose (mat4 m)
{
  return (mat4) {
    { m.x.x, m.y.x, m.z.x, m.t.x },
    { m.x.y, m.y.y, m.z.y, m.t.y },
    { m.x.z, m.y.z, m.z.z, m.t.z },
    { m.x.t, m.y.t, m.z.t, m.t.t }
  };
}

static inline mat4 mat4_invert_transpose (const mat4 m)
{
  return mat4_transpose (mat4_invert (m));
}

typedef struct { vec3 x, y; } mat23;

static inline int set_col23 (mat23 * m, int idx, vec2 v)
{
  assert (idx >= 0 && idx < 3);
  *(&m->x.x + idx) = v.x;
  *(&m->y.x + idx) = v.y;
  return 0; // fixme  
}

static inline vec2 mat23_mul (mat23 m, vec3 v)
{
  return (vec2) {
    m.x.x*v.x + m.x.y*v.y + m.x.z*v.z,
    m.y.x*v.x + m.y.y*v.y + m.y.z*v.z
  };
}

static inline vec3 cross(vec3 v1, vec3 v2)
{
  return (vec3){
    v1.y*v2.z - v1.z*v2.y,
    v1.z*v2.x - v1.x*v2.z,
    v1.x*v2.y - v1.y*v2.x
  };
}
