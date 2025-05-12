// Number of points to interpolate :
#define interpolation_data_N 5

// Values to interpolate :
double interpolation_data_values[interpolation_data_N] = {
					      -1.,
					      -1.,
					      0.,
					      -1.,
					      1.,
};


// We define the interpolation order and the type of interpolation.
// If the interpolation order is 0, it means nearest neighbor interpolation is used,
// 1 means linear interpolation,
// >1 means interpolation by part and
// <0 means interpolation over all points.
// The type of interpolation is 0 if Lagrange interpolation is used
// and otherwise Bézier approximation is used.
#define interpolation_order 2
#define interpolation_type 1


