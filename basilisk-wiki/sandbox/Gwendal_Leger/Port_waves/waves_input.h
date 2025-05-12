// Number of waves
#define waves_N 10

// Normalize amplitude ? (0 for no; else for yes)
#define waves_normalize 0

// Amplification of waves
#define waves_amplification 7.

// Dilatation of wavelength
#define waves_dilatation 1.

// Amplitude
double waves_A[waves_N] = {
			   0.5,
			   1.,
			   0.5,
			   0.25,
			   1.,
			   0.5,
			   0.75,
			   0.5,
			   0.4,
			   0.33			   
};

// Wavelength
double waves_lambda[waves_N] = {
				1.,
				2.,
				10.,
				11.,
				7.,
				3.,
				4.,
				25.,
				13.,
				12.
};

// Angle of diffusion (0 means toward increasing x)
double waves_theta[waves_N] = {
			       0.,
			       0.,
			       0.,
			       0.,
			       0.,
			       0.,
			       0.,
			       0.,
			       0.,
			       0.
};
