// Number of waves
#define waves_N 4

// Normalize amplitude ? (0 for no; else for yes)
#define waves_normalize 1

// Amplification of waves
#define waves_amp 0.4

// Amplitude
double waves_A[waves_N] = {
			   0.5,
			   1.,
			   0.5,
			   0.25
};

// Wavelength
double waves_lambda[waves_N] = {
				10.,
				20.,
				100.,
				0.5
};

// Angle of diffusion (0 means toward increasing x)
double waves_theta[waves_N] = {
			       0.,
			       0.,
			       0.,
			       0.
};
