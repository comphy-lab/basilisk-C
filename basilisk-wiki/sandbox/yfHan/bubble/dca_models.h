/** This headfile includes different dynamic contact angle models. */
#ifndef DYNAMIC_CONTACT_ANGLE_H
#define DYNAMIC_CONTACT_ANGLE_H
#endif



/** Structure to store model parameters*/
/** Ratio:  ratio of macro to micro length scales L/l */
/** xi_mkt:  MKT friction coefficient */
/** Cst01: fitting constant for pinning effects */
/** mu_liquid: liquid viscosity */

typedef struct {
    double Ratio;
    double xi_mkt;
    double Cst01;
    double mu_liquid;
    // Add more parameters as needed
} ModelParams;


double dynamic_contact_angle(double Ca, double theta0, double theta_adv, double theta_rec, ModelParams params, const char* model_name) {

    //double Temp = 293.15, kB_cst = 1.38e-23; //room temperature and Boltzman constant
    double theta_d = theta0; // Default to equilibrium contact angle

  /** Equation 1: Hydrodynamic Theory(HT) model */
    if (strcmp(model_name, "HT") == 0) {
        theta_d = pow(theta0/180.*pi,3.)+9.0*Ca*log(params.Ratio);
        theta_d = pow(theta_d,1.0/3)*180./pi;
    }

/** Equation 2: Molecular Kinetic Theory (MKT) model */
    else if (strcmp(model_name, "MKT") == 0) {
        theta_d = acos(cos(theta0*pi/180.) - Ca*params.xi_mkt/params.mu_liquid)*180./pi;
    }

/** Equation 3: Combined HT and MKT model */
    else if (strcmp(model_name, "COMBINED") == 0) {
        double theta_local;     

        theta_local = acos(cos(theta0*pi/180.) - Ca*params.xi_mkt/params.mu_liquid)*180./pi; // classical MKT
        
        theta_d = pow(theta_local / 180. * pi, 3.) + 9.0 * Ca * log(params.Ratio);
        theta_d = pow(theta_d, 1.0 / 3) * 180. / pi;
    }

/** Equation 4: Combined HT and MKT model with pinning effects */
    else if (strcmp(model_name, "COMBINED_PIN") == 0) {
        
double C_pin;
double theta_local;
if (Ca>=0.0)
  C_pin = f.sigma*(cos(theta0 * pi / 180.)-cos(theta_adv * pi / 180.));
else
  C_pin = f.sigma*(cos(theta_rec * pi / 180.)-cos(theta0 * pi / 180.));
        theta_local = acos(cos(theta0*pi/180.) - params.xi_mkt*Ca/params.mu_liquid - C_pin*tanh(params.Cst01*Ca)/f.sigma )*180./pi;
        theta_d = pow(theta_local / 180. * pi, 3) + 9.0 * Ca * log(params.Ratio);
        theta_d = pow(theta_d, 1.0 / 3) * 180. / pi;
    }

    return theta_d;
}