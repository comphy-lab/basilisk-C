#ifndef PID_H
#define PID_H

/**
This header enables the enforcement of PID dampening on a part of
the domain in the multilayer solver. 
Currently this implementation only supports the PID dampening of a single
scalar field with one set of dampening coefficients. 

To use this header, call on the init_corrector function in the initialising event.
The function takes three arguments, all double type, namely the dampening 
coefficients Kp (Proportional), Ki (Integral) and Kd (Derivative).
A general rule of advice is to let all these coefficients be in the interval
[0,1).

This first function call will ensure that all variables and 
fields are initialized correctly. Then, to enforce the PID correction,
simply call on the pid_scalar function at all desired timepoints. This function
takes four arguments, namely

scalar SP  - Scalar field specifying the desired value for the computed field

scalar PV  - The computed field

coord x0   - coordinates of the lower left corner of the desired dampening area 

coord x1   - coordinates of the upper right corner of the desired dampening area 

*/
struct PIDcorrector{
    double Kp;
    double Ki;
    double Kd;
};


struct PIDcorrector corrector;
scalar  _Pf; 
scalar  _If;

void init_corrector(double Kp, double Ki, double Kd){
    corrector.Kp = Kp;
    corrector.Ki = Ki;
    corrector.Kd = Kd;
    _Pf = new scalar[nl];
    _If = new scalar[nl];
    foreach(){
        foreach_layer(){
            _Pf[] = 0;
            _If[] = 0;
        }
    }
}

void pid_scalar(scalar SP, scalar PV, coord x0={X0,Y0,0}, coord x1={L0+X0,L0+Y0,0}){
    double P,I,D, cval;
    foreach(){
        foreach_layer(){
            int apply_PID = 0;  
            if (x > x0.x && x < x1.x){
                apply_PID +=1;
            }
            #if dimension > 1
                if (y > x0.y && y < x1.y){
                    apply_PID += 1;
                }
            #endif

            if (apply_PID == dimension){
                P = (SP[] - PV[]);
                I = 0.5*(P + _Pf[])*dt; //Trapezoidal rule
                D = (P-_Pf[])/dt;
                cval = corrector.Kp*P + corrector.Ki*(_If[] + I) + corrector.Kd*D;
                PV[] += cval;
                _Pf[] = P;
                _If[] += I;
            }
        }
    }
}

event cleanup_PID(t=end){
    delete((scalar *){_Pf, _If});
}
#endif