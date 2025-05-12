/**
In this file, we implement a routine to find the smallest real root of a
third-order polynomial, in absolute value.
*/

typedef struct realQuadraticRoots {
    bool real;
    double roots[2];
} realQuadraticRoots;

realQuadraticRoots real_quadratic_roots(double* p) {
    realQuadraticRoots result;
    double a, b, c;
    a = p[2]; b = p[1]; c = p[0];
    double epsilon = 1.e-16;
    if (fabs(a) < epsilon) {
        result.real = true;
        result.roots[0] = (fabs(b) > epsilon) ? -c/b : HUGE;
        result.roots[1] = HUGE;
        return result;
    }

    double delta = sq(b) - 4*a*c;
    if (fabs(delta) > epsilon && delta < 0) {
        result.real = false;
        result.roots[0] = HUGE;
        result.roots[1] = HUGE;
    }
    if (fabs(delta) > epsilon && delta > 0) {
        result.real = true;
        double r1, r2;
        r1 = (-b - sqrt(delta))/(2*a);
        r2 = (-b + sqrt(delta))/(2*a);
        result.roots[0] = fabs(r1) > fabs(r2) ? r2 : r1;
        result.roots[1] = fabs(r1) > fabs(r2) ? r1 : r2;
    }
    if (fabs(delta) < epsilon) {
        result.real = true;
        result.roots[0] = -b/(2*a);
        result.roots[1] = -b/(2*a);
    }
    return result;
}

trace
double find_smallest_real_root(double* a) {
    double epsilon = 1.e-10;
        
    if (fabs(a[3]) > epsilon) {
        /** Solve 3rd-order polynomial: */
        double alpha = a[2]/a[3];
        double beta = a[1]/a[3];
        double gamma = a[0]/a[3];
        double Q = (sq(alpha) - 3*beta)/9;
        double R = (2*cube(alpha) - 9*alpha*beta + 27*gamma)/54;

        if (sq(R) < cube(Q)) {
            double theta = acos(R/sqrt(cube(Q)));
            double r1 = -2*sqrt(Q)*cos(theta/3) - alpha/3;
            double r2 = -2*sqrt(Q)*cos((theta + 2*pi)/3) - alpha/3;
            double r3 = -2*sqrt(Q)*cos((theta - 2*pi)/3) - alpha/3;
            double min_root = (fabs(r1) < fabs(r2)) ? (fabs(r1) < fabs(r3) ?
                r1 : r3) : (fabs(r2) < fabs(r3) ? r2 : r3);
            return min_root;
        }
        else {
            double A = -sign(R)*cbrt(fabs(R) + sqrt(sq(R) - cube(Q)));
            double B = (fabs(A) > epsilon) ? Q/A : 0;
            double root = (A + B) - alpha/3;
            return root;
        }
    }

    if (fabs(a[2]) > epsilon) {
        /** Solve 2nd-order polynomial */
        realQuadraticRoots my_roots = real_quadratic_roots(a);
        if (my_roots.real == false) {
            fprintf(stderr, "Warning: no real roots found in find_smallest_real_root\n");
            return 0.;
        }
        double root1, root2;
        root1 = my_roots.roots[0];
        root2 = my_roots.roots[1];
        double smallest_real_root = (fabs(root1) > fabs(root2)) ? root2 : root1;
        return smallest_real_root;
    }

    if (fabs(a[1]) > epsilon) {
        /** Solve linear function */
        return -a[0]/a[1];
    }

    fprintf(stderr, "Warning: giving constant function instead of third-order \
        polynomial in find_smallest_real_root.");
    return 0.;
}