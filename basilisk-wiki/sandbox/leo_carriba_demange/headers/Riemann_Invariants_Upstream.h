// 16/04/2021
// Stage Injection-Basilisk
// Implementation des invariants de Riemann amont

// resout le trinome : ax^2 + bx + c = 0 et affiche une erreur s'il n'y a pas de racine reelle

void Quadratic_Equation(double a, double b, double c, double * x_plus, double * x_minus)
{
    double Delta = b * b - 4 * a * c;
    if (Delta < 0.)
    {
        printf("Error : this quadratic equation has no real roots ! We have Delta = %g\n", Delta);
    }
    else if (Delta == 0.)
    {
        printf("Beware : both roots have the same value ! We have Delta = %g\n", Delta);
        *x_minus = -b / (2 * a);
        *x_plus = -b / (2 * a);
    }
    else
    {
        double sqrt_Delta = sqrt(Delta);
        *x_minus = (-b - sqrt_Delta) / (2 * a);
        *x_plus = (-b + sqrt_Delta) / (2 * a);
    }
    return;
}


// Calcule la vitesse et la hauteur d'eau à injecter en fonction du debit et des conditions a l'interieur du domaine

void Solve_Riemann_Invariants_Upstream(double q0, double g, double delta_t,
                                  double u_int, double h_int, double Sf_int, double S0_int,
                                   double * u_ghost, double * h_ghost)
//the subscript "_int" indicates the interior of the domain, whereas "_ghost" refers to the upstream ghost cell
{
    double c_int = sqrt(g * h_int);
    double AA = g / c_int ;
    double BB = u_int - c_int + g * (S0_int - Sf_int) * delta_t;
    double h_plus, h_minus;
    double CC = -q0;
    Quadratic_Equation(AA, BB, CC, &h_plus, &h_minus);

    *h_ghost = h_plus;
    *u_ghost = q0 / *h_ghost;

    return;

}



