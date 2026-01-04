/**
We define a simple macro that returns the levelset function of a superquadric.
This family of shapes are defined by the following equation:
$$
\left(\frac{x}{a}\right)^m + \left(\frac{y}{b}\right)^m = 1
$$
where 'm' is the shape parameter, 'a' is the radius in the x direction, and 'b' is the radius in the y direction.
The levelset function is defined as the distance from the surface of the superquadric.
*/

macro double superquadric (double x, double y, 
                            double m, double a, double b,
                            double xc = 0., double yc = 0.)
{
    return 1. - pow((x-xc*L0)/a, m) - pow((y-yc*L0)/b, m);
}