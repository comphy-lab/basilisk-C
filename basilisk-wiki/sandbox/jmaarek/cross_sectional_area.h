#define MAX_POINTS 50
#define TOL 1e-8

#if dimension == 2


double distance(coord p1, coord p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return sqrt(dx * dx + dy * dy);
}

int is_duplicate(coord p, coord *points, int count) {
    for (int i = 0; i < count; i++) {
        if (fabs(p.x - points[i].x) < TOL && fabs(p.y - points[i].y) < TOL)
            return 1;
    }
    return 0;
}

/* 
 * Computes the intersection line segment length between the line
 * a*x + b*y + c = 0 and the rectangle defined by [xmin, xmax] and [ymin, ymax].
 */
double compute_cross_section_area(coord n, double alpha, coord aa, coord bb) {
    
    double a = n.x;
    double b = n.y;
    double c = alpha;
    double xmin = aa.x;
    double xmax = bb.x;
    double ymin = aa.y;
    double ymax = bb.y;



    coord intersections[4];
    int count = 0;

    // Bottom edge: y = ymin
    if (fabs(a) > TOL) {
        double x = -(b * ymin + c) / a;
        if (x >= xmin - TOL && x <= xmax + TOL) {
            coord p = {x, ymin};
            if (!is_duplicate(p, intersections, count))
                intersections[count++] = p;
        }
    } else {
        // a is 0 -> horizontal line: check if it lies on y = ymin.
        if (fabs(b * ymin + c) < TOL) {
            coord p1 = {xmin, ymin};
            coord p2 = {xmax, ymin};
            if (!is_duplicate(p1, intersections, count))
                intersections[count++] = p1;
            if (!is_duplicate(p2, intersections, count))
                intersections[count++] = p2;
        }
    }

    // Top edge: y = ymax
    if (fabs(a) > TOL) {
        double x = -(b * ymax + c) / a;
        if (x >= xmin - TOL && x <= xmax + TOL) {
            coord p = {x, ymax};
            if (!is_duplicate(p, intersections, count))
                intersections[count++] = p;
        }
    } else {
        if (fabs(b * ymax + c) < TOL) {
            coord p1 = {xmin, ymax};
            coord p2 = {xmax, ymax};
            if (!is_duplicate(p1, intersections, count))
                intersections[count++] = p1;
            if (!is_duplicate(p2, intersections, count))
                intersections[count++] = p2;
        }
    }

    // Left edge: x = xmin
    if (fabs(b) > TOL) {
        double y = -(a * xmin + c) / b;
        if (y >= ymin - TOL && y <= ymax + TOL) {
            coord p = {xmin, y};
            if (!is_duplicate(p, intersections, count))
                intersections[count++] = p;
        }
    } else {
        // b is 0 -> vertical line: check if it lies on x = xmin.
        if (fabs(a * xmin + c) < TOL) {
            coord p1 = {xmin, ymin};
            coord p2 = {xmin, ymax};
            if (!is_duplicate(p1, intersections, count))
                intersections[count++] = p1;
            if (!is_duplicate(p2, intersections, count))
                intersections[count++] = p2;
        }
    }

    // Right edge: x = xmax
    if (fabs(b) > TOL) {
        double y = -(a * xmax + c) / b;
        if (y >= ymin - TOL && y <= ymax + TOL) {
            coord p = {xmax, y};
            if (!is_duplicate(p, intersections, count))
                intersections[count++] = p;
        }
    } else {
        if (fabs(a * xmax + c) < TOL) {
            coord p1 = {xmax, ymin};
            coord p2 = {xmax, ymax};
            if (!is_duplicate(p1, intersections, count))
                intersections[count++] = p1;
            if (!is_duplicate(p2, intersections, count))
                intersections[count++] = p2;
        }
    }

    if (count < 2)
        return 0.0;

    // For a convex intersection (a segment), the farthest two points are the endpoints.
    double maxDist = 0.0;
    for (int i = 0; i < count; i++) {
        for (int j = i + 1; j < count; j++) {
            double d = distance(intersections[i], intersections[j]);
            if (d > maxDist)
                maxDist = d;
        }
    }
    return maxDist;
}



#else 

typedef struct {
    coord p;
    double angle;
} PointAngle;

/* Basic vector operations */
coord vec_add(coord a, coord b) {
    coord res = { a.x + b.x, a.y + b.y, a.z + b.z };
    return res;
}

coord vec_sub(coord a, coord b) {
    coord res = { a.x - b.x, a.y - b.y, a.z - b.z };
    return res;
}

coord vec_scale(coord a, double s) {
    coord res = { a.x * s, a.y * s, a.z * s };
    return res;
}

double vec_dot(coord a, coord b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

coord vec_cross(coord a, coord b) {
    coord res = { a.y * b.z - a.z * b.y,
                 a.z * b.x - a.x * b.z,
                 a.x * b.y - a.y * b.x };
    return res;
}

double vec_norm(coord a) {
    return sqrt(vec_dot(a, a));
}

coord vec_normalize(coord a) {
    double norm = vec_norm(a);
    if (norm < TOL)
        return a;
    return vec_scale(a, 1.0 / norm);
}

/* Evaluate plane equation at a point */
double plane_eval(double a, double b, double c, double d, coord p) {
    return a * p.x + b * p.y + c * p.z + d;
}

/* Check if two points are equal within tolerance */
int is_duplicate(coord a, coord b) {
    return (vec_norm(vec_sub(a, b)) < TOL);
}

/* Compute intersection points between the plane and the box edges.
   plane: [a, b, c, d]
   box: [xmin, xmax, ymin, ymax, zmin, zmax]
   points: array to store intersection points; returns the count.
*/
int compute_intersection_points(coord n, double alpha, coord aa, coord bb, coord *points) {
    double a = n.x, b = n.y, c = n.z, d = alpha;
    double xmin = aa.x, xmax = bb.x,
           ymin = aa.y, ymax = bb.y,
           zmin = aa.z, zmax = bb.z;

    /* Define the 8 vertices of the box */
    coord verts[8] = {
        {xmin, ymin, zmin},
        {xmax, ymin, zmin},
        {xmin, ymax, zmin},
        {xmax, ymax, zmin},
        {xmin, ymin, zmax},
        {xmax, ymin, zmax},
        {xmin, ymax, zmax},
        {xmax, ymax, zmax}
    };

    /* Define the 12 edges by pairs of vertex indices */
    int edges[12][2] = {
        {0, 1}, {0, 2}, {0, 4},
        {1, 3}, {1, 5},
        {2, 3}, {2, 6},
        {3, 7},
        {4, 5}, {4, 6},
        {5, 7},
        {6, 7}
    };

    int count = 0;
    for (int i = 0; i < 12; i++) {
        coord p1 = verts[edges[i][0]];
        coord p2 = verts[edges[i][1]];
        double f1 = plane_eval(a, b, c, d, p1);
        double f2 = plane_eval(a, b, c, d, p2);

        if (fabs(f1) < TOL && fabs(f2) < TOL) {
            /* Both endpoints lie on the plane */
            if (count == 0 || !is_duplicate(p1, points[count - 1])) {
                points[count++] = p1;
            }
            if (count == 0 || !is_duplicate(p2, points[count - 1])) {
                points[count++] = p2;
            }
        } else if (fabs(f1) < TOL) {
            /* p1 lies on the plane */
            int duplicate = 0;
            for (int j = 0; j < count; j++) {
                if (is_duplicate(p1, points[j])) {
                    duplicate = 1;
                    break;
                }
            }
            if (!duplicate) points[count++] = p1;
        } else if (fabs(f2) < TOL) {
            /* p2 lies on the plane */
            int duplicate = 0;
            for (int j = 0; j < count; j++) {
                if (is_duplicate(p2, points[j])) {
                    duplicate = 1;
                    break;
                }
            }
            if (!duplicate) points[count++] = p2;
        } else if (f1 * f2 < 0) {
            /* The edge straddles the plane */
            double t = -f1 / (f2 - f1);
            coord intersect = vec_add(p1, vec_scale(vec_sub(p2, p1), t));
            int duplicate = 0;
            for (int j = 0; j < count; j++) {
                if (is_duplicate(intersect, points[j])) {
                    duplicate = 1;
                    break;
                }
            }
            if (!duplicate) points[count++] = intersect;
        }
    }
    return count;
}

/* Compare function for qsort (sorting by angle) */
int compare_point_angle(const void *a, const void *b) {
    double angleA = ((PointAngle *)a)->angle;
    double angleB = ((PointAngle *)b)->angle;
    if (angleA < angleB) return -1;
    if (angleA > angleB) return 1;
    return 0;
}

/* Order the polygon points in a consistent order.
   points: array of points (of length count)
   count: number of points
   plane_normal: the [a, b, c] part of the plane.
   The function reorders the points in place.
*/
void order_polygon_points(coord *points, int count, coord plane_normal, coord *ordered) {
    /* Compute centroid */
    coord centroid = {0, 0, 0};
    for (int i = 0; i < count; i++) {
        centroid = vec_add(centroid, points[i]);
    }
    centroid = vec_scale(centroid, 1.0 / count);

    /* Choose a reference vector in the plane */
    coord ref;
    if (fabs(plane_normal.x) < fabs(plane_normal.y))
        ref = (coord){1, 0, 0};
    else
        ref = (coord){0, 1, 0};
    double norm_n = vec_norm(plane_normal);
    /* Project ref onto the plane */
    double dot_rn = vec_dot(ref, plane_normal);
    coord proj = vec_sub(ref, vec_scale(plane_normal, dot_rn / (norm_n * norm_n)));
    ref = vec_normalize(proj);

    /* Compute angles for each point relative to the centroid */
    PointAngle pa[MAX_POINTS];
    for (int i = 0; i < count; i++) {
        coord v = vec_sub(points[i], centroid);
        /* Project v onto the plane */
        double dot_vn = vec_dot(v, plane_normal);
        coord v_proj = vec_sub(v, vec_scale(plane_normal, dot_vn / (norm_n * norm_n)));
        double dot_rv = vec_dot(ref, v_proj);
        coord cross_rv = vec_cross(ref, v_proj);
        double cp = vec_dot(cross_rv, plane_normal);
        double angle = atan2(cp, dot_rv);
        pa[i].p = points[i];
        pa[i].angle = angle;
    }

    qsort(pa, count, sizeof(PointAngle), compare_point_angle);

    for (int i = 0; i < count; i++) {
        ordered[i] = pa[i].p;
    }
}

/* Compute the area of a planar polygon in 3D by triangulating from the centroid */
double polygon_area_3d(coord *points, int count) {
    if (count < 3)
        return 0.0;

    coord centroid = {0, 0, 0};
    for (int i = 0; i < count; i++) {
        centroid = vec_add(centroid, points[i]);
    }
    centroid = vec_scale(centroid, 1.0 / count);

    double area = 0.0;
    for (int i = 0; i < count; i++) {
        coord v1 = vec_sub(points[i], centroid);
        coord v2 = vec_sub(points[(i + 1) % count], centroid);
        coord cross = vec_cross(v1, v2);
        area += 0.5 * vec_norm(cross);
    }
    return area;
}

/* Compute the cross-sectional area of the intersection */
double compute_cross_section_area(coord n, double alpha, coord aa, coord bb) {
    coord points[MAX_POINTS];
    int count = compute_intersection_points(n, alpha, aa, bb, points);
    if (count < 3)
        return 0.0; // Not enough points for a polygon

    coord ordered[MAX_POINTS];
    order_polygon_points(points, count, n, ordered);
    return polygon_area_3d(ordered, count);
}

#endif