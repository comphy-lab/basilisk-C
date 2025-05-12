/**
## Introduction
In addition to the main recipe for the immersed boundary approach, such as force distribution and velocity interpolation, this file contains the auxiliary components including the constructions of markers and particles, supports to the AMR, and other useful functions.

## Variable list
* `p_n`: number of particles
* `tn_node`: maximum number of markers on each particle
* `S_SIZE`: number of cells within a stencil
* `stiffness`: IBM stiffness
* `motion_type`: 0 - fixed, 1 - free motion, 2 - precribed motion
* `particle_volume`, `moment_inertia`: as named
* `pr`, `pc`, `pv`, `pa`: radius, center coordinates, prescibed velocity and acceleration of particles
* `is_solid`: to obtain fluid cells. 
* `force`: force on fluid cells
*/

#include "lag_stencil.h"
#include "embed.h"

#ifndef p_n
#define p_n (1)
#endif

#if dimension == 2
#define tn_node (1200)
#elif dimension == 3
#define tn_node (33000)
#endif

#if dimension == 2
#if IBM_stencil == 1 || IBM_stencil == 11
#define S_SIZE (9)
#elif IBM_stencil == 4 || IBM_stencil == 14
#define S_SIZE (25)
#endif
#elif dimension == 3
#if IBM_stencil == 1 || IBM_stencil == 11
#define S_SIZE (27)
#elif IBM_stencil == 4 || IBM_stencil == 14
#define S_SIZE (125)
#endif
#endif

#ifndef particle_stiffness
#if dimension == 2
#define particle_stiffness (0.1)
#elif dimension == 3
#define particle_stiffness (0.8)
#endif
#endif

#ifndef motion_type
#define motion_type (0)
#endif

#if dimension == 2
#define particle_volume (M_PI * sq(pr[k]))
#define moment_inertia  (M_PI * sq(pr[k]) * sq(pr[k]) / 2.)
#elif dimension == 3
#define particle_volume (4./3.*M_PI*cube(pr[k]))
#define moment_inertia  (2. * particle_volume * sq(pr[k]) / 5.)
#endif 

double pr[p_n];
coord pc[p_n];
coord pv[p_n];
coord pa[p_n];
scalar is_solid[];
vector force[];

/**
We define a struct for each Lagrangian marker, in which contains:

* `pos`: current position
* `pos_ref`: reference position
* `vel`: velocity
* `agl_vel`: angular velocity
* `acl`: acceleration
* `agl_acl`: angular acceleration
* `lag_F` & `lag_T`: force and torque on the marker
* `my_stencil`: interpolation stencil 
*/

typedef struct
{
    coord pos;     
    coord pos_ref;
    coord vel;
    #if motion_type == 1
    coord vel_pre;
    coord acl;
    coord agl_vel;
    coord agl_vel_pre;
    coord agl_acl;
    #endif
    coord lag_F;
    coord lag_T;
    Cache my_stencil;
}node_struct;

/**
We define a struct for each particle, in which contains:

* `num_nodes`: number of Lagrangian markers
* `radius`: calibrated radius
* `rho_ratio`: density (ratio to the fluid density)
* `F_tot` & `T_tot`: Total force and torque on the particle
* `phi`: for visualization purpose
* `lag_F` & `lag_T`: force and torque on the marker

Each particle has a center and multiple markers.
*/

typedef struct 
{
    int num_nodes;    
    double radius;    
    double rho_ratio;
    coord F_tot;
    coord T_tot;
    coord angle;
    #if dimension == 2
        double (* phi) (double x, double y, double radius);
    #else
        double (* phi) (double x, double y, double z, double radius);
    #endif
    node_struct center;
    node_struct* node;  
}particle_struct;

particle_struct particle[p_n] = {0};

#if dimension == 2
double p_phi (double x, double y, double radius)
{
    return sq (x) + sq (y) - sq (radius);
}
#else // dimension == 3
double p_phi (double x, double y, double z, double radius)
{
    return sq (x) + sq (y) + sq(z) - sq (radius);
}
#endif

#define ADD_PARTICLE_PERIODICITY_X for (double xp = -(L0); xp <= (L0); xp += (L0))
#define ADD_PARTICLE_PERIODICITY_Y for (double yp = -(L0); yp <= (L0); yp += (L0))
#define ADD_PARTICLE_PERIODICITY_Z for (double zp = -(L0); zp <= (L0); zp += (L0))

/**
## Functions
To visualize the particles (copy from [Arthur's sandbox](http://basilisk.fr/sandbox/ghigo/src/myembed-moving-bicolor.h)):
*/

void p_shape(scalar c, face vector f, particle_struct *p)
{
    vertex scalar phi[];
    foreach_vertex()
    {
        phi[] = HUGE;
        for (int k = 0; k < p_n; k++)
        {
            
            ADD_PARTICLE_PERIODICITY_X
                ADD_PARTICLE_PERIODICITY_Y
                    #if dimension == 2
                    phi[] = intersection(phi[], p[k].phi((x + xp - p[k].center.pos_ref.x), (y + yp - p[k].center.pos_ref.y), p[k].radius));
                    #else // dimension == 3
                    ADD_PARTICLE_PERIODICITY_Z
                        phi[] = intersection(phi[], p[k].phi((x + xp - p[k].center.pos_ref.x), (y + yp - p[k].center.pos_ref.y), (z + zp - p[k].center.pos_ref.z), p[k].radius));
                    #endif
        }
    }
    fractions(phi, c, f);
}

/**
To compute the total force and torque on particles:
*/
void write_force(particle_struct* p, double uc)
{
    for (int k = 0; k < p_n; k++)
    {
        foreach_dimension()
        {
            particle[k].F_tot.x = 0.;
            particle[k].T_tot.z = 0.;
        }

        for (int n = 0; n < particle[k].num_nodes; n++)
        {
            foreach_dimension()
            {
                particle[k].F_tot.x += particle[k].node[n].lag_F.x;
                particle[k].T_tot.z += particle[k].node[n].lag_T.z;
            }
        }

/*
        foreach_dimension()
        {
            #if dimension == 2
            particle[k].F_tot.x /= -(uc * uc * particle[k].radius);
            particle[k].T_tot.z /= -(uc * uc * sq(particle[k].radius));
            #elif dimension == 3
            particle[k].F_tot.x /= -(0.5 * sq(uc) * M_PI * sq(particle[k].radius));
            particle[k].T_tot.z /= -(0.5 * sq(uc) * M_PI * cube(particle[k].radius));
            #endif
        }
*/
    }
}

/**
This function is used at the initialization stage when create the interpolation stencil for each Lagrangian marker. Cells in certain range will be included in the stencil. Users will see warnings if any selected cells are not at the finest grid level.
*/

void create_stencil_cache(particle_struct *p)
{
    double lx, ly;
    #if dimension == 3
    double lz;
    #endif
    Point lpoint;
    #if (IBM_stencil == 1 || IBM_stencil == 11)
    {
        for (int k = 0; k < p_n; k++)
        {
            for (int n = 0; n < p[k].num_nodes; n++)
            {
                bool stencil_check = true;
                for (int ix = -1; ix <= 1; ix++)
                {
                    lx = floor(p[k].node[n].pos.x) + 0.5 + ix;
                    for (int iy = -1; iy <= 1; iy++)
                    {
                        ly = floor(p[k].node[n].pos.y) + 0.5 + iy;
                        #if dimension == 2
                            lpoint = locate(lx, ly);
                            cache_append(&(p[k].node[n].my_stencil), lpoint, 0);
                            if (lpoint.level > -1 && lpoint.level < grid->maxdepth && stencil_check == true)
                            {
                                fprintf(stderr,"Fault stencil: No.%d node of No.%d particle!\n", n, k);
                                stencil_check = false;
                            }
                        #else // dimension == 3
                        for (int iz = -1; iz <= 1; iz++)
                        {
                            lz = floor(p[k].node[n].pos.z) + 0.5 + iz;
                            lpoint = locate(lx, ly, lz);
                            cache_append(&(p[k].node[n].my_stencil), lpoint, 0);
                            if (lpoint.level > -1 && lpoint.level < grid->maxdepth && stencil_check == true)
                            {
                                fprintf(stderr,"Fault stencil: No.%d node of No.%d particle!\n", n, k);
                                stencil_check = false;
                            }
                        }
                        #endif
                    }
                }
                cache_shrink(&(p[k].node[n].my_stencil));
            }
        }
    }
    #else // (IBM_stencil == 4 || IBM_stencil == 14)
    {
        for (int k = 0; k < p_n; k++)
        {
            for (int n = 0; n < p[k].num_nodes; n++)
            {
                bool stencil_check = true;
                for (int ix = -2; ix <= 2; ix++)
                {
                    lx = floor(p[k].node[n].pos.x) + 0.5 + ix;
                    for (int iy = -2; iy <= 2; iy++)
                    {
                        ly = floor(p[k].node[n].pos.y) + 0.5 + iy;
                        #if dimension == 2
                            lpoint = locate(lx, ly);
                            cache_append(&(p[k].node[n].my_stencil), lpoint, 0);
                            if (lpoint.level > -1 && lpoint.level < grid->maxdepth && stencil_check == true)
                            {
                                fprintf(stderr,"Fault stencil: No.%d node of No.%d particle!\n", n, k);
                                stencil_check = false;
                            }
                        #else // dimension == 3
                        for (int iz = -2; iz <= 2; iz++)
                        {
                            lz = floor(p[k].node[n].pos.z) + 0.5 + iz;
                            lpoint = locate(lx, ly, lz);
                            cache_append(&(p[k].node[n].my_stencil), lpoint, 0);
                            if (lpoint.level > -1 && lpoint.level < grid->maxdepth && stencil_check == true)
                            {
                                fprintf(stderr,"Fault stencil: No.%d node of No.%d particle!\n", n, k);
                                stencil_check = false;
                            }
                        }
                        #endif
                    }
                }
                cache_shrink(&(p[k].node[n].my_stencil));
            }
        }
    }
    #endif
}

/**
On an octree-grid, cells are identified by its global indices, not the coordinates.
*/
static void change_cache_entry(Cache *s, int i, Point pt, int flag)
{
    if (i > s->n)
        fprintf(stderr, "Error: Cache index out of range.\n");
    s->p[i].i = pt.i;
    s->p[i].j = pt.j;
    #if dimension == 3
    s->p[i].k = pt.k;
    #endif
    s->p[i].level = pt.level;
    s->p[i].flags = flag;
}

/**
When the AMR is employed, we need to give the correct indices to stencil cells after each adaptation.
*/
void Rearrange_indices(particle_struct *p)
{
    double lx, ly;
    #if dimension == 3
    double lz;
    #endif
    #if (IBM_stencil == 1 || IBM_stencil == 11)
    {
        for (int k = 0; k < p_n; k++)
        {
            for (int n = 0; n < p[k].num_nodes; n++)
            {
                int c = 0;
                bool stencil_check = true;
                for (int ix = -1; ix <= 1; ix++)
                {
                    lx = floor(p[k].node[n].pos.x) + 0.5 + ix;
                    for (int iy = -1; iy <= 1; iy++)
                    {
                        ly = floor(p[k].node[n].pos.y) + 0.5 + iy;
                        #if dimension == 2
                        Point lpoint = locate(lx, ly);
                        if (lpoint.level > -1 && lpoint.level < grid->maxdepth  && stencil_check == true)
                        {
                            fprintf(stderr, "==Fault stencil: No.%d node of No.%d particle!\n", n, k);
                            stencil_check = false;
                        }
                        change_cache_entry(&(p[k].node[n].my_stencil), c, lpoint, 0);
                        c++;
                        #else // dimension == 3
                        for (int iz = -1; iz <= 1; iz++)
                        {
                            lz = floor(p[k].node[n].pos.z) + 0.5 + iz;
                            Point lpoint = locate(lx, ly, lz);
                            if (lpoint.level > -1 && lpoint.level < grid->maxdepth && stencil_check == true)
                            {
                                fprintf(stderr,"==Fault stencil: No.%d node of No.%d particle! Cell level is %d\n", n, k, lpoint.level);
                                stencil_check = false;
                            }
                            change_cache_entry(&(p[k].node[n].my_stencil), c, lpoint, 0);
                            c++;
                        }
                        #endif
                    }
                }
            }
        }
    }
    #else // (IBM_stencil == 4 || IBM_stencil == 14)
    {
        for (int k = 0; k < p_n; k++)
        {
            for (int n = 0; n < p[k].num_nodes; n++)
            {
                int c = 0;
                bool stencil_check = true;
                for (int ix = -2; ix <= 2; ix++)
                {
                    lx = floor(p[k].node[n].pos.x) + 0.5 + ix;
                    for (int iy = -2; iy <= 2; iy++)
                    {
                        ly = floor(p[k].node[n].pos.y) + 0.5 + iy;
                        #if dimension == 2
                        Point lpoint = locate(lx, ly);
                        if (lpoint.level > -1 && lpoint.level < grid->maxdepth && stencil_check == true)
                        {
                            fprintf(stderr, "==Fault stencil: No.%d node of No.%d particle!\n", n, k);
                            stencil_check = false;
                        }
                        change_cache_entry(&(p[k].node[n].my_stencil), c, lpoint, 0);
                        c++;
                        #else // dimension == 3
                        for (int iz = -2; iz <= 2; iz++)
                        {
                            lz = floor(p[k].node[n].pos.z) + 0.5 + iz;
                            Point lpoint = locate(lx, ly, lz);
                            if (lpoint.level > -1 && lpoint.level < grid->maxdepth && stencil_check == true)
                            {
                                fprintf(stderr,"==Fault stencil: No.%d node of No.%d particle!\n", n, k);
                                stencil_check = false;
                            }
                            change_cache_entry(&(p[k].node[n].my_stencil), c, lpoint, 0);
                            c++;
                        }
                        #endif
                    }
                }
            }
        }
    }
    #endif
}

/**
We certainly want to store the information of particles for restarting simulations.
*/
void dump_particle(particle_struct *p, int time)
{
    FILE *fp;
    char file[80];
    sprintf(file,"dump_particle_%d",time);

    if ((fp = fopen(file, "w")) == NULL)
    {
        perror(file);
        exit(1);
    }
    assert(fp);

    for(int k = 0; k < p_n; k++)
    {
        fwrite(&(p[k].num_nodes), sizeof(int), 1, fp);
        fwrite(&(p[k].radius), sizeof(double), 1, fp);
        fwrite(&(p[k].stiffness), sizeof(double), 1, fp);
        foreach_dimension()
        {
            fwrite(&(p[k].F_tot.x), sizeof(double), 1, fp);
            fwrite(&(p[k].T_tot.x), sizeof(double), 1, fp);
            fwrite(&(p[k].angle.z), sizeof(double), 1, fp);
            fwrite(&(p[k].center.pos.x), sizeof(double), 1, fp);
            fwrite(&(p[k].center.pos_ref.x), sizeof(double), 1, fp);
            fwrite(&(p[k].center.vel.x), sizeof(double), 1, fp);
            fwrite(&(p[k].center.lag_F.x), sizeof(double), 1, fp);
            fwrite(&(p[k].center.lag_T.x), sizeof(double), 1, fp);
            #if motion_type == 1
            fwrite(&(p[k].center.vel_pre.x), sizeof(double), 1, fp);
            fwrite(&(p[k].center.acl.x), sizeof(double), 1, fp);
            fwrite(&(p[k].center.agl_vel.x), sizeof(double), 1, fp);
            fwrite(&(p[k].center.agl_vel_pre.x), sizeof(double), 1, fp);
            fwrite(&(p[k].center.agl_acl.x), sizeof(double), 1, fp);
            #endif
        }
        for (int n = 0; n < p[k].num_nodes; n++)
        {
            foreach_dimension()
            {
                fwrite(&(p[k].node[n].pos.x), sizeof(double), 1, fp);
                fwrite(&(p[k].node[n].pos_ref.x), sizeof(double), 1, fp);
                fwrite(&(p[k].node[n].vel.x), sizeof(double), 1, fp);
                fwrite(&(p[k].node[n].lag_F.x), sizeof(double), 1, fp);
                fwrite(&(p[k].node[n].lag_T.x), sizeof(double), 1, fp);
            }
        }
    }
    fclose(fp);
}

/**
Correspondingly, we need to restore simulations from the dump files.
*/
bool restore_particle(particle_struct *p, int time)
{
    FILE *fp;
    char file[80];
    sprintf(file,"dump_particle_%d",time);
    if ((fp = fopen(file, "r")) == NULL)
        return false;
    assert(fp);

    for(int k = 0; k < p_n; k++)
    {
        fread(&(p[k].num_nodes), sizeof(int), 1, fp);
        fread(&(p[k].radius), sizeof(double), 1, fp);
        fread(&(p[k].stiffness), sizeof(double), 1, fp);
        foreach_dimension()
        {
            fread(&(p[k].F_tot.x), sizeof(double), 1, fp);
            fread(&(p[k].T_tot.x), sizeof(double), 1, fp);
            fread(&(p[k].angle.z), sizeof(double), 1, fp);
            fread(&(p[k].center.pos.x), sizeof(double), 1, fp);
            fread(&(p[k].center.pos_ref.x), sizeof(double), 1, fp);
            fread(&(p[k].center.vel.x), sizeof(double), 1, fp);
            fread(&(p[k].center.lag_F.x), sizeof(double), 1, fp);
            fread(&(p[k].center.lag_T.x), sizeof(double), 1, fp);
            #if motion_type == 1
            fread(&(p[k].center.vel_pre.x), sizeof(double), 1, fp);
            fread(&(p[k].center.acl.x), sizeof(double), 1, fp);
            fread(&(p[k].center.agl_vel.x), sizeof(double), 1, fp);
            fread(&(p[k].center.agl_vel_pre.x), sizeof(double), 1, fp);
            fread(&(p[k].center.agl_acl.x), sizeof(double), 1, fp);
            #endif
        }
        for (int n = 0; n < p[k].num_nodes; n++)
        {
            foreach_dimension()
            {
                fread(&(p[k].node[n].pos.x), sizeof(double), 1, fp);
                fread(&(p[k].node[n].pos_ref.x), sizeof(double), 1, fp);
                fread(&(p[k].node[n].vel.x), sizeof(double), 1, fp);
                fread(&(p[k].node[n].lag_F.x), sizeof(double), 1, fp);
                fread(&(p[k].node[n].lag_T.x), sizeof(double), 1, fp);
            }
        }
        p[k].phi = p_phi;
    }
        
    create_stencil_cache(p);
    fclose(fp);
    return true;
}

/**
Don't forget to return anything you borrowed.
*/
void free_stencil_cache(particle_struct *p)
{
    for (int k = 0; k < p_n; k++)
    {
        for (int n = 0; n < p[k].num_nodes; n++)
        {
            Cache * c = &(p[k].node[n].my_stencil);
            free(c->p);
        }
        free(p[k].node);
    }
}

/**
Here we initialize the particle system, including giving the proper values of parameters, assigning memory for markers, distribute markers, defining periodicity, and generating cache for interpolation stencils.
*/
event Init_Particle(i = 0)
{
    for (int k = 0; k < p_n; k++)
    {
        particle[k].phi = p_phi;
        particle[k].stiffness = particle_stiffness;
        particle[k].rho_ratio = 1.;
        foreach_dimension()
        {
            particle[k].F_tot.x = 0.;
            particle[k].T_tot.x = 0.;
            particle[k].angle.z = 0.;
            particle[k].center.pos.x = fmod(pc[k].x,L0);
            particle[k].center.pos_ref.x = particle[k].center.pos.x;
        }
        #if dimension == 2
        particle[k].radius = pr[k];
        particle[k].num_nodes = 64 * (int)(pr[k] / 10);
        for (int n = 0; n < particle[k].num_nodes; n++)
        {
            particle[k].node[n].pos.x = particle[k].center.pos.x + particle[k].radius * cos(2. * M_PI * (double)n / particle[k].num_nodes);
            particle[k].node[n].pos.y = particle[k].center.pos.y + particle[k].radius * sin(2. * M_PI * (double)n / particle[k].num_nodes);
            foreach_dimension()
                particle[k].node[n].pos_ref.x = particle[k].node[n].pos.x;
        }
        #else // dimension == 3
        particle[k].radius = pow((cube(pr[k]) + cube(pr[k] - 1.))/2.,1./3.);
        particle[k].num_nodes = ((int)(M_PI * sq(particle[k].radius)));
        particle[k].node = malloc(particle[k].num_nodes*sizeof(node_struct));
        double indice = 0.;
        double phi    = 0.;
        double theta  = 0.;
        for (int n = 0; n < particle[k].num_nodes; n++)
        {
            indice = (double)n + 0.5;
            phi = acos(1. - 2. * indice / particle[k].num_nodes);
            theta = M_PI * (1. + sqrt(5.)) * indice;
            particle[k].node[n].pos.x = particle[k].center.pos.x + particle[k].radius * cos(theta) * sin(phi);
            if (particle[k].node[n].pos.x < -L0 / 2)
                particle[k].node[n].pos.x += L0;
            else if (particle[k].node[n].pos.x > L0 / 2)
                particle[k].node[n].pos.x -= L0;
            particle[k].node[n].pos.y = particle[k].center.pos.y + particle[k].radius * sin(theta) * sin(phi);
            if (particle[k].node[n].pos.y < -L0 / 2)
                particle[k].node[n].pos.y += L0;
            else if (particle[k].node[n].pos.y > L0 / 2)
                particle[k].node[n].pos.y -= L0;
            particle[k].node[n].pos.z = particle[k].center.pos.z + particle[k].radius * cos(phi);
            if (particle[k].node[n].pos.z < -L0 / 2)
                particle[k].node[n].pos.z += L0;
            else if (particle[k].node[n].pos.z > L0 / 2)
                particle[k].node[n].pos.z -= L0;
            foreach_dimension()
                particle[k].node[n].pos_ref.x = particle[k].node[n].pos.x;

            particle[k].node[n].my_stencil.n = 0;
            particle[k].node[n].my_stencil.nm = S_SIZE;
            particle[k].node[n].my_stencil.p = (Index*) malloc(S_SIZE*sizeof(Index));
         
        }
        #endif
    }

    FILE *fp;
    if ((fp = fopen("restart", "r")) == NULL)
        create_stencil_cache(particle);
}