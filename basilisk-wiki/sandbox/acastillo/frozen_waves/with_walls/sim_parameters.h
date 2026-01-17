/**
sim_parameters.h

Clean struct-based parameter reading for main_json.c.
Replaces preprocessor defines with runtime JSON configuration.

*/

/** Feature detection for cJSON library */
#if __has_include(<cjson/cJSON.h>)
  #define HAVE_CJSON 1
  #pragma autolink -lcjson
  #include <cjson/cJSON.h>
  
  /** Helper macro to safely extract double values from JSON */
  #define GET_DOUBLE(parent, name, dest) do { \
      cJSON *item = cJSON_GetObjectItem(parent, #name); \
      if (item && cJSON_IsNumber(item)) { \
          dest = item->valuedouble; \
      } \
  } while(0)
  
  /** Helper macro to safely extract integer values from JSON */
  #define GET_INT(parent, name, dest) do { \
      cJSON *item = cJSON_GetObjectItem(parent, #name); \
      if (item && cJSON_IsNumber(item)) { \
          dest = (int)item->valuedouble; \
      } \
  } while(0)
  
  /** Helper macro to safely extract boolean values from JSON */
  #define GET_BOOL(parent, name, dest) do { \
      cJSON *item = cJSON_GetObjectItem(parent, #name); \
      if (item && cJSON_IsBool(item)) { \
          dest = cJSON_IsTrue(item); \
      } \
  } while(0)
#else
  #define HAVE_CJSON 0
  #warning "cJSON library not found. Using hardcoded default parameters."
#endif

/**
Simulation parameters structure
*/
typedef struct {
    /** Domain parameters */
    double width;           // Domain width
    double height;          // Domain height
    double depth;           // Domain depth
    double aspectratio_x;   // Aspect ratio (nx dimension multiplier)
    double aspectratio_y;   // Aspect ratio (ny dimension multiplier)
    double aspectratio_z;   // Aspect ratio (nz dimension multiplier)
    
    /** Dimensionless parameters */
    double atwood;          // Atwood number
    double reynolds;        // Reynolds number
    double froude;          // Froude number
    double weber;           // Weber number
    double viscosity;       // Viscosity ratio (upsilon)
    
    /** Forcing parameters */
    double ramp_start;      // Time when ramp starts
    double ramp_slope;      // Slope of sigmoid ramp
    
    /** Initial conditions */
    double amplitude;       // RMS amplitude of initial perturbation
    double kmin;            // Minimum wavenumber
    double kmax;            // Maximum wavenumber
    
    /** Contact angle */
    double theta0;          // Contact angle in degrees
    
    /** Numerical parameters */
    int level;              // Grid refinement level
    double cfl;             // CFL number
    double dt;              // Time step
    double tolerance;       // Solver tolerance
    int nitermin;           // Minimum iterations
    
    /** Runtime parameters */
    double walltime;        // Maximum wall time (seconds)
    double endtime;         // Simulation end time
} SimParams;

/**
Initialize parameters with default values
*/
void init_default_params(SimParams *p) {
    /** Domain parameters */
    p->width = 300.0;
    p->height = 150.0;
    p->depth = 75.0;
    p->aspectratio_x = 2.0;
    p->aspectratio_y = 1.0;
    p->aspectratio_z = 1.0;
    
    /** Dimensionless parameters */
    p->atwood = 0.5;
    p->reynolds = 43.62;
    p->froude = 64.05;
    p->weber = 604.80;
    p->viscosity = 11.5;
    
    /** Forcing parameters */
    p->ramp_start = 100.0;
    p->ramp_slope = 0.04;
    
    /** Initial conditions */
    p->amplitude = 0.05;
    p->kmin = 15.0;
    p->kmax = 30.0;
    
    /** Contact angle */
    p->theta0 = 90.0;
    
    /** Numerical parameters */
    p->level = 8;
    p->cfl = 0.25;
    p->dt = 0.1;
    p->tolerance = 1e-4;
    p->nitermin = 1;
    
    /** Runtime parameters */
    p->walltime = 3600.0;
    p->endtime = 500.0;
}

#if HAVE_CJSON
/**
Parse domain parameters from JSON object
*/
void parse_domain_section(cJSON *domain, SimParams *p) {
    if (!domain) return;
    
    GET_DOUBLE(domain, width, p->width);
    GET_DOUBLE(domain, height, p->height);
    GET_DOUBLE(domain, depth, p->depth);
    
    double min_dim = p->width < p->height ? p->width : p->height;
    p->aspectratio_x = p->width / min_dim;
    p->aspectratio_y = p->height / min_dim;
    p->aspectratio_z = p->depth / min_dim;
}

/**
Parse dimensionless parameters from JSON object
*/
void parse_dimensionless_section(cJSON *dimless, SimParams *p) {
    if (!dimless) return;
    
    GET_DOUBLE(dimless, atwood, p->atwood);
    GET_DOUBLE(dimless, reynolds, p->reynolds);
    GET_DOUBLE(dimless, froude, p->froude);
    GET_DOUBLE(dimless, weber, p->weber);
    GET_DOUBLE(dimless, viscosity, p->viscosity);
}

/**
Parse forcing parameters from JSON object
*/
void parse_forcing_section(cJSON *forcing, SimParams *p) {
    if (!forcing) return;
    
    GET_DOUBLE(forcing, ramp_start, p->ramp_start);
    GET_DOUBLE(forcing, ramp_slope, p->ramp_slope);
}

/**
Parse initial condition parameters from JSON object
*/
void parse_initial_conditions_section(cJSON *ic, SimParams *p) {
    if (!ic) return;
    
    GET_DOUBLE(ic, amplitude, p->amplitude);
    GET_DOUBLE(ic, kmin, p->kmin);
    GET_DOUBLE(ic, kmax, p->kmax);
}

/**
Parse contact angle from JSON object
*/
void parse_contact_section(cJSON *contact, SimParams *p) {
    if (!contact) return;
    
    GET_DOUBLE(contact, theta0, p->theta0);
}

/**
Parse runtime parameters from JSON object
*/
void parse_runtime_section(cJSON *runtime, SimParams *p) {
    if (!runtime) return;
    
    GET_DOUBLE(runtime, walltime, p->walltime);
    GET_DOUBLE(runtime, endtime, p->endtime);
}

/**
Parse numerical parameters from JSON object
*/
void parse_numerics_section(cJSON *numerics, SimParams *p) {
    if (!numerics) return;
    
    GET_INT(numerics, level, p->level);
    GET_DOUBLE(numerics, cfl, p->cfl);
    GET_DOUBLE(numerics, dt, p->dt);
    GET_DOUBLE(numerics, tolerance, p->tolerance);
    GET_INT(numerics, nitermin, p->nitermin);
}
#endif



/**
Read and parse JSON file into SimParams structure

Returns 0 on success, -1 on error. If the file cannot be read or parsed,
default parameters are used and a warning is printed.
*/
int read_sim_params(const char *filename, SimParams *params) {
    /** Initialize with defaults first */
    init_default_params(params);
    
#if !HAVE_CJSON
    fprintf(stderr, "Warning: cJSON not available, using default parameters\n");
    return -1;
#else
    /** Read file */
    FILE *file = fopen(filename, "rb");
    if (!file) {
        fprintf(stderr, "Warning: Could not open file '%s'\n", filename);
        fprintf(stderr, "Using default parameters\n");
        return -1;
    }
    
    /** Get file size */
    fseek(file, 0, SEEK_END);
    long length = ftell(file);
    fseek(file, 0, SEEK_SET);
    
    /** Read content */
    char *content = (char *)malloc(length + 1);
    if (!content) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        fclose(file);
        return -1;
    }
    
    fread(content, 1, length, file);
    content[length] = '\0';
    fclose(file);
    
    /** Parse JSON */
    cJSON *json = cJSON_Parse(content);
    free(content);
    
    if (!json) {
        const char *error_ptr = cJSON_GetErrorPtr();
        if (error_ptr) {
            fprintf(stderr, "Error: JSON parse error before: %s\n", error_ptr);
        }
        fprintf(stderr, "Using default parameters\n");
        return -1;
    }
    
    /** Parse sections */
    parse_domain_section(cJSON_GetObjectItem(json, "domain"), params);
    parse_dimensionless_section(cJSON_GetObjectItem(json, "dimensionless"), params);
    parse_forcing_section(cJSON_GetObjectItem(json, "forcing"), params);
    parse_initial_conditions_section(cJSON_GetObjectItem(json, "initial_conditions"), params);
    parse_contact_section(cJSON_GetObjectItem(json, "contact"), params);
    parse_numerics_section(cJSON_GetObjectItem(json, "numerics"), params);
    parse_runtime_section(cJSON_GetObjectItem(json, "runtime"), params);
    
    cJSON_Delete(json);
    return 0;
#endif
}

/**
MPI-aware version: only rank 0 reads, then broadcasts
*/
void read_sim_params_mpi(const char *filename, SimParams *params) {
    if (pid() == 0) {
        read_sim_params(filename, params);
    }
    
@if _MPI
    /** Broadcast the entire struct */
    MPI_Bcast(params, sizeof(SimParams), MPI_BYTE, 0, MPI_COMM_WORLD);
@endif
}

/**
Print parameters for verification
*/
void print_sim_params(const SimParams *p) {
    fprintf(stderr, "\n=== SIMULATION PARAMETERS ===\n");
    fprintf(stderr, "Domain:\n");
    fprintf(stderr, "  Width       : %g\n", p->width);
    fprintf(stderr, "  Height      : %g\n", p->height);
    fprintf(stderr, "  Depth       : %g\n", p->depth);
    fprintf(stderr, "  Aspect Ratio: %g\n", p->aspectratio_x);
    fprintf(stderr, "  Aspect Ratio: %g\n", p->aspectratio_y);
    fprintf(stderr, "  Aspect Ratio: %g\n", p->aspectratio_z);
    fprintf(stderr, "\nDimensionless:\n");
    fprintf(stderr, "  Atwood      : %g\n", p->atwood);
    fprintf(stderr, "  Reynolds    : %g\n", p->reynolds);
    fprintf(stderr, "  Froude      : %g\n", p->froude);
    fprintf(stderr, "  Weber       : %g\n", p->weber);
    fprintf(stderr, "  Viscosity   : %g\n", p->viscosity);
    fprintf(stderr, "\nForcing:\n");
    fprintf(stderr, "  Ramp Start  : %g\n", p->ramp_start);
    fprintf(stderr, "  Ramp Slope  : %g\n", p->ramp_slope);
    fprintf(stderr, "\nInitial Conditions:\n");
    fprintf(stderr, "  Amplitude   : %g\n", p->amplitude);
    fprintf(stderr, "  kmin        : %g\n", p->kmin);
    fprintf(stderr, "  kmax        : %g\n", p->kmax);
    fprintf(stderr, "\nContact:\n");
    fprintf(stderr, "  Theta0      : %g°\n", p->theta0);
    fprintf(stderr, "\nNumerics:\n");
    fprintf(stderr, "  Level       : %d\n", p->level);
    fprintf(stderr, "  CFL         : %g\n", p->cfl);
    fprintf(stderr, "  DT          : %g\n", p->dt);
    fprintf(stderr, "  Tolerance   : %g\n", p->tolerance);
    fprintf(stderr, "  NITERMIN    : %d\n", p->nitermin);
    fprintf(stderr, "\nRuntime:\n");
    fprintf(stderr, "  Wall Time   : %g s\n", p->walltime);
    fprintf(stderr, "  End Time    : %g\n", p->endtime);
    fprintf(stderr, "=============================\n\n");
}
