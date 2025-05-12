#pragma autolink -L${CJSON_LIBDIR} -lcjson
#include <cjson/cJSON.h>

int file_exists(const char *filename)
{
	FILE *file = fopen(filename, "r");
	if (file)
	{
		fclose(file);
		return 1; // File exists
	}
	return 0; // File does not exist
}

// Reads the entire content of a file and returns it as a null-terminated string.
char *read_file(const char *filename)
{
	int read_size;
	NOT_UNUSED(read_size);

	FILE *file = fopen(filename, "rb"); // Open file for reading in binary mode
	if (file == NULL)
	{
		printf("Could not open file %s\n", filename);
		return NULL;
	}

	// Seek to the end of the file to determine its length
	fseek(file, 0, SEEK_END);
	long length = ftell(file); // Get the file length
	fseek(file, 0, SEEK_SET);  // Reset file pointer to the beginning

	// Allocate memory to hold the file contents plus a null terminator
	char *content = (char *)malloc(length + 1);
	if (content == NULL)
	{
		printf("Memory allocation error\n");
		fclose(file);
		return NULL;
	}

	// Read the file contents into the allocated memory
	read_size = fread(content, 1, length, file);
	content[length] = '\0'; // Null-terminate the string

	fclose(file); // Close the file
	return content;
}

// Parses the JSON string and assigns the values to the provided array.
void parse_and_assign(const char *json_string, double *values)
{
	cJSON *json;
	cJSON *fluids;
	cJSON *geometry;
	cJSON *experiment;
	cJSON *perturbation;
	cJSON *numerics;

	// Parse the JSON string
	json = cJSON_Parse(json_string);
	if (json == NULL)
	{
		printf("Error parsing JSON\n");
		return;
	}

	// Parse fluids section
	fluids = cJSON_GetObjectItem(json, "fluids");
	if (fluids != NULL)
	{
		values[0] = cJSON_GetObjectItem(fluids, "atwood")->valuedouble;
		values[1] = cJSON_GetObjectItem(fluids, "density1")->valuedouble;
		values[2] = cJSON_GetObjectItem(fluids, "density2")->valuedouble;
		values[3] = cJSON_GetObjectItem(fluids, "viscosity1")->valuedouble;
		values[4] = cJSON_GetObjectItem(fluids, "viscosity2")->valuedouble;
		values[5] = cJSON_GetObjectItem(fluids, "tension")->valuedouble;
	}

	// Parse geometry section
	geometry = cJSON_GetObjectItem(json, "geometry");
	if (geometry != NULL)
	{
		values[6] = cJSON_GetObjectItem(geometry, "aspect_ratio_x")->valuedouble;
		values[7] = cJSON_GetObjectItem(geometry, "aspect_ratio_y")->valuedouble;
		values[8] = cJSON_GetObjectItem(geometry, "aspect_ratio_z")->valuedouble;
		values[9] = cJSON_GetObjectItem(geometry, "width")->valuedouble;
	}

	// Parse experiment section
	experiment = cJSON_GetObjectItem(json, "experiment");
	if (experiment != NULL)
	{
		values[10] = cJSON_GetObjectItem(experiment, "gravity")->valuedouble;
		values[11] = cJSON_GetObjectItem(experiment, "velocity")->valuedouble;
		values[12] = cJSON_GetObjectItem(experiment, "frequency")->valuedouble;
		values[13] = cJSON_GetObjectItem(experiment, "velocity_prev")->valuedouble;
		values[14] = cJSON_GetObjectItem(experiment, "ramp_type")->valuedouble;
		values[15] = cJSON_GetObjectItem(experiment, "time_start")->valuedouble;
		values[16] = cJSON_GetObjectItem(experiment, "shape_m")->valuedouble;
		values[17] = cJSON_GetObjectItem(experiment, "shape_k")->valuedouble;
		values[18] = cJSON_GetObjectItem(experiment, "time_end")->valuedouble;		
	}

	// Parse perturbation section
	perturbation = cJSON_GetObjectItem(json, "perturbation");
	if (perturbation != NULL)
	{
		values[19] = cJSON_GetObjectItem(perturbation, "eta0")->valuedouble;
		values[20] = cJSON_GetObjectItem(perturbation, "k0")->valuedouble;
		values[21] = cJSON_GetObjectItem(perturbation, "delta_k")->valuedouble;
	}

	// Parse numerics section
	numerics = cJSON_GetObjectItem(json, "numerics");
	if (numerics != NULL)
	{
		values[22] = cJSON_GetObjectItem(numerics, "CFL")->valuedouble;
		values[23] = cJSON_GetObjectItem(numerics, "DT")->valuedouble;
		values[24] = cJSON_GetObjectItem(numerics, "TOLERANCE")->valuedouble;
		values[25] = cJSON_GetObjectItem(numerics, "NITERMIN")->valuedouble;		
	}

	// Cleanup the cJSON object
	cJSON_Delete(json);
}

void read_and_parse_json_mpi(const char *filename, double *values)
{
	char *json_string = NULL; // Pointer to hold JSON string

	if (pid() == 0)
	{
		// Root process reads the file
		json_string = read_file(filename);

		if (json_string != NULL)
		{
			parse_and_assign(json_string, values); // Parse JSON and assign values
			free(json_string);					   // Free the allocated memory for the JSON string
		}
	}


@if _MPI
	// Broadcast the parsed values to all processes
	MPI_Bcast(values, 26, MPI_DOUBLE, 0, MPI_COMM_WORLD);
@endif
}

// Define structures to hold reference and simulation parameters
struct ReferenceValues
{
	double m;
	double Lx;
	double G0;
	double T;
	double Omega;
	double rho;
	double sigma;
	double mu;
	double Re;
	double Fr;
	double We;
	double At;
	double Bo;
	double Oh;
	double B;
	double S;
	double D;
};

struct PerturbationDefect
{
	double eta0;
	double k0;
	double dk;
	double kmin;
	double kmax;
};

struct PeriodicForcing
{
	double eta0;
	double kmin;
	double kmax;
	double G0;
	double V0;
	double V0prev;
	double F0;
	double F0prev;
	double freq0;
	double omega0;
	double period0;
	double Gn;
	double Gnm1;
	double ramp;
};

struct ReferenceValues ref = {0., 1., 1., 1., 6.28, 1., 1., 1., 1., 1., 1., 1., 1., 1., 0., 0., 0.};
struct ReferenceValues sim = {0., 1., 1., 1., 6.28, 1., 1., 1., 1., 1., 1., 1., 1., 1., 0., 0., 0.};
struct PeriodicForcing force;
struct PerturbationDefect defect;

void print_parameters(struct ReferenceValues ref, struct PeriodicForcing force, struct PerturbationDefect defect, struct ReferenceValues sim)
{
	fputs("REFERENCE VALUES\n", stderr);
	fprintf(stderr, " Mode       : %g \n", ref.m);
	fprintf(stderr, " Length     : %g \t\t[m]    \n", ref.Lx);
	fprintf(stderr, " Gravity    : %g \t\t[m/s2] \n", ref.G0);
	fprintf(stderr, " Time       : %g \t\t[s]\n", ref.T);
	fprintf(stderr, " Frequency  : %g \t\t[rad/s]\n", ref.Omega);
	fprintf(stderr, " Density    : %g \t\t[kg/m3]\n", ref.rho);
	fprintf(stderr, " Viscosity  : %g \t\t[Pa s] \n", ref.mu);
	fprintf(stderr, " Tension C. : %g \t\t[N/m]  \n", ref.sigma);
	
	fputs("SIMULATION PARAMETERS\n", stderr);	
	fprintf(stderr, " Length     : %g \n", L0);
	fprintf(stderr, " eta0       : %g \n", defect.eta0);
	fprintf(stderr, " k0         : %g \n", defect.k0);
	fprintf(stderr, " delta k    : %g \n", defect.dk);
	fprintf(stderr, " Wav. min.  : %g \n", defect.kmin);
	fprintf(stderr, " Wav. max.  : %g \n", defect.kmax);
	fprintf(stderr, " Gravity    : %g \n", force.G0);
	fprintf(stderr, " Forcing    : %g \n", force.F0);
	fprintf(stderr, " Velocity   : %g \n", force.V0);
	fprintf(stderr, " Frequency  : %g \n", force.omega0);
	fprintf(stderr, " Period     : %g \n", force.period0);
	fprintf(stderr, " rho1       : %g \n", rho1);
	fprintf(stderr, " rho2       : %g \n", rho2);
	fprintf(stderr, " rho1/rho2  : %g \n", rho1 / rho2);
	fprintf(stderr, " mu1        : %g \n", mu1);
	fprintf(stderr, " mu2        : %g \n", mu2);
	fprintf(stderr, " mu1/mu2    : %g \n", mu1 / mu2);
	fprintf(stderr, " nu1        : %g \n", mu1 / rho1);
	fprintf(stderr, " nu2        : %g \n", mu2 / rho2);
	fprintf(stderr, " nu1/nu2    : %g \n", (mu1 / rho1) / (mu2 / rho2));
	fprintf(stderr, " sigma      : %g \n\n", f.sigma);
	fprintf(stderr, " Atwood     : %g \n", sim.At);
	fprintf(stderr, " Weber      : %g \n", sim.We);
	fprintf(stderr, " R          : %g \n", sim.Re);
	fprintf(stderr, " B          : %g \n", sim.B);
	fprintf(stderr, " S          : %g \n", sim.S);	
	fprintf(stderr, " D          : %g \n\n", sim.D);
}