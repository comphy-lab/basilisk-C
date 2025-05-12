/**
# Convert dump file to binary/ascii format

Convert dump file produced by dump() in Basilisk to binary/ascii format.
Use %.12le format for human readable floating point expression or
use %la format for exact floating point conversion so that the reconstructed
binary file is bit-wise correct.
Author: Yuxuan Liu
Email: yuxuan.liu@eng.ox.ac.uk
Date: 2023/05/21
*/

#include "common.h"

struct DumpHeader {
  double t;
  long len;
  int i, depth, npe, version;
  coord n;
};

void bin_to_ascii(const char *binFile, const char *asciiFile) {
    FILE *fp_bin = fopen(binFile, "rb");
    FILE *fp_ascii = fopen(asciiFile, "w");

    if (fp_bin == NULL) {
        perror(binFile);
        exit(1);
    }
  
    if (fp_ascii == NULL) {
        perror(asciiFile);
        exit(1);
    }
  
    struct DumpHeader header;
    fread(&header, sizeof(header), 1, fp_bin);

    // Read and write the header data
    fprintf(fp_ascii, "%.12le %ld %d %d %d %d\n", header.t, header.len, header.i, header.depth, header.npe, header.version);

    // Read and write the MPI dimensions
    fprintf(fp_ascii, "%.12le %.12le %.12le\n", header.n.x, header.n.y, header.n.z);

    // Read and write the scalar names
    fprintf(fp_ascii, "isleaf ");

    for (int i = 0; i < header.len; i++) {
        unsigned len;
        fread(&len, sizeof(len), 1, fp_bin);
        char *name = malloc((len + 1) * sizeof(char));
        fread(name, sizeof(char), len, fp_bin);
        name[len] = '\0';
        fprintf(fp_ascii, "%s ", name);
        free(name);
    }

    // Read and write the coordinates
    double o[4];
    fread(o, sizeof(double), 4, fp_bin);
    fprintf(fp_ascii, "%.12le %.12le %.12le %.12le\n", o[0], o[1], o[2], o[3]);
  
    unsigned flags;
    int cell_size = sizeof(unsigned) + header.len*sizeof(double);
    long pos = ftell(fp_bin);
    int n;
  
    while ((n = fread(&flags, sizeof(flags), 1, fp_bin)) == 1) {
        fseek(fp_bin, -sizeof(unsigned), SEEK_CUR);
        unsigned flags;
        fread(&flags, sizeof(flags), 1, fp_bin);
        fprintf(fp_ascii, "%u ", flags);
      
        double val;
        for (int i = 0; i < header.len; i++) {
            fread(&val, sizeof(val), 1, fp_bin);
            fprintf(fp_ascii, "%.12le ", val);
        }
        fprintf(fp_ascii, "\n");

        pos += cell_size;
        fseek(fp_bin, pos, SEEK_SET);
    }

    fclose(fp_bin);
    fclose(fp_ascii);
}

void ascii_to_bin(const char *asciiFile, const char *binFile) {
    FILE *fp_ascii = fopen(asciiFile, "r");
    FILE *fp_bin = fopen(binFile, "wb");

    if (fp_ascii == NULL) {
        perror(asciiFile);
        exit(1);
    }
  
    if (fp_bin == NULL) {
        perror(binFile);
        exit(1);
    }
  
    struct DumpHeader header;

    // Read the header data
    fscanf(fp_ascii, "%le %ld %d %d %d %d", &header.t, &header.len, &header.i, &header.depth, &header.npe, &header.version);

    // Read the MPI dimensions
    fscanf(fp_ascii, "%le %le %le", &header.n.x, &header.n.y, &header.n.z);

    // Write the header data
    fwrite(&header, sizeof(struct DumpHeader), 1, fp_bin);

    // Read and write the scalar names
    char name[256];
    fscanf(fp_ascii, "%s", name);

    for (int i = 0; i < header.len; i++) {
        char name[256]; // assuming a maximum length for the scalar names
        fscanf(fp_ascii, "%s", name);
        unsigned len = strlen(name);
        fwrite(&len, sizeof(len), 1, fp_bin);
        fwrite(name, sizeof(char), len, fp_bin);
        fprintf(stdout, "Read scalar name %s %d\n", name, len);
    }

    // Read and write the coordinates
    double o[4];
    fscanf(fp_ascii, "%le %le %le %le", &o[0], &o[1], &o[2], &o[3]);
    fwrite(o, sizeof(double), 4, fp_bin);
  
    unsigned flags;
    double val;

    // Notice the change in the loop condition here
    while (fscanf(fp_ascii, "%u", &flags) == 1) {
        fwrite(&flags, sizeof(flags), 1, fp_bin);
      
        for (int i = 0; i < header.len; i++) {
            fscanf(fp_ascii, "%le", &val);
            fwrite(&val, sizeof(val), 1, fp_bin);
        }
    }

    fclose(fp_ascii);
    fclose(fp_bin);
}

bool compare_bin_files(const char *file1, const char *file2) {
    FILE *fp1 = fopen(file1, "rb");
    FILE *fp2 = fopen(file2, "rb");

    if (fp1 == NULL) {
        perror(file1);
        return false;
    }

    if (fp2 == NULL) {
        fclose(fp1);
        perror(file2);
        return false;
    }

    bool match = true;
    unsigned char byte1, byte2;
    size_t bytesRead1, bytesRead2;

    do {
        bytesRead1 = fread(&byte1, sizeof(byte1), 1, fp1);
        bytesRead2 = fread(&byte2, sizeof(byte2), 1, fp2);

        if (bytesRead1 != bytesRead2 || byte1 != byte2) {
            match = false;
            break;
        }
    } while (bytesRead1 > 0 && bytesRead2 > 0);

    fclose(fp1);
    fclose(fp2);

    return match;
}

int main()
{
    bin_to_ascii("restart", "restart.dat");
    ascii_to_bin("restart.dat", "restart.o");
    bin_to_ascii("restart.o", "restart.o.dat");
    fprintf(stdout, compare_bin_files("restart.o", "restart.o") ? "Read binary file success!\n":"Read binary file failed\n");
    fprintf(stdout, compare_bin_files("restart", "restart.o") ? "Reconstruct binary file success!\n":"Reconstruct binary file failed\n");
    exit(0);
}