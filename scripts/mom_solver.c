#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.141592653589793

/* ========================================================================
   PROBLEM PARAMETERS
   ======================================================================== */
typedef struct {
    double eps1;        /* Permittivity above interface (z > 0) */
    double eps2;        /* Permittivity below interface (z < 0) */
    double V0;          /* Applied voltage on conductor */
    
    /* Conductor geometry (rectangular patch at z=0) */
    double patch_width;   /* Width in x-direction */
    double patch_height;  /* Height in y-direction */
    
    /* MoM discretization */
    int N_x;            /* Number of cells in x */
    int N_y;            /* Number of cells in y */
    int N_total;        /* Total cells = N_x * N_y */
    
    /* Field evaluation domain */
    double x_min, x_max;
    double y_min, y_max;
    double z_min, z_max;
    int N_field_x, N_field_y, N_field_z;
} Params;

/* ========================================================================
   CELL STRUCTURE
   ======================================================================== */
typedef struct {
    double x_center, y_center;  /* Cell center coordinates */
    double dx, dy;              /* Cell dimensions */
    double charge_density;      /* Solution: sigma */
} Cell;

/* ========================================================================
   GAUSSIAN QUADRATURE DATA (3-point)
   ======================================================================== */
#define NGAUSS 3
static const double gauss_points[NGAUSS] = {
    -0.7745966692414834,
     0.0,
     0.7745966692414834
};
static const double gauss_weights[NGAUSS] = {
    0.5555555555555556,
    0.8888888888888888,
    0.5555555555555556
};

/* ========================================================================
   MACRO FOR 2D ARRAY INDEXING
   ======================================================================== */
#define IDX(i, j, cols) ((i) * (cols) + (j))

/* ========================================================================
   FUNCTION PROTOTYPES
   ======================================================================== */
void initialize_cells(Cell *cells, const Params *p);
double green_function_interface(double x, double y, double z, 
                                 double xp, double yp, double zp, 
                                 const Params *p);
double integrate_green_over_cell(double x_obs, double y_obs, double z_obs,
                                 const Cell *cell, const Params *p);
double self_term_potential(double L, const Params *p);
void build_mom_matrix(double **Z, Cell *cells, const Params *p);
void build_rhs_vector(double *b, const Params *p);
int gaussian_elimination(double **A, double *b, double *x, int n);
double calculate_potential(double x, double y, double z, 
                          Cell *cells, const Params *p);
void calculate_electric_field(double x, double y, double z,
                             Cell *cells, const Params *p,
                             double *Ex, double *Ey, double *Ez);
void save_results(Cell *cells, const Params *p);
void save_field_data(Cell *cells, const Params *p);
double** allocate_matrix(int rows, int cols);
void free_matrix(double **matrix, int rows);
double compute_residual(double **Z, double *sigma, double *b, int n);

/* ========================================================================
   GREEN'S FUNCTION FOR DIELECTRIC INTERFACE (IMAGE CHARGE METHOD)
   ======================================================================== */
double green_function_interface(double x, double y, double z, 
                                 double xp, double yp, double zp, 
                                 const Params *p) {
    /* Point source at (xp, yp, zp)
       Observer at (x, y, z)
       
       For z > 0 (region 1):
       G = 1/(4*pi*eps1) * [1/R_direct + q_image/R_image]
       
       For z < 0 (region 2):
       G = 1/(4*pi*eps2) * [q_trans/R_direct]
       
       where:
       q_image = (eps1 - eps2)/(eps1 + eps2)
       q_trans = 2*eps1/(eps1 + eps2)
    */
    
    double dx = x - xp;
    double dy = y - yp;
    double dz = z - zp;
    double R_direct = sqrt(dx*dx + dy*dy + dz*dz);
    
    /* Avoid singularity */
    if (R_direct < 1e-12) R_direct = 1e-12;
    
    double G;
    
    if (z >= 0.0) {
        /* Region 1: above interface */
        /* Image source at (xp, yp, -zp) */
        double dz_image = z + zp;
        double R_image = sqrt(dx*dx + dy*dy + dz_image*dz_image);
        if (R_image < 1e-12) R_image = 1e-12;
        
        double q_image = (p->eps1 - p->eps2) / (p->eps1 + p->eps2);
        G = (1.0 / R_direct + q_image / R_image) / (4.0 * PI * p->eps1);
    } else {
        /* Region 2: below interface */
        double q_trans = 2.0 * p->eps1 / (p->eps1 + p->eps2);
        G = q_trans / R_direct / (4.0 * PI * p->eps2);
    }
    
    return G;
}

/* ========================================================================
   INTEGRATE GREEN'S FUNCTION OVER SOURCE CELL (GAUSSIAN QUADRATURE)
   ======================================================================== */
double integrate_green_over_cell(double x_obs, double y_obs, double z_obs,
                                 const Cell *cell, const Params *p) {
    /* Integrate G(x_obs, y_obs, z_obs | x', y', 0) over source cell
       using 3x3 Gaussian quadrature */
    
    double x_c = cell->x_center;
    double y_c = cell->y_center;
    double dx_half = cell->dx / 2.0;
    double dy_half = cell->dy / 2.0;
    
    double sum = 0.0;
    
    for (int i = 0; i < NGAUSS; i++) {
        for (int j = 0; j < NGAUSS; j++) {
            /* Map quadrature points from [-1,1] to cell domain */
            double xp = x_c + dx_half * gauss_points[i];
            double yp = y_c + dy_half * gauss_points[j];
            double zp = 0.0;  /* Source on interface */
            
            double G = green_function_interface(x_obs, y_obs, z_obs, 
                                               xp, yp, zp, p);
            sum += gauss_weights[i] * gauss_weights[j] * G;
        }
    }
    
    /* Jacobian for transformation: dx_half * dy_half */
    return sum * dx_half * dy_half;
}

/* ========================================================================
   ANALYTICAL SELF-TERM FOR MOM MATRIX
   ======================================================================== */
double self_term_potential(double L, const Params *p) {
    /* For a uniform charge density sigma on a square patch with side L
       at the interface (eps1, eps2), the potential at the center is:
       
       V(0,0,0) = 1/(eps1 + eps2) * (L/pi) * [ln(1+sqrt(2)) + sqrt(2)]
       
       In the MoM formulation with pulse basis functions:
       - We expand: sigma(r) = sum_j sigma_j * f_j(r)
       - Where f_j(r) = 1 inside cell j, 0 outside
       - The potential is: V(r) = sum_j sigma_j * integral[G(r,r') dS']
       
       The matrix element Z_ij represents:
       Z_ij = integral_over_cell_j [G(r_i, r') dS']
       
       For the self-term (i=j), this is the potential at the cell center
       due to unit charge density over that cell.
       
       The analytical formula above gives exactly this quantity.
       We do NOT multiply by L^2 again because:
       - The formula already accounts for integration over the L×L cell
       - Multiplying by L^2 would be double-counting the area
    */
    
    double sqrt2 = sqrt(2.0);
    double factor = log(1.0 + sqrt2) + sqrt2;
    double eps_sum = p->eps1 + p->eps2;
    
    /* This is Z_ii: potential at center due to unit charge density on L×L cell */
    double Z_ii = (L / PI) * factor / eps_sum;
    
    return Z_ii;
}

/* ========================================================================
   INITIALIZE MOM CELLS ON CONDUCTING PATCH
   ======================================================================== */
void initialize_cells(Cell *cells, const Params *p) {
    double dx = p->patch_width / p->N_x;
    double dy = p->patch_height / p->N_y;
    double x_start = -p->patch_width / 2.0;
    double y_start = -p->patch_height / 2.0;
    
    for (int i = 0; i < p->N_x; i++) {
        for (int j = 0; j < p->N_y; j++) {
            int idx = IDX(i, j, p->N_y);
            cells[idx].x_center = x_start + (i + 0.5) * dx;
            cells[idx].y_center = y_start + (j + 0.5) * dy;
            cells[idx].dx = dx;
            cells[idx].dy = dy;
            cells[idx].charge_density = 0.0;
        }
    }
}

/* ========================================================================
   BUILD MOM IMPEDANCE MATRIX (PULSE BASIS + POINT MATCHING)
   ======================================================================== */
void build_mom_matrix(double **Z, Cell *cells, const Params *p) {
    printf("Building MoM matrix...\n");
    
    for (int i = 0; i < p->N_total; i++) {
        if (i % 10 == 0 || i == p->N_total - 1) {
            printf("  Row %d/%d (%.1f%%)\r", i+1, p->N_total, 
                   100.0 * (i+1) / p->N_total);
            fflush(stdout);
        }
        
        /* Observation point at cell center */
        double xi = cells[i].x_center;
        double yi = cells[i].y_center;
        double zi = 0.0;  /* On conductor surface */
        
        for (int j = 0; j < p->N_total; j++) {
            if (i == j) {
                /* Self-term: use analytical formula */
                double L = cells[i].dx;  /* Assume square cells (dx = dy) */
                Z[i][j] = self_term_potential(L, p);
            } else {
                /* Mutual term: integrate Green's function over source cell j */
                Z[i][j] = integrate_green_over_cell(xi, yi, zi, &cells[j], p);
            }
        }
    }
    printf("\nMoM matrix built successfully.\n");
}

/* ========================================================================
   BUILD RIGHT-HAND SIDE VECTOR (APPLIED VOLTAGE)
   ======================================================================== */
void build_rhs_vector(double *b, const Params *p) {
    for (int i = 0; i < p->N_total; i++) {
        b[i] = p->V0;
    }
}

/* ========================================================================
   GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
   ======================================================================== */
int gaussian_elimination(double **A, double *b, double *x, int n) {
    printf("Solving linear system (%dx%d)...\n", n, n);
    
    /* Create augmented matrix */
    double **aug = allocate_matrix(n, n + 1);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aug[i][j] = A[i][j];
        }
        aug[i][n] = b[i];
    }
    
    /* Forward elimination with partial pivoting */
    for (int k = 0; k < n; k++) {
        if (k % 20 == 0 || k == n - 1) {
            printf("  Elimination step %d/%d (%.1f%%)\r", k+1, n, 
                   100.0 * (k+1) / n);
            fflush(stdout);
        }
        
        /* Find pivot */
        int max_row = k;
        double max_val = fabs(aug[k][k]);
        for (int i = k + 1; i < n; i++) {
            if (fabs(aug[i][k]) > max_val) {
                max_val = fabs(aug[i][k]);
                max_row = i;
            }
        }
        
        /* Swap rows */
        if (max_row != k) {
            double *temp = aug[k];
            aug[k] = aug[max_row];
            aug[max_row] = temp;
        }
        
        /* Check for singular matrix */
        if (fabs(aug[k][k]) < 1e-15) {
            fprintf(stderr, "\nError: Matrix is singular at row %d\n", k);
            free_matrix(aug, n);
            return -1;
        }
        
        /* Eliminate column */
        for (int i = k + 1; i < n; i++) {
            double factor = aug[i][k] / aug[k][k];
            for (int j = k; j <= n; j++) {
                aug[i][j] -= factor * aug[k][j];
            }
        }
    }
    
    /* Back substitution */
    for (int i = n - 1; i >= 0; i--) {
        x[i] = aug[i][n];
        for (int j = i + 1; j < n; j++) {
            x[i] -= aug[i][j] * x[j];
        }
        x[i] /= aug[i][i];
    }
    
    printf("\nLinear system solved.\n");
    free_matrix(aug, n);
    return 0;
}

/* ========================================================================
   COMPUTE RESIDUAL NORM
   ======================================================================== */
double compute_residual(double **Z, double *sigma, double *b, int n) {
    double residual = 0.0;
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += Z[i][j] * sigma[j];
        }
        double diff = sum - b[i];
        residual += diff * diff;
    }
    return sqrt(residual / n);
}

/* ========================================================================
   CALCULATE POTENTIAL AT ANY POINT
   ======================================================================== */
double calculate_potential(double x, double y, double z, 
                          Cell *cells, const Params *p) {
    double V = 0.0;
    
    for (int i = 0; i < p->N_total; i++) {
        double xp = cells[i].x_center;
        double yp = cells[i].y_center;
        double zp = 0.0;
        double sigma = cells[i].charge_density;
        double dS = cells[i].dx * cells[i].dy;
        
        /* V = sum of sigma * G * dS */
        double G = green_function_interface(x, y, z, xp, yp, zp, p);
        V += sigma * G * dS;
    }
    
    return V;
}

/* ========================================================================
   CALCULATE ELECTRIC FIELD E = -grad(V)
   ======================================================================== */
void calculate_electric_field(double x, double y, double z,
                             Cell *cells, const Params *p,
                             double *Ex, double *Ey, double *Ez) {
    double h = 1e-6;  /* Small step for numerical derivative */
    
    double V_xp = calculate_potential(x + h, y, z, cells, p);
    double V_xm = calculate_potential(x - h, y, z, cells, p);
    double V_yp = calculate_potential(x, y + h, z, cells, p);
    double V_ym = calculate_potential(x, y - h, z, cells, p);
    double V_zp = calculate_potential(x, y, z + h, cells, p);
    double V_zm = calculate_potential(x, y, z - h, cells, p);
    
    *Ex = -(V_xp - V_xm) / (2.0 * h);
    *Ey = -(V_yp - V_ym) / (2.0 * h);
    *Ez = -(V_zp - V_zm) / (2.0 * h);
}

/* ========================================================================
   SAVE CHARGE DISTRIBUTION AND STATISTICS
   ======================================================================== */
void save_results(Cell *cells, const Params *p) {
    FILE *fp = fopen("charge_distribution.csv", "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open charge_distribution.csv\n");
        return;
    }
    
    fprintf(fp, "x,y,sigma\n");
    
    double total_charge = 0.0;
    double max_sigma = -1e100;
    double min_sigma = 1e100;
    
    for (int i = 0; i < p->N_total; i++) {
        double x = cells[i].x_center;
        double y = cells[i].y_center;
        double sigma = cells[i].charge_density;
        double dS = cells[i].dx * cells[i].dy;
        
        total_charge += sigma * dS;
        
        if (sigma > max_sigma) max_sigma = sigma;
        if (sigma < min_sigma) min_sigma = sigma;
        
        fprintf(fp, "%.8f,%.8f,%.8e\n", x, y, sigma);
    }
    
    fclose(fp);
    
    /* Calculate capacitance */
    double capacitance = total_charge / p->V0;
    
    printf("\n=== RESULTS ===\n");
    printf("Total charge: %.6e C\n", total_charge);
    printf("Capacitance (MoM): %.6e F\n", capacitance);
    printf("Max charge density: %.6e C/m^2\n", max_sigma);
    printf("Min charge density: %.6e C/m^2\n", min_sigma);
    printf("Charge density ratio (max/min): %.2f\n", max_sigma / min_sigma);
    printf("\nCharge distribution saved to: charge_distribution.csv\n");
}

/* ========================================================================
   SAVE POTENTIAL AND FIELD DATA ON GRID
   ======================================================================== */
void save_field_data(Cell *cells, const Params *p) {
    printf("\nCalculating field data...\n");
    
    FILE *fp = fopen("field_data.csv", "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open field_data.csv\n");
        return;
    }
    
    fprintf(fp, "x,y,z,V,Ex,Ey,Ez,E_mag\n");
    
    double dx = (p->x_max - p->x_min) / (p->N_field_x - 1);
    double dy = (p->y_max - p->y_min) / (p->N_field_y - 1);
    double dz = (p->z_max - p->z_min) / (p->N_field_z - 1);
    
    int count = 0;
    int total = p->N_field_x * p->N_field_y * p->N_field_z;
    
    for (int i = 0; i < p->N_field_x; i++) {
        for (int j = 0; j < p->N_field_y; j++) {
            for (int k = 0; k < p->N_field_z; k++) {
                double x = p->x_min + i * dx;
                double y = p->y_min + j * dy;
                double z = p->z_min + k * dz;
                
                /* Skip points on the conductor */
                if (fabs(z) < 1e-8 && fabs(x) <= p->patch_width/2.0 && 
                    fabs(y) <= p->patch_height/2.0) {
                    continue;
                }
                
                double V = calculate_potential(x, y, z, cells, p);
                double Ex, Ey, Ez;
                calculate_electric_field(x, y, z, cells, p, &Ex, &Ey, &Ez);
                double E_mag = sqrt(Ex*Ex + Ey*Ey + Ez*Ez);
                
                fprintf(fp, "%.6f,%.6f,%.6f,%.8e,%.8e,%.8e,%.8e,%.8e\n",
                       x, y, z, V, Ex, Ey, Ez, E_mag);
                
                count++;
                if (count % 100 == 0 || count == total) {
                    printf("  Progress: %d/%d points (%.1f%%)\r", 
                           count, total, 100.0 * count / total);
                    fflush(stdout);
                }
            }
        }
    }
    
    fclose(fp);
    printf("\nField data saved to: field_data.csv\n");
}

/* ========================================================================
   MEMORY ALLOCATION HELPERS
   ======================================================================== */
double** allocate_matrix(int rows, int cols) {
    double **matrix = (double**)malloc(rows * sizeof(double*));
    if (!matrix) {
        fprintf(stderr, "Error: Cannot allocate matrix rows\n");
        exit(1);
    }
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double*)malloc(cols * sizeof(double));
        if (!matrix[i]) {
            fprintf(stderr, "Error: Cannot allocate matrix row %d\n", i);
            exit(1);
        }
    }
    return matrix;
}

void free_matrix(double **matrix, int rows) {
    for (int i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

/* ========================================================================
   MAIN PROGRAM
   ======================================================================== */
int main(void) {
    Params p;
    
    /* Set problem parameters */
    p.eps1 = 8.854e-12;        /* Permittivity of air (F/m) */
    p.eps2 = 4.0 * 8.854e-12;  /* Permittivity of dielectric (4*eps0) */
    p.V0 = 1.0;                /* Applied voltage (V) */
    
    /* Conductor geometry */
    p.patch_width = 0.1;       /* 10 cm square patch */
    p.patch_height = 0.1;
    
    /* MoM discretization */
    p.N_x = 20;                /* 20x20 = 400 unknowns */
    p.N_y = 20;
    p.N_total = p.N_x * p.N_y;
    
    /* Field evaluation domain */
    p.x_min = -0.15;
    p.x_max = 0.15;
    p.y_min = -0.15;
    p.y_max = 0.15;
    p.z_min = -0.1;
    p.z_max = 0.1;
    p.N_field_x = 31;
    p.N_field_y = 31;
    p.N_field_z = 21;
    
    printf("=== Method of Moments: Conducting Patch on Dielectric Interface ===\n\n");
    printf("Problem setup:\n");
    printf("  Permittivity (z>0): %.3e F/m\n", p.eps1);
    printf("  Permittivity (z<0): %.3e F/m\n", p.eps2);
    printf("  Dielectric contrast: eps2/eps1 = %.2f\n", p.eps2/p.eps1);
    printf("  Patch size: %.3f x %.3f m\n", p.patch_width, p.patch_height);
    printf("  Applied voltage: %.2f V\n", p.V0);
    printf("  MoM cells: %d x %d = %d unknowns\n\n", p.N_x, p.N_y, p.N_total);
    
    /* Allocate cells */
    Cell *cells = (Cell*)malloc(p.N_total * sizeof(Cell));
    if (!cells) {
        fprintf(stderr, "Error: Cannot allocate cells\n");
        return 1;
    }
    initialize_cells(cells, &p);
    
    /* Allocate MoM matrix and vectors */
    double **Z = allocate_matrix(p.N_total, p.N_total);
    double *b = (double*)malloc(p.N_total * sizeof(double));
    double *sigma = (double*)malloc(p.N_total * sizeof(double));
    
    if (!b || !sigma) {
        fprintf(stderr, "Error: Cannot allocate vectors\n");
        return 1;
    }
    
    /* Build and solve MoM system */
    build_mom_matrix(Z, cells, &p);
    build_rhs_vector(b, &p);
    
    if (gaussian_elimination(Z, b, sigma, p.N_total) != 0) {
        fprintf(stderr, "Error: Failed to solve linear system\n");
        return 1;
    }
    
    /* Store solution in cells */
    for (int i = 0; i < p.N_total; i++) {
        cells[i].charge_density = sigma[i];
    }
    
    /* Compute residual */
    double residual = compute_residual(Z, sigma, b, p.N_total);
    printf("Residual norm: %.6e\n", residual);
    
    /* Save results */
    save_results(cells, &p);
    save_field_data(cells, &p);
    
    /* Cleanup */
    free_matrix(Z, p.N_total);
    free(b);
    free(sigma);
    free(cells);
    
    printf("\n=== Simulation complete! ===\n");
    printf("Generated files:\n");
    printf("  - charge_distribution.csv\n");
    printf("  - field_data.csv\n\n");
    
    return 0;
}
