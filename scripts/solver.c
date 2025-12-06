#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.141592653589793

/* Parámetros del problema */
typedef struct {
    double eps1;        /* Permitividad sobre la interfaz (z > 0) */
    double eps2;        /* Permitividad bajo la interfaz (z < 0) */
    double V0;          /* Voltaje aplicado en el conductor */
    
    double ancho_parche;   /* Ancho en dirección x */
    double alto_parche;    /* Alto en dirección y */
    
    int N_x;            /* Número de celdas en x */
    int N_y;            /* Número de celdas en y */
    int N_total;        /* Celdas totales = N_x * N_y */
    
    double x_min, x_max;
    double y_min, y_max;
    double z_min, z_max;
    int N_campo_x, N_campo_y, N_campo_z;
} Parametros;

/* Estructura de celda */
typedef struct {
    double x_centro, y_centro;
    double dx, dy;
    double densidad_carga;  /* sigma */
} Celda;

/* Cuadratura Gaussiana (3 puntos) */
#define NGAUSS 3
static const double puntos_gauss[NGAUSS] = {
    -0.7745966692414834, 0.0, 0.7745966692414834
};
static const double pesos_gauss[NGAUSS] = {
    0.5555555555555556, 0.8888888888888888, 0.5555555555555556
};

#define IDX(i, j, cols) ((i) * (cols) + (j))

/* Prototipos */
void inicializar_celdas(Celda *celdas, const Parametros *p);
double funcion_green_interfaz(double x, double y, double z, 
                              double xp, double yp, double zp, 
                              const Parametros *p);
double integrar_green_sobre_celda(double x_obs, double y_obs, double z_obs,
                                  const Celda *celda, const Parametros *p);
double termino_auto_potencial(double L, const Parametros *p);
void construir_matriz_mom(double **Z, Celda *celdas, const Parametros *p);
void construir_vector_rhs(double *b, const Parametros *p);
int eliminacion_gaussiana(double **A, double *b, double *x, int n);
double calcular_potencial(double x, double y, double z, 
                         Celda *celdas, const Parametros *p);
void calcular_campo_electrico(double x, double y, double z,
                             Celda *celdas, const Parametros *p,
                             double *Ex, double *Ey, double *Ez);
void guardar_resultados(Celda *celdas, const Parametros *p);
void guardar_datos_campo(Celda *celdas, const Parametros *p);
double** asignar_matriz(int filas, int cols);
void liberar_matriz(double **matriz, int filas);
double calcular_residual(double **Z, double *sigma, double *b, int n);

/* Función de Green para interfaz dieléctrica (método de carga imagen) */
double funcion_green_interfaz(double x, double y, double z, 
                              double xp, double yp, double zp, 
                              const Parametros *p) {
    double dx = x - xp;
    double dy = y - yp;
    double dz = z - zp;
    double R_directo = sqrt(dx*dx + dy*dy + dz*dz);
    
    if (R_directo < 1e-12) R_directo = 1e-12;
    
    double G;
    
    if (z >= 0.0) {
        /* Región 1: sobre la interfaz */
        double dz_imagen = z + zp;
        double R_imagen = sqrt(dx*dx + dy*dy + dz_imagen*dz_imagen);
        if (R_imagen < 1e-12) R_imagen = 1e-12;
        
        double q_imagen = (p->eps1 - p->eps2) / (p->eps1 + p->eps2);
        G = (1.0 / R_directo + q_imagen / R_imagen) / (4.0 * PI * p->eps1);
    } else {
        /* Región 2: bajo la interfaz */
        double q_trans = 2.0 * p->eps1 / (p->eps1 + p->eps2);
        G = q_trans / R_directo / (4.0 * PI * p->eps2);
    }
    
    return G;
}

/* Integrar función de Green sobre celda fuente */
double integrar_green_sobre_celda(double x_obs, double y_obs, double z_obs,
                                  const Celda *celda, const Parametros *p) {
    double x_c = celda->x_centro;
    double y_c = celda->y_centro;
    double dx_mitad = celda->dx / 2.0;
    double dy_mitad = celda->dy / 2.0;
    
    double suma = 0.0;
    
    for (int i = 0; i < NGAUSS; i++) {
        for (int j = 0; j < NGAUSS; j++) {
            double xp = x_c + dx_mitad * puntos_gauss[i];
            double yp = y_c + dy_mitad * puntos_gauss[j];
            double zp = 0.0;
            
            double G = funcion_green_interfaz(x_obs, y_obs, z_obs, 
                                             xp, yp, zp, p);
            suma += pesos_gauss[i] * pesos_gauss[j] * G;
        }
    }
    
    return suma * dx_mitad * dy_mitad;
}

/* Término auto-analítico para matriz MoM */
double termino_auto_potencial(double L, const Parametros *p) {
    double sqrt2 = sqrt(2.0);
    double factor = log(1.0 + sqrt2) + sqrt2;
    double eps_suma = p->eps1 + p->eps2;
    
    double Z_ii = (L / PI) * factor / eps_suma;
    
    return Z_ii;
}

/* Inicializar celdas MoM en el parche conductor */
void inicializar_celdas(Celda *celdas, const Parametros *p) {
    double dx = p->ancho_parche / p->N_x;
    double dy = p->alto_parche / p->N_y;
    double x_inicio = -p->ancho_parche / 2.0;
    double y_inicio = -p->alto_parche / 2.0;
    
    for (int i = 0; i < p->N_x; i++) {
        for (int j = 0; j < p->N_y; j++) {
            int idx = IDX(i, j, p->N_y);
            celdas[idx].x_centro = x_inicio + (i + 0.5) * dx;
            celdas[idx].y_centro = y_inicio + (j + 0.5) * dy;
            celdas[idx].dx = dx;
            celdas[idx].dy = dy;
            celdas[idx].densidad_carga = 0.0;
        }
    }
}

/* Construir matriz de impedancia MoM */
void construir_matriz_mom(double **Z, Celda *celdas, const Parametros *p) {
    for (int i = 0; i < p->N_total; i++) {
        double xi = celdas[i].x_centro;
        double yi = celdas[i].y_centro;
        double zi = 0.0;
        
        for (int j = 0; j < p->N_total; j++) {
            if (i == j) {
                double L = celdas[i].dx;
                Z[i][j] = termino_auto_potencial(L, p);
            } else {
                Z[i][j] = integrar_green_sobre_celda(xi, yi, zi, &celdas[j], p);
            }
        }
    }
}

/* Construir vector del lado derecho */
void construir_vector_rhs(double *b, const Parametros *p) {
    for (int i = 0; i < p->N_total; i++) {
        b[i] = p->V0;
    }
}

/* Eliminación Gaussiana con pivoteo parcial */
int eliminacion_gaussiana(double **A, double *b, double *x, int n) {
    double **aug = asignar_matriz(n, n + 1);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aug[i][j] = A[i][j];
        }
        aug[i][n] = b[i];
    }
    
    /* Eliminación hacia adelante con pivoteo parcial */
    for (int k = 0; k < n; k++) {
        int fila_max = k;
        double val_max = fabs(aug[k][k]);
        for (int i = k + 1; i < n; i++) {
            if (fabs(aug[i][k]) > val_max) {
                val_max = fabs(aug[i][k]);
                fila_max = i;
            }
        }
        
        if (fila_max != k) {
            double *temp = aug[k];
            aug[k] = aug[fila_max];
            aug[fila_max] = temp;
        }
        
        if (fabs(aug[k][k]) < 1e-15) {
            liberar_matriz(aug, n);
            return -1;
        }
        
        for (int i = k + 1; i < n; i++) {
            double factor = aug[i][k] / aug[k][k];
            for (int j = k; j <= n; j++) {
                aug[i][j] -= factor * aug[k][j];
            }
        }
    }
    
    /* Sustitución hacia atrás */
    for (int i = n - 1; i >= 0; i--) {
        x[i] = aug[i][n];
        for (int j = i + 1; j < n; j++) {
            x[i] -= aug[i][j] * x[j];
        }
        x[i] /= aug[i][i];
    }
    
    liberar_matriz(aug, n);
    return 0;
}

/* Calcular norma del residual */
double calcular_residual(double **Z, double *sigma, double *b, int n) {
    double residual = 0.0;
    for (int i = 0; i < n; i++) {
        double suma = 0.0;
        for (int j = 0; j < n; j++) {
            suma += Z[i][j] * sigma[j];
        }
        double dif = suma - b[i];
        residual += dif * dif;
    }
    return sqrt(residual / n);
}

/* Calcular potencial en cualquier punto */
double calcular_potencial(double x, double y, double z, 
                         Celda *celdas, const Parametros *p) {
    double V = 0.0;
    
    for (int i = 0; i < p->N_total; i++) {
        double xp = celdas[i].x_centro;
        double yp = celdas[i].y_centro;
        double zp = 0.0;
        double sigma = celdas[i].densidad_carga;
        double dS = celdas[i].dx * celdas[i].dy;
        
        double G = funcion_green_interfaz(x, y, z, xp, yp, zp, p);
        V += sigma * G * dS;
    }
    
    return V;
}

/* Calcular campo eléctrico E = -grad(V) */
void calcular_campo_electrico(double x, double y, double z,
                             Celda *celdas, const Parametros *p,
                             double *Ex, double *Ey, double *Ez) {
    double h = 1e-6;
    
    double V_xp = calcular_potencial(x + h, y, z, celdas, p);
    double V_xm = calcular_potencial(x - h, y, z, celdas, p);
    double V_yp = calcular_potencial(x, y + h, z, celdas, p);
    double V_ym = calcular_potencial(x, y - h, z, celdas, p);
    double V_zp = calcular_potencial(x, y, z + h, celdas, p);
    double V_zm = calcular_potencial(x, y, z - h, celdas, p);
    
    *Ex = -(V_xp - V_xm) / (2.0 * h);
    *Ey = -(V_yp - V_ym) / (2.0 * h);
    *Ez = -(V_zp - V_zm) / (2.0 * h);
}

/* Guardar distribución de carga y estadísticas */
void guardar_resultados(Celda *celdas, const Parametros *p) {
    FILE *fp = fopen("distribucion_carga.csv", "w");
    if (!fp) return;
    
    fprintf(fp, "x,y,sigma\n");
    
    double carga_total = 0.0;
    double sigma_max = -1e100;
    double sigma_min = 1e100;
    
    for (int i = 0; i < p->N_total; i++) {
        double x = celdas[i].x_centro;
        double y = celdas[i].y_centro;
        double sigma = celdas[i].densidad_carga;
        double dS = celdas[i].dx * celdas[i].dy;
        
        carga_total += sigma * dS;
        
        if (sigma > sigma_max) sigma_max = sigma;
        if (sigma < sigma_min) sigma_min = sigma;
        
        fprintf(fp, "%.8f,%.8f,%.8e\n", x, y, sigma);
    }
    
    fclose(fp);
    
    double capacitancia = carga_total / p->V0;
    
    printf("Carga total: %.6e C\n", carga_total);
    printf("Capacitancia: %.6e F\n", capacitancia);
    printf("Densidad max: %.6e C/m^2\n", sigma_max);
    printf("Densidad min: %.6e C/m^2\n", sigma_min);
}

/* Guardar datos de potencial y campo en malla */
void guardar_datos_campo(Celda *celdas, const Parametros *p) {
    FILE *fp = fopen("datos_campo.csv", "w");
    if (!fp) return;
    
    fprintf(fp, "x,y,z,V,Ex,Ey,Ez,E_mag\n");
    
    double dx = (p->x_max - p->x_min) / (p->N_campo_x - 1);
    double dy = (p->y_max - p->y_min) / (p->N_campo_y - 1);
    double dz = (p->z_max - p->z_min) / (p->N_campo_z - 1);
    
    for (int i = 0; i < p->N_campo_x; i++) {
        for (int j = 0; j < p->N_campo_y; j++) {
            for (int k = 0; k < p->N_campo_z; k++) {
                double x = p->x_min + i * dx;
                double y = p->y_min + j * dy;
                double z = p->z_min + k * dz;
                
                if (fabs(z) < 1e-8 && fabs(x) <= p->ancho_parche/2.0 && 
                    fabs(y) <= p->alto_parche/2.0) {
                    continue;
                }
                
                double V = calcular_potencial(x, y, z, celdas, p);
                double Ex, Ey, Ez;
                calcular_campo_electrico(x, y, z, celdas, p, &Ex, &Ey, &Ez);
                double E_mag = sqrt(Ex*Ex + Ey*Ey + Ez*Ez);
                
                fprintf(fp, "%.6f,%.6f,%.6f,%.8e,%.8e,%.8e,%.8e,%.8e\n",
                       x, y, z, V, Ex, Ey, Ez, E_mag);
            }
        }
    }
    
    fclose(fp);
}

/* Asignación de memoria */
double** asignar_matriz(int filas, int cols) {
    double **matriz = (double**)malloc(filas * sizeof(double*));
    if (!matriz) exit(1);
    for (int i = 0; i < filas; i++) {
        matriz[i] = (double*)malloc(cols * sizeof(double));
        if (!matriz[i]) exit(1);
    }
    return matriz;
}

void liberar_matriz(double **matriz, int filas) {
    for (int i = 0; i < filas; i++) {
        free(matriz[i]);
    }
    free(matriz);
}

/* Programa principal */
int main(void) {
    Parametros p;
    
    p.eps1 = 8.854e-12;
    p.eps2 = 4.0 * 8.854e-12;
    p.V0 = 1.0;
    
    p.ancho_parche = 0.1;
    p.alto_parche = 0.1;
    
    p.N_x = 20;
    p.N_y = 20;
    p.N_total = p.N_x * p.N_y;
    
    p.x_min = -0.15;
    p.x_max = 0.15;
    p.y_min = -0.15;
    p.y_max = 0.15;
    p.z_min = -0.1;
    p.z_max = 0.1;
    p.N_campo_x = 31;
    p.N_campo_y = 31;
    p.N_campo_z = 21;
    
    Celda *celdas = (Celda*)malloc(p.N_total * sizeof(Celda));
    if (!celdas) return 1;
    inicializar_celdas(celdas, &p);
    
    double **Z = asignar_matriz(p.N_total, p.N_total);
    double *b = (double*)malloc(p.N_total * sizeof(double));
    double *sigma = (double*)malloc(p.N_total * sizeof(double));
    
    if (!b || !sigma) return 1;
    
    construir_matriz_mom(Z, celdas, &p);
    construir_vector_rhs(b, &p);
    
    if (eliminacion_gaussiana(Z, b, sigma, p.N_total) != 0) {
        return 1;
    }
    
    for (int i = 0; i < p.N_total; i++) {
        celdas[i].densidad_carga = sigma[i];
    }
    
    double residual = calcular_residual(Z, sigma, b, p.N_total);
    printf("Residual: %.6e\n", residual);
    
    guardar_resultados(celdas, &p);
    guardar_datos_campo(celdas, &p);
    
    liberar_matriz(Z, p.N_total);
    free(b);
    free(sigma);
    free(celdas);
    
    return 0;
}
