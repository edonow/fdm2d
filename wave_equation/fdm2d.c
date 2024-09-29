#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void set_initialize(double *x, double *y, double *u, double *v, double *du, int nx, int ny, double width, double height);
void update_wave(double *u, double *v, double *du, double c, double dt, double dx, double dy, int nx, int ny);

int main() {
    int nx = 50;
    int ny = 50;
    double width = 100.0;
    double height = 100.0;
    double dx = width / (nx - 1);
    double dy = height / (ny - 1);
    double c = 1.0;
    double dt = 0.1;
    int time_steps = 1000;

    double *x = (double *)malloc(nx * ny * sizeof(double));
    double *y = (double *)malloc(nx * ny * sizeof(double));
    double *u = (double *)malloc(nx * ny * sizeof(double));
    double *v = (double *)malloc(nx * ny * sizeof(double));
    double *du = (double *)malloc(nx * ny * sizeof(double));

    FILE *fp;
    char *fname;

    set_initialize(x, y, u, v, du, nx, ny, width, height);

    
    system("rm -rf ./output/*.csv");

    for (int t = 0; t <= time_steps; t++) {
        update_wave(u, v, du, c, dt, dx, dy, nx, ny);

        // -- output --
        if(t%50==0){
            fname = malloc(255*sizeof(char));
            sprintf(fname,"./output/%d.csv",t);
            fp = fopen(fname,"w");

            for(int i=0;i<nx;i++){
                for(int j=0;j<ny;j++){
                    fprintf(fp,"%f,%f,%f\n",x[ny*i+j],y[ny*i+j],u[ny*i+j]);
                }
            }
            fclose(fp);
        }
    }

    free(x);
    free(y);
    free(u);
    free(v);
    free(du);

    printf("\n    loading output file... \n\n");
    printf("\n    > moving on to python.\n\n");
    system("python view.py");

    printf("\n     > success!\n\n");
    
    return 0;
}

// --------------------------------------------------------------
double calc_distance(double x1, double y1, double x2, double y2) {
    return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}
void set_initialize(double *x, double *y, double *u, double *v, double *du, int nx, int ny, double width, double height) {
    double center_x = width / 2.0;
    double center_y = height / 2.0;
    double radius = 5.0;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            x[ny * i + j] = i * (width / (nx - 1));
            y[ny * i + j] = j * (height / (ny - 1));
            if (calc_distance(center_x,center_y,x[ny * i + j],y[ny * i + j]) < radius) {
                u[ny * i + j] = 5.0;
            } else {
                u[ny * i + j] = 0.0;
            }
            v[ny * i + j] = 0.0;
            du[ny * i + j] = 0.0;            
        }
    }
}
// --------------------------------------------------------------
void update_wave(double *u, double *v, double *du, double c, double dt, double dx, double dy, int nx, int ny) {
    // -- update velocity --
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            double dudx2 = (u[ny * (i + 1) + j] - 2 * u[ny * i + j] + u[ny * (i - 1) + j]) / (dx * dx);
            double dudy2 = (u[ny * i + (j + 1)] - 2 * u[ny * i + j] + u[ny * i + (j - 1)]) / (dy * dy);
            v[ny * i + j] += dt * c * c * (dudx2 + dudy2); //update velocity
            du[ny * i + j] = u[ny * i + j] + dt * v[ny * i + j]; //update displacement
        }
    }
    // -- set boundary condition --
    for (int i = 0; i < nx; i++) {
        du[ny * i + 0] = 0.0;
        du[ny * i + (ny - 1)] = 0.0;
    }
    for (int j = 0; j < ny; j++) {
        du[ny * 0 + j] = 0.0;
        du[ny * (nx - 1) + j] = 0.0;
    }

    //  -- Neumann boundary condition --
    // for(int i = 0; i < nx; i++) {
    //     u[ny*i + 0     ] = u[ny*i + 1];
    //     u[ny*i + ny - 1] = u[ny*i + ny - 2];
    // }
    // for(int j = 0; j < ny; j++) {
    //     u[ny*0        + j] = u[ny*1 + j];
    //     u[ny*(nx - 1) + j] = u[ny*(nx - 2) + j];
    // }


    // -- update displacement --
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            u[ny * i + j] = du[ny * i + j];
        }
    }
}
