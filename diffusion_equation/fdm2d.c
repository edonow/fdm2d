# include  <stdio.h>
# include  <string.h>
# include  <stdlib.h>
# include  <math.h>
# include  <time.h>

double cal_distance(double x1, double y1, double x2, double y2);
void initial_condition(double *x, double *y, double *u, int nx, int ny, double width, double height);
void set_boundary(double *x, double *y, double *u, int nx, int ny, double width, double height);
void set_initialize(double *x, double *y, double *u, double *du, int nx, int ny,double width, double height, int div_width, int div_height);
void finalize(double *x, double *y, double *u, double *du);

// --------------------------------------------------------------
int main(int argc, char *argv[]){

    int i,j;
    double width = 100;
    double height = 100;
    int div_width = 50;
    int div_height = 50;
    int nx = div_width+1;
    int ny = div_height+1;
    double dx = width/div_width;
    double dy = height/div_height;

    // --------------------
    // > Diffusion coefficient
    double D = 1.0e-1;
    // > Time loop
    int loop = 10000;
    // CFL condition
    double cfl = 0.1;
    double dt = cfl * (width / div_width) / D;
    printf(" >> dt = %f\n",dt);
    // --------------------
    // file pointer 
    FILE *fp;
    char *fname;


    // Allocate memory for the grid
    double *x, *y, *u, *du;
    x = malloc(nx*ny*sizeof(double));
    y = malloc(nx*ny*sizeof(double));
    u = malloc(nx*ny*sizeof(double));
    du = malloc(nx*ny*sizeof(double));


    // setting x, y corrdinates and initial value of u
    set_initialize(x, y, u, du, nx, ny, width, height, div_width, div_height);
    initial_condition(x, y, u, nx, ny, width, height);

    system("rm -f ./output/*.csv");

    for(int l=1;l<=loop;l++){

        for(i=1; i<nx-1;i++){
            for(j=1;j<ny-1;j++){
                du[ny*i+j] = D * (u[ny*(i+1) + j] + u[ny*(i-1) + j] - 2 * u[ny*i + j]) / (dx * dx)
                    + D * (u[ny*i + j + 1] + u[ny*i + j - 1] - 2 * u[ny*i + j]) / (dy * dy);

            }
        }
        for(i=0;i<nx;i++){
            for(j=0;j<ny;j++){
                u[ny*i + j] = u[ny*i + j] + dt*du[ny*i + j];
            }
        }

        set_boundary(x, y, u, nx, ny, width, height);

        if(l==1 || l%100==0){
            fname = malloc(255*sizeof(char));
            sprintf(fname,"./output/%d.csv",l);
            fp = fopen(fname,"w");

            for(i=0;i<nx;i++){
                for(j=0;j<ny;j++){
                    fprintf(fp,"%f,%f,%f\n",x[ny*i+j],y[ny*i+j],u[ny*i+j]);
                }
            }
            fclose(fp);
        }
    }

    finalize(x, y, u, du);


    printf("\n    loading output file... \n\n");
    printf("\n    > moving on to python.\n\n");
    system("python view.py");

    printf("\n     > success!\n\n");

    return 0;
}
// --------------------------------------------------------------
void set_initialize(double *x, double *y, double *u, double *du, int nx, int ny,double width, double height, int div_width, int div_height){
    int i,j;
    double dx = width/div_width;
    double dy = height/div_height;
    for(i=0;i<nx;i++){
        for(j=0;j<ny+1;j++){
            x[ny*i+j] = i*dx;
            y[ny*i+j] = j*dy;
            u[ny*i+j] = 0.0;
            du[ny*i+j] = 0.0;
        }
    }
    return;
}
// --------------------------------------------------------------
void set_boundary(double *x, double *y, double *u, int nx, int ny, double width, double height){
    int i,j;
    // -- Dirichlet boundary condition --
    // for(i=0;i<nx;i++){
    //     u[ny*i + 0] = 0.0;
    //     u[ny*i + ny-1] = 0.0;
    // }
    // for(j=0;j<ny;j++){
    //     u[ny*0      + j] = 0.0;
    //     u[ny*(nx-1) + j] = 0.0;
    // }
    // return;

    //  -- Neumann boundary condition --
    for(i = 0; i < nx; i++) {
        u[ny*i + 0     ] = u[ny*i + 1];
        u[ny*i + ny - 1] = u[ny*i + ny - 2];
    }
    for(j = 0; j < ny; j++) {
        u[ny*0        + j] = u[ny*1 + j];
        u[ny*(nx - 1) + j] = u[ny*(nx - 2) + j];
    }

    return;
}
// --------------------------------------------------------------
double cal_distance(double x1, double y1, double x2, double y2){
    return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}
// --------------------------------------------------------------
void initial_condition(double *x, double *y, double *u, int nx, int ny, double width, double height) {
    double initial_val = 10.0;
    int i, j;
    double r = 5.0;
    double center_x = width / 2.0;
    double center_y = height / 2.0;
    for(i=0;i<nx;i++){
        for(j=0;j<ny;j++){
            if(cal_distance(center_x, center_y, x[ny*i+j], y[ny*i+j]) < r){
                u[ny*i+j] = initial_val;
            }
        }
    }
    return;

}
// --------------------------------------------------------------
void finalize(double *x, double *y, double *u, double *du){
    free(x);
    free(y);
    free(u);
    free(du);

    return;
}
