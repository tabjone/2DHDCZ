#ifndef STRUCTS_H__
#define STRUCTS_H__

struct BackgroundVariables
{
    double *r_over_R;
    double *c_s;
    double *p0;
    double *T0;
    double *rho0;
    double *Gamma_1;
    double *H;
    double *superad_param;
    double *grad_s0;
    double *g;
};

struct ForegroundVariables
{
    double **p1;
    double **rho1;
    double **T1;
    double **s1;
    double **vx;
    double **vz;
};

void allocate_background_struct(int nz, struct BackgroundVariables *background_variables);
void allocate_foreground_struct(int nz, int nx, struct ForegroundVariables *foreground_variables);

void deallocate_background_struct(struct BackgroundVariables *background_variables);
void deallocate_foreground_struct(struct ForegroundVariables *foreground_variables);


#endif // STRUCTS_H__