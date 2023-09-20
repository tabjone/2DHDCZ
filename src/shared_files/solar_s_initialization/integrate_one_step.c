#include "solar_s_initialization.h"

void integrate_one_step(struct IntegrationVariables *bg, int i, bool updown)
{
    // false for upward integration, true for downward integration
    double k = get_k_value(bg->r[i]);
    double nabla_star = NABLA_AD + k;

    double dm_dr = 4 * M_PI * pow(bg->r[i],2) * bg->rho0[i];
    double dp_dr = - G * bg->m[i] /pow(bg->r[i],2) * bg->rho0[i];
    double dT_dr = nabla_star * bg->T0[i]/bg->p0[i] * dp_dr;
    double ds_dr = (nabla_star-NABLA_AD) / (bg->rho0[i] * bg->T0[i]) * dp_dr;

    if ((i+1) % bg->N == 0 && i != 0 && false)
    {
        bg->r = realloc(bg->r, sizeof(double)*(bg->N*bg->N_increment));
        bg->p0 = realloc(bg->p0, sizeof(double)*(bg->N*bg->N_increment));
        bg->T0 = realloc(bg->T0, sizeof(double)*(bg->N*bg->N_increment));
        bg->rho0 = realloc(bg->rho0, sizeof(double)*(bg->N*bg->N_increment));
        bg->s0 = realloc(bg->s0, sizeof(double)*(bg->N*bg->N_increment));
        bg->grad_s0 = realloc(bg->grad_s0, sizeof(double)*(bg->N*bg->N_increment));
        bg->m = realloc(bg->m, sizeof(double)*(bg->N*bg->N_increment));

        bg->N_increment = bg->N_increment + 1;
    }

    double dr1 = fabs(p_step * bg->m[i]/dm_dr);
    double dr2 = fabs(p_step * bg->p0[i]/dp_dr);
    double dr3 = fabs(p_step * bg->T0[i]/dT_dr);
    double dr4 = fabs(p_step * bg->s0[i]/ds_dr);

    double dr;
    if (dr1 < dr2 && dr1 < dr3 && dr1 < dr4) {
        dr = dr1;
    } else if (dr2 < dr1 && dr2 < dr3 && dr2 < dr4) {
        dr = dr2;
    } else if (dr3 < dr1 && dr3 < dr2 && dr3 < dr4) {
        dr = dr3;
    } else {
        dr = dr4;
    }
    if (updown)
    {
        dr = -dr;
    }
    
    bg->grad_s0[i] = ds_dr;
    
    bg->r[i+1] = bg->r[i] + dr;
    bg->m[i+1] = bg->m[i] + dm_dr * dr;
    bg->p0[i+1] = bg->p0[i] + dp_dr * dr;
    bg->T0[i+1] = bg->T0[i] + dT_dr * dr;
    bg->s0[i+1] = bg->s0[i] + ds_dr * dr;
    bg->rho0[i+1] = M_U * MU / K_B * bg->p0[i+1]/bg->T0[i+1]; //ideal gas law
}