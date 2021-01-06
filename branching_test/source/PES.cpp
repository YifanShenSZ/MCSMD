#include <torch/torch.h>

double V(const double & q) {
    const double k_peak = 10.0, width = 1.0;
    double q_dimless = q / width;
    double V = k_peak * k_peak / 2.0 / 2000.0 
             * exp(-0.5 * q_dimless * q_dimless);
}