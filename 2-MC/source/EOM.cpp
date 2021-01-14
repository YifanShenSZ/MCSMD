#include <torch/torch.h>

#include <tchem/gaussian.hpp>

#include "../include/basic.hpp"

// Return <▽V>, <q▽V>, <p▽V>
std::tuple<at::Tensor, at::Tensor, at::Tensor> compute_dV_expectations(tchem::gaussian::Gaussian & gaussian);

// Based on current expectations,
// fit a gaussian phase space distribution to compute gradient
std::vector<at::Tensor> compute_gradients(const std::vector<at::Tensor> & expectations) {
    // Fit a gaussian from expectations
    tchem::gaussian::Gaussian g(expectations[1], expectations[2] - expectations[1].outer(expectations[1]));
    // Get ▽V related integrals
    at::Tensor dV, qdV, pdV;
    std::tie(dV, qdV, pdV) = compute_dV_expectations(g);
    // Moyal equation of motion
    std::vector<at::Tensor> gradients(3);
    gradients[1] = expectations[1].new_empty(expectations[1].sizes());
    at::Tensor  q = expectations[1].slice(0, 0, dimension),
                p = expectations[1].slice(0, dimension),
               dq =    gradients[1].slice(0, 0, dimension),
               dp =    gradients[1].slice(0, dimension);
    dq.copy_(p / mass);
    dp.copy_(-dV     );
    gradients[2] = expectations[2].new_empty(expectations[2].sizes());
    at::Tensor  qq = expectations[2].slice(0, 0, dimension).slice(1, 0, dimension),
                qp = expectations[2].slice(0, 0, dimension).slice(1, dimension),
                pp = expectations[2].slice(0, dimension).slice(1, dimension),
               dqq =    gradients[2].slice(0, 0, dimension).slice(1, 0, dimension),
               dqp =    gradients[2].slice(0, 0, dimension).slice(1, dimension),
               dpp =    gradients[2].slice(0, dimension).slice(1, dimension);
    dqq.copy_((qp + qp.transpose(0, 1)) / mass);
    dqp.copy_(pp / mass - qdV                 );
    dpp.copy_(-(pdV + pdV.transpose(0, 1))    );
    return gradients;
}

// Update expectations with RK4
void update_expectations(std::vector<at::Tensor> & expectations, const double & dt) {
    double dtd2 = dt / 2.0, dtd6 = dt / 6.0;
    std::vector<at::Tensor> intermediates(expectations.size());
    // gradient 1
    std::vector<at::Tensor> gradients1 = compute_gradients(expectations);
    // gradient 2
    for (size_t i = 1; i < expectations.size(); i++)
    intermediates[i] = expectations[i] + gradients1[i] * dtd2;
    std::vector<at::Tensor> gradients2 = compute_gradients(intermediates);
    // gradient 3
    for (size_t i = 1; i < expectations.size(); i++)
    intermediates[i] = expectations[i] + gradients2[i] * dtd2;
    std::vector<at::Tensor> gradients3 = compute_gradients(intermediates);
    // gradient 4
    for (size_t i = 1; i < expectations.size(); i++)
    intermediates[i] = expectations[i] + gradients3[i] * dt;
    std::vector<at::Tensor> gradients4 = compute_gradients(intermediates);
    // extrapolate
    for (size_t i = 1; i < expectations.size(); i++)
    expectations[i] += dtd6 * (gradients1[i] + 2.0 * gradients2[i] + 2.0 * gradients3[i] + gradients4[i]);
}