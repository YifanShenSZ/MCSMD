#include <torch/torch.h>

#include <tchem/gaussian.hpp>

#include "../include/basic.hpp"

// Return <▽V>, <q▽V>, <p▽V>
std::tuple<at::Tensor, at::Tensor, at::Tensor> compute_dV_expectations(tchem::gaussian::Gaussian & gaussian);

// Based on current expectations,
// fit a gaussian phase space distribution to compute gradient
// Update expectations with simple Euler method
void update_expectations(std::vector<at::Tensor> & expectations, const double & dt) {
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
    // Propagate by Euler method
    expectations[1] += dt * gradients[1];
    expectations[2] += dt * gradients[2];
}