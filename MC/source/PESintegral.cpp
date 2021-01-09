#include <tchem/gaussian.hpp>

#include "../include/basic.hpp"

at::Tensor compute_dV(const at::Tensor & q);

// Return <▽V>, <q▽V>, <p▽V>
std::tuple<at::Tensor, at::Tensor, at::Tensor> compute_dV_expectations(tchem::gaussian::Gaussian & gaussian) {
    at::Tensor miu = gaussian.miu();
    at::Tensor average_dV = miu.new_zeros(dimension);
    at::Tensor average_qdV = miu.new_zeros({dimension, dimension});
    at::Tensor average_pdV = miu.new_zeros({dimension, dimension});
    // Monte Carlo integration
    gaussian.rand_init();
    const size_t NSamples = 10000;
    for (size_t i = 0; i < NSamples; i++) {
        // First half is q, latter half is p
        at::Tensor rand_vec = gaussian.rand(generator);
        at::Tensor q = rand_vec.slice(0, 0, dimension),
                   p = rand_vec.slice(0, dimension);
        // ▽V depends only on q
        at::Tensor dV = compute_dV(q);
        at::Tensor qdV = q.outer(dV),
                   pdV = p.outer(dV);
        // accumulate
        average_dV  +=  dV;
        average_qdV += qdV;
        average_pdV += pdV;
    }
    average_dV  /= (double)NSamples;
    average_qdV /= (double)NSamples;
    average_pdV /= (double)NSamples;
    return std::make_tuple(average_dV, average_qdV, average_pdV);
}