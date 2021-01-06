#include <tchem/gaussian.hpp>

#include "../include/basic.hpp"

// Return <▽V>, <q▽V>, <p▽V>
std::tuple<at::Tensor, at::Tensor, at::Tensor> compute_dV_expectations(tchem::gaussian::Gaussian & gaussian) {
    at::Tensor miu = gaussian.miu();
    at::Tensor average_dV = miu.new_zeros(dimension);
    at::Tensor average_qdV = miu.new_zeros({dimension, dimension});
    at::Tensor average_pdV = miu.new_zeros({dimension, dimension});

    tchem::polynomial::PolynomialSet var_set(2, 3);
    // 0, <1>
    // 1, <q>; 2, <p>
    // 3, <qq>; 4, <pq>; 5, <pp>
    // 6, <qqq>; 7, <pqq>
    at::Tensor integrals = gaussian.integral(var_set);
    // <▽V> = <q> + 0.5 <qq>
    average_dV [0]    = integrals[1] + 0.5 * integrals[3];
    // <q▽V> = <qq> + 0.5 <qqq>
    average_qdV[0][0] = integrals[3] + 0.5 * integrals[6];
    // <p▽V> = <pq> + 0.5 <pqq>
    average_pdV[0][0] = integrals[4] + 0.5 * integrals[7];

    return std::make_tuple(average_dV, average_qdV, average_pdV);
}