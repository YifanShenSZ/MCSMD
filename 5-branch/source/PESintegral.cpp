#include <tchem/linalg.hpp>
#include <tchem/gaussian.hpp>

#include "../include/basic.hpp"

at::Tensor compute_dV(const at::Tensor & q);

// Return <▽V>, <q▽V>, <p▽V>
std::tuple<at::Tensor, at::Tensor, at::Tensor> compute_dV_2expectations(tchem::gaussian::Gaussian & gaussian) {
    at::Tensor miu = gaussian.miu();
    at::Tensor average_dV  = miu.new_zeros(dimension);
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

// Return <▽V>, <q▽V>, <p▽V>,
//        <qq▽V>, <qp▽V>, <pp▽V>,
//        <qqq▽V>, <qqp▽V>, <qpp▽V>, <ppp▽V>
std::tuple<at::Tensor, at::Tensor, at::Tensor,
           at::Tensor, at::Tensor, at::Tensor,
           at::Tensor, at::Tensor, at::Tensor, at::Tensor>
compute_dV_4expectations(tchem::gaussian::Gaussian & gaussian) {
    at::Tensor miu = gaussian.miu();
    // 1st order
    at::Tensor average_dV = miu.new_zeros(dimension);
    // 2nd order
    at::Tensor average_qdV = miu.new_zeros({dimension, dimension});
    at::Tensor average_pdV = miu.new_zeros({dimension, dimension});
    // 3rd order
    at::Tensor average_qqdV = miu.new_zeros({dimension, dimension, dimension});
    at::Tensor average_qpdV = miu.new_zeros({dimension, dimension, dimension});
    at::Tensor average_ppdV = miu.new_zeros({dimension, dimension, dimension});
    // 4th order
    at::Tensor average_qqqdV = miu.new_zeros({dimension, dimension, dimension, dimension});
    at::Tensor average_qqpdV = miu.new_zeros({dimension, dimension, dimension, dimension});
    at::Tensor average_qppdV = miu.new_zeros({dimension, dimension, dimension, dimension});
    at::Tensor average_pppdV = miu.new_zeros({dimension, dimension, dimension, dimension});
    // Monte Carlo integration
    gaussian.rand_init();
    const size_t NSamples = 10000;
    for (size_t i = 0; i < NSamples; i++) {
        // First half is q, latter half is p
        at::Tensor rand_vec = gaussian.rand(generator);
        at::Tensor q = rand_vec.slice(0, 0, dimension),
                   p = rand_vec.slice(0, dimension);
        // 1st order
        at::Tensor dV = compute_dV(q);
        average_dV += dV;
        // 2nd order
        at::Tensor qdV = q.outer(dV),
                   pdV = p.outer(dV);
        average_qdV += qdV;
        average_pdV += pdV;
        // 3rd order
        at::Tensor qqdV = tchem::LA::outer_product(q, qdV),
                   qpdV = tchem::LA::outer_product(q, pdV),
                   ppdV = tchem::LA::outer_product(p, pdV);
        average_qqdV += qqdV;
        average_qpdV += qpdV;
        average_ppdV += ppdV;
        // 4th order
        at::Tensor qqqdV = tchem::LA::outer_product(q, qqdV),
                   qqpdV = tchem::LA::outer_product(q, qpdV),
                   qppdV = tchem::LA::outer_product(q, ppdV),
                   pppdV = tchem::LA::outer_product(p, ppdV);
        average_qqqdV += qqqdV;
        average_qqpdV += qqpdV;
        average_qppdV += qppdV;
        average_pppdV += pppdV;
    }
    average_dV    /= (double)NSamples; // 1st order
    average_qdV   /= (double)NSamples; // 2nd order
    average_pdV   /= (double)NSamples;
    average_qqdV  /= (double)NSamples; // 3rd order
    average_qpdV  /= (double)NSamples;
    average_ppdV  /= (double)NSamples;
    average_qqqdV /= (double)NSamples; // 4th order
    average_qqpdV /= (double)NSamples;
    average_qppdV /= (double)NSamples;
    average_pppdV /= (double)NSamples;
    return std::make_tuple(average_dV, average_qdV, average_pdV,
                           average_qqdV, average_qpdV, average_ppdV,
                           average_qqqdV, average_qqpdV, average_qppdV, average_pppdV);
}