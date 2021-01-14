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
    // Flat surface
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
    // Flat surface
    return std::make_tuple(average_dV, average_qdV, average_pdV,
                           average_qqdV, average_qpdV, average_ppdV,
                           average_qqqdV, average_qqpdV, average_qppdV, average_pppdV);
}