#include <torch/torch.h>

#include <tchem/linalg.hpp>
#include <tchem/gaussian.hpp>

#include "../include/basic.hpp"

// Return <▽V>, <q▽V>, <p▽V>
std::tuple<at::Tensor, at::Tensor, at::Tensor> compute_dV_2expectations(tchem::gaussian::Gaussian & gaussian);

// Return <▽V>, <q▽V>, <p▽V>,
//        <qq▽V>, <qp▽V>, <pp▽V>,
//        <qqq▽V>, <qqp▽V>, <qpp▽V>, <ppp▽V>
std::tuple<at::Tensor, at::Tensor, at::Tensor,
           at::Tensor, at::Tensor, at::Tensor,
           at::Tensor, at::Tensor, at::Tensor, at::Tensor>
compute_dV_4expectations(tchem::gaussian::Gaussian & gaussian);

// Based on current expectations,
// fit a gaussian phase space distribution to compute gradient
std::tuple<std::vector<at::Tensor>, std::vector<at::Tensor>> compute_gradients(
const std::vector<at::Tensor> & expectations,
const tchem::gaussian::Gaussian & window,
const std::vector<at::Tensor> & windowed_expectations) {
    // Fit a gaussian from expectations
    tchem::gaussian::Gaussian g(expectations[1], expectations[2] - expectations[1].outer(expectations[1]));
    // Get ▽V related integrals
    at::Tensor dV, qdV, pdV;
    std::tie(dV, qdV, pdV) = compute_dV_2expectations(g);
    // Moyal equation of motion
    std::vector<at::Tensor> gradients(3);
    gradients[1] = expectations[1].new_empty(expectations[1].sizes());
    gradients[2] = expectations[2].new_empty(expectations[2].sizes());
    at::Tensor q  = expectations[1].slice(0, 0, dimension),
               p  = expectations[1].slice(0, dimension),
               qq = expectations[2].slice(0, 0, dimension).slice(1, 0, dimension),
               qp = expectations[2].slice(0, 0, dimension).slice(1, dimension),
               pp = expectations[2].slice(0, dimension).slice(1, dimension),
              dq  = gradients   [1].slice(0, 0, dimension),
              dp  = gradients   [1].slice(0, dimension),
              dqq = gradients   [2].slice(0, 0, dimension).slice(1, 0, dimension),
              dqp = gradients   [2].slice(0, 0, dimension).slice(1, dimension),
              dpp = gradients   [2].slice(0, dimension).slice(1, dimension);
    dq.copy_(p / mass);
    dp.copy_(-dV     );
    dqq.copy_((qp + qp.transpose(0, 1)) / mass);
    dqp.copy_(pp / mass - qdV                 );
    dpp.copy_(-(pdV + pdV.transpose(0, 1))    );

    // Get the product between distribution and window
    at::Tensor coeff;
    tchem::gaussian::Gaussian product;
    std::tie(coeff, product) = g * window;
    // Get higher order windowed expectations
    tchem::polynomial::PolynomialSet variable_set(2 * dimension, 4);
    at::Tensor integrals = coeff * product.integral(variable_set);
    std::vector<at::Tensor> views = variable_set.views(integrals);
    at::Tensor windowed_3 = tchem::LA::vec2sytensor(views[3], {2 * dimension, 2 * dimension, 2 * dimension});
    at::Tensor wqqq = windowed_3.slice(0, 0, dimension).slice(1, 0, dimension).slice(2, 0, dimension),
               wqqp = windowed_3.slice(0, 0, dimension).slice(1, 0, dimension).slice(2, dimension),
               wqpp = windowed_3.slice(0, 0, dimension).slice(1, dimension).slice(2, dimension),
               wppp = windowed_3.slice(0, dimension).slice(1, dimension).slice(2, dimension);
    at::Tensor windowed_4 = tchem::LA::vec2sytensor(views[4], {2 * dimension, 2 * dimension, 2 * dimension, 2 * dimension});
    at::Tensor wqqqq = windowed_4.slice(0, 0, dimension).slice(1, 0, dimension).slice(2, 0, dimension).slice(3, 0, dimension),
               wqqqp = windowed_4.slice(0, 0, dimension).slice(1, 0, dimension).slice(2, 0, dimension).slice(3, dimension),
               wqqpp = windowed_4.slice(0, 0, dimension).slice(1, 0, dimension).slice(2, dimension).slice(3, dimension),
               wqppp = windowed_4.slice(0, 0, dimension).slice(1, dimension).slice(2, dimension).slice(3, dimension),
               wpppp = windowed_4.slice(0, dimension).slice(1, dimension).slice(2, dimension).slice(3, dimension);
    // Get ▽V related integrals
    at::Tensor wdV, wqdV, wpdV,
               wqqdV, wqpdV, wppdV,
               wqqqdV, wqqpdV, wqppdV, wpppdV;
    std::tie(wdV, wqdV, wpdV,
             wqqdV, wqpdV, wppdV,
             wqqqdV, wqqpdV, wqppdV, wpppdV) = compute_dV_4expectations(product);
    wdV    *= coeff; // 1st order
    wqdV   *= coeff; // 2nd order
    wpdV   *= coeff;
    wqqdV  *= coeff; // 3rd order
    wqpdV  *= coeff;
    wppdV  *= coeff;
    wqqqdV *= coeff; // 4th order
    wqqpdV *= coeff;
    wqppdV *= coeff;
    wpppdV *= coeff;
    // Moyal equation of motion with window
    at::Tensor miu = window.miu(), var = window.var();
    at::Tensor varcholesky = var.cholesky(true);
    at::Tensor varinv = at::cholesky_inverse(varcholesky, true);
    at::Tensor varinv_qq = varinv.slice(0, 0, dimension).slice(1, 0, dimension),
               varinv_qp = varinv.slice(0, 0, dimension).slice(1, dimension),
               varinv_pq = varinv.slice(0, dimension).slice(1, 0, dimension),
               varinv_pp = varinv.slice(0, dimension).slice(1, dimension);
    at::Tensor varinv_miu   = varinv.mv(miu);
    at::Tensor varinv_miu_q = varinv_miu.slice(0, 0, dimension),
               varinv_miu_p = varinv_miu.slice(0, dimension);
    std::vector<at::Tensor> windowed_gradients(3);
    windowed_gradients[0] = (varinv_miu_q.dot(wp ) - at::trace(varinv_qq.mm(wqp )) - at::trace(varinv_qp.mm(wpp ))) / mass       
                          - (varinv_miu_p.dot(wdV) - at::trace(varinv_pq.mm(wqdV)) - at::trace(varinv_pp.mm(wpdV)));
    windowed_gradients[1] = windowed_expectations[1].new_empty(windowed_expectations[1].sizes());
    windowed_gradients[2] = windowed_expectations[2].new_empty(windowed_expectations[2].sizes());
    at::Tensor  wq  = windowed_expectations[1].slice(0, 0, dimension),
                wp  = windowed_expectations[1].slice(0, dimension),
                wqq = windowed_expectations[2].slice(0, 0, dimension).slice(1, 0, dimension),
                wqp = windowed_expectations[2].slice(0, 0, dimension).slice(1, dimension),
                wpp = windowed_expectations[2].slice(0, dimension).slice(1, dimension),
               dwq  = windowed_gradients   [1].slice(0, 0, dimension),
               dwp  = windowed_gradients   [1].slice(0, dimension),
               dwqq = windowed_gradients   [2].slice(0, 0, dimension).slice(1, 0, dimension),
               dwqp = windowed_gradients   [2].slice(0, 0, dimension).slice(1, dimension),
               dwpp = windowed_gradients   [2].slice(0, dimension).slice(1, dimension);
    

    return std::make_tuple(gradients, windowed_gradients);
}

// Update expectations with RK4
void update_expectations(
std::vector<at::Tensor> & expectations,
const tchem::gaussian::Gaussian & window,
const at::Tensor & windowed_expectations, const double & dt) {
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
