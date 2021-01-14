#include <torch/torch.h>

#include <tchem/linalg.hpp>
#include <tchem/gaussian.hpp>

#include "../include/basic.hpp"

// Return <▽V>, <q▽V>, <p▽V>
std::tuple<at::Tensor, at::Tensor, at::Tensor> compute_dV_2expectations(tchem::gaussian::Gaussian & gaussian);

// Based on current expectations,
// fit a gaussian phase space distribution to compute gradient
std::vector<at::Tensor> compute_gradients(const std::vector<at::Tensor> & expectations) {
    std::vector<at::Tensor> gradients(3);
    gradients[1] = expectations[1].new_empty(expectations[1].sizes());
    gradients[2] = expectations[2].new_empty(expectations[2].sizes());
    // Creat views according to the form in Moyal equation of motion
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
    // Fit a distribution from expectations
    tchem::gaussian::Gaussian g(expectations[1], expectations[2] - expectations[1].outer(expectations[1]));
    // Get ▽V related integrals
    at::Tensor dV, qdV, pdV;
    std::tie(dV, qdV, pdV) = compute_dV_2expectations(g);
    // Compute gradient from Moyal equation of motion
    dq .copy_(p / mass                        );
    dp .copy_(-dV                             );
    dqq.copy_((qp + qp.transpose(0, 1)) / mass);
    dqp.copy_(pp / mass - qdV                 );
    dpp.copy_(-(pdV + pdV.transpose(0, 1))    );
    return gradients;
}

// Return <▽V>, <q▽V>, <p▽V>,
//        <qq▽V>, <qp▽V>, <pp▽V>,
//        <qqq▽V>, <qqp▽V>, <qpp▽V>, <ppp▽V>
std::tuple<at::Tensor, at::Tensor, at::Tensor,
           at::Tensor, at::Tensor, at::Tensor,
           at::Tensor, at::Tensor, at::Tensor, at::Tensor>
compute_dV_4expectations(tchem::gaussian::Gaussian & gaussian);

// Based on current expectations,
// fit a gaussian phase space distribution to compute gradient
std::vector<at::Tensor> compute_windowed_gradients(
const std::vector<at::Tensor> & expectations,
const tchem::gaussian::Gaussian & window,
const std::vector<at::Tensor> & windowed_expectations) {
    std::vector<at::Tensor> windowed_gradients(3);
    // windowed_gradients[0] will be created as the output of an expression
    windowed_gradients[1] = windowed_expectations[1].new_empty(windowed_expectations[1].sizes());
    windowed_gradients[2] = windowed_expectations[2].new_empty(windowed_expectations[2].sizes());
    // Creat views according to the form in Moyal equation of motion
    at::Tensor wq  = windowed_expectations[1].slice(0, 0, dimension),
               wp  = windowed_expectations[1].slice(0, dimension),
               wqq = windowed_expectations[2].slice(0, 0, dimension).slice(1, 0, dimension),
               wqp = windowed_expectations[2].slice(0, 0, dimension).slice(1, dimension),
               // wpp is special since it will perform matrix multiplication
              dwq  = windowed_gradients   [1].slice(0, 0, dimension),
              dwp  = windowed_gradients   [1].slice(0, dimension),
              dwqq = windowed_gradients   [2].slice(0, 0, dimension).slice(1, 0, dimension),
              dwqp = windowed_gradients   [2].slice(0, 0, dimension).slice(1, dimension),
              dwpp = windowed_gradients   [2].slice(0, dimension)   .slice(1, dimension);
    at::Tensor wpp = wqq.new_empty(wqq.sizes());
    wpp.copy_(windowed_expectations[2].slice(0, dimension).slice(1, dimension));
    for (size_t i = 0    ; i < dimension; i++)
    for (size_t j = i + 1; j < dimension; j++) wpp[j][i] = wpp[i][j];
    // Fit a distribution from expectations
    tchem::gaussian::Gaussian g(expectations[1], expectations[2] - expectations[1].outer(expectations[1]));
    // Get the product between distribution and window
    at::Tensor coeff;
    tchem::gaussian::Gaussian product;
    std::tie(coeff, product) = g * window;
    // Get higher order windowed expectations
    tchem::polynomial::PolynomialSet variable_set(2 * dimension, 4);
    at::Tensor integrals = coeff * product.integral(variable_set);
    std::vector<at::Tensor> views = variable_set.views(integrals);
    at::Tensor windowed_3 = tchem::LA::vec2sytensor(views[3], 3, 2 * dimension);
    at::Tensor wqqq = windowed_3.slice(0, 0, dimension).slice(1, 0, dimension).slice(2, 0, dimension),
               wqqp = windowed_3.slice(0, 0, dimension).slice(1, 0, dimension).slice(2, dimension),
               wqpp = windowed_3.slice(0, 0, dimension).slice(1, dimension)   .slice(2, dimension),
               wppp = windowed_3.slice(0, dimension)   .slice(1, dimension)   .slice(2, dimension);
    at::Tensor windowed_4 = tchem::LA::vec2sytensor(views[4], 4, 2 * dimension);
    at::Tensor wqqqq = windowed_4.slice(0, 0, dimension).slice(1, 0, dimension).slice(2, 0, dimension).slice(3, 0, dimension),
               wqqqp = windowed_4.slice(0, 0, dimension).slice(1, 0, dimension).slice(2, 0, dimension).slice(3, dimension),
               wqqpp = windowed_4.slice(0, 0, dimension).slice(1, 0, dimension).slice(2, dimension)   .slice(3, dimension),
               wqppp = windowed_4.slice(0, 0, dimension).slice(1, dimension)   .slice(2, dimension)   .slice(3, dimension),
               wpppp = windowed_4.slice(0, dimension)   .slice(1, dimension)   .slice(2, dimension)   .slice(3, dimension);
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
    // Prepare terms required in differentiating the window
    at::Tensor miu = window.miu(), var = window.var();
    at::Tensor varcholesky = var.cholesky(true);
    at::Tensor varinv = at::cholesky_inverse(varcholesky, true);
    at::Tensor varinv_qq = varinv.slice(0, 0, dimension).slice(1, 0, dimension),
               varinv_qp = varinv.slice(0, 0, dimension).slice(1, dimension),
               varinv_pq = varinv.slice(0, dimension)   .slice(1, 0, dimension),
               varinv_pp = varinv.slice(0, dimension)   .slice(1, dimension);
    at::Tensor varinv_miu   = varinv.mv(miu);
    at::Tensor varinv_miu_q = varinv_miu.slice(0, 0, dimension),
               varinv_miu_p = varinv_miu.slice(0, dimension);
    // Compute gradient from Moyal equation of motion with window
    windowed_gradients[0] = (varinv_miu_q.dot(wp ) - at::trace(varinv_qq.mm(wqp )) - at::trace(varinv_qp.mm(wpp ))) / mass       
                          - (varinv_miu_p.dot(wdV) - at::trace(varinv_pq.mm(wqdV)) - at::trace(varinv_pp.mm(wpdV)));
    // dwq
    at::Tensor sumvec1 = dwq.new_zeros(dwq.sizes()),
               sumvec2 = dwq.new_zeros(dwq.sizes());
    for (size_t a = 0; a < dimension; a++)
    for (size_t i = 0; i < dimension; i++)
    for (size_t j = 0; j < dimension; j++) {
        // 1
        std::vector<size_t> indices = {a, j};
        std::sort(indices.begin(), indices.end());
        sumvec1[a] += varinv_qq[i][j] * wqqp[indices[0]][indices[1]][i];
        indices = {i, j};
        std::sort(indices.begin(), indices.end());
        sumvec1[a] += varinv_qp[i][j] * wqpp[a][indices[0]][indices[1]];
        // 2
        indices = {a, j};
        std::sort(indices.begin(), indices.end());
        sumvec2[a] += varinv_pq[i][j] * wqqdV[indices[0]][indices[1]][i];
        sumvec2[a] += varinv_pp[i][j] * wqpdV[a][j][i];
    }
    dwq.copy_(wp / mass
            + (wqp .mv(varinv_miu_q) - sumvec1) / mass
            - (wqdV.mv(varinv_miu_p) - sumvec2));
    // dwp
    sumvec1.zero_();
    sumvec2.zero_();
    for (size_t a = 0; a < dimension; a++)
    for (size_t i = 0; i < dimension; i++)
    for (size_t j = 0; j < dimension; j++) {
        // 1
        std::vector<size_t> indices = {a, i};
        std::sort(indices.begin(), indices.end());
        sumvec1[a] += varinv_qq[i][j] * wqpp[j][indices[0]][indices[1]];
        indices = {a, i, j};
        std::sort(indices.begin(), indices.end());
        sumvec1[a] += varinv_qp[i][j] * wppp[indices[0]][indices[1]][indices[2]];
        // 2
        sumvec2[a] += varinv_pq[i][j] * wqpdV[j][a][i];
        indices = {a, j};
        std::sort(indices.begin(), indices.end());
        sumvec2[a] += varinv_pp[i][j] * wppdV[indices[0]][indices[1]][i];
    }
    dwp.copy_(-wdV
            + (wpp .mv(varinv_miu_q) - sumvec1) / mass
            - (wpdV.mv(varinv_miu_p) - sumvec2));
    // dwqq
    at::Tensor summat11 = dwqq.new_zeros(dwqq.sizes()),
               summat12 = dwqq.new_zeros(dwqq.sizes()),
               summat21 = dwqq.new_zeros(dwqq.sizes()),
               summat22 = dwqq.new_zeros(dwqq.sizes());
    for (size_t a = 0; a < dimension; a++)
    for (size_t b = a; b < dimension; b++)
    for (size_t i = 0; i < dimension; i++) {
        // 11
        summat11[a][b] += varinv_miu_q[i] * wqqp [a][b][i];
        // 21
        summat21[a][b] += varinv_miu_p[i] * wqqdV[a][b][i];
        for (size_t j = 0; j < dimension; j++) {
            // 12
            std::vector<size_t> indices = {a, b, j};
            std::sort(indices.begin(), indices.end());
            summat12[a][b] += varinv_qq[i][j] * wqqqp[indices[0]][indices[1]][indices[2]][i];
            indices = {i, j};
            std::sort(indices.begin(), indices.end());
            summat12[a][b] += varinv_qp[i][j] * wqqpp[a][b][indices[0]][indices[1]];
            // 22
            indices = {a, b, j};
            std::sort(indices.begin(), indices.end());
            summat22[a][b] += varinv_pq[i][j] * wqqqdV[indices[0]][indices[1]][indices[2]][i];
            summat22[a][b] += varinv_pp[i][j] * wqqpdV[a][b][j][i];
        }
    }
    dwqq.copy_((wqp + wqp.transpose(0, 1)) / mass
             + (summat11 - summat12) / mass
             - (summat21 - summat22));
    // dwqp
    summat11.zero_();
    summat12.zero_();
    summat21.zero_();
    summat22.zero_();
    for (size_t a = 0; a < dimension; a++)
    for (size_t b = a; b < dimension; b++)
    for (size_t i = 0; i < dimension; i++) {
        // 11
        std::vector<size_t> indices = {b, i};
        std::sort(indices.begin(), indices.end());
        summat11[a][b] += varinv_miu_q[i] * wqpp[a][indices[0]][indices[1]];
        // 21
        summat21[a][b] += varinv_miu_p[i] * wqpdV[a][b][i];
        for (size_t j = 0; j < dimension; j++) {
            // 12
            std::vector<size_t> indicesq = {a, j},
                                indicesp = {b, i};
            std::sort(indicesq.begin(), indicesq.end());
            std::sort(indicesp.begin(), indicesp.end());
            summat12[a][b] += varinv_qq[i][j] * wqqpp[indicesq[0]][indicesq[1]][indicesp[0]][indicesp[1]];
            indices = {b, i ,j};
            std::sort(indices.begin(), indices.end());
            summat12[a][b] += varinv_qp[i][j] * wqppp[a][indices[0]][indices[1]][indices[2]];
            // 22
            indices = {a, j};
            std::sort(indices.begin(), indices.end());
            summat22[a][b] += varinv_pq[i][j] * wqqpdV[indices[0]][indices[1]][b][i];
            indices = {b, j};
            std::sort(indices.begin(), indices.end());
            summat22[a][b] += varinv_pp[i][j] * wqppdV[a][indices[0]][indices[1]][i];
        }
    }
    dwqp.copy_(wpp / mass - wqdV
             + (summat11 - summat12) / mass
             - (summat21 - summat22));
    // dwpp
    summat11.zero_();
    summat12.zero_();
    summat21.zero_();
    summat22.zero_();
    for (size_t a = 0; a < dimension; a++)
    for (size_t b = a; b < dimension; b++)
    for (size_t i = 0; i < dimension; i++) {
        // 11
        std::vector<size_t> indices = {a, b, i};
        std::sort(indices.begin(), indices.end());
        summat11[a][b] += varinv_miu_q[i] * wppp[indices[0]][indices[1]][indices[2]];
        // 21
        summat12[a][b] += varinv_miu_p[i] * wppdV[a][b][i];
        for (size_t j = 0; j < dimension; j++) {
            // 12
            summat12[a][b] += varinv_qq[i][j] * wqppp[j][indices[0]][indices[1]][indices[2]];
            indices = {a, b, i, j};
            std::sort(indices.begin(), indices.end());
            summat12[a][b] += varinv_qp[i][j] * wpppp[indices[0]][indices[1]][indices[2]][indices[3]];
            // 22
            summat22[a][b] += varinv_pq[i][j] * wqppdV[j][a][b][i];
            indices = {a, b, j};
            std::sort(indices.begin(), indices.end());
            summat22[a][b] += varinv_pp[i][j] * wpppdV[indices[0]][indices[1]][indices[2]][i];
        }
    }
    dwpp.copy_(-(wpdV + wpdV.transpose(0, 1))
             + (summat11 - summat12) / mass
             - (summat21 - summat22));
    return windowed_gradients;
}

// Update expectations with RK4
void update_expectations(
std::vector<at::Tensor> & expectations,
std::vector<at::Tensor> & windowed_expectations,
const tchem::gaussian::Gaussian & window,
const double & dt) {
    double dtd2 = dt / 2.0, dtd6 = dt / 6.0;
    std::vector<at::Tensor> intermediates(expectations.size()),
                            windowed_intermediates(windowed_expectations.size());
    // gradient 1
    std::vector<at::Tensor> gradients1 = compute_gradients(expectations),
                            windowed_gradients1 = compute_windowed_gradients(expectations, window, windowed_expectations);
    // gradient 2
    windowed_intermediates[0] = windowed_expectations[0] + dtd2 * windowed_gradients1[0];
    for (size_t i = 1; i < expectations.size(); i++) {
                 intermediates[i] =          expectations[i] + dtd2 *          gradients1[i];
        windowed_intermediates[i] = windowed_expectations[i] + dtd2 * windowed_gradients1[i];
    }
    std::vector<at::Tensor> gradients2 = compute_gradients(intermediates),
                            windowed_gradients2 = compute_windowed_gradients(intermediates, window, windowed_intermediates);
    // gradient 3
    windowed_intermediates[0] = windowed_expectations[0] + dtd2 * windowed_gradients2[0];
    for (size_t i = 1; i < expectations.size(); i++) {
                 intermediates[i] =          expectations[i] + dtd2 *          gradients2[i];
        windowed_intermediates[i] = windowed_expectations[i] + dtd2 * windowed_gradients2[i];
    }
    std::vector<at::Tensor> gradients3 = compute_gradients(intermediates),
                            windowed_gradients3 = compute_windowed_gradients(intermediates, window, windowed_intermediates);
    // gradient 4
    windowed_intermediates[0] = windowed_expectations[0] + dt * windowed_gradients3[0];
    for (size_t i = 1; i < expectations.size(); i++) {
                 intermediates[i] =          expectations[i] + dt *          gradients3[i];
        windowed_intermediates[i] = windowed_expectations[i] + dt * windowed_gradients3[i];
    }
    std::vector<at::Tensor> gradients4 = compute_gradients(intermediates),
                            windowed_gradients4 = compute_windowed_gradients(intermediates, window, windowed_intermediates);
    // extrapolate
    windowed_expectations[0] += dtd6 * (windowed_gradients1[0] + 2.0 * windowed_gradients2[0] + 2.0 * windowed_gradients3[0] + windowed_gradients4[0]);
    for (size_t i = 1; i < expectations.size(); i++) {
                 expectations[i] += dtd6 * (         gradients1[i] + 2.0 *          gradients2[i] + 2.0 *          gradients3[i] +          gradients4[i]);
        windowed_expectations[i] += dtd6 * (windowed_gradients1[i] + 2.0 * windowed_gradients2[i] + 2.0 * windowed_gradients3[i] + windowed_gradients4[i]);
    }
}
