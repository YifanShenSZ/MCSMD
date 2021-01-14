#include <torch/torch.h>

#include <tchem/linalg.hpp>
#include <tchem/gaussian.hpp>

#include "../include/basic.hpp"

// Fit a gaussian phase space distribution based on current expectations
// If the current window is too far away from the distribution,
// then move the window to the current distribution
//
// We define the distance by:
//     the product of gaussian distribution and window is
//     some other coefficient * a gaussian coefficient * a gaussian function
//     so if the exponent of that gaussian coefficient > n^2, means it exceeds n sigma
void update_windows(
const std::vector<at::Tensor> & expectations,
std::vector<at::Tensor> & windowed_expectations,
tchem::gaussian::Gaussian & window) {
    // Fit a distribution from expectations
    tchem::gaussian::Gaussian distribution(expectations[1], expectations[2] - expectations[1].outer(expectations[1]));
    // Compute the distance
    at::Tensor miu1 = distribution.miu(), var1 = distribution.var(),
               miu2 = window      .miu(), var2 = window      .var();
    at::Tensor cholesky_var1 = var1.cholesky(true),
               cholesky_var2 = var2.cholesky(true);
    at::Tensor inv_var1 = at::cholesky_inverse(cholesky_var1, true),
               inv_var2 = at::cholesky_inverse(cholesky_var2, true);
    at::Tensor inv_var3 = inv_var1 + inv_var2;
    at::Tensor cholesky_inv_var3 = inv_var3.cholesky(true);
    at::Tensor var3 = at::cholesky_inverse(cholesky_inv_var3, true);
    at::Tensor temp = inv_var1.mv(miu1) + inv_var2.mv(miu2);
    at::Tensor miu3 = var3.mv(temp);
    double distance = (miu1.dot(inv_var1.mv(miu1)) + miu2.dot(inv_var2.mv(miu2)) - temp.dot(miu3)).item<double>();

    // Move the window if too far away
    if (distance > 1.0) {
        window = distribution;
        at::Tensor coeff;
        tchem::gaussian::Gaussian product;
        std::tie(coeff, product) = distribution * window;
        tchem::polynomial::PolynomialSet variable_set(2 * dimension, 2);
        at::Tensor integrals = coeff * product.integral(variable_set);
        windowed_expectations = variable_set.views(integrals);
        windowed_expectations[0] = windowed_expectations[0][0];
        // 1 is a vector, so no need to transform to scalor nor tensor
        windowed_expectations[2] = tchem::LA::vec2sytensor(windowed_expectations[2], 2, 2 * dimension);
    }
}