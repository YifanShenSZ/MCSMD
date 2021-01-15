#include <torch/torch.h>

#include <tchem/linalg.hpp>
#include <tchem/gaussian.hpp>

#include "../include/basic.hpp"
#include "../include/expectation.hpp"

// Fit a gaussian phase space distribution based on current expectations
// If the current window is too far away from the distribution,
// then move the window to the current distribution
//
// We define the distance by:
//     the product of gaussian distribution and window is
//     some other coefficient * a gaussian coefficient * a gaussian function
//     so if the exponent of that gaussian coefficient > n^2, means it exceeds n sigma
void update_windows(
ExpectationSet & expectation_set,
WindowedSet & windowed_set) {
    // Fit a distribution from expectations
    tchem::gaussian::Gaussian distribution(expectation_set[1], expectation_set[2] - expectation_set[1].outer(expectation_set[1]));
    // Compute the distance
    at::Tensor miu1 = distribution.miu(), miu2 = windowed_set.window().miu(),
               var1 = distribution.var(), var2 = windowed_set.window().var();
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
        tchem::gaussian::Gaussian window = distribution.clone();
        at::Tensor coeff;
        tchem::gaussian::Gaussian product;
        std::tie(coeff, product) = distribution * window;
        tchem::polynomial::PolynomialSet variable_set(2 * dimension, 2);
        at::Tensor integrals = coeff * product.integral(variable_set);
        std::vector<at::Tensor> expectations = variable_set.views(integrals);
        expectations[0] = expectations[0][0];
        // 1 is a vector, so no need to transform to scalor nor tensor
        expectations[2] = tchem::LA::vec2sytensor(expectations[2], 2, 2 * dimension);
        windowed_set = WindowedSet(expectations, window);
    }
}