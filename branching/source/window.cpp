#include <torch/torch.h>

#include <tchem/linalg.hpp>
#include <tchem/gaussian.hpp>

#include "../include/basic.hpp"

// Based on current expectations, fit a gaussian phase space distribution
// If the predictions differ from windowed expectations too much,
// then fit the expectations with/without window with 2 gaussians
// else if <w> is too small,
// then move the windows to the current gaussians
void update_windows(
const std::vector<at::Tensor> & expectations,
const std::vector<at::Tensor> & windowed_expectations,
tchem::gaussian::Gaussian & window) {
    // Fit a distribution from expectations
    tchem::gaussian::Gaussian distribution(expectations[1], expectations[2] - expectations[1].outer(expectations[1]));
    // Get the product between distribution and window
    at::Tensor coeff;
    tchem::gaussian::Gaussian product;
    std::tie(coeff, product) = distribution * window;
    // Predict 2nd and lower order windowed expectations
    tchem::polynomial::PolynomialSet variable_set(2 * dimension, 2);
    at::Tensor integrals = coeff * product.integral(variable_set);
    std::vector<at::Tensor> predictions = variable_set.views(integrals);
    predictions[0] = predictions[0][0];
    // 1 is a vector, so no need to transform to scalor nor tensor
    predictions[2] = tchem::LA::vec2sytensor(predictions[2], 2, 2 * dimension);
    return predictions;
}