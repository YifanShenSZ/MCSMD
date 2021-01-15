#include <torch/torch.h>

#include <tchem/linalg.hpp>
#include <tchem/gaussian.hpp>

#include "../include/basic.hpp"

// Fit a gaussian phase space distribution based on current expectations
// In this version we only compute the predictions to compare with the windowed expectations
std::vector<at::Tensor> update_windows(
const std::vector<at::Tensor> & expectations,
const tchem::gaussian::Gaussian & window) {
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