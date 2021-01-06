#include <torch/torch.h>

#include <tchem/gaussian.hpp>

// Single centre case:
//     1. fit a gaussian from expectations
//     2. compare to the windowed expectations
at::Tensor single_difference(
const at::Tensor & expectations,
const tchem::gaussian::Gaussian & window,
const at::Tensor & windowed_expectations,
const tchem::polynomial::PolynomialSet & variable_set) {
    // Fit a gaussian from expectations
    at::Tensor miu = expectations.new_empty(2);
    miu[0] = expectations[0];
    miu[1] = expectations[1];
    at::Tensor var = expectations.new_empty({2, 2});
    var[0][0] = expectations[2];
    var[0][1] = expectations[3];
    var[1][1] = expectations[4];
    tchem::gaussian::Gaussian g(miu, var);
    // Calculate windowed expectations from this gaussian
    at::Tensor coeff;
    tchem::gaussian::Gaussian g_product;
    std::tie(coeff, g_product) = g * window;
    at::Tensor predictions = coeff * g_product.integral(variable_set);
    // Return difference
    return predictions - windowed_expectations;
}