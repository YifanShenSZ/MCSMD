#include <torch/torch.h>

#include <tchem/gaussian.hpp>

// Based on current expectations and window and windowed expectations,
// 
std::tuple<at::Tensor, at::Tensor, at::Tensor> update_variables(
const at::Tensor & expectations,
const tchem::gaussian::Gaussian & window,
const at::Tensor & windowed_expectations,
const tchem::polynomial::PolynomialSet & variable_set,
const double & dt
)