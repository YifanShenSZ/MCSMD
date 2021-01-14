#include <torch/torch.h>

at::Tensor compute_dV(const at::Tensor & q) {
    // Gaussian barrier
    const double k_peak = 10.0, width = 2.0;
    at::Tensor x = q.new_empty(q.sizes());
    x.copy_(q);
    x.set_requires_grad(true);
    if (x.grad().defined()) {
        x.grad().detach_();
        x.grad().zero_();
    }
    at::Tensor x_dimless = x / (6.283185307179586 / k_peak * width);
    at::Tensor V = k_peak * k_peak / 2.0 / 2000.0 
                 * exp(-0.5 * x_dimless * x_dimless);
    V.backward();
    return x.grad();
}
