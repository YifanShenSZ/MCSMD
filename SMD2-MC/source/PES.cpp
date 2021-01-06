// Oleg Prezhdo 2002 JCP 3rd order potential

#include <torch/torch.h>

at::Tensor compute_dV(const at::Tensor & q) {
    at::Tensor dV = q.new_empty(q.sizes());
    dV[0] = q[0] + 0.5 * q[0] * q[0];
    return dV;
}