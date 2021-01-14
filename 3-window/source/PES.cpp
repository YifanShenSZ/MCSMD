#include <torch/torch.h>

at::Tensor compute_dV(const at::Tensor & q) {
    // Flat surface
    return q.new_zeros(q.sizes());
}
