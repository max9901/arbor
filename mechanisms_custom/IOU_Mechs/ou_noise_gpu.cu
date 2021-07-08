//
// Created by max on 08-07-21.
//

#include <arbor/gpu/gpu_common.hpp>
#include <arbor/gpu/math_cu.hpp>
#include <arbor/gpu/reduce_by_key.hpp>
#include <arbor/mechanism_abi.h>


#include "../tests/util_cuda.h"
#include <Random123/philox.h>

namespace arb {
    namespace IOU_catalogue {

#define PPACK_IFACE_BLOCK \
auto  _pp_var_width             __attribute__((unused)) = params_.width;\
auto  _pp_var_n_detectors       __attribute__((unused)) = params_.n_detectors;\
auto* _pp_var_vec_ci            __attribute__((unused)) = params_.vec_ci;\
auto* _pp_var_vec_di            __attribute__((unused)) = params_.vec_di;\
auto* _pp_var_vec_t             __attribute__((unused)) = params_.vec_t;\
auto* _pp_var_vec_dt            __attribute__((unused)) = params_.vec_dt;\
auto* _pp_var_vec_v             __attribute__((unused)) = params_.vec_v;\
auto* _pp_var_vec_i             __attribute__((unused)) = params_.vec_i;\
auto* _pp_var_vec_g             __attribute__((unused)) = params_.vec_g;\
auto* _pp_var_temperature_degC  __attribute__((unused)) = params_.temperature_degC;\
auto* _pp_var_diam_um           __attribute__((unused)) = params_.diam_um;\
auto* _pp_var_time_since_spike  __attribute__((unused)) = params_.time_since_spike;\
auto* _pp_var_node_index        __attribute__((unused)) = params_.node_index;\
auto* _pp_var_multiplicity      __attribute__((unused)) = params_.multiplicity;\
auto* _pp_var_state_vars        __attribute__((unused)) = params_.state_vars;\
auto* _pp_var_weight            __attribute__((unused)) = params_.weight;\
auto& _pp_var_events            __attribute__((unused)) = params_.events;\
auto& _pp_var_mechanism_id      __attribute__((unused)) = params_.mechanism_id;\
auto& _pp_var_index_constraints __attribute__((unused)) = params_.index_constraints; \
                          \
auto _pp_var_theta  __attribute__((unused)) = pp->globals[0];      \
auto _pp_var_sigma  __attribute__((unused)) = pp->globals[1];      \
auto _pp_var_mu     __attribute__((unused)) = pp->globals[2];      \
auto _pp_var_alpha  __attribute__((unused)) = pp->globals[3];       \
auto _pp_var_seed   __attribute__((unused)) = pp->globals[4];      \
auto* _pp_var_ouNoise __attribute__((unused)) = pp->state_vars[0];\
//End of IFACEBLOCK

namespace {

using ::arb::gpu::exprelr;
using ::arb::gpu::safeinv;
using ::arb::gpu::min;
using ::arb::gpu::max;

__global__
void init(arb_mechanism_ppack params_) {
    int n_ = params_.width;
    int tid_ = threadIdx.x + blockDim.x*blockIdx.x;

    PPACK_IFACE_BLOCK;

    if (tid_<n_) {
        _pp_var_ouNoise[tid_] = 0;
    }
}

__global__
void multiply(arb_mechanism_ppack params_) {
    PPACK_IFACE_BLOCK;

    auto tid_ = threadIdx.x + blockDim.x*blockIdx.x;
    auto idx_ = blockIdx.y;
    if(tid_<_pp_var_width) {
        _pp_var_state_vars[idx_][tid_] *= _pp_var_multiplicity[tid_];
    }

}

__global__
void advance_state(arb_mechanism_ppack params_) {}

__global__
void compute_currents(arb_mechanism_ppack params_) {
    int n_ = params_.width;
    int tid_ = threadIdx.x + blockDim.x*blockIdx.x;
    PPACK_IFACE_BLOCK;
    if (tid_<n_) {

    }

}

__global__
void write_ions(arb_mechanism_ppack params_) {}

} // namespace


void mechanism_ca_conc_gpu_init_(arb_mechanism_ppack* p) {
    unsigned block_dim = 128;
    unsigned grid_dim = ::arb::gpu::impl::block_count(p->width, block_dim);

    init<<<grid_dim, block_dim>>>(*p);
}

void mechanism_ca_conc_gpu_compute_currents_(arb_mechanism_ppack* p) {
    unsigned block_dim = 128;
    unsigned grid_dim = ::arb::gpu::impl::block_count(p->width, block_dim);

    compute_currents<<<grid_dim, block_dim>>>(*p);
}

void mechanism_ca_conc_gpu_advance_state_(arb_mechanism_ppack* p) {}
void mechanism_ou_noise_gpu_write_ions_(arb_mechanism_ppack* p) {}
void mechanism_ca_conc_gpu_post_event_(arb_mechanism_ppack* p) {}
void mechanism_ca_conc_gpu_apply_events_(arb_mechanism_ppack* p) {}

    } // namespace smol_catalogue
} // namespace arb
