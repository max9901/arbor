#include <arbor/gpu/gpu_common.hpp>
#include <arbor/gpu/math_cu.hpp>
#include <arbor/gpu/reduce_by_key.hpp>
#include <arbor/mechanism_abi.h>

namespace arb {
namespace new_default_catalogue {

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
auto* _pp_var_peer_index        __attribute__((unused)) = params_.peer_index;\
auto* _pp_var_multiplicity      __attribute__((unused)) = params_.multiplicity;\
auto* _pp_var_state_vars        __attribute__((unused)) = params_.state_vars;\
auto* _pp_var_weight            __attribute__((unused)) = params_.weight;\
auto& _pp_var_events            __attribute__((unused)) = params_.events;\
auto& _pp_var_mechanism_id      __attribute__((unused)) = params_.mechanism_id;\
auto& _pp_var_index_constraints __attribute__((unused)) = params_.index_constraints;\
auto* _pp_var_g __attribute__((unused)) = params_.parameters[0];\
//End of IFACEBLOCK

namespace {

using ::arb::gpu::exprelr;
using ::arb::gpu::safeinv;
using ::arb::gpu::min;
using ::arb::gpu::max;

__global__
void compute_currents(arb_mechanism_ppack params_) {
    int n_ = params_.width;
    int tid_ = threadIdx.x + blockDim.x*blockIdx.x;
    unsigned lane_mask_ = arb::gpu::ballot(0xffffffff, tid_<n_);
    PPACK_IFACE_BLOCK;
    if (tid_<n_) {
        auto peer_indexi_ = _pp_var_peer_index[tid_];
        auto node_indexi_ = _pp_var_node_index[tid_];
        arb_value_type conductivity_ = 0;
        arb_value_type current_ = 0;
        arb_value_type v_peer = _pp_var_vec_v[peer_indexi_];
        arb_value_type v = _pp_var_vec_v[node_indexi_];
        arb_value_type i = 0;
        i = _pp_var_g[tid_]*(v-v_peer);
        current_ = i;
        conductivity_ = _pp_var_g[tid_];
        ::arb::gpu::reduce_by_key(_pp_var_weight[tid_]*conductivity_,_pp_var_vec_g, node_indexi_, lane_mask_);
        ::arb::gpu::reduce_by_key(_pp_var_weight[tid_]*current_,_pp_var_vec_i, node_indexi_, lane_mask_);
    }
}

} // namespace

void mechanism_gj_gpu_init_(arb_mechanism_ppack* p) {}

void mechanism_gj_gpu_compute_currents_(arb_mechanism_ppack* p) {
    auto n = p->width;
    unsigned block_dim = 128;
    unsigned grid_dim = ::arb::gpu::impl::block_count(n, block_dim);
    compute_currents<<<grid_dim, block_dim>>>(*p);
}

void mechanism_gj_gpu_advance_state_(arb_mechanism_ppack* p) {}

void mechanism_gj_gpu_write_ions_(arb_mechanism_ppack* p) {}

void mechanism_gj_gpu_post_event_(arb_mechanism_ppack* p) {}
void mechanism_gj_gpu_apply_events_(arb_mechanism_ppack* p, arb_deliverable_event_stream* events) {}

} // namespace new_default_catalogue
} // namespace arb
