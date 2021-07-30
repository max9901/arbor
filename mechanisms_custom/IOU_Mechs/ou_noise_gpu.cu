//
// Created by max on 08-07-21.
//

#include <arbor/gpu/gpu_common.hpp>
#include <arbor/mechanism_abi.h>
#include <Random123/philox.h>
#include <Random123/boxmuller.hpp>

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
auto _pp_var_theta              __attribute__((unused)) = params_.globals[0];        \
auto _pp_var_sigma              __attribute__((unused)) = params_.globals[1];        \
auto _pp_var_mu                 __attribute__((unused)) = params_.globals[2];        \
auto _pp_var_alpha              __attribute__((unused)) = params_.globals[3];        \
auto _pp_var_seed               __attribute__((unused)) = params_.globals[4];        \
auto& _pp_var_cnt               __attribute__((unused)) = params_.globals[5];        \
auto* _pp_var_ouNoise           __attribute__((unused)) = params_.state_vars[0];     \
//End of IFACEBLOCK

namespace {

__global__
void init(arb_mechanism_ppack params_) {  //done
    const size_t tid_ = threadIdx.x + blockDim.x*blockIdx.x;
    if(!tid_) params_.globals[5] = 0;
    if (tid_<params_.width) {
        params_.state_vars[0][tid_] = 0;
    }
}

__global__ void compute_currents(arb_mechanism_ppack params_) {
    PPACK_IFACE_BLOCK;
    const size_t tid_ = threadIdx.x + blockDim.x*blockIdx.x;
    philox2x32_key_t k = {{(uint32_t)_pp_var_seed}};
    philox2x32_ctr_t c = {{(uint32_t)_pp_var_cnt}};
    philox2x32_ctr_t cresult = philox2x32(c, k);
    const float  rand_global = r123::boxmuller(cresult.v[0],cresult.v[1]).x;
    if (tid_<params_.width) {
        c.v[0] = _pp_var_cnt+tid_+1;
        const auto node_indexi_          = _pp_var_node_index[tid_];
        const arb_value_type dt          = _pp_var_vec_dt[node_indexi_];
        const arb_value_type sqrt_dt     = std::sqrt(dt);
        const arb_value_type Iapp_global = _pp_var_sigma * sqrt_dt * rand_global;
        cresult = philox2x32(c, k);
        const float rand_local = r123::boxmuller(cresult.v[0],cresult.v[1]).x;
        _pp_var_ouNoise[tid_] +=
                _pp_var_theta * (_pp_var_mu - _pp_var_ouNoise[tid_]) * dt +
                (1-_pp_var_alpha) * _pp_var_sigma * sqrt_dt * rand_local
                + _pp_var_alpha * Iapp_global;
        _pp_var_vec_i[node_indexi_] -= _pp_var_ouNoise[tid_];
    }
    if(!tid_) _pp_var_cnt += params_.width+1;
}
} // namespace

void mechanism_ou_noise_gpu_init_(arb_mechanism_ppack* p) {
    unsigned block_dim = 128;
    unsigned grid_dim = ::arb::gpu::impl::block_count(p->width, block_dim);
    init<<<grid_dim, block_dim>>>(*p);
}
void mechanism_ou_noise_gpu_compute_currents_(arb_mechanism_ppack* p) {
    unsigned block_dim = 128;
    unsigned grid_dim = ::arb::gpu::impl::block_count(p->width, block_dim);
    compute_currents<<<grid_dim, block_dim>>>(*p);
}
void mechanism_ou_noise_gpu_advance_state_(arb_mechanism_ppack* p) {}
void mechanism_ou_noise_gpu_write_ions_(arb_mechanism_ppack* p) {}
void mechanism_ou_noise_gpu_post_event_(arb_mechanism_ppack* p) {}
void mechanism_ou_noise_gpu_apply_events_(arb_mechanism_ppack* p, arb_deliverable_event_stream*) {}

} // namespace IOU_catalogue
} // namespace arb
