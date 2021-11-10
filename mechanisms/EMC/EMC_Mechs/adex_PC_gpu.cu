//
// Created by max on 08-07-21.
//

#include <arbor/gpu/gpu_common.hpp>
#include <arbor/mechanism_abi.h>

namespace arb {
    namespace EMC_catalogue {

#define PPACK_IFACE_BLOCK \
auto  _pp_var_width             __attribute__((unused)) = params_.width; \
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
auto* _pp_var_weight            __attribute__((unused)) = params_.weight;\
auto& _pp_var_events            __attribute__((unused)) = params_.events;\
auto& _pp_var_mechanism_id      __attribute__((unused)) = params_.mechanism_id;\
auto& _pp_var_index_constraints __attribute__((unused)) = params_.index_constraints; \
auto* _pp_var_w                 __attribute__((unused)) = params_.state_vars[0];\
auto* _pp_var_C                 __attribute__((unused)) = params_.parameters[0];\
auto* _pp_var_gL                __attribute__((unused)) = params_.parameters[1];\
auto* _pp_var_EL                __attribute__((unused)) = params_.parameters[2];\
auto* _pp_var_VT                __attribute__((unused)) = params_.parameters[3];\
auto* _pp_var_DeltaT            __attribute__((unused)) = params_.parameters[4];\
auto* _pp_var_Vr                __attribute__((unused)) = params_.parameters[5];\
auto* _pp_var_Vcut              __attribute__((unused)) = params_.parameters[6];\
auto* _pp_var_tauw              __attribute__((unused)) = params_.parameters[7];\
auto* _pp_var_a                 __attribute__((unused)) = params_.parameters[8];\
auto* _pp_var_b                 __attribute__((unused)) = params_.parameters[9];\
auto* _pp_var_I                 __attribute__((unused)) = params_.parameters[10];\
//End of IFACEBLOCK

        namespace {
            //init fixed
            __global__ void init(arb_mechanism_ppack params_) {
                PPACK_IFACE_BLOCK;
                const size_t i = threadIdx.x + blockDim.x*blockIdx.x;
                if (i<params_.width) {
                    const auto nidx = _pp_var_node_index[i];
                    _pp_var_Vcut[i] = _pp_var_VT[i] + 5. *_pp_var_DeltaT[i];
                    _pp_var_vec_i[nidx] = 0;
                }
            }

            __global__ void advance_state(arb_mechanism_ppack params_) {
                PPACK_IFACE_BLOCK;
                const size_t i = threadIdx.x + blockDim.x * blockIdx.x;
                if (i < params_.width) {
                    const auto nidx = _pp_var_node_index[i];
                    const arb_value_type vm = _pp_var_vec_v[nidx] * 1e-3  /*mV -> V*/;
                    const arb_value_type dt = _pp_var_vec_dt[nidx] * 1e-3 /*ms -> s*/;
                    const arb_value_type grad_w = (_pp_var_a[i]*(vm-_pp_var_EL[i]) - _pp_var_w[i])/_pp_var_tauw[i];
                    arb_value_type il = _pp_var_gL[i]*(_pp_var_EL[i]-vm);
                    arb_value_type id = _pp_var_gL[i]*_pp_var_DeltaT[i]*exp((vm-_pp_var_VT[i])/_pp_var_DeltaT[i]);
                    const arb_value_type grad_vm = (il + id + _pp_var_I[i] -_pp_var_w[i])/_pp_var_C[i];
                    _pp_var_w[i] += 1e3 * dt * grad_w;
                    _pp_var_vec_v[nidx] += 1e3 * dt*grad_vm;
                    if (vm > _pp_var_Vcut[i]) {
                        _pp_var_vec_v[i] = 1e3 * _pp_var_Vr[i];
                        _pp_var_w[i] += _pp_var_b[i];
                    }
                }
            }

            __global__ void apply_events(arb_mechanism_ppack params_, arb_deliverable_event_stream stream){
                PPACK_IFACE_BLOCK;
                auto tid_ = threadIdx.x + blockDim.x*blockIdx.x;
                if(tid_<stream.n_streams) {
                    auto begin  = stream.events + stream.begin[tid_];
                    auto end    = stream.events + stream.end[tid_];
                    for (auto p = begin; p<end; ++p) {
                        auto i_     = p->mech_index;
                        auto weight = p->weight;
                        if (p->mech_id==_pp_var_mechanism_id) {
                            _pp_var_w[i_] = _pp_var_w[i_] + weight;
                        }
                    }
                }
            }
        } // empty namespace to guard gpu kernels

        void mechanism_adex_PC_gpu_init_(arb_mechanism_ppack* p) {
            unsigned block_dim = 128;
            unsigned grid_dim = ::arb::gpu::impl::block_count(p->width, block_dim);
            init<<<grid_dim, block_dim>>>(*p);
        }
        void mechanism_adex_PC_gpu_compute_currents_(arb_mechanism_ppack* p) {}
        void mechanism_adex_PC_gpu_advance_state_(arb_mechanism_ppack* p) {
            unsigned block_dim = 128;
            unsigned grid_dim = ::arb::gpu::impl::block_count(p->width, block_dim);
            advance_state<<<grid_dim, block_dim>>>(*p);
        }
        void mechanism_adex_PC_gpu_write_ions_(arb_mechanism_ppack* p) {}
        void mechanism_adex_PC_gpu_post_event_(arb_mechanism_ppack* p) {}
        void mechanism_adex_PC_gpu_apply_events_(arb_mechanism_ppack* p, arb_deliverable_event_stream* stream_ptr) {
            auto n = stream_ptr->n_streams;
            unsigned block_dim = 128;
            unsigned grid_dim = ::arb::gpu::impl::block_count(n, block_dim);
            apply_events<<<grid_dim, block_dim>>>(*p, *stream_ptr);
        }

    } // namespace EMC_catalogue
} // namespace arb
