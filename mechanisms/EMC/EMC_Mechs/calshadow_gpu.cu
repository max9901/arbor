#include <arbor/gpu/gpu_common.hpp>
#include <arbor/gpu/math_cu.hpp>
#include <arbor/gpu/reduce_by_key.hpp>
#include <arbor/mechanism_abi.h>
#include <stdio.h>

namespace arb {
namespace EMC_catalogue {

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
auto* _pp_var_k_q __attribute__((unused)) = params_.state_vars[0];\
auto* _pp_var_l_q __attribute__((unused)) = params_.state_vars[1];\
auto* _pp_var_gmax_actual __attribute__((unused)) = params_.state_vars[2];\
auto* _pp_var_vmin __attribute__((unused)) = params_.state_vars[3];\
auto* _pp_var_vmax __attribute__((unused)) = params_.state_vars[4];\
auto* _pp_var_timer __attribute__((unused)) = params_.state_vars[5];\
auto* _pp_var_vamp __attribute__((unused)) = params_.state_vars[6];\
auto* _pp_var_deltagmax __attribute__((unused)) = params_.state_vars[7];\
auto* _pp_var_k_q_lp __attribute__((unused)) = params_.state_vars[8];\
auto* _pp_var_totalTime __attribute__((unused)) = params_.state_vars[9];\
auto* _pp_var_stopped __attribute__((unused)) = params_.state_vars[10];\
auto* _pp_var_gion __attribute__((unused)) = params_.state_vars[11];\
auto* _pp_var_v __attribute__((unused)) = params_.state_vars[12];\
auto* _pp_var_celsius __attribute__((unused)) = params_.state_vars[13];\
auto* _pp_var_temperature __attribute__((unused)) = params_.state_vars[14];\
auto* _pp_var_eca __attribute__((unused)) = params_.state_vars[15];\
auto* _pp_var_k_steadyState_x __attribute__((unused)) = params_.state_vars[16];\
auto* _pp_var_k_timeCourse_t __attribute__((unused)) = params_.state_vars[17];\
auto* _pp_var_k_rateScale __attribute__((unused)) = params_.state_vars[18];\
auto* _pp_var_k_fcond __attribute__((unused)) = params_.state_vars[19];\
auto* _pp_var_k_inf __attribute__((unused)) = params_.state_vars[20];\
auto* _pp_var_k_tauUnscaled __attribute__((unused)) = params_.state_vars[21];\
auto* _pp_var_k_tau __attribute__((unused)) = params_.state_vars[22];\
auto* _pp_var_l_steadyState_x __attribute__((unused)) = params_.state_vars[23];\
auto* _pp_var_l_timeCourse_V __attribute__((unused)) = params_.state_vars[24];\
auto* _pp_var_l_timeCourse_t __attribute__((unused)) = params_.state_vars[25];\
auto* _pp_var_l_rateScale __attribute__((unused)) = params_.state_vars[26];\
auto* _pp_var_l_fcond __attribute__((unused)) = params_.state_vars[27];\
auto* _pp_var_l_inf __attribute__((unused)) = params_.state_vars[28];\
auto* _pp_var_l_tauUnscaled __attribute__((unused)) = params_.state_vars[29];\
auto* _pp_var_l_tau __attribute__((unused)) = params_.state_vars[30];\
auto* _pp_var_conductanceScale __attribute__((unused)) = params_.state_vars[31];\
auto* _pp_var_fopen0 __attribute__((unused)) = params_.state_vars[32];\
auto* _pp_var_fopen __attribute__((unused)) = params_.state_vars[33];\
auto* _pp_var_g __attribute__((unused)) = params_.state_vars[34];\
auto* _pp_var_rate_k_q __attribute__((unused)) = params_.state_vars[35];\
auto* _pp_var_rate_l_q __attribute__((unused)) = params_.state_vars[36];\
auto* _pp_var_rate_vstd __attribute__((unused)) = params_.state_vars[37];\
auto* _pp_var_rate_gmax __attribute__((unused)) = params_.state_vars[38];\
auto* _pp_var_rate_k_q_lp __attribute__((unused)) = params_.state_vars[39];\
auto* _pp_var_rate_totalTime __attribute__((unused)) = params_.state_vars[40];\
auto* _pp_var_mirror_cv __attribute__((unused)) = params_.state_vars[41];\
auto* _pp_var_gmax __attribute__((unused)) = params_.parameters[0];\
auto* _pp_var_conductance __attribute__((unused)) = params_.parameters[1];\
auto* _pp_var_k_instances __attribute__((unused)) = params_.parameters[2];\
auto* _pp_var_k_steadyState_rate __attribute__((unused)) = params_.parameters[3];\
auto* _pp_var_k_steadyState_midpoint __attribute__((unused)) = params_.parameters[4];\
auto* _pp_var_k_steadyState_scale __attribute__((unused)) = params_.parameters[5];\
auto* _pp_var_k_timeCourse_tau __attribute__((unused)) = params_.parameters[6];\
auto* _pp_var_l_instances __attribute__((unused)) = params_.parameters[7];\
auto* _pp_var_l_steadyState_rate __attribute__((unused)) = params_.parameters[8];\
auto* _pp_var_l_steadyState_midpoint __attribute__((unused)) = params_.parameters[9];\
auto* _pp_var_l_steadyState_scale __attribute__((unused)) = params_.parameters[10];\
auto* _pp_var_l_timeCourse_TIME_SCALE __attribute__((unused)) = params_.parameters[11];\
auto* _pp_var_l_timeCourse_VOLT_SCALE __attribute__((unused)) = params_.parameters[12];\
auto* _pp_var_stopAfter __attribute__((unused)) = params_.parameters[13];\
auto* _pp_var_silenceAfter __attribute__((unused)) = params_.parameters[14];\
auto* _pp_var_gmaxScaleAfterStop __attribute__((unused)) = params_.parameters[15];\
auto* _pp_var_target_vamp_min __attribute__((unused)) = params_.parameters[16];\
auto* _pp_var_target_vamp_max __attribute__((unused)) = params_.parameters[17];\
auto* _pp_var_target_vamp_silent __attribute__((unused)) = params_.parameters[18];\
auto* _pp_var_neuron_gid __attribute__((unused)) = params_.parameters[19];\
auto* _pp_var_other_cell __attribute__((unused)) = params_.parameters[20];\
auto* _pp_var_on_shadow_network __attribute__((unused)) = params_.parameters[21];\
auto& _pp_var_ion_ca __attribute__((unused)) = params_.ion_states[0];\
auto* _pp_var_ion_ca_index __attribute__((unused)) = params_.ion_states[0].index;\
//End of IFACEBLOCK

namespace {

using ::arb::gpu::exprelr;
using ::arb::gpu::safeinv;
using ::arb::gpu::min;
using ::arb::gpu::max;

__device__
void rates(arb_mechanism_ppack params_, int tid_, arb_value_type v) {
    PPACK_IFACE_BLOCK;
    _pp_var_k_steadyState_x[tid_] = _pp_var_k_steadyState_rate[tid_]/( 1.0+exp( -((v-_pp_var_k_steadyState_midpoint[tid_])/_pp_var_k_steadyState_scale[tid_])));
    _pp_var_k_timeCourse_t[tid_] = _pp_var_k_timeCourse_tau[tid_];
    _pp_var_k_rateScale[tid_] =  1.0;
    _pp_var_k_fcond[tid_] = pow(_pp_var_k_q[tid_], _pp_var_k_instances[tid_]);
    _pp_var_k_inf[tid_] = _pp_var_k_steadyState_x[tid_];
    _pp_var_k_tauUnscaled[tid_] = _pp_var_k_timeCourse_t[tid_];
    _pp_var_k_tau[tid_] = _pp_var_k_tauUnscaled[tid_]/_pp_var_k_rateScale[tid_];
    _pp_var_l_steadyState_x[tid_] = _pp_var_l_steadyState_rate[tid_]/( 1.0+exp( -((v-_pp_var_l_steadyState_midpoint[tid_])/_pp_var_l_steadyState_scale[tid_])));
    _pp_var_l_timeCourse_V[tid_] = v/_pp_var_l_timeCourse_VOLT_SCALE[tid_];
    _pp_var_l_timeCourse_t[tid_] = _pp_var_l_timeCourse_TIME_SCALE[tid_]*( 20.0*exp((_pp_var_l_timeCourse_V[tid_]+ 160.0)* 0.033333333333333333)/( 1.0+exp((_pp_var_l_timeCourse_V[tid_]+ 84.0)* 0.13698630136986301))+ 35.0);
    _pp_var_l_rateScale[tid_] =  1.0;
    _pp_var_l_fcond[tid_] = pow(_pp_var_l_q[tid_], _pp_var_l_instances[tid_]);
    _pp_var_l_inf[tid_] = _pp_var_l_steadyState_x[tid_];
    _pp_var_l_tauUnscaled[tid_] = _pp_var_l_timeCourse_t[tid_];
    _pp_var_l_tau[tid_] = _pp_var_l_tauUnscaled[tid_]/_pp_var_l_rateScale[tid_];
    _pp_var_rate_k_q[tid_] = (_pp_var_k_inf[tid_]-_pp_var_k_q[tid_])/_pp_var_k_tau[tid_];
    _pp_var_rate_l_q[tid_] = (_pp_var_l_inf[tid_]-_pp_var_l_q[tid_])/_pp_var_l_tau[tid_];
    if (_pp_var_vamp[tid_]< 3.0&&_pp_var_vmax[tid_]-_pp_var_vmin[tid_]< 3.0&&_pp_var_stopped[tid_]< 0.5) {
        _pp_var_rate_k_q_lp[tid_] =  0.01*(_pp_var_k_q[tid_]-_pp_var_k_q_lp[tid_]);
    }
    else {
        _pp_var_rate_k_q_lp[tid_] =  0.;
    }
    _pp_var_rate_totalTime[tid_] =  1.0;
}

__global__
void init(arb_mechanism_ppack params_) {
    int n_ = params_.width;
    int tid_ = threadIdx.x + blockDim.x*blockIdx.x;
    PPACK_IFACE_BLOCK;
    if (tid_<n_) {
        arb_size_type me = tid_;
        _pp_var_mirror_cv[me] = -1;
        for (int other = 0; other < n_; ++other) {
            // Sorry for O(n^2), but this is portable across CPU & GPU
            if (_pp_var_neuron_gid[other] == _pp_var_other_cell[me]) {
                _pp_var_mirror_cv[me] = other;
                break;
            }
        }
        if (_pp_var_mirror_cv[me] == -1) {
            printf("ERROR calpid: could not match neuron %d : gid=%d other=%d with another neuron\n",
                    (int)tid_,
                    (int)_pp_var_neuron_gid[me],
                    (int)_pp_var_other_cell[me]);
            printf("Force cellgroups?\n");
            _pp_var_mirror_cv[me] = _pp_var_mirror_cv[me*12313123131]; // just crash
        }
        auto node_indexi_ = _pp_var_node_index[tid_];
        arb_value_type v = _pp_var_vec_v[node_indexi_];
        _pp_var_eca[tid_] =  120.0;
        _pp_var_temperature[tid_] = _pp_var_celsius[tid_]+ 273.14999999999998;
        rates(params_, tid_, v);
        _pp_var_k_q[tid_] = _pp_var_k_inf[tid_];
        _pp_var_k_q_lp[tid_] = _pp_var_k_inf[tid_];
        _pp_var_l_q[tid_] = _pp_var_l_inf[tid_];
        _pp_var_gmax_actual[tid_] = _pp_var_gmax[tid_];
        _pp_var_timer[tid_] =  2.0 + tid_ / 1000;
        _pp_var_vamp[tid_] =  0.;
        _pp_var_totalTime[tid_] =  0.;
        _pp_var_stopped[tid_] =  0.;
    }
}

__global__
void multiply(arb_mechanism_ppack params_) {
    PPACK_IFACE_BLOCK;
    auto tid_ = threadIdx.x + blockDim.x*blockIdx.x;
    auto idx_ = blockIdx.y;    if(tid_<_pp_var_width) {
        _pp_var_state_vars[idx_][tid_] *= _pp_var_multiplicity[tid_];
    }
}

__global__
void advance_state(arb_mechanism_ppack params_) {
    int n_ = params_.width;
    int tid_ = threadIdx.x + blockDim.x*blockIdx.x;
    PPACK_IFACE_BLOCK;
    if (tid_<n_) {
        auto node_indexi_ = _pp_var_node_index[tid_];
        arb_value_type dt = _pp_var_vec_dt[node_indexi_];
        arb_value_type v = _pp_var_vec_v[node_indexi_];
        arb_value_type b_4_, b_3_, b_2_, b_0_, b_1_;
        rates(params_, tid_, v);
        b_0_ = _pp_var_rate_k_q[tid_];
        _pp_var_k_q[tid_] = _pp_var_k_q[tid_]+b_0_*dt;
        b_1_ = _pp_var_rate_k_q_lp[tid_];
        _pp_var_k_q_lp[tid_] = _pp_var_k_q_lp[tid_]+b_1_*dt;
        b_2_ = _pp_var_rate_l_q[tid_];
        _pp_var_l_q[tid_] = _pp_var_l_q[tid_]+b_2_*dt;
        b_3_ =  -0.0030000000000000001;
        _pp_var_timer[tid_] = _pp_var_timer[tid_]+b_3_*dt;
        b_4_ = _pp_var_rate_totalTime[tid_];
        _pp_var_totalTime[tid_] = _pp_var_totalTime[tid_]+b_4_*dt;
    }
}

__global__
void compute_currents(arb_mechanism_ppack params_) {
    int n_ = params_.width;
    int tid_ = threadIdx.x + blockDim.x*blockIdx.x;
    PPACK_IFACE_BLOCK;
    if (tid_<n_) {
        auto ion_ca_indexi_ = _pp_var_ion_ca_index[tid_];
        arb_size_type i_ = tid_;
        arb_size_type j_ = (arb_size_type)_pp_var_mirror_cv[i_];
        auto node_indexi_ = _pp_var_node_index[tid_];
        arb_value_type conductivity_ = 0;
        arb_value_type current_ = 0;
        arb_value_type v = _pp_var_vec_v[node_indexi_];
        arb_value_type ica = 0;
        arb_value_type dt = _pp_var_vec_dt[node_indexi_];
        _pp_var_conductanceScale[tid_] =  1.0;
        _pp_var_fopen0[tid_] = _pp_var_k_fcond[tid_]*_pp_var_l_fcond[tid_];
        _pp_var_fopen[tid_] = _pp_var_conductanceScale[tid_]*_pp_var_fopen0[tid_];
        _pp_var_g[tid_] = _pp_var_conductance[tid_]*_pp_var_fopen[tid_];
        _pp_var_gion[tid_] = _pp_var_gmax_actual[tid_]*_pp_var_fopen[tid_];
        ica = _pp_var_gion[tid_]*(v-_pp_var_eca[tid_]);
        if (v>_pp_var_vmax[tid_]) {
            _pp_var_vmax[tid_] = v;
        }
        if (v<_pp_var_vmin[tid_]) {
            _pp_var_vmin[tid_] = v;
        }
        _pp_var_vamp[tid_] = _pp_var_vmax[tid_]-_pp_var_vmin[tid_];
        if (_pp_var_stopped[tid_]< 0.5&&_pp_var_totalTime[tid_]>_pp_var_stopAfter[tid_]) {
            _pp_var_stopped[tid_] =  1.0;
            _pp_var_gmax_actual[tid_] = _pp_var_gmax_actual[tid_]*_pp_var_gmaxScaleAfterStop[tid_];
            _pp_var_vmin[tid_] = v;
            _pp_var_vmax[tid_] = v;
        }

        if (_pp_var_on_shadow_network[i_] > 0.5) {
            _pp_var_gmax_actual[i_] = _pp_var_gmax_actual[j_];
        }

        if (_pp_var_stopped[tid_]< 0.5&&_pp_var_timer[tid_]< 0.) {
            _pp_var_timer[tid_] =  1.0 + (tid_ / 10000.);
            if (_pp_var_on_shadow_network[i_] < 0.5) {
                if (_pp_var_silenceAfter[tid_]> 0.&&_pp_var_totalTime[tid_]>_pp_var_silenceAfter[tid_]) {
                    if (_pp_var_vamp[tid_]>_pp_var_target_vamp_silent[tid_]) {
                        _pp_var_gmax_actual[tid_] = _pp_var_gmax_actual[tid_]- 0.002 * (dt / 0.025);
                    }
                }
                else {
                    // Current cell (non-shadow)
                    bool ok = true;
                    if (_pp_var_target_vamp_max[tid_]> 3.0) {
                        if (_pp_var_vamp[tid_]<_pp_var_target_vamp_min[tid_]) {
                            if (_pp_var_k_q_lp[tid_]> 0.77000000000000002) {
                                _pp_var_gmax_actual[tid_] = _pp_var_gmax_actual[tid_] - 0.003 * (dt / 0.025);
                                ok = false;
                            }
                            else {
                                _pp_var_gmax_actual[tid_] = _pp_var_gmax_actual[tid_] + 0.003 * (dt / 0.025);
                                ok = false;
                            }
                        }
                        if (_pp_var_vamp[tid_]>_pp_var_target_vamp_max[tid_]&&_pp_var_vmax[tid_]< 0.) {
                            _pp_var_gmax_actual[tid_] = _pp_var_gmax_actual[tid_]- 0.001 * (dt / 0.025);
                            ok = false;
                        }
                    }
                    else {
                        if (_pp_var_vamp[tid_]>_pp_var_target_vamp_max[tid_]) {
                            _pp_var_gmax_actual[tid_] = _pp_var_gmax_actual[tid_]- 0.001 * (dt / 0.025);
                            ok = false;
                        }
                    }

                    if (ok) {
                        // Other cell (shadow)
                        if (_pp_var_target_vamp_max[j_] > 3.0) {
                            if (_pp_var_vamp[j_] < _pp_var_target_vamp_min[j_]) {
                                if (_pp_var_k_q_lp[j_] > 0.77) {
                                    _pp_var_gmax_actual[i_] = _pp_var_gmax_actual[i_] - 0.003/2.0 * (dt / 0.025);
                                }
                                else {
                                    _pp_var_gmax_actual[i_] = _pp_var_gmax_actual[i_] + 0.003/2.0 * (dt / 0.025);
                                }
                            }
                            if (_pp_var_vamp[j_] > _pp_var_target_vamp_max[j_]&&_pp_var_vmax[j_] < 0.) {
                                _pp_var_gmax_actual[i_] = _pp_var_gmax_actual[i_] - 0.001/2.0 * (dt / 0.025);
                            }
                        }
                        else {
                            if (_pp_var_vamp[j_] > _pp_var_target_vamp_max[j_]) {
                                _pp_var_gmax_actual[i_] = _pp_var_gmax_actual[i_] - 0.001/2.0 * (dt / 0.025);
                            }
                        }
                    }

                }
                if (_pp_var_gmax_actual[tid_] < 0.) {
                    _pp_var_gmax_actual[tid_] =  0.;
                }
                _pp_var_gmax_actual[j_] = _pp_var_gmax_actual[i_];
                _pp_var_vmin[tid_] = v;
                _pp_var_vmax[tid_] = v;
                _pp_var_vmin[j_] = v;
                _pp_var_vmax[j_] = v;
            }
        }
        if (_pp_var_stopped[tid_]> 0.5) {
            _pp_var_vamp[tid_] = _pp_var_vmax[tid_]-_pp_var_vmin[tid_];
        }
        current_ = ica;
        conductivity_ = _pp_var_gion[tid_];
        _pp_var_vec_g[node_indexi_] = fma(10.0*_pp_var_weight[tid_], conductivity_, _pp_var_vec_g[node_indexi_]);
        _pp_var_vec_i[node_indexi_] = fma(10.0*_pp_var_weight[tid_], current_, _pp_var_vec_i[node_indexi_]);
        _pp_var_ion_ca.current_density[ion_ca_indexi_] = fma(10.0*_pp_var_weight[tid_], ica, _pp_var_ion_ca.current_density[ion_ca_indexi_]);
    }
}

} // namespace

void mechanism_calshadow_gpu_init_(arb_mechanism_ppack* p) {
    auto n = p->width;
    unsigned block_dim = 128;
    unsigned grid_dim = ::arb::gpu::impl::block_count(n, block_dim);
    init<<<grid_dim, block_dim>>>(*p);
    if (!p->multiplicity) return;
    multiply<<<dim3{grid_dim, 11}, block_dim>>>(*p);
}

void mechanism_calshadow_gpu_compute_currents_(arb_mechanism_ppack* p) {
    auto n = p->width;
    unsigned block_dim = 128;
    unsigned grid_dim = ::arb::gpu::impl::block_count(n, block_dim);
    compute_currents<<<grid_dim, block_dim>>>(*p);
}

void mechanism_calshadow_gpu_advance_state_(arb_mechanism_ppack* p) {
    auto n = p->width;
    unsigned block_dim = 128;
    unsigned grid_dim = ::arb::gpu::impl::block_count(n, block_dim);
    advance_state<<<grid_dim, block_dim>>>(*p);
}

void mechanism_calshadow_gpu_write_ions_(arb_mechanism_ppack* p) {}

void mechanism_calshadow_gpu_post_event_(arb_mechanism_ppack* p) {}
void mechanism_calshadow_gpu_apply_events_(arb_mechanism_ppack* p, arb_deliverable_event_stream* events) {}

} // namespace EMC_catalogue
} // namespace arb
