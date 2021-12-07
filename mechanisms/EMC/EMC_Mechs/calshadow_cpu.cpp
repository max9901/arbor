#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <arbor/mechanism_abi.h>
#include <arbor/math.hpp>
#include <iostream>

namespace arb {
namespace EMC_catalogue {
namespace kernel_calshadow {

using ::arb::math::exprelr;
using ::arb::math::safeinv;
using ::std::abs;
using ::std::cos;
using ::std::exp;
using ::std::log;
using ::std::max;
using ::std::min;
using ::std::pow;
using ::std::sin;

static constexpr unsigned simd_width_ = 1;
static constexpr unsigned min_align_ = std::max(alignof(arb_value_type), alignof(arb_index_type));

#define PPACK_IFACE_BLOCK \
[[maybe_unused]] auto  _pp_var_width             = pp->width;\
[[maybe_unused]] auto  _pp_var_n_detectors       = pp->n_detectors;\
[[maybe_unused]] auto* _pp_var_vec_ci            = pp->vec_ci;\
[[maybe_unused]] auto* _pp_var_vec_di            = pp->vec_di;\
[[maybe_unused]] auto* _pp_var_vec_t             = pp->vec_t;\
[[maybe_unused]] auto* _pp_var_vec_dt            = pp->vec_dt;\
[[maybe_unused]] auto* _pp_var_vec_v             = pp->vec_v;\
[[maybe_unused]] auto* _pp_var_vec_i             = pp->vec_i;\
[[maybe_unused]] auto* _pp_var_vec_g             = pp->vec_g;\
[[maybe_unused]] auto* _pp_var_temperature_degC  = pp->temperature_degC;\
[[maybe_unused]] auto* _pp_var_diam_um           = pp->diam_um;\
[[maybe_unused]] auto* _pp_var_time_since_spike  = pp->time_since_spike;\
[[maybe_unused]] auto* _pp_var_node_index        = pp->node_index;\
[[maybe_unused]] auto* _pp_var_peer_index        = pp->peer_index;\
[[maybe_unused]] auto* _pp_var_multiplicity      = pp->multiplicity;\
[[maybe_unused]] auto* _pp_var_weight            = pp->weight;\
[[maybe_unused]] auto& _pp_var_events            = pp->events;\
[[maybe_unused]] auto& _pp_var_mechanism_id      = pp->mechanism_id;\
[[maybe_unused]] auto& _pp_var_index_constraints = pp->index_constraints;\
[[maybe_unused]] auto* _pp_var_k_q = pp->state_vars[0];\
[[maybe_unused]] auto* _pp_var_l_q = pp->state_vars[1];\
[[maybe_unused]] auto* _pp_var_gmax_actual = pp->state_vars[2];\
[[maybe_unused]] auto* _pp_var_vmin = pp->state_vars[3];\
[[maybe_unused]] auto* _pp_var_vmax = pp->state_vars[4];\
[[maybe_unused]] auto* _pp_var_timer = pp->state_vars[5];\
[[maybe_unused]] auto* _pp_var_vamp = pp->state_vars[6];\
[[maybe_unused]] auto* _pp_var_deltagmax = pp->state_vars[7];\
[[maybe_unused]] auto* _pp_var_k_q_lp = pp->state_vars[8];\
[[maybe_unused]] auto* _pp_var_totalTime = pp->state_vars[9];\
[[maybe_unused]] auto* _pp_var_stopped = pp->state_vars[10];\
[[maybe_unused]] auto* _pp_var_gion = pp->state_vars[11];\
[[maybe_unused]] auto* _pp_var_v = pp->state_vars[12];\
[[maybe_unused]] auto* _pp_var_celsius = pp->state_vars[13];\
[[maybe_unused]] auto* _pp_var_temperature = pp->state_vars[14];\
[[maybe_unused]] auto* _pp_var_eca = pp->state_vars[15];\
[[maybe_unused]] auto* _pp_var_k_steadyState_x = pp->state_vars[16];\
[[maybe_unused]] auto* _pp_var_k_timeCourse_t = pp->state_vars[17];\
[[maybe_unused]] auto* _pp_var_k_rateScale = pp->state_vars[18];\
[[maybe_unused]] auto* _pp_var_k_fcond = pp->state_vars[19];\
[[maybe_unused]] auto* _pp_var_k_inf = pp->state_vars[20];\
[[maybe_unused]] auto* _pp_var_k_tauUnscaled = pp->state_vars[21];\
[[maybe_unused]] auto* _pp_var_k_tau = pp->state_vars[22];\
[[maybe_unused]] auto* _pp_var_l_steadyState_x = pp->state_vars[23];\
[[maybe_unused]] auto* _pp_var_l_timeCourse_V = pp->state_vars[24];\
[[maybe_unused]] auto* _pp_var_l_timeCourse_t = pp->state_vars[25];\
[[maybe_unused]] auto* _pp_var_l_rateScale = pp->state_vars[26];\
[[maybe_unused]] auto* _pp_var_l_fcond = pp->state_vars[27];\
[[maybe_unused]] auto* _pp_var_l_inf = pp->state_vars[28];\
[[maybe_unused]] auto* _pp_var_l_tauUnscaled = pp->state_vars[29];\
[[maybe_unused]] auto* _pp_var_l_tau = pp->state_vars[30];\
[[maybe_unused]] auto* _pp_var_conductanceScale = pp->state_vars[31];\
[[maybe_unused]] auto* _pp_var_fopen0 = pp->state_vars[32];\
[[maybe_unused]] auto* _pp_var_fopen = pp->state_vars[33];\
[[maybe_unused]] auto* _pp_var_g = pp->state_vars[34];\
[[maybe_unused]] auto* _pp_var_rate_k_q = pp->state_vars[35];\
[[maybe_unused]] auto* _pp_var_rate_l_q = pp->state_vars[36];\
[[maybe_unused]] auto* _pp_var_rate_vstd = pp->state_vars[37];\
[[maybe_unused]] auto* _pp_var_rate_gmax = pp->state_vars[38];\
[[maybe_unused]] auto* _pp_var_rate_k_q_lp = pp->state_vars[39];\
[[maybe_unused]] auto* _pp_var_rate_totalTime = pp->state_vars[40];\
[[maybe_unused]] auto* _pp_var_mirror_cv = pp->state_vars[41];\
[[maybe_unused]] auto* _pp_var_gmax = pp->parameters[0];\
[[maybe_unused]] auto* _pp_var_conductance = pp->parameters[1];\
[[maybe_unused]] auto* _pp_var_k_instances = pp->parameters[2];\
[[maybe_unused]] auto* _pp_var_k_steadyState_rate = pp->parameters[3];\
[[maybe_unused]] auto* _pp_var_k_steadyState_midpoint = pp->parameters[4];\
[[maybe_unused]] auto* _pp_var_k_steadyState_scale = pp->parameters[5];\
[[maybe_unused]] auto* _pp_var_k_timeCourse_tau = pp->parameters[6];\
[[maybe_unused]] auto* _pp_var_l_instances = pp->parameters[7];\
[[maybe_unused]] auto* _pp_var_l_steadyState_rate = pp->parameters[8];\
[[maybe_unused]] auto* _pp_var_l_steadyState_midpoint = pp->parameters[9];\
[[maybe_unused]] auto* _pp_var_l_steadyState_scale = pp->parameters[10];\
[[maybe_unused]] auto* _pp_var_l_timeCourse_TIME_SCALE = pp->parameters[11];\
[[maybe_unused]] auto* _pp_var_l_timeCourse_VOLT_SCALE = pp->parameters[12];\
[[maybe_unused]] auto* _pp_var_stopAfter = pp->parameters[13];\
[[maybe_unused]] auto* _pp_var_silenceAfter = pp->parameters[14];\
[[maybe_unused]] auto* _pp_var_gmaxScaleAfterStop = pp->parameters[15];\
[[maybe_unused]] auto* _pp_var_target_vamp_min = pp->parameters[16];\
[[maybe_unused]] auto* _pp_var_target_vamp_max = pp->parameters[17];\
[[maybe_unused]] auto* _pp_var_target_vamp_silent = pp->parameters[18];\
[[maybe_unused]] auto* _pp_var_neuron_gid = pp->parameters[19];\
[[maybe_unused]] auto* _pp_var_other_cell = pp->parameters[20];\
[[maybe_unused]] auto* _pp_var_on_shadow_network = pp->parameters[21];\
[[maybe_unused]] auto& _pp_var_ion_ca = pp->ion_states[0];\
[[maybe_unused]] auto* _pp_var_ion_ca_index = pp->ion_states[0].index;\
//End of IFACEBLOCK

// procedure prototypes
[[maybe_unused]] static void rates(arb_mechanism_ppack* pp, int i_, arb_value_type v);

        //// node_indexi only for current and voltage and other globals like dt etc
        //

// interface methods
static void init(arb_mechanism_ppack* pp) {
    PPACK_IFACE_BLOCK;
    std::cout << "INIT --------------- " << std::endl;
    for (arb_size_type i_ = 0; i_ < _pp_var_width; ++i_) {
        for (arb_size_type j_ = 0; j_ < _pp_var_width; ++j_) {
            // Sorry for O(n^2), but this is portably across CPU & GPU
            _pp_var_mirror_cv[i_] = -1;
            if (_pp_var_neuron_gid[j_] == _pp_var_other_cell[i_]) {
                _pp_var_mirror_cv[i_] = j_;
                break;
            }
        }
        if (_pp_var_mirror_cv[i_] == -1) {
            std::cout << "No correspondence between neurons and mirror neurons" << std::endl;
            exit(42);
        }
        std::cout <<
            " cv: " << i_ <<
            " gid: " << _pp_var_neuron_gid[i_]<<
            " oth: " << _pp_var_other_cell[i_] <<
            " othcv: " << _pp_var_mirror_cv[i_] <<
            " osn: " << _pp_var_on_shadow_network[i_] <<
            " min/max:" << _pp_var_target_vamp_min[i_] << "/" << _pp_var_target_vamp_max[i_] <<
            std::endl;
    }
    for (arb_size_type i_ = 0; i_ < _pp_var_width; ++i_) {
        auto node_indexi_ = _pp_var_node_index[i_];
        arb_value_type v = _pp_var_vec_v[node_indexi_];
        _pp_var_eca[i_] =  120.0;
        _pp_var_temperature[i_] = _pp_var_celsius[i_]+ 273.14999999999998;
        rates(pp, i_, v);
        _pp_var_k_q[i_] = _pp_var_k_inf[i_];
        _pp_var_k_q_lp[i_] = _pp_var_k_inf[i_];
        _pp_var_l_q[i_] = _pp_var_l_inf[i_];
        _pp_var_gmax_actual[i_] = _pp_var_gmax[i_];
        _pp_var_timer[i_] =  2.0;
        _pp_var_vamp[i_] =  0.;
        _pp_var_totalTime[i_] =  0.;
        _pp_var_stopped[i_] =  0.;
    }
    if (!_pp_var_multiplicity) return;
    for (arb_size_type ix = 0; ix < 11; ++ix) {
        for (arb_size_type iy = 0; iy < _pp_var_width; ++iy) {
            pp->state_vars[ix][iy] *= _pp_var_multiplicity[iy];
        }
    }
}

static void advance_state(arb_mechanism_ppack* pp) {
    PPACK_IFACE_BLOCK;
    for (arb_size_type i_ = 0; i_ < _pp_var_width; ++i_) {
        auto node_indexi_ = _pp_var_node_index[i_];
        arb_value_type dt = _pp_var_vec_dt[node_indexi_];
        arb_value_type v = _pp_var_vec_v[node_indexi_];
        arb_value_type b_4_, b_3_, b_2_, b_0_, b_1_;
        rates(pp, i_, v);
        b_0_ = _pp_var_rate_k_q[i_];
        _pp_var_k_q[i_] = _pp_var_k_q[i_]+b_0_*dt;
        b_1_ = _pp_var_rate_k_q_lp[i_];
        _pp_var_k_q_lp[i_] = _pp_var_k_q_lp[i_]+b_1_*dt;
        b_2_ = _pp_var_rate_l_q[i_];
        _pp_var_l_q[i_] = _pp_var_l_q[i_]+b_2_*dt;
        b_3_ =  -0.0030000000000000001;
        _pp_var_timer[i_] = _pp_var_timer[i_]+b_3_*dt;
        b_4_ = _pp_var_rate_totalTime[i_];
        _pp_var_totalTime[i_] = _pp_var_totalTime[i_]+b_4_*dt;
    }
}

static void compute_currents(arb_mechanism_ppack* pp) {
    PPACK_IFACE_BLOCK;
    for (arb_size_type i_ = 0; i_ < _pp_var_width; ++i_) {
        auto ion_ca_indexi_ = _pp_var_ion_ca_index[i_];
        auto node_indexi_ = _pp_var_node_index[i_];
        arb_size_type j_ = (arb_size_type)_pp_var_mirror_cv[i_];
        arb_value_type dt = _pp_var_vec_dt[node_indexi_];
        arb_value_type conductivity_ = 0;
        arb_value_type current_ = 0;
        arb_value_type v = _pp_var_vec_v[node_indexi_];
        arb_value_type ica = 0;
        _pp_var_conductanceScale[i_] =  1.0;
        _pp_var_fopen0[i_] = _pp_var_k_fcond[i_]*_pp_var_l_fcond[i_];
        _pp_var_fopen[i_] = _pp_var_conductanceScale[i_]*_pp_var_fopen0[i_];
        _pp_var_g[i_] = _pp_var_conductance[i_]*_pp_var_fopen[i_];
        _pp_var_gion[i_] = _pp_var_gmax_actual[i_]*_pp_var_fopen[i_];
        ica = _pp_var_gion[i_]*(v-_pp_var_eca[i_]);
        if (_pp_var_timer[i_] > 0.75) {
            // discard transient
            _pp_var_vmin[i_] = v;
            _pp_var_vmax[i_] = v;
        }
        if (v > _pp_var_vmax[i_]) {
            _pp_var_vmax[i_] = v;
        }
        if (v < _pp_var_vmin[i_]) {
            _pp_var_vmin[i_] = v;
        }
        _pp_var_vamp[i_] = _pp_var_vmax[i_]-_pp_var_vmin[i_];
        if (_pp_var_stopped[i_] < 0.5&&_pp_var_totalTime[i_] > _pp_var_stopAfter[i_]) {
            _pp_var_stopped[i_] =  1.0;
            _pp_var_gmax_actual[i_] = _pp_var_gmax_actual[i_]*_pp_var_gmaxScaleAfterStop[i_];
            _pp_var_vmin[i_] = v;
            _pp_var_vmax[i_] = v;
        }
        if (_pp_var_stopped[i_]< 0.5&&_pp_var_timer[i_]< 0. && _pp_var_on_shadow_network[i_] < 0.5) {
            _pp_var_timer[j_] = _pp_var_timer[i_] = 1.5 + (fmod(i_*3.4231423, 0.455));
            _pp_var_vamp[i_] = _pp_var_vmax[i_]-_pp_var_vmin[i_];
            float f = 1 * (1 - _pp_var_totalTime[i_] / _pp_var_stopAfter[i_]);
            // Current cell (non-shadow)
            bool ok = true;

            if (_pp_var_vmax[i_] > 0) {
                // we don't want spikes
                _pp_var_gmax_actual[i_] = _pp_var_gmax_actual[i_] - f * 0.003 * 2 * (dt / 0.025);
            } else if (_pp_var_vamp[i_] < _pp_var_target_vamp_min[i_]) {
                if (_pp_var_k_q_lp[i_] > 0.77) {
                    _pp_var_gmax_actual[i_] = _pp_var_gmax_actual[i_] - f * 0.003 * 2 * (dt / 0.025);
                    ok = false;
                }
                else {
                    _pp_var_gmax_actual[i_] = _pp_var_gmax_actual[i_] + f * 0.003 * 2 * (dt / 0.025);
                    ok = false;
                }
            } else if (_pp_var_vamp[i_] > _pp_var_target_vamp_max[i_]&&_pp_var_vmax[i_] < 0.) {
                _pp_var_gmax_actual[i_] = _pp_var_gmax_actual[i_] - f * 0.001 * 2 * (dt / 0.025);
                ok = false;
            }

            float f2 = f * 2;
            if (ok) {
                // Other cell (shadow)
                if (_pp_var_target_vamp_max[j_] > 3.0) {
                    if (_pp_var_vamp[j_] < _pp_var_target_vamp_min[j_]) {
                        if (_pp_var_k_q_lp[j_] > 0.77) {
                            _pp_var_gmax_actual[i_] = _pp_var_gmax_actual[i_] - f2 * 0.003/2.0 * (dt / 0.025);
                        }
                        else {
                            _pp_var_gmax_actual[i_] = _pp_var_gmax_actual[i_] + f2 * 0.003/2.0 * (dt / 0.025);
                        }
                    }
                    if (_pp_var_vamp[j_] > _pp_var_target_vamp_max[j_]&&_pp_var_vmax[j_] < 0.) {
                        _pp_var_gmax_actual[i_] = _pp_var_gmax_actual[i_] - f2 * 0.001/2.0 * (dt / 0.025);
                    }
                }
                else {
                    if (_pp_var_vamp[j_] > _pp_var_target_vamp_max[j_]) {
                        _pp_var_gmax_actual[i_] = _pp_var_gmax_actual[i_] - f2 * 0.001/2.0 * (dt / 0.025);
                    }
                }
            }

            if (_pp_var_gmax_actual[i_] < 0.) {
                _pp_var_gmax_actual[i_] =  0.;
            }
            _pp_var_gmax_actual[j_] = _pp_var_gmax_actual[i_];
        }
        if (_pp_var_stopped[i_]> 0.5) {
            _pp_var_vamp[i_] = _pp_var_vmax[i_]-_pp_var_vmin[i_];
        }
        current_ = ica;
        conductivity_ = _pp_var_gion[i_];
        _pp_var_vec_g[node_indexi_] = fma(10.0*_pp_var_weight[i_], conductivity_, _pp_var_vec_g[node_indexi_]);
        _pp_var_vec_i[node_indexi_] = fma(10.0*_pp_var_weight[i_], current_, _pp_var_vec_i[node_indexi_]);
        _pp_var_ion_ca.current_density[ion_ca_indexi_] = fma(10.0*_pp_var_weight[i_], ica, _pp_var_ion_ca.current_density[ion_ca_indexi_]);
    }
}

static void write_ions(arb_mechanism_ppack* pp) {
}

static void apply_events(arb_mechanism_ppack*, arb_deliverable_event_stream*) {}

static void post_event(arb_mechanism_ppack*) {}

// Procedure definitions
[[maybe_unused]] static void rates(arb_mechanism_ppack* pp, int i_, arb_value_type v) {
    PPACK_IFACE_BLOCK;
    _pp_var_k_steadyState_x[i_] = _pp_var_k_steadyState_rate[i_]/( 1.0+exp( -((v-_pp_var_k_steadyState_midpoint[i_])/_pp_var_k_steadyState_scale[i_])));
    _pp_var_k_timeCourse_t[i_] = _pp_var_k_timeCourse_tau[i_];
    _pp_var_k_rateScale[i_] =  1.0;
    _pp_var_k_fcond[i_] = pow(_pp_var_k_q[i_], _pp_var_k_instances[i_]);
    _pp_var_k_inf[i_] = _pp_var_k_steadyState_x[i_];
    _pp_var_k_tauUnscaled[i_] = _pp_var_k_timeCourse_t[i_];
    _pp_var_k_tau[i_] = _pp_var_k_tauUnscaled[i_]/_pp_var_k_rateScale[i_];
    _pp_var_l_steadyState_x[i_] = _pp_var_l_steadyState_rate[i_]/( 1.0+exp( -((v-_pp_var_l_steadyState_midpoint[i_])/_pp_var_l_steadyState_scale[i_])));
    _pp_var_l_timeCourse_V[i_] = v/_pp_var_l_timeCourse_VOLT_SCALE[i_];
    _pp_var_l_timeCourse_t[i_] = _pp_var_l_timeCourse_TIME_SCALE[i_]*( 20.0*exp((_pp_var_l_timeCourse_V[i_]+ 160.0)* 0.033333333333333333)/( 1.0+exp((_pp_var_l_timeCourse_V[i_]+ 84.0)* 0.13698630136986301))+ 35.0);
    _pp_var_l_rateScale[i_] =  1.0;
    _pp_var_l_fcond[i_] = pow(_pp_var_l_q[i_], _pp_var_l_instances[i_]);
    _pp_var_l_inf[i_] = _pp_var_l_steadyState_x[i_];
    _pp_var_l_tauUnscaled[i_] = _pp_var_l_timeCourse_t[i_];
    _pp_var_l_tau[i_] = _pp_var_l_tauUnscaled[i_]/_pp_var_l_rateScale[i_];
    _pp_var_rate_k_q[i_] = (_pp_var_k_inf[i_]-_pp_var_k_q[i_])/_pp_var_k_tau[i_];
    _pp_var_rate_l_q[i_] = (_pp_var_l_inf[i_]-_pp_var_l_q[i_])/_pp_var_l_tau[i_];
    if (_pp_var_vamp[i_]< 3.0&&_pp_var_vmax[i_]-_pp_var_vmin[i_]< 3.0&&_pp_var_stopped[i_]< 0.5) {
        _pp_var_rate_k_q_lp[i_] =  0.01*(_pp_var_k_q[i_]-_pp_var_k_q_lp[i_]);
    }
    else {
        _pp_var_rate_k_q_lp[i_] =  0.;
    }
    _pp_var_rate_totalTime[i_] =  1.0;
}
#undef PPACK_IFACE_BLOCK
} // namespace kernel_calshadow
} // namespace EMC_catalogue
} // namespace arb

extern "C" {
  arb_mechanism_interface* make_arb_EMC_catalogue_calshadow_interface_multicore() {
    static arb_mechanism_interface result;
    result.partition_width = arb::EMC_catalogue::kernel_calshadow::simd_width_;
    result.backend = arb_backend_kind_cpu;
    result.alignment = arb::EMC_catalogue::kernel_calshadow::min_align_;
    result.init_mechanism = arb::EMC_catalogue::kernel_calshadow::init;
    result.compute_currents = arb::EMC_catalogue::kernel_calshadow::compute_currents;
    result.apply_events = arb::EMC_catalogue::kernel_calshadow::apply_events;
    result.advance_state = arb::EMC_catalogue::kernel_calshadow::advance_state;
    result.write_ions = arb::EMC_catalogue::kernel_calshadow::write_ions;
    result.post_event = arb::EMC_catalogue::kernel_calshadow::post_event;
    return &result;
  }}

