#include <iostream>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <arbor/mechanism_abi.h>
#include <arbor/math.hpp>
#include <arbor/profile/profiler.hpp>

namespace arb {
namespace IOU_catalogue {
namespace kernel_smol_dend {

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

static constexpr unsigned simd_width_ = 0;

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
[[maybe_unused]] auto* _pp_var_multiplicity      = pp->multiplicity;\
[[maybe_unused]] auto* _pp_var_weight            = pp->weight;\
[[maybe_unused]] auto& _pp_var_events            = pp->events;\
[[maybe_unused]] auto& _pp_var_mechanism_id      = pp->mechanism_id;\
[[maybe_unused]] auto& _pp_var_index_constraints = pp->index_constraints; \
                       /* cah */   \
[[maybe_unused]] auto* _pp_var_cah_r_q                           = pp->state_vars[0];\
[[maybe_unused]] auto* _pp_var_cah_rate_r_q                      = pp->parameters[0];\
[[maybe_unused]] auto* _pp_var_cah_fopen                         = pp->parameters[1];\
[[maybe_unused]] auto* _pp_var_cah_fopen0                        = pp->parameters[2];\
[[maybe_unused]] auto* _pp_var_cah_conductanceScale              = pp->parameters[3];\
[[maybe_unused]] auto* _pp_var_cah_r_reverseRate_rate            = pp->parameters[4];\
[[maybe_unused]] auto* _pp_var_cah_gion                          = pp->parameters[5];\
[[maybe_unused]] auto* _pp_var_cah_r_forwardRate_scale           = pp->parameters[6];\
[[maybe_unused]] auto* _pp_var_cah_r_forwardRate_midpoint        = pp->parameters[7];\
[[maybe_unused]] auto* _pp_var_cah_conductance                   = pp->parameters[8];\
[[maybe_unused]] auto* _pp_var_cah_r_beta                        = pp->parameters[9];\
[[maybe_unused]] auto* _pp_var_cah_r_fcond                       = pp->parameters[10];\
[[maybe_unused]] auto* _pp_var_cah_r_tau                         = pp->parameters[11];\
[[maybe_unused]] auto* _pp_var_cah_r_instances                   = pp->parameters[12];\
[[maybe_unused]] auto* _pp_var_cah_gmax                          = pp->parameters[13];\
[[maybe_unused]] auto* _pp_var_cah_temperature                   = pp->parameters[14];\
[[maybe_unused]] auto* _pp_var_cah_g                             = pp->parameters[15];\
[[maybe_unused]] auto* _pp_var_cah_r_reverseRate_scale           = pp->parameters[16];\
[[maybe_unused]] auto* _pp_var_cah_r_reverseRate_midpoint        = pp->parameters[17];\
[[maybe_unused]] auto* _pp_var_cah_r_forwardRate_rate            = pp->parameters[18];\
[[maybe_unused]] auto* _pp_var_cah_r_q10Settings_q10             = pp->parameters[19];\
[[maybe_unused]] auto* _pp_var_cah_eca                           = pp->parameters[20];\
[[maybe_unused]] auto* _pp_var_cah_r_q10Settings_fixedQ10        = pp->parameters[21];\
[[maybe_unused]] auto* _pp_var_cah_r_forwardRate_r               = pp->parameters[22];\
[[maybe_unused]] auto* _pp_var_cah_r_reverseRate_r               = pp->parameters[23];\
[[maybe_unused]] auto* _pp_var_cah_r_reverseRate_x               = pp->parameters[24];\
[[maybe_unused]] auto* _pp_var_cah_r_rateScale                   = pp->parameters[25];\
[[maybe_unused]] auto* _pp_var_cah_r_inf                         = pp->parameters[26];\
[[maybe_unused]] auto* _pp_var_cah_celsius                       = pp->parameters[27];\
[[maybe_unused]] auto* _pp_var_cah_r_alpha                       = pp->parameters[28];\
[[maybe_unused]] auto& _pp_var_cah_ion_ca                        = pp->ion_states[0];\
[[maybe_unused]] auto* _pp_var_cah_ion_ca_index                  = pp->ion_states[0].index; \
                        /* kca */   \
[[maybe_unused]] auto* _pp_var_kca_z_q                       = pp->state_vars[0+1];\
[[maybe_unused]] auto* _pp_var_kca_celsius                   = pp->parameters[0+29];\
[[maybe_unused]] auto* _pp_var_kca_z_instances               = pp->parameters[1+29];\
[[maybe_unused]] auto* _pp_var_kca_z_alpha                   = pp->parameters[2+29];\
[[maybe_unused]] auto* _pp_var_kca_conductance               = pp->parameters[3+29];\
[[maybe_unused]] auto* _pp_var_kca_z_forwardRate_r           = pp->parameters[4+29];\
[[maybe_unused]] auto* _pp_var_kca_conductanceScale          = pp->parameters[5+29];\
[[maybe_unused]] auto* _pp_var_kca_z_forwardRate_TIME_SCALE  = pp->parameters[6+29];\
[[maybe_unused]] auto* _pp_var_kca_z_reverseRate_TIME_SCALE  = pp->parameters[7+29];\
[[maybe_unused]] auto* _pp_var_kca_z_reverseRate_r           = pp->parameters[8+29];\
[[maybe_unused]] auto* _pp_var_kca_gmax                      = pp->parameters[9+29];\
[[maybe_unused]] auto* _pp_var_kca_temperature               = pp->parameters[10+29];\
[[maybe_unused]] auto* _pp_var_kca_z_beta                    = pp->parameters[11+29];\
[[maybe_unused]] auto* _pp_var_kca_fopen                     = pp->parameters[12+29];\
[[maybe_unused]] auto* _pp_var_kca_z_forwardRate_ca_norm     = pp->parameters[13+29];\
[[maybe_unused]] auto* _pp_var_kca_z_rateScale               = pp->parameters[14+29];\
[[maybe_unused]] auto* _pp_var_kca_gion                      = pp->parameters[15+29];\
[[maybe_unused]] auto* _pp_var_kca_z_forwardRate_CONC_SCALE  = pp->parameters[16+29];\
[[maybe_unused]] auto* _pp_var_kca_z_fcond                   = pp->parameters[17+29];\
[[maybe_unused]] auto* _pp_var_kca_z_tau                     = pp->parameters[18+29];\
[[maybe_unused]] auto* _pp_var_kca_fopen0                    = pp->parameters[19+29];\
[[maybe_unused]] auto* _pp_var_kca_z_inf                     = pp->parameters[20+29];\
[[maybe_unused]] auto* _pp_var_kca_ek                        = pp->parameters[21+29];\
[[maybe_unused]] auto* _pp_var_kca_g                         = pp->parameters[22+29];\
[[maybe_unused]] auto* _pp_var_kca_rate_z_q                  = pp->parameters[23+29];\
[[maybe_unused]] auto& _pp_var_kca_ion_ca                    = pp->ion_states[0+1];\
[[maybe_unused]] auto* _pp_var_kca_ion_ca_index              = pp->ion_states[0+1].index;\
                           /* h */   \
[[maybe_unused]] auto* _pp_var_h_n_q                         = pp->state_vars[0  + 1  + 1];\
[[maybe_unused]] auto* _pp_var_h_rate_n_q                    = pp->parameters[0  + 29 + 24];\
[[maybe_unused]] auto* _pp_var_h_n_tauUnscaled               = pp->parameters[1  + 29 + 24];\
[[maybe_unused]] auto* _pp_var_h_n_timeCourse_TIME_SCALE     = pp->parameters[2  + 29 + 24];\
[[maybe_unused]] auto* _pp_var_h_n_instances                 = pp->parameters[3  + 29 + 24];\
[[maybe_unused]] auto* _pp_var_h_n_timeCourse_VOLT_SCALE     = pp->parameters[4  + 29 + 24];\
[[maybe_unused]] auto* _pp_var_h_celsius                     = pp->parameters[5  + 29 + 24];\
[[maybe_unused]] auto* _pp_var_h_n_steadyState_rate          = pp->parameters[6  + 29 + 24];\
[[maybe_unused]] auto* _pp_var_h_n_tau                       = pp->parameters[7  + 29 + 24];\
[[maybe_unused]] auto* _pp_var_h_gmax                        = pp->parameters[8  + 29 + 24];\
[[maybe_unused]] auto* _pp_var_h_temperature                 = pp->parameters[9  + 29 + 24];\
[[maybe_unused]] auto* _pp_var_h_fopen                       = pp->parameters[10 + 29 + 24];\
[[maybe_unused]] auto* _pp_var_h_gion                        = pp->parameters[11 + 29 + 24];\
[[maybe_unused]] auto* _pp_var_h_eh                          = pp->parameters[12 + 29 + 24];\
[[maybe_unused]] auto* _pp_var_h_n_steadyState_x             = pp->parameters[13 + 29 + 24];\
[[maybe_unused]] auto* _pp_var_h_n_timeCourse_V              = pp->parameters[14 + 29 + 24];\
[[maybe_unused]] auto* _pp_var_h_n_timeCourse_t              = pp->parameters[15 + 29 + 24];\
[[maybe_unused]] auto* _pp_var_h_conductance                 = pp->parameters[16 + 29 + 24];\
[[maybe_unused]] auto* _pp_var_h_n_rateScale                 = pp->parameters[17 + 29 + 24];\
[[maybe_unused]] auto* _pp_var_h_n_fcond                     = pp->parameters[18 + 29 + 24];\
[[maybe_unused]] auto* _pp_var_h_n_inf                       = pp->parameters[19 + 29 + 24];\
[[maybe_unused]] auto* _pp_var_h_n_steadyState_midpoint      = pp->parameters[20 + 29 + 24];\
[[maybe_unused]] auto* _pp_var_h_conductanceScale            = pp->parameters[21 + 29 + 24];\
[[maybe_unused]] auto* _pp_var_h_n_steadyState_scale         = pp->parameters[22 + 29 + 24];\
[[maybe_unused]] auto* _pp_var_h_fopen0                      = pp->parameters[23 + 29 + 24];\
[[maybe_unused]] auto* _pp_var_h_g                           = pp->parameters[24 + 29 + 24];\
                    /* cacc */   \
[[maybe_unused]] auto* _pp_var_cacc_celsius                  = pp->parameters[0  + 29 + 24 + 25];\
[[maybe_unused]] auto* _pp_var_cacc_m_q                      = pp->parameters[1  + 29 + 24 + 25];\
[[maybe_unused]] auto* _pp_var_cacc_m_steadyState_x          = pp->parameters[2  + 29 + 24 + 25];\
[[maybe_unused]] auto* _pp_var_cacc_m_steadyState_VOLT_SCALE = pp->parameters[3  + 29 + 24 + 25];\
[[maybe_unused]] auto* _pp_var_cacc_m_SEC                    = pp->parameters[4  + 29 + 24 + 25];\
[[maybe_unused]] auto* _pp_var_cacc_ecl                      = pp->parameters[5  + 29 + 24 + 25];\
[[maybe_unused]] auto* _pp_var_cacc_conductance              = pp->parameters[6  + 29 + 24 + 25];\
[[maybe_unused]] auto* _pp_var_cacc_gmax                     = pp->parameters[7  + 29 + 24 + 25];\
[[maybe_unused]] auto* _pp_var_cacc_fopen                    = pp->parameters[8  + 29 + 24 + 25];\
[[maybe_unused]] auto* _pp_var_cacc_m_instances              = pp->parameters[9  + 29 + 24 + 25];\
[[maybe_unused]] auto* _pp_var_cacc_m_steadyState_CONC_SCALE = pp->parameters[10 + 29 + 24 + 25];\
[[maybe_unused]] auto* _pp_var_cacc_gion                     = pp->parameters[11 + 29 + 24 + 25];\
[[maybe_unused]] auto* _pp_var_cacc_m_inf                    = pp->parameters[12 + 29 + 24 + 25];\
[[maybe_unused]] auto* _pp_var_cacc_temperature              = pp->parameters[13 + 29 + 24 + 25];\
[[maybe_unused]] auto* _pp_var_cacc_m_tau                    = pp->parameters[14 + 29 + 24 + 25];\
[[maybe_unused]] auto* _pp_var_cacc_fopen0                   = pp->parameters[15 + 29 + 24 + 25];\
[[maybe_unused]] auto* _pp_var_cacc_m_steadyState_V          = pp->parameters[16 + 29 + 24 + 25];\
[[maybe_unused]] auto* _pp_var_cacc_m_steadyState_ca_conc    = pp->parameters[17 + 29 + 24 + 25];\
[[maybe_unused]] auto* _pp_var_cacc_m_fcond                  = pp->parameters[18 + 29 + 24 + 25];\
[[maybe_unused]] auto* _pp_var_cacc_conductanceScale         = pp->parameters[19 + 29 + 24 + 25];\
[[maybe_unused]] auto* _pp_var_cacc_g                        = pp->parameters[20 + 29 + 24 + 25];\
[[maybe_unused]] auto& _pp_var_cacc_ion_ca                   = pp->ion_states[0 + 1 + 0 + 1];\
[[maybe_unused]] auto* _pp_var_cacc_ion_ca_index             = pp->ion_states[0 + 1 + 0 + 1].index;\
 //End of IFACEBLOCK

// procedure prototypes
static void rates_cah(arb_mechanism_ppack* pp, int i_, arb_value_type v);
static void rates_kca(arb_mechanism_ppack* pp, int i_, arb_value_type v, arb_value_type cai);
static void rates_h(arb_mechanism_ppack* pp, int i_, arb_value_type v);
static void rates_cacc(arb_mechanism_ppack* pp, int i_, arb_value_type v, arb_value_type cai);

// interface methods


static void init(arb_mechanism_ppack* pp) {
    std::cout << "init dend CPU" << std::endl;
    PPACK_IFACE_BLOCK;
    for (arb_size_type i_ = 0; i_ < _pp_var_width; ++i_) {
//globals
        auto node_indexi_ = _pp_var_node_index[i_];
        arb_value_type v = _pp_var_vec_v[node_indexi_];
//CAH
        _pp_var_cah_eca[i_] =  120.0;
        _pp_var_cah_temperature[i_] = _pp_var_cah_celsius[i_]+ 273.14999999999998;
        rates_cah(pp, i_, v);
        rates_cah(pp, i_, v);
        _pp_var_cah_r_q[i_] = _pp_var_cah_r_inf[i_];
//KCA
        const arb_value_type cai_kca = _pp_var_kca_ion_ca.internal_concentration[_pp_var_kca_ion_ca_index[i_]];
        _pp_var_kca_ek[i_] =  -75.0;
        _pp_var_kca_temperature[i_] = _pp_var_kca_celsius[i_]+ 273.14999999999998;
        rates_kca(pp, i_, v, cai_kca);
        rates_kca(pp, i_, v, cai_kca);
        _pp_var_kca_z_q[i_] = _pp_var_kca_z_inf[i_];
//H
        _pp_var_h_eh[i_] =  -43.0;
        _pp_var_h_temperature[i_] = _pp_var_h_celsius[i_]+ 273.14999999999998;
        rates_h(pp, i_, v);
        rates_h(pp, i_, v);
        _pp_var_h_n_q[i_] = _pp_var_h_n_inf[i_];
//CACC
        const auto ion_ca_cacc_indexi_ = _pp_var_cacc_ion_ca_index[i_];
        const arb_value_type cai_cacc = _pp_var_cacc_ion_ca.internal_concentration[ion_ca_cacc_indexi_];
        _pp_var_cacc_ecl[i_] =  -45.0;
        _pp_var_cacc_temperature[i_] = _pp_var_cacc_celsius[i_]+ 273.14999999999998;
        rates_cacc(pp, i_, v, cai_cacc);
        rates_cacc(pp, i_, v, cai_cacc);
    }
}

static void compute_currents(arb_mechanism_ppack* pp) {
            PPACK_IFACE_BLOCK;
            for (arb_size_type i_ = 0; i_ < _pp_var_width; ++i_) {
                const auto node_indexi_ = _pp_var_node_index[i_];
                const arb_value_type v = _pp_var_vec_v[node_indexi_];
                const auto ion_ca_cah_indexi_ = _pp_var_cah_ion_ca_index[i_];
                const auto ion_ca_cacc_indexi_ = _pp_var_cacc_ion_ca_index[i_];
//CAH
                {
                    _pp_var_cah_conductanceScale[i_] = 1.0;
                    _pp_var_cah_fopen0[i_] = _pp_var_cah_r_fcond[i_];
                    _pp_var_cah_fopen[i_] = _pp_var_cah_conductanceScale[i_] * _pp_var_cah_fopen0[i_];
                    _pp_var_cah_g[i_] = _pp_var_cah_conductance[i_] * _pp_var_cah_fopen[i_];
                    _pp_var_cah_gion[i_] = _pp_var_cah_gmax[i_] * _pp_var_cah_fopen[i_];
                    const arb_value_type current_ = _pp_var_cah_gion[i_] * (v - _pp_var_cah_eca[i_]);
                    const arb_value_type conductivity_ = _pp_var_cah_gion[i_];
                    _pp_var_vec_g[node_indexi_] = fma(10.0 * _pp_var_weight[i_], conductivity_, _pp_var_vec_g[node_indexi_]);
                    _pp_var_vec_i[node_indexi_] = fma(10.0 * _pp_var_weight[i_], current_, _pp_var_vec_i[node_indexi_]);
                    _pp_var_cah_ion_ca.current_density[ion_ca_cah_indexi_] = fma(10.0 * _pp_var_weight[i_], current_, _pp_var_cah_ion_ca.current_density[ion_ca_cah_indexi_]);
                }
//KCA
                {
                    _pp_var_kca_conductanceScale[i_] = 1.0;
                    _pp_var_kca_fopen0[i_] = _pp_var_kca_z_fcond[i_];
                    _pp_var_kca_fopen[i_] = _pp_var_kca_conductanceScale[i_] * _pp_var_kca_fopen0[i_];
                    _pp_var_kca_g[i_] = _pp_var_kca_conductance[i_] * _pp_var_kca_fopen[i_];
                    _pp_var_kca_gion[i_] = _pp_var_kca_gmax[i_] * _pp_var_kca_fopen[i_];
                    const arb_value_type current_ = _pp_var_kca_gion[i_] * (v - _pp_var_kca_ek[i_]);
                    const arb_value_type conductivity_ = _pp_var_kca_gion[i_];
                    _pp_var_vec_g[node_indexi_] = fma(10.0 * _pp_var_weight[i_], conductivity_, _pp_var_vec_g[node_indexi_]);
                    _pp_var_vec_i[node_indexi_] = fma(10.0 * _pp_var_weight[i_], current_, _pp_var_vec_i[node_indexi_]);
                }
//H
                {
                    _pp_var_h_conductanceScale[i_] = 1.0;
                    _pp_var_h_fopen0[i_] = _pp_var_h_n_fcond[i_];
                    _pp_var_h_fopen[i_] = _pp_var_h_conductanceScale[i_] * _pp_var_h_fopen0[i_];
                    _pp_var_h_g[i_] = _pp_var_h_conductance[i_] * _pp_var_h_fopen[i_];
                    _pp_var_h_gion[i_] = _pp_var_h_gmax[i_] * _pp_var_h_fopen[i_];
                    const arb_value_type current_ = _pp_var_h_gion[i_] * (v - _pp_var_h_eh[i_]);
                    const arb_value_type conductivity_ = _pp_var_h_gion[i_];
                    _pp_var_vec_g[node_indexi_] = fma(10.0 * _pp_var_weight[i_], conductivity_, _pp_var_vec_g[node_indexi_]);
                    _pp_var_vec_i[node_indexi_] = fma(10.0 * _pp_var_weight[i_], current_, _pp_var_vec_i[node_indexi_]);
                }
//CACC
                {
                    rates_cacc(pp, i_, v, _pp_var_cacc_ion_ca.internal_concentration[ion_ca_cacc_indexi_]);
                    _pp_var_cacc_conductanceScale[i_] = 1.0;
                    _pp_var_cacc_fopen0[i_] = _pp_var_cacc_m_fcond[i_];
                    _pp_var_cacc_fopen[i_] = _pp_var_cacc_conductanceScale[i_] * _pp_var_cacc_fopen0[i_];
                    _pp_var_cacc_g[i_] = _pp_var_cacc_conductance[i_] * _pp_var_cacc_fopen[i_];
                    _pp_var_cacc_gion[i_] = _pp_var_cacc_gmax[i_] * _pp_var_cacc_fopen[i_];
                    const arb_value_type current_ = _pp_var_cacc_gion[i_] * (v - _pp_var_cacc_ecl[i_]);
                    const arb_value_type conductivity_ = _pp_var_cacc_gion[i_];
                    _pp_var_vec_g[node_indexi_] = fma(10.0 * _pp_var_weight[i_], conductivity_, _pp_var_vec_g[node_indexi_]);
                    _pp_var_vec_i[node_indexi_] = fma(10.0 * _pp_var_weight[i_], current_, _pp_var_vec_i[node_indexi_]);
                }
            }
        }

static void advance_state(arb_mechanism_ppack* pp) {
    PPACK_IFACE_BLOCK;
    for (arb_size_type i_ = 0; i_ < _pp_var_width; ++i_) {
//globals
        auto node_indexi_ = _pp_var_node_index[i_];
        arb_value_type dt = _pp_var_vec_dt[node_indexi_];
        arb_value_type v = _pp_var_vec_v[node_indexi_];
//CAH
        rates_cah(pp, i_, v);
        _pp_var_cah_r_q[i_] = _pp_var_cah_r_q[i_] + _pp_var_cah_rate_r_q[i_] * dt;
//KCA
        rates_kca(pp, i_, v, _pp_var_kca_ion_ca.internal_concentration[_pp_var_kca_ion_ca_index[i_]]);
        _pp_var_kca_z_q[i_] = _pp_var_kca_z_q[i_] + _pp_var_kca_rate_z_q[i_] * dt;
//H
        rates_h(pp, i_, v);
        _pp_var_h_n_q[i_] = _pp_var_h_n_q[i_] + _pp_var_h_rate_n_q[i_] * dt;
//CACC
        //empty
    }
}

static void write_ions(arb_mechanism_ppack* pp) {}

static void apply_events(arb_mechanism_ppack*) {}

static void post_event(arb_mechanism_ppack*) {}

// Procedure definitions
static void rates_cah(arb_mechanism_ppack* pp, int i_, arb_value_type v) {
    PPACK_IFACE_BLOCK;
    _pp_var_cah_r_forwardRate_r[i_] = _pp_var_cah_r_forwardRate_rate[i_]/( 1.0+exp( -((v-_pp_var_cah_r_forwardRate_midpoint[i_])/_pp_var_cah_r_forwardRate_scale[i_])));
    _pp_var_cah_r_reverseRate_x[i_] = (v-_pp_var_cah_r_reverseRate_midpoint[i_])/_pp_var_cah_r_reverseRate_scale[i_];
    if (_pp_var_cah_r_reverseRate_x[i_]!= 0.) {
        _pp_var_cah_r_reverseRate_r[i_] = _pp_var_cah_r_reverseRate_rate[i_]*_pp_var_cah_r_reverseRate_x[i_]/( 1.0-exp( -_pp_var_cah_r_reverseRate_x[i_]));
    }
    else {
        if (_pp_var_cah_r_reverseRate_x[i_]== 0.) {
            _pp_var_cah_r_reverseRate_r[i_] = _pp_var_cah_r_reverseRate_rate[i_];
        }
    }
    _pp_var_cah_r_q10Settings_q10[i_] = _pp_var_cah_r_q10Settings_fixedQ10[i_];
    _pp_var_cah_r_rateScale[i_] = _pp_var_cah_r_q10Settings_q10[i_];
    _pp_var_cah_r_alpha[i_] = _pp_var_cah_r_forwardRate_r[i_];
    _pp_var_cah_r_beta[i_] = _pp_var_cah_r_reverseRate_r[i_];
    _pp_var_cah_r_fcond[i_] = pow(_pp_var_cah_r_q[i_], _pp_var_cah_r_instances[i_]);
    _pp_var_cah_r_inf[i_] = _pp_var_cah_r_alpha[i_]/(_pp_var_cah_r_alpha[i_]+_pp_var_cah_r_beta[i_]);
    _pp_var_cah_r_tau[i_] =  1.0/((_pp_var_cah_r_alpha[i_]+_pp_var_cah_r_beta[i_])*_pp_var_cah_r_rateScale[i_]);
    _pp_var_cah_rate_r_q[i_] = (_pp_var_cah_r_inf[i_]-_pp_var_cah_r_q[i_])/_pp_var_cah_r_tau[i_];
}
static void rates_kca(arb_mechanism_ppack* pp, int i_, arb_value_type v, arb_value_type cai) {
            PPACK_IFACE_BLOCK;
            arb_value_type caConc;
            caConc = cai;
            _pp_var_kca_z_forwardRate_ca_norm[i_] =  2.0000000000000002e-05*caConc/_pp_var_kca_z_forwardRate_CONC_SCALE[i_];
            if (_pp_var_kca_z_forwardRate_ca_norm[i_]> 0.01) {
                _pp_var_kca_z_forwardRate_r[i_] =  0.01/_pp_var_kca_z_forwardRate_TIME_SCALE[i_];
            }
            else {
                if (_pp_var_kca_z_forwardRate_ca_norm[i_]<= 0.01) {
                    _pp_var_kca_z_forwardRate_r[i_] = _pp_var_kca_z_forwardRate_ca_norm[i_]/_pp_var_kca_z_forwardRate_TIME_SCALE[i_];
                }
            }
            _pp_var_kca_z_reverseRate_r[i_] =  0.014999999999999999/_pp_var_kca_z_reverseRate_TIME_SCALE[i_];
            _pp_var_kca_z_rateScale[i_] =  1.0;
            _pp_var_kca_z_alpha[i_] = _pp_var_kca_z_forwardRate_r[i_];
            _pp_var_kca_z_beta[i_] = _pp_var_kca_z_reverseRate_r[i_];
            _pp_var_kca_z_fcond[i_] = pow(_pp_var_kca_z_q[i_], _pp_var_kca_z_instances[i_]);
            _pp_var_kca_z_inf[i_] = _pp_var_kca_z_alpha[i_]/(_pp_var_kca_z_alpha[i_]+_pp_var_kca_z_beta[i_]);
            _pp_var_kca_z_tau[i_] =  1.0/((_pp_var_kca_z_alpha[i_]+_pp_var_kca_z_beta[i_])*_pp_var_kca_z_rateScale[i_]);
            _pp_var_kca_rate_z_q[i_] = (_pp_var_kca_z_inf[i_]-_pp_var_kca_z_q[i_])/_pp_var_kca_z_tau[i_];
        }
static void rates_h(arb_mechanism_ppack* pp, int i_, arb_value_type v) {
    PPACK_IFACE_BLOCK;
    _pp_var_h_n_steadyState_x[i_] = _pp_var_h_n_steadyState_rate[i_]/( 1.0+exp( -((v-_pp_var_h_n_steadyState_midpoint[i_])/_pp_var_h_n_steadyState_scale[i_])));
    _pp_var_h_n_timeCourse_V[i_] = v/_pp_var_h_n_timeCourse_VOLT_SCALE[i_];
    _pp_var_h_n_timeCourse_t[i_] = _pp_var_h_n_timeCourse_TIME_SCALE[i_]/(exp( -0.085999999999999993*_pp_var_h_n_timeCourse_V[i_]- 14.6)+exp( 0.070000000000000007*_pp_var_h_n_timeCourse_V[i_]- 1.8700000000000001));
    _pp_var_h_n_rateScale[i_] =  1.0;
    _pp_var_h_n_fcond[i_] = pow(_pp_var_h_n_q[i_], _pp_var_h_n_instances[i_]);
    _pp_var_h_n_inf[i_] = _pp_var_h_n_steadyState_x[i_];
    _pp_var_h_n_tauUnscaled[i_] = _pp_var_h_n_timeCourse_t[i_];
    _pp_var_h_n_tau[i_] = _pp_var_h_n_tauUnscaled[i_]/_pp_var_h_n_rateScale[i_];
    _pp_var_h_rate_n_q[i_] = (_pp_var_h_n_inf[i_]-_pp_var_h_n_q[i_])/_pp_var_h_n_tau[i_];
}
static void rates_cacc(arb_mechanism_ppack* pp, int i_, arb_value_type v, arb_value_type cai) {
    PPACK_IFACE_BLOCK;
    arb_value_type caConc;
    caConc = cai;
    _pp_var_cacc_m_steadyState_V[i_] = v/_pp_var_cacc_m_steadyState_VOLT_SCALE[i_];
    _pp_var_cacc_m_steadyState_ca_conc[i_] = caConc/_pp_var_cacc_m_steadyState_CONC_SCALE[i_];
    _pp_var_cacc_m_steadyState_x[i_] =  1.0/( 1.0+exp(( 0.00036999999999999999-_pp_var_cacc_m_steadyState_ca_conc[i_])* 11.111111111111111));
    _pp_var_cacc_m_inf[i_] = _pp_var_cacc_m_steadyState_x[i_];
    _pp_var_cacc_m_tau[i_] =  0.;
    _pp_var_cacc_m_q[i_] = _pp_var_cacc_m_inf[i_];
    _pp_var_cacc_m_fcond[i_] = pow(_pp_var_cacc_m_q[i_], _pp_var_cacc_m_instances[i_]);
}

#undef PPACK_IFACE_BLOCK
} // namespace kernel_cah
} // namespace smol_catalogue
} // namespace arb

extern "C" {
  arb_mechanism_interface* make_arb_IOU_catalogue_smol_dend_interface_multicore() {
    static arb_mechanism_interface result;
    result.partition_width = arb::IOU_catalogue::kernel_smol_dend::simd_width_;
    result.backend=arb_backend_kind_cpu;
    result.alignment=1;
    result.init_mechanism=(arb_mechanism_method)arb::IOU_catalogue::kernel_smol_dend::init;
    result.compute_currents=(arb_mechanism_method)arb::IOU_catalogue::kernel_smol_dend::compute_currents;
    result.apply_events=(arb_mechanism_method)arb::IOU_catalogue::kernel_smol_dend::apply_events;
    result.advance_state=(arb_mechanism_method)arb::IOU_catalogue::kernel_smol_dend::advance_state;
    result.write_ions=(arb_mechanism_method)arb::IOU_catalogue::kernel_smol_dend::write_ions;
    result.post_event=(arb_mechanism_method)arb::IOU_catalogue::kernel_smol_dend::post_event;
    return &result;
  }}

