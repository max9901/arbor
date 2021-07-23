//
// Created by max on 08-07-21.
//
#include <iostream>

#include <arbor/gpu/gpu_common.hpp>
#include <arbor/gpu/math_cu.hpp>
#include <arbor/gpu/reduce_by_key.hpp>
#include <arbor/mechanism_abi.h>

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
                       /* cah */   \
auto* _pp_var_cah_r_q                        __attribute__((unused))  = params_.state_vars[0];\
auto* _pp_var_cah_rate_r_q                   __attribute__((unused))  = params_.parameters[0];\
auto* _pp_var_cah_fopen                      __attribute__((unused))  = params_.parameters[1];\
auto* _pp_var_cah_fopen0                     __attribute__((unused))  = params_.parameters[2];\
auto* _pp_var_cah_conductanceScale           __attribute__((unused))  = params_.parameters[3];\
auto* _pp_var_cah_r_reverseRate_rate         __attribute__((unused))  = params_.parameters[4];\
auto* _pp_var_cah_gion                       __attribute__((unused))  = params_.parameters[5];\
auto* _pp_var_cah_r_forwardRate_scale        __attribute__((unused))  = params_.parameters[6];\
auto* _pp_var_cah_r_forwardRate_midpoint     __attribute__((unused))  = params_.parameters[7];\
auto* _pp_var_cah_conductance                __attribute__((unused))  = params_.parameters[8];\
auto* _pp_var_cah_r_beta                     __attribute__((unused))  = params_.parameters[9];\
auto* _pp_var_cah_r_fcond                    __attribute__((unused))  = params_.parameters[10];\
auto* _pp_var_cah_r_tau                      __attribute__((unused))  = params_.parameters[11];\
auto* _pp_var_cah_r_instances                __attribute__((unused))  = params_.parameters[12];\
auto* _pp_var_cah_gmax                       __attribute__((unused))  = params_.parameters[13];\
auto* _pp_var_cah_temperature                __attribute__((unused))  = params_.parameters[14];\
auto* _pp_var_cah_g                          __attribute__((unused))  = params_.parameters[15];\
auto* _pp_var_cah_r_reverseRate_scale        __attribute__((unused))  = params_.parameters[16];\
auto* _pp_var_cah_r_reverseRate_midpoint     __attribute__((unused))  = params_.parameters[17];\
auto* _pp_var_cah_r_forwardRate_rate         __attribute__((unused))  = params_.parameters[18];\
auto* _pp_var_cah_r_q10Settings_q10          __attribute__((unused))  = params_.parameters[19];\
auto* _pp_var_cah_eca                        __attribute__((unused))  = params_.parameters[20];\
auto* _pp_var_cah_r_q10Settings_fixedQ10     __attribute__((unused))  = params_.parameters[21];\
auto* _pp_var_cah_r_forwardRate_r            __attribute__((unused))  = params_.parameters[22];\
auto* _pp_var_cah_r_reverseRate_r            __attribute__((unused))  = params_.parameters[23];\
auto* _pp_var_cah_r_reverseRate_x            __attribute__((unused))  = params_.parameters[24];\
auto* _pp_var_cah_r_rateScale                __attribute__((unused))  = params_.parameters[25];\
auto* _pp_var_cah_r_inf                      __attribute__((unused))  = params_.parameters[26];\
auto* _pp_var_cah_celsius                    __attribute__((unused))  = params_.parameters[27];\
auto* _pp_var_cah_r_alpha                    __attribute__((unused))  = params_.parameters[28];\
auto& _pp_var_cah_ion_ca                     __attribute__((unused))  = params_.ion_states[0];\
auto* _pp_var_cah_ion_ca_index               __attribute__((unused))  = params_.ion_states[0].index; \
                        /* kca */   \
auto* _pp_var_kca_z_q                      __attribute__((unused))  = params_.state_vars[0+1];\
auto* _pp_var_kca_celsius                  __attribute__((unused))  = params_.parameters[0+29];\
auto* _pp_var_kca_z_instances              __attribute__((unused))  = params_.parameters[1+29];\
auto* _pp_var_kca_z_alpha                  __attribute__((unused))  = params_.parameters[2+29];\
auto* _pp_var_kca_conductance              __attribute__((unused))  = params_.parameters[3+29];\
auto* _pp_var_kca_z_forwardRate_r          __attribute__((unused))  = params_.parameters[4+29];\
auto* _pp_var_kca_conductanceScale         __attribute__((unused))  = params_.parameters[5+29];\
auto* _pp_var_kca_z_forwardRate_TIME_SCALE __attribute__((unused))  = params_.parameters[6+29];\
auto* _pp_var_kca_z_reverseRate_TIME_SCALE __attribute__((unused))  = params_.parameters[7+29];\
auto* _pp_var_kca_z_reverseRate_r          __attribute__((unused))  = params_.parameters[8+29];\
auto* _pp_var_kca_gmax                     __attribute__((unused))  = params_.parameters[9+29];\
auto* _pp_var_kca_temperature              __attribute__((unused))  = params_.parameters[10+29];\
auto* _pp_var_kca_z_beta                   __attribute__((unused))  = params_.parameters[11+29];\
auto* _pp_var_kca_fopen                    __attribute__((unused))  = params_.parameters[12+29];\
auto* _pp_var_kca_z_forwardRate_ca_norm    __attribute__((unused))  = params_.parameters[13+29];\
auto* _pp_var_kca_z_rateScale              __attribute__((unused))  = params_.parameters[14+29];\
auto* _pp_var_kca_gion                     __attribute__((unused))  = params_.parameters[15+29];\
auto* _pp_var_kca_z_forwardRate_CONC_SCALE __attribute__((unused))  = params_.parameters[16+29];\
auto* _pp_var_kca_z_fcond                  __attribute__((unused))  = params_.parameters[17+29];\
auto* _pp_var_kca_z_tau                    __attribute__((unused))  = params_.parameters[18+29];\
auto* _pp_var_kca_fopen0                   __attribute__((unused))  = params_.parameters[19+29];\
auto* _pp_var_kca_z_inf                    __attribute__((unused))  = params_.parameters[20+29];\
auto* _pp_var_kca_ek                       __attribute__((unused))  = params_.parameters[21+29];\
auto* _pp_var_kca_g                        __attribute__((unused))  = params_.parameters[22+29];\
auto* _pp_var_kca_rate_z_q                 __attribute__((unused))  = params_.parameters[23+29];\
auto& _pp_var_kca_ion_ca                   __attribute__((unused))  = params_.ion_states[0+1];\
auto* _pp_var_kca_ion_ca_index             __attribute__((unused))  = params_.ion_states[0+1].index;\
                           /* h */   \
auto* _pp_var_h_n_q                     __attribute__((unused))  = params_.state_vars[0  + 1  + 1];\
auto* _pp_var_h_rate_n_q                __attribute__((unused))  = params_.parameters[0  + 29 + 24];\
auto* _pp_var_h_n_tauUnscaled           __attribute__((unused))  = params_.parameters[1  + 29 + 24];\
auto* _pp_var_h_n_timeCourse_TIME_SCALE __attribute__((unused))  = params_.parameters[2  + 29 + 24];\
auto* _pp_var_h_n_instances             __attribute__((unused))  = params_.parameters[3  + 29 + 24];\
auto* _pp_var_h_n_timeCourse_VOLT_SCALE __attribute__((unused))  = params_.parameters[4  + 29 + 24];\
auto* _pp_var_h_celsius                 __attribute__((unused))  = params_.parameters[5  + 29 + 24];\
auto* _pp_var_h_n_steadyState_rate      __attribute__((unused))  = params_.parameters[6  + 29 + 24];\
auto* _pp_var_h_n_tau                   __attribute__((unused))  = params_.parameters[7  + 29 + 24];\
auto* _pp_var_h_gmax                    __attribute__((unused))  = params_.parameters[8  + 29 + 24];\
auto* _pp_var_h_temperature             __attribute__((unused))  = params_.parameters[9  + 29 + 24];\
auto* _pp_var_h_fopen                   __attribute__((unused))  = params_.parameters[10 + 29 + 24];\
auto* _pp_var_h_gion                    __attribute__((unused))  = params_.parameters[11 + 29 + 24];\
auto* _pp_var_h_eh                      __attribute__((unused))  = params_.parameters[12 + 29 + 24];\
auto* _pp_var_h_n_steadyState_x         __attribute__((unused))  = params_.parameters[13 + 29 + 24];\
auto* _pp_var_h_n_timeCourse_V          __attribute__((unused))  = params_.parameters[14 + 29 + 24];\
auto* _pp_var_h_n_timeCourse_t          __attribute__((unused))  = params_.parameters[15 + 29 + 24];\
auto* _pp_var_h_conductance             __attribute__((unused))  = params_.parameters[16 + 29 + 24];\
auto* _pp_var_h_n_rateScale             __attribute__((unused))  = params_.parameters[17 + 29 + 24];\
auto* _pp_var_h_n_fcond                 __attribute__((unused))  = params_.parameters[18 + 29 + 24];\
auto* _pp_var_h_n_inf                   __attribute__((unused))  = params_.parameters[19 + 29 + 24];\
auto* _pp_var_h_n_steadyState_midpoint  __attribute__((unused))  = params_.parameters[20 + 29 + 24];\
auto* _pp_var_h_conductanceScale        __attribute__((unused))  = params_.parameters[21 + 29 + 24];\
auto* _pp_var_h_n_steadyState_scale     __attribute__((unused))  = params_.parameters[22 + 29 + 24];\
auto* _pp_var_h_fopen0                  __attribute__((unused))  = params_.parameters[23 + 29 + 24];\
auto* _pp_var_h_g                       __attribute__((unused))  = params_.parameters[24 + 29 + 24];\
                    /* cacc */   \
auto* _pp_var_cacc_celsius                  __attribute__((unused)) = params_.parameters[0  + 29 + 24 + 25];\
auto* _pp_var_cacc_m_q                      __attribute__((unused)) = params_.parameters[1  + 29 + 24 + 25];\
auto* _pp_var_cacc_m_steadyState_x          __attribute__((unused)) = params_.parameters[2  + 29 + 24 + 25];\
auto* _pp_var_cacc_m_steadyState_VOLT_SCALE __attribute__((unused)) = params_.parameters[3  + 29 + 24 + 25];\
auto* _pp_var_cacc_m_SEC                    __attribute__((unused)) = params_.parameters[4  + 29 + 24 + 25];\
auto* _pp_var_cacc_ecl                      __attribute__((unused)) = params_.parameters[5  + 29 + 24 + 25];\
auto* _pp_var_cacc_conductance              __attribute__((unused)) = params_.parameters[6  + 29 + 24 + 25];\
auto* _pp_var_cacc_gmax                     __attribute__((unused)) = params_.parameters[7  + 29 + 24 + 25];\
auto* _pp_var_cacc_fopen                    __attribute__((unused)) = params_.parameters[8  + 29 + 24 + 25];\
auto* _pp_var_cacc_m_instances              __attribute__((unused)) = params_.parameters[9  + 29 + 24 + 25];\
auto* _pp_var_cacc_m_steadyState_CONC_SCALE __attribute__((unused)) = params_.parameters[10 + 29 + 24 + 25];\
auto* _pp_var_cacc_gion                     __attribute__((unused)) = params_.parameters[11 + 29 + 24 + 25];\
auto* _pp_var_cacc_m_inf                    __attribute__((unused)) = params_.parameters[12 + 29 + 24 + 25];\
auto* _pp_var_cacc_temperature              __attribute__((unused)) = params_.parameters[13 + 29 + 24 + 25];\
auto* _pp_var_cacc_m_tau                    __attribute__((unused)) = params_.parameters[14 + 29 + 24 + 25];\
auto* _pp_var_cacc_fopen0                   __attribute__((unused)) = params_.parameters[15 + 29 + 24 + 25];\
auto* _pp_var_cacc_m_steadyState_V          __attribute__((unused)) = params_.parameters[16 + 29 + 24 + 25];\
auto* _pp_var_cacc_m_steadyState_ca_conc    __attribute__((unused)) = params_.parameters[17 + 29 + 24 + 25];\
auto* _pp_var_cacc_m_fcond                  __attribute__((unused)) = params_.parameters[18 + 29 + 24 + 25];\
auto* _pp_var_cacc_conductanceScale         __attribute__((unused)) = params_.parameters[19 + 29 + 24 + 25];\
auto* _pp_var_cacc_g                        __attribute__((unused)) = params_.parameters[20 + 29 + 24 + 25];\
auto& _pp_var_cacc_ion_ca                   __attribute__((unused)) = params_.ion_states[0 + 1 + 0 + 1];\
auto* _pp_var_cacc_ion_ca_index             __attribute__((unused)) = params_.ion_states[0 + 1 + 0 + 1].index;\
//End of IFACEBLOCK

namespace {

using ::arb::gpu::exprelr;
using ::arb::gpu::safeinv;
using ::arb::gpu::min;
using ::arb::gpu::max;

__device__ void rates_cah(arb_mechanism_ppack params_, int i_, arb_value_type v) {
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
__device__ void rates_kca(arb_mechanism_ppack params_, int i_, arb_value_type v, arb_value_type cai) {
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
__device__ void rates_h(arb_mechanism_ppack params_, int i_, arb_value_type v) {
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
__device__ void rates_cacc(arb_mechanism_ppack params_, int i_, arb_value_type v, arb_value_type cai) {
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

__global__ void init(arb_mechanism_ppack params_) {
    PPACK_IFACE_BLOCK;
    const size_t i_ = threadIdx.x + blockDim.x*blockIdx.x;
    if (i_<params_.width) {
//globals
        auto node_indexi_ = _pp_var_node_index[i_];
        arb_value_type v = _pp_var_vec_v[node_indexi_];
//CAH
        _pp_var_cah_eca[i_] =  120.0;
        _pp_var_cah_temperature[i_] = _pp_var_cah_celsius[i_]+ 273.14999999999998;
        rates_cah(params_, i_, v);
        rates_cah(params_, i_, v);
        _pp_var_cah_r_q[i_] = _pp_var_cah_r_inf[i_];
//KCA
        const arb_value_type cai_kca = _pp_var_kca_ion_ca.internal_concentration[_pp_var_kca_ion_ca_index[i_]];
        _pp_var_kca_ek[i_] =  -75.0;
        _pp_var_kca_temperature[i_] = _pp_var_kca_celsius[i_]+ 273.14999999999998;
        rates_kca(params_, i_, v, cai_kca);
        rates_kca(params_, i_, v, cai_kca);
        _pp_var_kca_z_q[i_] = _pp_var_kca_z_inf[i_];
//H
        _pp_var_h_eh[i_] =  -43.0;
        _pp_var_h_temperature[i_] = _pp_var_h_celsius[i_]+ 273.14999999999998;
        rates_h(params_, i_, v);
        rates_h(params_, i_, v);
        _pp_var_h_n_q[i_] = _pp_var_h_n_inf[i_];
//CACC
        const auto ion_ca_cacc_indexi_ = _pp_var_cacc_ion_ca_index[i_];
        const arb_value_type cai_cacc = _pp_var_cacc_ion_ca.internal_concentration[ion_ca_cacc_indexi_];
        _pp_var_cacc_ecl[i_] =  -45.0;
        _pp_var_cacc_temperature[i_] = _pp_var_cacc_celsius[i_]+ 273.14999999999998;
        rates_cacc(params_, i_, v, cai_cacc);
        rates_cacc(params_, i_, v, cai_cacc);
    }
}
__global__ void compute_currents(arb_mechanism_ppack params_) {
    PPACK_IFACE_BLOCK;
    const size_t i_ = threadIdx.x + blockDim.x*blockIdx.x;
    if (i_<params_.width) {
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
            rates_cacc(params_, i_, v, _pp_var_cacc_ion_ca.internal_concentration[ion_ca_cacc_indexi_]);
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
__global__ void advance_state(arb_mechanism_ppack params_) {
    PPACK_IFACE_BLOCK;
    const size_t i_ = threadIdx.x + blockDim.x*blockIdx.x;
    if (i_<params_.width) {
//globals
        auto node_indexi_ = _pp_var_node_index[i_];
        arb_value_type dt = _pp_var_vec_dt[node_indexi_];
        arb_value_type v = _pp_var_vec_v[node_indexi_];
//CAH
        rates_cah(params_, i_, v);
        _pp_var_cah_r_q[i_] = _pp_var_cah_r_q[i_] + _pp_var_cah_rate_r_q[i_] * dt;
//KCA
        rates_kca(params_, i_, v, _pp_var_kca_ion_ca.internal_concentration[_pp_var_kca_ion_ca_index[i_]]);
        _pp_var_kca_z_q[i_] = _pp_var_kca_z_q[i_] + _pp_var_kca_rate_z_q[i_] * dt;
//H
        rates_h(params_, i_, v);
        _pp_var_h_n_q[i_] = _pp_var_h_n_q[i_] + _pp_var_h_rate_n_q[i_] * dt;
//CACC
        //empty
    }
}

} // namespace

void mechanism_smol_dend_gpu_init_(arb_mechanism_ppack* p) {
    std::cout << "init dend GPU" << std::endl;
    unsigned block_dim = 128;
    unsigned grid_dim = ::arb::gpu::impl::block_count(p->width, block_dim);
    init<<<grid_dim, block_dim>>>(*p);
}
void mechanism_smol_dend_gpu_compute_currents_(arb_mechanism_ppack* p) {
    unsigned block_dim = 128;
    unsigned grid_dim = ::arb::gpu::impl::block_count(p->width, block_dim);
    compute_currents<<<grid_dim, block_dim>>>(*p);
}
void mechanism_smol_dend_gpu_advance_state_(arb_mechanism_ppack* p) {
    unsigned block_dim = 128;
    unsigned grid_dim = ::arb::gpu::impl::block_count(p->width, block_dim);
    advance_state<<<grid_dim, block_dim>>>(*p);
}

void mechanism_smol_dend_gpu_write_ions_(arb_mechanism_ppack* p) {}
void mechanism_smol_dend_gpu_post_event_(arb_mechanism_ppack* p) {}
void mechanism_smol_dend_gpu_apply_events_(arb_mechanism_ppack* p) {}

} // namespace IOU_catalogue
} // namespace arb
