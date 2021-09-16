//
// Created by max on 08-07-21.
//
#include <iostream>

#include <arbor/gpu/gpu_common.hpp>
#include <arbor/mechanism_abi.h>

namespace arb {
    namespace IOU_catalogue {
        
#define PPACK_IFACE_BLOCK \
auto* _pp_var_vec_i               __attribute__((unused)) = params_.vec_i;\
auto* _pp_var_vec_v               __attribute__((unused)) = params_.vec_v;\
auto* _pp_var_vec_dt              __attribute__((unused)) = params_.vec_dt;\
auto& _pp_var_gap_junctions       __attribute__((unused)) = params_.gap_junctions;\
auto& _pp_var_gap_junction_width  __attribute__((unused)) = params_.gap_junction_width; \
auto& _pp_var_Cm                  __attribute__((unused)) = params_.globals[0];\
auto& _pp_var_D_CA                __attribute__((unused)) = params_.globals[1];\
auto& _pp_var_D_CL                __attribute__((unused)) = params_.globals[2];\
auto& _pp_var_F                   __attribute__((unused)) = params_.globals[3];\
auto& _pp_var_R                   __attribute__((unused)) = params_.globals[4];\
auto& _pp_var_T                   __attribute__((unused)) = params_.globals[5];\
auto& _pp_var_dend1_diam          __attribute__((unused)) = params_.globals[6];\
auto& _pp_var_dend1_half_G_CA     __attribute__((unused)) = params_.globals[7];\
auto& _pp_var_dend1_half_G_CL     __attribute__((unused)) = params_.globals[8];\
auto& _pp_var_dend1_length        __attribute__((unused)) = params_.globals[9];\
auto& _pp_var_dend1_n_CA          __attribute__((unused)) = params_.globals[10];\
auto& _pp_var_dend1_n_CL          __attribute__((unused)) = params_.globals[11];\
auto& _pp_var_dend1_n_K           __attribute__((unused)) = params_.globals[12];\
auto& _pp_var_dend2_diam          __attribute__((unused)) = params_.globals[13];\
auto& _pp_var_dend2_half_G_CA     __attribute__((unused)) = params_.globals[14];\
auto& _pp_var_dend2_half_G_CL     __attribute__((unused)) = params_.globals[15];\
auto& _pp_var_dend2_length        __attribute__((unused)) = params_.globals[16];\
auto& _pp_var_dend2_n_CA          __attribute__((unused)) = params_.globals[17];\
auto& _pp_var_dend2_n_CL          __attribute__((unused)) = params_.globals[18];\
auto& _pp_var_dend2_n_K           __attribute__((unused)) = params_.globals[19];\
auto& _pp_var_head1_E_K_NA        __attribute__((unused)) = params_.globals[20];\
auto& _pp_var_head1_E_head_ca     __attribute__((unused)) = params_.globals[21];\
auto& _pp_var_head1_Gleak         __attribute__((unused)) = params_.globals[22];\
auto& _pp_var_head1_I_max         __attribute__((unused)) = params_.globals[23];\
auto& _pp_var_head1_Ip            __attribute__((unused)) = params_.globals[24];\
auto& _pp_var_head1_Kp            __attribute__((unused)) = params_.globals[25];\
auto& _pp_var_head1_beta          __attribute__((unused)) = params_.globals[26];\
auto& _pp_var_head1_diam          __attribute__((unused)) = params_.globals[27];\
auto& _pp_var_head1_eta_gaba      __attribute__((unused)) = params_.globals[28];\
auto& _pp_var_head1_gamma         __attribute__((unused)) = params_.globals[29];\
auto& _pp_var_head1_gbar_gaba     __attribute__((unused)) = params_.globals[30];\
auto& _pp_var_head1_gbar_nmda     __attribute__((unused)) = params_.globals[31];\
auto& _pp_var_head1_ggap          __attribute__((unused)) = params_.globals[32];\
auto& _pp_var_head1_kCL           __attribute__((unused)) = params_.globals[33];\
auto& _pp_var_head1_kK            __attribute__((unused)) = params_.globals[34];\
auto& _pp_var_head1_kb            __attribute__((unused)) = params_.globals[35];\
auto& _pp_var_head1_kf            __attribute__((unused)) = params_.globals[36];\
auto& _pp_var_head1_length        __attribute__((unused)) = params_.globals[37];\
auto& _pp_var_head1_nTB           __attribute__((unused)) = params_.globals[38];\
auto& _pp_var_head1_n_K           __attribute__((unused)) = params_.globals[39];\
auto& _pp_var_head1_n_Mg_out      __attribute__((unused)) = params_.globals[40];\
auto& _pp_var_head1_sigma_gaba    __attribute__((unused)) = params_.globals[41];\
auto& _pp_var_head1_tau_d         __attribute__((unused)) = params_.globals[42];\
auto& _pp_var_head1_tau_decay_nmda __attribute__((unused)) = params_.globals[43];\
auto& _pp_var_head1_tau_gaba0     __attribute__((unused)) = params_.globals[44];\
auto& _pp_var_head1_tau_r         __attribute__((unused)) = params_.globals[45];\
auto& _pp_var_head1_tau_rise_nmda __attribute__((unused)) = params_.globals[46];\
auto& _pp_var_head1_theta_gaba    __attribute__((unused)) = params_.globals[47];\
auto& _pp_var_head2_E_K_NA        __attribute__((unused)) = params_.globals[48];\
auto& _pp_var_head2_E_head_ca     __attribute__((unused)) = params_.globals[49];\
auto& _pp_var_head2_Gleak         __attribute__((unused)) = params_.globals[50];\
auto& _pp_var_head2_I_max         __attribute__((unused)) = params_.globals[51];\
auto& _pp_var_head2_Ip            __attribute__((unused)) = params_.globals[52];\
auto& _pp_var_head2_Kp            __attribute__((unused)) = params_.globals[53];\
auto& _pp_var_head2_beta          __attribute__((unused)) = params_.globals[54];\
auto& _pp_var_head2_diam          __attribute__((unused)) = params_.globals[55];\
auto& _pp_var_head2_eta_gaba      __attribute__((unused)) = params_.globals[56];\
auto& _pp_var_head2_gamma         __attribute__((unused)) = params_.globals[57];\
auto& _pp_var_head2_gbar_gaba     __attribute__((unused)) = params_.globals[58];\
auto& _pp_var_head2_gbar_nmda     __attribute__((unused)) = params_.globals[59];\
auto& _pp_var_head2_ggap          __attribute__((unused)) = params_.globals[60];\
auto& _pp_var_head2_kCL           __attribute__((unused)) = params_.globals[61];\
auto& _pp_var_head2_kK            __attribute__((unused)) = params_.globals[62];\
auto& _pp_var_head2_kb            __attribute__((unused)) = params_.globals[63];\
auto& _pp_var_head2_kf            __attribute__((unused)) = params_.globals[64];\
auto& _pp_var_head2_length        __attribute__((unused)) = params_.globals[65];\
auto& _pp_var_head2_nTB           __attribute__((unused)) = params_.globals[66];\
auto& _pp_var_head2_n_K           __attribute__((unused)) = params_.globals[67];\
auto& _pp_var_head2_n_Mg_out      __attribute__((unused)) = params_.globals[68];\
auto& _pp_var_head2_sigma_gaba    __attribute__((unused)) = params_.globals[69];\
auto& _pp_var_head2_tau_d         __attribute__((unused)) = params_.globals[70];\
auto& _pp_var_head2_tau_decay_nmda __attribute__((unused)) = params_.globals[71];\
auto& _pp_var_head2_tau_gaba0     __attribute__((unused)) = params_.globals[72];\
auto& _pp_var_head2_tau_r         __attribute__((unused)) = params_.globals[73];\
auto& _pp_var_head2_tau_rise_nmda __attribute__((unused)) = params_.globals[74];\
auto& _pp_var_head2_theta_gaba    __attribute__((unused)) = params_.globals[75];\
auto& _pp_var_neck1_Gleak         __attribute__((unused)) = params_.globals[76];\
auto& _pp_var_neck1_diam          __attribute__((unused)) = params_.globals[77];\
auto& _pp_var_neck1_length        __attribute__((unused)) = params_.globals[78];\
auto& _pp_var_neck1_n_K           __attribute__((unused)) = params_.globals[79];\
auto& _pp_var_neck2_Gleak         __attribute__((unused)) = params_.globals[80];\
auto& _pp_var_neck2_diam          __attribute__((unused)) = params_.globals[81];\
auto& _pp_var_neck2_length        __attribute__((unused)) = params_.globals[82];\
auto& _pp_var_neck2_n_K           __attribute__((unused)) = params_.globals[83];\
auto& _pp_var_out_n_CA            __attribute__((unused)) = params_.globals[84];\
auto& _pp_var_out_n_CL            __attribute__((unused)) = params_.globals[85];\
auto& _pp_var_out_n_K             __attribute__((unused)) = params_.globals[86];\
auto* _pp_var_neck1_V             __attribute__((unused)) = params_.state_vars[0];\
auto* _pp_var_neck1_n_CA          __attribute__((unused)) = params_.state_vars[1];\
auto* _pp_var_neck1_n_CL          __attribute__((unused)) = params_.state_vars[2];\
auto* _pp_var_head1_V             __attribute__((unused)) = params_.state_vars[3];\
auto* _pp_var_head1_n_CL          __attribute__((unused)) = params_.state_vars[4];\
auto* _pp_var_head1_n_CA          __attribute__((unused)) = params_.state_vars[5];\
auto* _pp_var_head1_n_B           __attribute__((unused)) = params_.state_vars[6];\
auto* _pp_var_head1_X             __attribute__((unused)) = params_.state_vars[7];\
auto* _pp_var_head1_Y             __attribute__((unused)) = params_.state_vars[8];\
auto* _pp_var_head1_ggaba         __attribute__((unused)) = params_.state_vars[9];\
auto* _pp_var_head1_ca_presyn     __attribute__((unused)) = params_.state_vars[10];\
auto* _pp_var_neck2_V             __attribute__((unused)) = params_.state_vars[11];\
auto* _pp_var_neck2_n_CA          __attribute__((unused)) = params_.state_vars[12];\
auto* _pp_var_neck2_n_CL          __attribute__((unused)) = params_.state_vars[13];\
auto* _pp_var_head2_V             __attribute__((unused)) = params_.state_vars[14];\
auto* _pp_var_head2_n_CL          __attribute__((unused)) = params_.state_vars[15];\
auto* _pp_var_head2_n_CA          __attribute__((unused)) = params_.state_vars[16];\
auto* _pp_var_head2_n_B           __attribute__((unused)) = params_.state_vars[17];\
auto* _pp_var_head2_X             __attribute__((unused)) = params_.state_vars[18];\
auto* _pp_var_head2_Y             __attribute__((unused)) = params_.state_vars[19];\
auto* _pp_var_head2_ggaba         __attribute__((unused)) = params_.state_vars[20];\
auto* _pp_var_head2_ca_presyn     __attribute__((unused)) = params_.state_vars[21];\
//End of IFACEBLOCK

namespace {
    __global__ void compute_currents(arb_mechanism_ppack params_) {
        PPACK_IFACE_BLOCK;
        const size_t i = threadIdx.x + blockDim.x*blockIdx.x;
        if(i < _pp_var_gap_junction_width){
            const arb_size_type cv1 = _pp_var_gap_junctions[i].loc.first;
            const arb_size_type cv2 = _pp_var_gap_junctions[i].loc.second;

            const arb_value_type dend1_V = _pp_var_vec_v[cv1]/1000;   //[mv] -> [V]
            const arb_value_type dend2_V = _pp_var_vec_v[cv2]/1000;   //[mv] -> [V]
            const arb_value_type dt      = _pp_var_vec_dt[cv1];

            const arb_value_type x0 = 1.0/_pp_var_neck1_length;
            const arb_value_type x1 = 1.0/_pp_var_neck1_n_CA[i];
            const arb_value_type x2 = 1.0/_pp_var_F;
            const arb_value_type x3 = _pp_var_R*_pp_var_T;
            const arb_value_type x4 = x2*x3;
            const arb_value_type x5 = (1.0/2.0)*x4;
            const arb_value_type x6 = -dend1_V + _pp_var_neck1_V[i];
            const arb_value_type x7 = pow(_pp_var_neck1_diam, -2);
            const arb_value_type x8 = M_1_PI;
            const arb_value_type x9 = x3*x8/pow(_pp_var_F, 2);
            const arb_value_type x10 = _pp_var_neck1_length*x7*x9;
            const arb_value_type x11 = 0.5/_pp_var_D_CA;
            const arb_value_type x12 = x1*x10*x11;
            const arb_value_type x13 = -(x5*log(_pp_var_dend1_n_CA*x1) + x6)/(x12 + 1.0/_pp_var_dend1_half_G_CA);
            const arb_value_type x14 = 1.0/_pp_var_neck1_n_CL[i];
            const arb_value_type x15 = 2.0/_pp_var_D_CL;
            const arb_value_type x16 = x10*x14*x15;
            const arb_value_type x17 = -(-x4*log(_pp_var_dend1_n_CL*x14) + x6)/(x16 + 1.0/_pp_var_dend1_half_G_CL);
            const arb_value_type x18 = x13 + x17;
            const arb_value_type x19 = -_pp_var_neck1_V[i];
            const arb_value_type x20 = -_pp_var_head1_V[i];
            const arb_value_type x21 = _pp_var_neck1_V[i] + x20;
            const arb_value_type x22 = 1.0/_pp_var_head1_n_CL[i];
            const arb_value_type x23 = pow(_pp_var_head1_diam, 2);
            const arb_value_type x24 = 1.0/x23;
            const arb_value_type x25 = _pp_var_head1_length*x24;
            const arb_value_type x26 = x15*x9;
            const arb_value_type x27 = x22*x25*x26;
            const arb_value_type x28 = 1.0/(x16 + x27);
            const arb_value_type x29 = _pp_var_neck1_Gleak*(x19 - x4*log(_pp_var_out_n_CL*x14)) - x28*(x21 - x4*log(_pp_var_head1_n_CL[i]*x14));
            const arb_value_type x30 = 1.0/_pp_var_head1_n_CA[i];
            const arb_value_type x31 = x11*x9;
            const arb_value_type x32 = x25*x30*x31;
            const arb_value_type x33 = 1.0/(x12 + x32);
            const arb_value_type x34 = 2.5000000000000002e-20*_pp_var_F;
            const arb_value_type x35 = -x33*(x21 + x5*log(_pp_var_head1_n_CA[i]*x1)) - x34*(_pp_var_neck1_n_CA[i] - 5.0000000000000002e-5);
            const arb_value_type x36 = x8/_pp_var_Cm;
            const arb_value_type x37 = x2*x8;
            const arb_value_type x38 = 4*x37;
            const arb_value_type x39 = x0*x7;
            const arb_value_type x40 = 2*x37;
            const arb_value_type x41 = 1.0/_pp_var_head1_length;
            const arb_value_type x42 = _pp_var_head2_V[i] + x20;
            const arb_value_type x43 = 1.0/_pp_var_head2_n_CL[i];
            const arb_value_type x44 = pow(_pp_var_head2_diam, 2);
            const arb_value_type x45 = 1.0/x44;
            const arb_value_type x46 = _pp_var_head2_length*x45;
            const arb_value_type x47 = x26*x43*x46;
            const arb_value_type x48 = 1.0/_pp_var_head2_ggap + 1.0/_pp_var_head1_ggap;
            const arb_value_type x49 = 1.0/(x27 + x47 + x48);
            const arb_value_type x50 = _pp_var_head1_V[i] + x19;
            const arb_value_type x51 = 1.0/_pp_var_head1_kK;
            const arb_value_type x52 = 1.0/_pp_var_head1_kCL;
            const arb_value_type x53 = _pp_var_head1_n_CL[i]*x52;
            const arb_value_type x54 = _pp_var_head1_n_K*x51;
            const arb_value_type x55 = _pp_var_out_n_CL*x52;
            const arb_value_type x56 = _pp_var_out_n_K*x51;
            const arb_value_type x57 = _pp_var_out_n_CL*_pp_var_out_n_K;
            const arb_value_type x58 = M_PI*_pp_var_head1_length;
            const arb_value_type x59 = _pp_var_head1_Gleak*(x20 - x4*log(_pp_var_out_n_CL*x22)) - _pp_var_head1_I_max*_pp_var_head1_diam*_pp_var_head1_kCL*x51*x58*(-_pp_var_head1_n_CL[i]*_pp_var_head1_n_K + x57)/((x53 + 1)*(x54 + 1)*(x55*x56 + 1) + (x55 + 1)*(x56 + 1)*(x53*x54 + 1)) - x28*(-x4*log(_pp_var_neck1_n_CL[i]*x22) + x50) + x42*x49;
            const arb_value_type x60 = 1.0/_pp_var_head2_n_CA[i];
            const arb_value_type x61 = x31*x46*x60;
            const arb_value_type x62 = 1.0/(x32 + x48 + x61);
            const arb_value_type x63 = -x33*(x5*log(_pp_var_neck1_n_CA[i]*x30) + x50) - x34*(_pp_var_head1_n_CA[i] - 5.0000000000000002e-5) + x42*x62;
            const arb_value_type x64 = x24*x41;
            const arb_value_type x65 = _pp_var_head1_kf*_pp_var_head1_n_B[i]*_pp_var_head1_n_CA[i];
            const arb_value_type x66 = _pp_var_head1_kb*(_pp_var_head1_nTB - _pp_var_head1_n_B[i]);
            const arb_value_type x67 = (1.0/2.0)*_pp_var_F;
            const arb_value_type x68 = pow(_pp_var_head1_ca_presyn[0], 2);
            const arb_value_type x69 = 1.0/_pp_var_neck2_length;
            const arb_value_type x70 = 1.0/_pp_var_neck2_n_CA[i];
            const arb_value_type x71 = -dend2_V + _pp_var_neck2_V[0];
            const arb_value_type x72 = pow(_pp_var_neck2_diam, -2);
            const arb_value_type x73 = _pp_var_neck2_length*x72;
            const arb_value_type x74 = x31*x70*x73;
            const arb_value_type x75 = -(x5*log(_pp_var_dend2_n_CA*x70) + x71)/(x74 + 1.0/_pp_var_dend2_half_G_CA);
            const arb_value_type x76 = 1.0/_pp_var_neck2_n_CL[i];
            const arb_value_type x77 = x26*x73*x76;
            const arb_value_type x78 = -(-x4*log(_pp_var_dend2_n_CL*x76) + x71)/(x77 + 1.0/_pp_var_dend2_half_G_CL);
            const arb_value_type x79 = x75 + x78;
            const arb_value_type x80 = -_pp_var_neck2_V[i];
            const arb_value_type x81 = -_pp_var_head2_V[i];
            const arb_value_type x82 = _pp_var_neck2_V[i] + x81;
            const arb_value_type x83 = 1.0/(x47 + x77);
            const arb_value_type x84 = _pp_var_neck2_Gleak*(-x4*log(_pp_var_out_n_CL*x76) + x80) - x83*(-x4*log(_pp_var_head2_n_CL[i]*x76) + x82);
            const arb_value_type x85 = 1.0/(x61 + x74);
            const arb_value_type x86 = -x34*(_pp_var_neck2_n_CA[i] - 5.0000000000000002e-5) - x85*(x5*log(_pp_var_head2_n_CA[i]*x70) + x82);
            const arb_value_type x87 = x69*x72;
            const arb_value_type x88 = 1.0/_pp_var_head2_length;
            const arb_value_type x89 = _pp_var_head1_V[i] + x81;
            const arb_value_type x90 = _pp_var_head2_V[i] + x80;
            const arb_value_type x91 = 1.0/_pp_var_head2_kK;
            const arb_value_type x92 = 1.0/_pp_var_head2_kCL;
            const arb_value_type x93 = _pp_var_head2_n_CL[i]*x92;
            const arb_value_type x94 = _pp_var_head2_n_K*x91;
            const arb_value_type x95 = M_PI*_pp_var_head2_length;
            const arb_value_type x96 = _pp_var_head2_Gleak*(-x4*log(_pp_var_out_n_CL*x43) + x81) - _pp_var_head2_I_max*_pp_var_head2_diam*_pp_var_head2_kCL*x91*x95*(-_pp_var_head2_n_CL[i]*_pp_var_head2_n_K + x57)/((x93 + 1)*(x94 + 1)*(x57*x91*x92 + 1) + (_pp_var_out_n_CL*x92 + 1)*(_pp_var_out_n_K*x91 + 1)*(x93*x94 + 1)) + x49*x89 - x83*(-x4*log(_pp_var_neck2_n_CL[i]*x43) + x90);
            const arb_value_type x97 = -x34*(_pp_var_head2_n_CA[i] - 5.0000000000000002e-5) + x62*x89 - x85*(x5*log(_pp_var_neck2_n_CA[i]*x60) + x90);
            const arb_value_type x98 = x45*x88;
            const arb_value_type x99 = _pp_var_head2_kf*_pp_var_head2_n_B[i]*_pp_var_head2_n_CA[i];
            const arb_value_type x100 = _pp_var_head2_kb*(_pp_var_head2_nTB - _pp_var_head2_n_B[i]);
            const arb_value_type x101 = pow(_pp_var_head2_ca_presyn[i], 2);

            const arb_value_type grad_neck1_V            = x0*x36*(x18 + x29 + x35)/_pp_var_neck1_diam;
            const arb_value_type grad_neck1_n_CL         = -x38*x39*(x17 + x29);
            const arb_value_type grad_neck1_n_CA         = x39*x40*(x13 + x35);
            const arb_value_type grad_head1_V            = x36*x41*(x59 + x63)/_pp_var_head1_diam;
            const arb_value_type grad_head1_n_CL         = -x38*x59*x64;
            const arb_value_type grad_head1_n_CA         = x40*x64*(x23*x58*x67*(x65 - x66) + x63);
            const arb_value_type grad_head1_n_B          = -2*x65 + 2*x66;
            const arb_value_type grad_head1_Y            = -_pp_var_head1_Y[i]/_pp_var_head1_tau_d;
            const arb_value_type grad_head1_X            = (-_pp_var_head1_X[i] - _pp_var_head1_Y[i] + 1)/_pp_var_head1_tau_r;
            const arb_value_type grad_head1_ggaba        = _pp_var_head1_Y[i]*_pp_var_head1_gbar_gaba - _pp_var_head1_ggaba[i]/(_pp_var_head1_eta_gaba/(exp((-_pp_var_head1_ca_presyn[i] + _pp_var_head1_theta_gaba)/_pp_var_head1_sigma_gaba) + 1.0) + _pp_var_head1_tau_gaba0);
            const arb_value_type grad_head1_ca_presyn    = _pp_var_head1_Ip - _pp_var_head1_beta*x68/(pow(_pp_var_head1_Kp, 2) + x68) + _pp_var_head1_gamma*log(2/_pp_var_head1_ca_presyn[i]);
            const arb_value_type grad_neck2_V            = x36*x69*(x79 + x84 + x86)/_pp_var_neck2_diam;
            const arb_value_type grad_neck2_n_CL         = -x38*x87*(x78 + x84);
            const arb_value_type grad_neck2_n_CA         = x40*x87*(x75 + x86);
            const arb_value_type grad_head2_V            = x36*x88*(x96 + x97)/_pp_var_head2_diam;
            const arb_value_type grad_head2_n_CL         = -x38*x96*x98;
            const arb_value_type grad_head2_n_CA         = x40*x98*(x44*x67*x95*(-x100 + x99) + x97);
            const arb_value_type grad_head2_n_B          = 2*x100 - 2*x99;
            const arb_value_type grad_head2_Y            = -_pp_var_head2_Y[i]/_pp_var_head2_tau_d;
            const arb_value_type grad_head2_X            = (-_pp_var_head2_X[i] - _pp_var_head2_Y[i] + 1)/_pp_var_head2_tau_r;
            const arb_value_type grad_head2_ggaba        = _pp_var_head2_Y[i]*_pp_var_head2_gbar_gaba - _pp_var_head2_ggaba[i]/(_pp_var_head2_eta_gaba/(exp((-_pp_var_head2_ca_presyn[i] + _pp_var_head2_theta_gaba)/_pp_var_head2_sigma_gaba) + 1.0) + _pp_var_head2_tau_gaba0);
            const arb_value_type grad_head2_ca_presyn    = _pp_var_head2_Ip - _pp_var_head2_beta*x101/(pow(_pp_var_head2_Kp, 2) + x101) + _pp_var_head2_gamma*log(2/_pp_var_head2_ca_presyn[i]);

            //update states
            _pp_var_neck1_V[i]         += dt * grad_neck1_V;
            _pp_var_neck1_n_CL[i]      += dt * grad_neck1_n_CL;
            _pp_var_neck1_n_CA[i]      += dt * grad_neck1_n_CA;
            _pp_var_head1_V[i]         += dt * grad_head1_V;
            _pp_var_head1_n_CL[i]      += dt * grad_head1_n_CL;
            _pp_var_head1_n_CA[i]      += dt * grad_head1_n_CA;
            _pp_var_head1_n_B[i]       += dt * grad_head1_n_B;
            _pp_var_head1_Y[i]         += dt * grad_head1_Y;
            _pp_var_head1_X[i]         += dt * grad_head1_X;
            _pp_var_head1_ggaba[i]     += dt * grad_head1_ggaba;
            _pp_var_head1_ca_presyn[i] += dt * grad_head1_ca_presyn;
            _pp_var_neck2_V[i]         += dt * grad_neck2_V;
            _pp_var_neck2_n_CL[i]      += dt * grad_neck2_n_CL;
            _pp_var_neck2_n_CA[i]      += dt * grad_neck2_n_CA;
            _pp_var_head2_V[i]         += dt * grad_head2_V;
            _pp_var_head2_n_CL[i]      += dt * grad_head2_n_CL;
            _pp_var_head2_n_CA[i]      += dt * grad_head2_n_CA;
            _pp_var_head2_n_B[i]       += dt * grad_head2_n_B;
            _pp_var_head2_Y[i]         += dt * grad_head2_Y;
            _pp_var_head2_X[i]         += dt * grad_head2_X;
            _pp_var_head2_ggaba[i]     += dt * grad_head2_ggaba;
            _pp_var_head2_ca_presyn[i] += dt * grad_head2_ca_presyn;

            if (_pp_var_head1_n_CA[i] <= 0) _pp_var_head1_n_CA[i] = 1e-20;
            if (_pp_var_head2_n_CA[i] <= 0) _pp_var_head2_n_CA[i] = 1e-20;
            if (_pp_var_head1_n_CL[i] <= 0) _pp_var_head1_n_CL[i] = 1e-20;
            if (_pp_var_head2_n_CL[i] <= 0) _pp_var_head2_n_CL[i] = 1e-20;
            if (_pp_var_neck1_n_CA[i] <= 0) _pp_var_neck1_n_CA[i] = 1e-20;
            if (_pp_var_neck2_n_CA[i] <= 0) _pp_var_neck2_n_CA[i] = 1e-20;
            if (_pp_var_neck1_n_CL[i] <= 0) _pp_var_neck1_n_CL[i] = 1e-20;
            if (_pp_var_neck2_n_CL[i] <= 0) _pp_var_neck2_n_CL[i] = 1e-20;
            if (_pp_var_head1_n_B[i]  <= 0)  _pp_var_head1_n_B[i] = 1e-20;
            if (_pp_var_head2_n_B[i]  <= 0)  _pp_var_head2_n_B[i] = 1e-20;

            //TODO fixen van weights zodat het minder calculaties worden :D
            //write current attributions
            const arb_value_type I_inward_dend1 = x18; //[A]
//            const arb_value_type I_inward_dend2 = x79; //[A]

            _pp_var_vec_i[cv1] +=  fma(_pp_var_gap_junctions[i].weight , I_inward_dend1, _pp_var_vec_i[cv1]); //[A] -> [A/m2]
//            _pp_var_vec_i[cv2] +=  fma( _pp_var_gap_junctions[i].weight , I_inward_dend2, _pp_var_vec_i[cv2]); //[A] -> [A/m2] This can not be used right now..

        }
    }
} // namespace

void mechanism_glomerulus_gj_gpu_init_(arb_mechanism_ppack*) {
    std::cout << "init glomerulus GPU" << std::endl;
}
void mechanism_glomerulus_gj_gpu_compute_currents_(arb_mechanism_ppack* p) {
    unsigned block_dim = 128;
    unsigned grid_dim = ::arb::gpu::impl::block_count(p->gap_junction_width, block_dim);
    compute_currents<<<grid_dim, block_dim>>>(*p);
}
void mechanism_glomerulus_gj_gpu_advance_state_(arb_mechanism_ppack*) {}
void mechanism_glomerulus_gj_gpu_write_ions_(arb_mechanism_ppack*) {}
void mechanism_glomerulus_gj_gpu_post_event_(arb_mechanism_ppack*) {}
void mechanism_glomerulus_gj_gpu_apply_events_(arb_mechanism_ppack*) {}

} // namespace IOU_catalogue
} // namespace arb
