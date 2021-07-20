//
// Created by max on 09-07-21.
//
#include "iostream" //debugging
//#define S(x) std::cout << #x << "\t\t" << x << std::endl;

#include <arbor/arb_types.h>
#include <arbor/mechanism_abi.h>
#include <cmath>

namespace arb::IOU_catalogue::kernel_glomerulus {
    static constexpr unsigned simd_width_ = 0;

#define PPACK_IFACE_BLOCK \
[[maybe_unused]] auto  _pp_var_width                = pp->width;\
[[maybe_unused]] auto  _pp_var_n_detectors          = pp->n_detectors;\
[[maybe_unused]] auto* _pp_var_vec_ci               = pp->vec_ci;\
[[maybe_unused]] auto* _pp_var_vec_di               = pp->vec_di;\
[[maybe_unused]] auto* _pp_var_vec_t                = pp->vec_t;\
[[maybe_unused]] auto* _pp_var_vec_dt               = pp->vec_dt;\
[[maybe_unused]] auto* _pp_var_vec_v                = pp->vec_v;\
[[maybe_unused]] auto* _pp_var_vec_i                = pp->vec_i;\
[[maybe_unused]] auto* _pp_var_vec_g                = pp->vec_g;\
[[maybe_unused]] auto* _pp_var_temperature_degC     = pp->temperature_degC;\
[[maybe_unused]] auto* _pp_var_diam_um              = pp->diam_um;\
[[maybe_unused]] auto* _pp_var_time_since_spike     = pp->time_since_spike;\
[[maybe_unused]] auto* _pp_var_node_index           = pp->node_index;\
[[maybe_unused]] auto* _pp_var_multiplicity         = pp->multiplicity;\
[[maybe_unused]] auto* _pp_var_weight               = pp->weight;\
[[maybe_unused]] auto& _pp_var_events               = pp->events;\
[[maybe_unused]] auto& _pp_var_mechanism_id         = pp->mechanism_id;\
[[maybe_unused]] auto& _pp_var_index_constraints    = pp->index_constraints; \
[[maybe_unused]] auto* _pp_var_gap_junctions        = pp->gap_junctions; \
[[maybe_unused]] auto _pp_var_gap_junction_width    = pp->gap_junction_width;\
[[maybe_unused]] auto& _pp_var_Cm                   = pp->globals[0];\
[[maybe_unused]] auto& _pp_var_D_CA                 = pp->globals[1];\
[[maybe_unused]] auto& _pp_var_D_CL                 = pp->globals[2];\
[[maybe_unused]] auto& _pp_var_F                    = pp->globals[3];\
[[maybe_unused]] auto& _pp_var_R                    = pp->globals[4];\
[[maybe_unused]] auto& _pp_var_T                    = pp->globals[5];\
[[maybe_unused]] auto& _pp_var_dend1_diam           = pp->globals[6];\
[[maybe_unused]] auto& _pp_var_dend1_half_G_CA      = pp->globals[7];\
[[maybe_unused]] auto& _pp_var_dend1_half_G_CL      = pp->globals[8];\
[[maybe_unused]] auto& _pp_var_dend1_length         = pp->globals[9];\
[[maybe_unused]] auto& _pp_var_dend1_n_CA           = pp->globals[10];\
[[maybe_unused]] auto& _pp_var_dend1_n_CL           = pp->globals[11];\
[[maybe_unused]] auto& _pp_var_dend1_n_K            = pp->globals[12];\
[[maybe_unused]] auto& _pp_var_dend2_diam           = pp->globals[13];\
[[maybe_unused]] auto& _pp_var_dend2_half_G_CA      = pp->globals[14];\
[[maybe_unused]] auto& _pp_var_dend2_half_G_CL      = pp->globals[15];\
[[maybe_unused]] auto& _pp_var_dend2_length         = pp->globals[16];\
[[maybe_unused]] auto& _pp_var_dend2_n_CA           = pp->globals[17];\
[[maybe_unused]] auto& _pp_var_dend2_n_CL           = pp->globals[18];\
[[maybe_unused]] auto& _pp_var_dend2_n_K            = pp->globals[19];\
[[maybe_unused]] auto& _pp_var_head1_E_K_NA         = pp->globals[20];\
[[maybe_unused]] auto& _pp_var_head1_E_head_ca      = pp->globals[21];\
[[maybe_unused]] auto& _pp_var_head1_Gleak          = pp->globals[22];\
[[maybe_unused]] auto& _pp_var_head1_I_max          = pp->globals[23];\
[[maybe_unused]] auto& _pp_var_head1_Ip             = pp->globals[24];\
[[maybe_unused]] auto& _pp_var_head1_Kp             = pp->globals[25];\
[[maybe_unused]] auto& _pp_var_head1_beta           = pp->globals[26];\
[[maybe_unused]] auto& _pp_var_head1_diam           = pp->globals[27];\
[[maybe_unused]] auto& _pp_var_head1_eta_gaba       = pp->globals[28];\
[[maybe_unused]] auto& _pp_var_head1_gamma          = pp->globals[29];\
[[maybe_unused]] auto& _pp_var_head1_gbar_gaba      = pp->globals[30];\
[[maybe_unused]] auto& _pp_var_head1_gbar_nmda      = pp->globals[31];\
[[maybe_unused]] auto& _pp_var_head1_ggap           = pp->globals[32];\
[[maybe_unused]] auto& _pp_var_head1_kCL            = pp->globals[33];\
[[maybe_unused]] auto& _pp_var_head1_kK             = pp->globals[34];\
[[maybe_unused]] auto& _pp_var_head1_kb             = pp->globals[35];\
[[maybe_unused]] auto& _pp_var_head1_kf             = pp->globals[36];\
[[maybe_unused]] auto& _pp_var_head1_length         = pp->globals[37];\
[[maybe_unused]] auto& _pp_var_head1_nTB            = pp->globals[38];\
[[maybe_unused]] auto& _pp_var_head1_n_K            = pp->globals[39];\
[[maybe_unused]] auto& _pp_var_head1_n_Mg_out       = pp->globals[40];\
[[maybe_unused]] auto& _pp_var_head1_sigma_gaba     = pp->globals[41];\
[[maybe_unused]] auto& _pp_var_head1_tau_d          = pp->globals[42];\
[[maybe_unused]] auto& _pp_var_head1_tau_decay_nmda = pp->globals[43];\
[[maybe_unused]] auto& _pp_var_head1_tau_gaba0      = pp->globals[44];\
[[maybe_unused]] auto& _pp_var_head1_tau_r          = pp->globals[45];\
[[maybe_unused]] auto& _pp_var_head1_tau_rise_nmda  = pp->globals[46];\
[[maybe_unused]] auto& _pp_var_head1_theta_gaba     = pp->globals[47];\
[[maybe_unused]] auto& _pp_var_head2_E_K_NA         = pp->globals[48];\
[[maybe_unused]] auto& _pp_var_head2_E_head_ca      = pp->globals[49];\
[[maybe_unused]] auto& _pp_var_head2_Gleak          = pp->globals[50];\
[[maybe_unused]] auto& _pp_var_head2_I_max          = pp->globals[51];\
[[maybe_unused]] auto& _pp_var_head2_Ip             = pp->globals[52];\
[[maybe_unused]] auto& _pp_var_head2_Kp             = pp->globals[53];\
[[maybe_unused]] auto& _pp_var_head2_beta           = pp->globals[54];\
[[maybe_unused]] auto& _pp_var_head2_diam           = pp->globals[55];\
[[maybe_unused]] auto& _pp_var_head2_eta_gaba       = pp->globals[56];\
[[maybe_unused]] auto& _pp_var_head2_gamma          = pp->globals[57];\
[[maybe_unused]] auto& _pp_var_head2_gbar_gaba      = pp->globals[58];\
[[maybe_unused]] auto& _pp_var_head2_gbar_nmda      = pp->globals[59];\
[[maybe_unused]] auto& _pp_var_head2_ggap           = pp->globals[60];\
[[maybe_unused]] auto& _pp_var_head2_kCL            = pp->globals[61];\
[[maybe_unused]] auto& _pp_var_head2_kK             = pp->globals[62];\
[[maybe_unused]] auto& _pp_var_head2_kb             = pp->globals[63];\
[[maybe_unused]] auto& _pp_var_head2_kf             = pp->globals[64];\
[[maybe_unused]] auto& _pp_var_head2_length         = pp->globals[65];\
[[maybe_unused]] auto& _pp_var_head2_nTB            = pp->globals[66];\
[[maybe_unused]] auto& _pp_var_head2_n_K            = pp->globals[67];\
[[maybe_unused]] auto& _pp_var_head2_n_Mg_out       = pp->globals[68];\
[[maybe_unused]] auto& _pp_var_head2_sigma_gaba     = pp->globals[69];\
[[maybe_unused]] auto& _pp_var_head2_tau_d          = pp->globals[70];\
[[maybe_unused]] auto& _pp_var_head2_tau_decay_nmda = pp->globals[71];\
[[maybe_unused]] auto& _pp_var_head2_tau_gaba0      = pp->globals[72];\
[[maybe_unused]] auto& _pp_var_head2_tau_r          = pp->globals[73];\
[[maybe_unused]] auto& _pp_var_head2_tau_rise_nmda  = pp->globals[74];\
[[maybe_unused]] auto& _pp_var_head2_theta_gaba     = pp->globals[75];\
[[maybe_unused]] auto& _pp_var_neck1_Gleak          = pp->globals[76];\
[[maybe_unused]] auto& _pp_var_neck1_diam           = pp->globals[77];\
[[maybe_unused]] auto& _pp_var_neck1_length         = pp->globals[78];\
[[maybe_unused]] auto& _pp_var_neck1_n_K            = pp->globals[79];\
[[maybe_unused]] auto& _pp_var_neck2_Gleak          = pp->globals[80];\
[[maybe_unused]] auto& _pp_var_neck2_diam           = pp->globals[81];\
[[maybe_unused]] auto& _pp_var_neck2_length         = pp->globals[82];\
[[maybe_unused]] auto& _pp_var_neck2_n_K            = pp->globals[83];\
[[maybe_unused]] auto& _pp_var_out_n_CA             = pp->globals[84];\
[[maybe_unused]] auto& _pp_var_out_n_CL             = pp->globals[85];\
[[maybe_unused]] auto& _pp_var_out_n_K              = pp->globals[86];\
[[maybe_unused]] auto* _pp_var_neck1_V              = pp->state_vars[0];\
[[maybe_unused]] auto* _pp_var_neck1_n_CA           = pp->state_vars[1];\
[[maybe_unused]] auto* _pp_var_neck1_n_CL           = pp->state_vars[2];\
[[maybe_unused]] auto* _pp_var_head1_V              = pp->state_vars[3];\
[[maybe_unused]] auto* _pp_var_head1_n_CL           = pp->state_vars[4];\
[[maybe_unused]] auto* _pp_var_head1_n_CA           = pp->state_vars[5];\
[[maybe_unused]] auto* _pp_var_head1_n_B            = pp->state_vars[6];\
[[maybe_unused]] auto* _pp_var_head1_X              = pp->state_vars[7];\
[[maybe_unused]] auto* _pp_var_head1_Y              = pp->state_vars[8];\
[[maybe_unused]] auto* _pp_var_head1_ggaba          = pp->state_vars[9];\
[[maybe_unused]] auto* _pp_var_head1_ca_presyn      = pp->state_vars[10];\
[[maybe_unused]] auto* _pp_var_neck2_V              = pp->state_vars[11];\
[[maybe_unused]] auto* _pp_var_neck2_n_CA           = pp->state_vars[12];\
[[maybe_unused]] auto* _pp_var_neck2_n_CL           = pp->state_vars[13];\
[[maybe_unused]] auto* _pp_var_head2_V              = pp->state_vars[14];\
[[maybe_unused]] auto* _pp_var_head2_n_CL           = pp->state_vars[15];\
[[maybe_unused]] auto* _pp_var_head2_n_CA           = pp->state_vars[16];\
[[maybe_unused]] auto* _pp_var_head2_n_B            = pp->state_vars[17];\
[[maybe_unused]] auto* _pp_var_head2_X              = pp->state_vars[18];\
[[maybe_unused]] auto* _pp_var_head2_Y              = pp->state_vars[19];\
[[maybe_unused]] auto* _pp_var_head2_ggaba          = pp->state_vars[20];\
[[maybe_unused]] auto* _pp_var_head2_ca_presyn      = pp->state_vars[21];\
//End of IFACEBLOCK

// interface methods
    static void init(arb_mechanism_ppack*) {
        std::cout << "init glomerulus" <<std::endl;
    }
    static void compute_currents(arb_mechanism_ppack* pp) {
        PPACK_IFACE_BLOCK
        for(arb_size_type i = 0; i < _pp_var_gap_junction_width; i++){

            const arb_size_type cv1 = _pp_var_gap_junctions[i].loc.first;
            const arb_size_type cv2 = _pp_var_gap_junctions[i].loc.first;

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


             //write current attributions
            const arb_value_type I_inward_dend1 = x18; //[A]
//            const arb_value_type I_inward_dend2 = x79; //[A]

            _pp_var_vec_i[cv1] +=  fma( _pp_var_gap_junctions[i].weight , I_inward_dend1, _pp_var_vec_i[cv1]); //[A] -> [A/m2]
//            _pp_var_vec_i[cv2] +=  fma( _pp_var_gap_junctions[i].weight , I_inward_dend2, _pp_var_vec_i[cv2]); //[A] -> [A/m2] This can not be used right now..
        }
    }
    static void advance_state(arb_mechanism_ppack*) {}
    static void write_ions(arb_mechanism_ppack*) {}
    static void apply_events(arb_mechanism_ppack*) {}
    static void post_event(arb_mechanism_ppack*) {}

#undef PPACK_IFACE_BLOCK
}

extern "C" {
arb_mechanism_interface* make_arb_IOU_catalogue_glomerulus_interface_multicore() {
    static arb_mechanism_interface result;
    result.partition_width = arb::IOU_catalogue::kernel_glomerulus::simd_width_;
    result.backend=arb_backend_kind_cpu;
    result.alignment=1;
    result.init_mechanism  = (arb_mechanism_method)arb::IOU_catalogue::kernel_glomerulus::init;
    result.compute_currents= (arb_mechanism_method)arb::IOU_catalogue::kernel_glomerulus::compute_currents;
    result.apply_events    = (arb_mechanism_method)arb::IOU_catalogue::kernel_glomerulus::apply_events;
    result.advance_state   = (arb_mechanism_method)arb::IOU_catalogue::kernel_glomerulus::advance_state;
    result.write_ions      = (arb_mechanism_method)arb::IOU_catalogue::kernel_glomerulus::write_ions;
    result.post_event      = (arb_mechanism_method)arb::IOU_catalogue::kernel_glomerulus::post_event;
    return &result;
}}
