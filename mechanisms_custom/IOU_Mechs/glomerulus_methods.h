#include <arbor/arb_types.h>

struct glomerulus_methods {
//state->
    arb_value_type dend1_V;     //is dit de voltage ?
    arb_value_type dend2_V;     //is dit de voltage ?

    arb_value_type I_inward_dend1; //is dit de stroom ?
    arb_value_type I_inward_dend2; //is dit de stroom ?

    arb_value_type neck1_V;
    arb_value_type neck1_n_CA;
    arb_value_type neck1_n_CL;

    arb_value_type head1_V;
    arb_value_type head1_n_CL;
    arb_value_type head1_n_CA;
    arb_value_type head1_n_B;
    arb_value_type head1_X;
    arb_value_type head1_Y;
    arb_value_type head1_ggaba;
    arb_value_type head1_ca_presyn;

    arb_value_type neck2_V;
    arb_value_type neck2_n_CA;
    arb_value_type neck2_n_CL;

    arb_value_type head2_V;
    arb_value_type head2_n_CL;
    arb_value_type head2_n_CA;
    arb_value_type head2_n_B;
    arb_value_type head2_X;
    arb_value_type head2_Y;
    arb_value_type head2_ggaba;
    arb_value_type head2_ca_presyn;
//----<

//static
    arb_value_type Cm;
    arb_value_type D_CA;
    arb_value_type D_CL;
    arb_value_type F;
    arb_value_type R;
    arb_value_type T;

    arb_value_type dend1_diam;
    arb_value_type dend1_half_G_CA;
    arb_value_type dend1_half_G_CL;
    arb_value_type dend1_length;
    arb_value_type dend1_n_CA;
    arb_value_type dend1_n_CL;
    arb_value_type dend1_n_K;

    arb_value_type dend2_diam;
    arb_value_type dend2_half_G_CA;
    arb_value_type dend2_half_G_CL;
    arb_value_type dend2_length;
    arb_value_type dend2_n_CA;
    arb_value_type dend2_n_CL;
    arb_value_type dend2_n_K;

    arb_value_type head1_E_K_NA;
    arb_value_type head1_E_head_ca;
    arb_value_type head1_Gleak;
    arb_value_type head1_I_max;
    arb_value_type head1_Ip;
    arb_value_type head1_Kp;
    arb_value_type head1_beta;
    arb_value_type head1_diam;
    arb_value_type head1_eta_gaba;
    arb_value_type head1_gamma;
    arb_value_type head1_gbar_gaba;
    arb_value_type head1_gbar_nmda;
    arb_value_type head1_ggap;
    arb_value_type head1_kCL;
    arb_value_type head1_kK;
    arb_value_type head1_kb;
    arb_value_type head1_kf;
    arb_value_type head1_length;
    arb_value_type head1_nTB;
    arb_value_type head1_n_K;
    arb_value_type head1_n_Mg_out;
    arb_value_type head1_sigma_gaba;
    arb_value_type head1_tau_d;
    arb_value_type head1_tau_decay_nmda;
    arb_value_type head1_tau_gaba0;
    arb_value_type head1_tau_r;
    arb_value_type head1_tau_rise_nmda;
    arb_value_type head1_theta_gaba;

    arb_value_type head2_E_K_NA;
    arb_value_type head2_E_head_ca;
    arb_value_type head2_Gleak;
    arb_value_type head2_I_max;
    arb_value_type head2_Ip;
    arb_value_type head2_Kp;
    arb_value_type head2_beta;
    arb_value_type head2_diam;
    arb_value_type head2_eta_gaba;
    arb_value_type head2_gamma;
    arb_value_type head2_gbar_gaba;
    arb_value_type head2_gbar_nmda;
    arb_value_type head2_ggap;
    arb_value_type head2_kCL;
    arb_value_type head2_kK;
    arb_value_type head2_kb;
    arb_value_type head2_kf;
    arb_value_type head2_length;
    arb_value_type head2_nTB;
    arb_value_type head2_n_K;
    arb_value_type head2_n_Mg_out;
    arb_value_type head2_sigma_gaba;
    arb_value_type head2_tau_d;
    arb_value_type head2_tau_decay_nmda;
    arb_value_type head2_tau_gaba0;
    arb_value_type head2_tau_r;
    arb_value_type head2_tau_rise_nmda;
    arb_value_type head2_theta_gaba;

    arb_value_type neck1_Gleak;
    arb_value_type neck1_diam;
    arb_value_type neck1_length;
    arb_value_type neck1_n_K;

    arb_value_type neck2_Gleak;
    arb_value_type neck2_diam;
    arb_value_type neck2_length;
    arb_value_type neck2_n_K;

    arb_value_type out_n_CA;
    arb_value_type out_n_CL;
    arb_value_type out_n_K;
}; // glomerulus_methods

void ggj_init(struct glomerulus_methods * self) {
    self->Cm = 2e-6*1e4;
    self->D_CA = 0.79e-9;
    self->D_CL = 2.08e-9;
    self->F = 96485;
    self->R = 8.314;
    self->T = 293;
    self->T = 293;
    self->dend1_half_G_CL = 100;
    self->dend1_half_G_CA = 100;
    self->dend2_half_G_CL = 100;
    self->dend2_half_G_CA = 100;
    self->neck1_Gleak = 1.0e-8;
    self->neck2_Gleak = 1.0e-8;
    self->head1_tau_d = 0.255;
    self->head1_tau_r = 0.0050000000000000001;
    self->head1_gbar_gaba = 1.0000000000000001e-9;
    self->head1_tau_gaba0 = 0.029999999999999999;
    self->head1_eta_gaba = 0.034000000000000002;
    self->head1_theta_gaba = 0.0011999999999999999;
    self->head1_sigma_gaba = 6.0000000000000002e-5;
    self->head1_beta = 0.080000000000000002;
    self->head1_Kp = 0.00020000000000000001;
    self->head1_Ip = 0.00046999999999999999;
    self->head1_gamma = 0;
    self->head1_I_max = 30;
    self->head1_kK = 9;
    self->head1_kCL = 1;
    self->head1_Gleak = 1.0e-8;
    self->head1_n_Mg_out = 2;
    self->head1_gbar_nmda = 2.0000000000000001e-10;
    self->head1_tau_decay_nmda = 0.088999999999999996;
    self->head1_tau_rise_nmda = 0.0011999999999999999;
    self->head1_E_K_NA = 0.01;
    self->head1_E_head_ca = 0.01;
    self->head1_nTB = 0.10000000000000001;
    self->head1_kf = 1000000.0;
    self->head1_kb = 500;

    self->head2_tau_d = 0.255;
    self->head2_tau_r = 0.0050000000000000001;
    self->head2_gbar_gaba = 1.0000000000000001e-9;
    self->head2_tau_gaba0 = 0.029999999999999999;
    self->head2_eta_gaba = 0.034000000000000002;
    self->head2_theta_gaba = 0.0011999999999999999;
    self->head2_sigma_gaba = 6.0000000000000002e-5;
    self->head2_beta = 0.080000000000000002;
    self->head2_Kp = 0.00020000000000000001;
    self->head2_Ip = 0.00046999999999999999;
    self->head2_gamma = 0;
    self->head2_I_max = 30;
    self->head2_kK = 9;
    self->head2_kCL = 1;
    self->head2_Gleak = 1.0e-8;
    self->head2_n_Mg_out = 2;
    self->head2_gbar_nmda = 2.0000000000000001e-10;
    self->head2_tau_decay_nmda = 0.088999999999999996;
    self->head2_tau_rise_nmda = 0.0011999999999999999;
    self->head2_E_K_NA = 0.01;
    self->head2_E_head_ca = 0.01;
    self->head2_nTB = 0.10000000000000001;
    self->head2_kf = 1000000.0;
    self->head2_kb = 500;

    self->out_n_K = 3;
    self->out_n_CL = 134;
    self->out_n_CA = 2;

    self->neck1_V = -0.060999999999999999;
    self->neck1_diam = 9.9999999999999995e-8;
    self->neck1_length = 9.9999999999999995e-8;
    self->neck1_n_K = 85;
    self->neck1_n_CL = 5;
    self->neck1_n_CA = 10;

    self->head1_V = -0.060999999999999999;
    self->head1_diam = 9.9999999999999995e-8;
    self->head1_length = 9.9999999999999995e-8;
    self->head1_n_K = 85;
    self->head1_n_CL = 3.5;
    self->head1_n_CA = 10;
    self->head1_n_B = 0;
    self->head1_ggaba = 0;
    self->head1_ggap = 1.5e-11;
    self->head1_Y = 0;
    self->head1_X = 0;
    self->head1_ca_presyn = 0;

    self->neck2_V = -0.060999999999999999;
    self->neck2_diam = 9.9999999999999995e-8;
    self->neck2_length = 9.9999999999999995e-8;
    self->neck2_n_K = 85;
    self->neck2_n_CL = 5;
    self->neck2_n_CA = 10;

    self->head2_V = -0.060999999999999999;
    self->head2_n_CL = 3.5;
    self->head2_n_CA = 10;
    self->head2_n_B = 0;
    self->head2_Y = 0;
    self->head2_X = 0;
    self->head2_ggaba = 0;
    self->head2_ca_presyn = 0;

    self->head2_diam = 9.9999999999999995e-8;
    self->head2_length = 9.9999999999999995e-8;
    self->head2_n_K = 85;


    self->head2_ggap = 1.5e-11;

    self->dend1_V = -0.060999999999999999;
    self->dend1_diam = 9.9999999999999995e-7;
    self->dend1_length = 9.9999999999999995e-7;
    self->dend1_n_K = 85;
    self->dend1_n_CL = 15;
    self->dend1_n_CA = 15;
    self->dend2_V = -0.060999999999999999;
    self->dend2_diam = 9.9999999999999995e-7;
    self->dend2_length = 9.9999999999999995e-7;
    self->dend2_n_K = 85;
    self->dend2_n_CL = 15;
    self->dend2_n_CA = 15;
}
void ggj_timestep(struct glomerulus_methods * self, double dt) {
//do some intermediate stuff
    const arb_value_type x0 = 1.0/self->neck1_length;
    const arb_value_type x1 = 1.0/self->neck1_n_CA;
    const arb_value_type x2 = 1.0/self->F;
    const arb_value_type x3 = self->R*self->T;
    const arb_value_type x4 = x2*x3;
    const arb_value_type x5 = (1.0/2.0)*x4;
    const arb_value_type x6 = -self->dend1_V + self->neck1_V;
    const arb_value_type x7 = pow(self->neck1_diam, -2);
    const arb_value_type x8 = M_1_PI;
    const arb_value_type x9 = x3*x8/pow(self->F, 2);
    const arb_value_type x10 = self->neck1_length*x7*x9;
    const arb_value_type x11 = 0.5/self->D_CA;
    const arb_value_type x12 = x1*x10*x11;
    const arb_value_type x13 = -(x5*log(self->dend1_n_CA*x1) + x6)/(x12 + 1.0/self->dend1_half_G_CA);
    const arb_value_type x14 = 1.0/self->neck1_n_CL;
    const arb_value_type x15 = 2.0/self->D_CL;
    const arb_value_type x16 = x10*x14*x15;
    const arb_value_type x17 = -(-x4*log(self->dend1_n_CL*x14) + x6)/(x16 + 1.0/self->dend1_half_G_CL);
    const arb_value_type x18 = x13 + x17;
    const arb_value_type x19 = -self->neck1_V;
    const arb_value_type x20 = -self->head1_V;
    const arb_value_type x21 = self->neck1_V + x20;
    const arb_value_type x22 = 1.0/self->head1_n_CL;
    const arb_value_type x23 = pow(self->head1_diam, 2);
    const arb_value_type x24 = 1.0/x23;
    const arb_value_type x25 = self->head1_length*x24;
    const arb_value_type x26 = x15*x9;
    const arb_value_type x27 = x22*x25*x26;
    const arb_value_type x28 = 1.0/(x16 + x27);
    const arb_value_type x29 = self->neck1_Gleak*(x19 - x4*log(self->out_n_CL*x14)) - x28*(x21 - x4*log(self->head1_n_CL*x14));
    const arb_value_type x30 = 1.0/self->head1_n_CA;
    const arb_value_type x31 = x11*x9;
    const arb_value_type x32 = x25*x30*x31;
    const arb_value_type x33 = 1.0/(x12 + x32);
    const arb_value_type x34 = 2.5000000000000002e-20*self->F;
    const arb_value_type x35 = -x33*(x21 + x5*log(self->head1_n_CA*x1)) - x34*(self->neck1_n_CA - 5.0000000000000002e-5);
    const arb_value_type x36 = x8/self->Cm;
    const arb_value_type x37 = x2*x8;
    const arb_value_type x38 = 4*x37;
    const arb_value_type x39 = x0*x7;
    const arb_value_type x40 = 2*x37;
    const arb_value_type x41 = 1.0/self->head1_length;
    const arb_value_type x42 = self->head2_V + x20;
    const arb_value_type x43 = 1.0/self->head2_n_CL;
    const arb_value_type x44 = pow(self->head2_diam, 2);
    const arb_value_type x45 = 1.0/x44;
    const arb_value_type x46 = self->head2_length*x45;
    const arb_value_type x47 = x26*x43*x46;
    const arb_value_type x48 = 1.0/self->head2_ggap + 1.0/self->head1_ggap;
    const arb_value_type x49 = 1.0/(x27 + x47 + x48);
    const arb_value_type x50 = self->head1_V + x19;
    const arb_value_type x51 = 1.0/self->head1_kK;
    const arb_value_type x52 = 1.0/self->head1_kCL;
    const arb_value_type x53 = self->head1_n_CL*x52;
    const arb_value_type x54 = self->head1_n_K*x51;
    const arb_value_type x55 = self->out_n_CL*x52;
    const arb_value_type x56 = self->out_n_K*x51;
    const arb_value_type x57 = self->out_n_CL*self->out_n_K;
    const arb_value_type x58 = M_PI*self->head1_length;
    const arb_value_type x59 = self->head1_Gleak*(x20 - x4*log(self->out_n_CL*x22)) - self->head1_I_max*self->head1_diam*self->head1_kCL*x51*x58*(-self->head1_n_CL*self->head1_n_K + x57)/((x53 + 1)*(x54 + 1)*(x55*x56 + 1) + (x55 + 1)*(x56 + 1)*(x53*x54 + 1)) - x28*(-x4*log(self->neck1_n_CL*x22) + x50) + x42*x49;
    const arb_value_type x60 = 1.0/self->head2_n_CA;
    const arb_value_type x61 = x31*x46*x60;
    const arb_value_type x62 = 1.0/(x32 + x48 + x61);
    const arb_value_type x63 = -x33*(x5*log(self->neck1_n_CA*x30) + x50) - x34*(self->head1_n_CA - 5.0000000000000002e-5) + x42*x62;
    const arb_value_type x64 = x24*x41;
    const arb_value_type x65 = self->head1_kf*self->head1_n_B*self->head1_n_CA;
    const arb_value_type x66 = self->head1_kb*(self->head1_nTB - self->head1_n_B);
    const arb_value_type x67 = (1.0/2.0)*self->F;
    const arb_value_type x68 = pow(self->head1_ca_presyn, 2);
    const arb_value_type x69 = 1.0/self->neck2_length;
    const arb_value_type x70 = 1.0/self->neck2_n_CA;
    const arb_value_type x71 = -self->dend2_V + self->neck2_V;
    const arb_value_type x72 = pow(self->neck2_diam, -2);
    const arb_value_type x73 = self->neck2_length*x72;
    const arb_value_type x74 = x31*x70*x73;
    const arb_value_type x75 = -(x5*log(self->dend2_n_CA*x70) + x71)/(x74 + 1.0/self->dend2_half_G_CA);
    const arb_value_type x76 = 1.0/self->neck2_n_CL;
    const arb_value_type x77 = x26*x73*x76;
    const arb_value_type x78 = -(-x4*log(self->dend2_n_CL*x76) + x71)/(x77 + 1.0/self->dend2_half_G_CL);
    const arb_value_type x79 = x75 + x78;
    const arb_value_type x80 = -self->neck2_V;
    const arb_value_type x81 = -self->head2_V;
    const arb_value_type x82 = self->neck2_V + x81;
    const arb_value_type x83 = 1.0/(x47 + x77);
    const arb_value_type x84 = self->neck2_Gleak*(-x4*log(self->out_n_CL*x76) + x80) - x83*(-x4*log(self->head2_n_CL*x76) + x82);
    const arb_value_type x85 = 1.0/(x61 + x74);
    const arb_value_type x86 = -x34*(self->neck2_n_CA - 5.0000000000000002e-5) - x85*(x5*log(self->head2_n_CA*x70) + x82);
    const arb_value_type x87 = x69*x72;
    const arb_value_type x88 = 1.0/self->head2_length;
    const arb_value_type x89 = self->head1_V + x81;
    const arb_value_type x90 = self->head2_V + x80;
    const arb_value_type x91 = 1.0/self->head2_kK;
    const arb_value_type x92 = 1.0/self->head2_kCL;
    const arb_value_type x93 = self->head2_n_CL*x92;
    const arb_value_type x94 = self->head2_n_K*x91;
    const arb_value_type x95 = M_PI*self->head2_length;
    const arb_value_type x96 = self->head2_Gleak*(-x4*log(self->out_n_CL*x43) + x81) - self->head2_I_max*self->head2_diam*self->head2_kCL*x91*x95*(-self->head2_n_CL*self->head2_n_K + x57)/((x93 + 1)*(x94 + 1)*(x57*x91*x92 + 1) + (self->out_n_CL*x92 + 1)*(self->out_n_K*x91 + 1)*(x93*x94 + 1)) + x49*x89 - x83*(-x4*log(self->neck2_n_CL*x43) + x90);
    const arb_value_type x97 = -x34*(self->head2_n_CA - 5.0000000000000002e-5) + x62*x89 - x85*(x5*log(self->neck2_n_CA*x60) + x90);
    const arb_value_type x98 = x45*x88;
    const arb_value_type x99 = self->head2_kf*self->head2_n_B*self->head2_n_CA;
    const arb_value_type x100 = self->head2_kb*(self->head2_nTB - self->head2_n_B);
    const arb_value_type x101 = pow(self->head2_ca_presyn, 2);
    const arb_value_type grad_neck1_V = x0*x36*(x18 + x29 + x35)/self->neck1_diam;
    const arb_value_type grad_neck1_n_CL = -x38*x39*(x17 + x29);
    const arb_value_type grad_neck1_n_CA = x39*x40*(x13 + x35);
    const arb_value_type grad_head1_V = x36*x41*(x59 + x63)/self->head1_diam;
    const arb_value_type grad_head1_n_CL = -x38*x59*x64;
    const arb_value_type grad_head1_n_CA = x40*x64*(x23*x58*x67*(x65 - x66) + x63);
    const arb_value_type grad_head1_n_B = -2*x65 + 2*x66;
    const arb_value_type grad_head1_Y = -self->head1_Y/self->head1_tau_d;
    const arb_value_type grad_head1_X = (-self->head1_X - self->head1_Y + 1)/self->head1_tau_r;
    const arb_value_type grad_head1_ggaba = self->head1_Y*self->head1_gbar_gaba - self->head1_ggaba/(self->head1_eta_gaba/(exp((-self->head1_ca_presyn + self->head1_theta_gaba)/self->head1_sigma_gaba) + 1.0) + self->head1_tau_gaba0);
    const arb_value_type grad_head1_ca_presyn = self->head1_Ip - self->head1_beta*x68/(pow(self->head1_Kp, 2) + x68) + self->head1_gamma*log(2/self->head1_ca_presyn);
    const arb_value_type grad_neck2_V = x36*x69*(x79 + x84 + x86)/self->neck2_diam;
    const arb_value_type grad_neck2_n_CL = -x38*x87*(x78 + x84);
    const arb_value_type grad_neck2_n_CA = x40*x87*(x75 + x86);
    const arb_value_type grad_head2_V = x36*x88*(x96 + x97)/self->head2_diam;
    const arb_value_type grad_head2_n_CL = -x38*x96*x98;
    const arb_value_type grad_head2_n_CA = x40*x98*(x44*x67*x95*(-x100 + x99) + x97);
    const arb_value_type grad_head2_n_B = 2*x100 - 2*x99;
    const arb_value_type grad_head2_Y = -self->head2_Y/self->head2_tau_d;
    const arb_value_type grad_head2_X = (-self->head2_X - self->head2_Y + 1)/self->head2_tau_r;
    const arb_value_type grad_head2_ggaba = self->head2_Y*self->head2_gbar_gaba - self->head2_ggaba/(self->head2_eta_gaba/(exp((-self->head2_ca_presyn + self->head2_theta_gaba)/self->head2_sigma_gaba) + 1.0) + self->head2_tau_gaba0);
    const arb_value_type grad_head2_ca_presyn = self->head2_Ip - self->head2_beta*x101/(pow(self->head2_Kp, 2) + x101) + self->head2_gamma*log(2/self->head2_ca_presyn);


//update states
    self->I_inward_dend1 = x18;
    self->I_inward_dend2 = x79;

    self->neck1_V += dt * grad_neck1_V;
    self->neck1_n_CL += dt * grad_neck1_n_CL;
    self->neck1_n_CA += dt * grad_neck1_n_CA;

    self->head1_V += dt * grad_head1_V;
    self->head1_n_CL += dt * grad_head1_n_CL;
    self->head1_n_CA += dt * grad_head1_n_CA;
    self->head1_n_B += dt * grad_head1_n_B;
    self->head1_Y += dt * grad_head1_Y;
    self->head1_X += dt * grad_head1_X;
    self->head1_ggaba += dt * grad_head1_ggaba;
    self->head1_ca_presyn += dt * grad_head1_ca_presyn;

    self->neck2_V += dt * grad_neck2_V;
    self->neck2_n_CL += dt * grad_neck2_n_CL;
    self->neck2_n_CA += dt * grad_neck2_n_CA;

    self->head2_V += dt * grad_head2_V;
    self->head2_n_CL += dt * grad_head2_n_CL;
    self->head2_n_CA += dt * grad_head2_n_CA;
    self->head2_n_B += dt * grad_head2_n_B;
    self->head2_Y += dt * grad_head2_Y;
    self->head2_X += dt * grad_head2_X;
    self->head2_ggaba += dt * grad_head2_ggaba;
    self->head2_ca_presyn += dt * grad_head2_ca_presyn;
}