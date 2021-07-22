//
// Created by lennart on 22-07-21.
//

#include "iostream" //debugging
#include <cmath>
#include <arbor/mechanism_abi.h>

namespace arb::IOU_catalogue::kernel_adex {
    static constexpr unsigned simd_width_ = 0;

#define S(x) std::cout << #x << "\t\t" << x << std::endl;

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
[[maybe_unused]] auto* _pp_var_w                 = pp->state_vars[0];\
[[maybe_unused]] auto* _pp_var_C                 = pp->parameters[0];\
[[maybe_unused]] auto* _pp_var_gL                = pp->parameters[1];\
[[maybe_unused]] auto* _pp_var_EL                = pp->parameters[2];\
[[maybe_unused]] auto* _pp_var_VT                = pp->parameters[3];\
[[maybe_unused]] auto* _pp_var_DeltaT            = pp->parameters[4];\
[[maybe_unused]] auto* _pp_var_Vr                = pp->parameters[5];\
[[maybe_unused]] auto* _pp_var_Vcut              = pp->parameters[6];\
[[maybe_unused]] auto* _pp_var_tauw              = pp->parameters[7];\
[[maybe_unused]] auto* _pp_var_a                 = pp->parameters[8];\
[[maybe_unused]] auto* _pp_var_b                 = pp->parameters[9];\
[[maybe_unused]] auto* _pp_var_I                 = pp->parameters[10];\
//End of IFACEBLOCK

// interface methods
    static void init(arb_mechanism_ppack* pp) {
        PPACK_IFACE_BLOCK;
        std::cout << "init adex "<< std::endl;
        for (unsigned i = 0; i < _pp_var_width; i++) {
            const auto nidx = _pp_var_node_index[i];
            _pp_var_Vcut[i] = _pp_var_VT[i] + 5. *_pp_var_DeltaT[i];
            S(_pp_var_vec_v[nidx]);
            _pp_var_vec_i[nidx] = 0;
            //_pp_var_vec_v[nidx] = -50;
        }
    }

    static void compute_currents(arb_mechanism_ppack * pp) {
    }

    static void advance_state(arb_mechanism_ppack* pp) {
        /*
        // Mechanism weights are F·α where α ∈ [0, 1] is the proportional
        // contribution in the CV, and F is the scaling factor required
        // to convert from the mechanism current contribution units to A/m².
        switch (config.kind) {
        case arb_mechanism_kind_point:
            // Point mechanism contributions are in [nA]; CV area A in [µm^2].
            // F = 1/A * [nA/µm²] / [A/m²] = 1000/A.
            // */
        PPACK_IFACE_BLOCK;
        for (unsigned i = 0; i < _pp_var_width; i++) {
            const auto nidx = _pp_var_node_index[i];
            S(_pp_var_vec_v[nidx]);
            const arb_value_type vm = _pp_var_vec_v[nidx] * 1e-3 /*mV -> V*/;
            const arb_value_type dt = _pp_var_vec_dt[nidx] * 1e-3 /*ms -> s*/;
            /*  dvm/dt = (gL*(EL-vm)+gL*DeltaT*exp((vm-VT)/DeltaT) + I -w)/C : volt
             *  dw/dt = (a*(vm-EL)-w)/tauw : amp 
             * */
            const arb_value_type grad_w = (
                    _pp_var_a[i]*(vm-_pp_var_EL[i]) -
                    _pp_var_w[i]
                    )/_pp_var_tauw[i];
            arb_value_type il = _pp_var_gL[i]*(_pp_var_EL[i]-vm);
            arb_value_type id = _pp_var_gL[i]*_pp_var_DeltaT[i]*exp((vm-_pp_var_VT[i])/_pp_var_DeltaT[i]);
            const arb_value_type grad_vm = (
                    il +
                    id +
                    _pp_var_I[i] -_pp_var_w[i])/_pp_var_C[i];
            _pp_var_w[i] += 1e3 * dt * grad_w;
            _pp_var_vec_v[nidx] += 1e3 * dt*grad_vm;

            //plot -2.118e-09 5.27682 1.5e-09 4.90278e-14
            std::cout << "plot " << i << "=> " << nidx  << " / " << il << " " << id << " " << _pp_var_I[i] << " " << _pp_var_w[i] << std::endl;

            if (vm > _pp_var_Vcut[i]) {
                std::cout << "I'm triggered " << vm << " > " << _pp_var_Vcut[i] << std::endl;
                _pp_var_vec_v[i] = 1e3 * _pp_var_Vr[i];
                _pp_var_w[i] += _pp_var_b[i];
            }
        }
    }
    static void write_ions(arb_mechanism_ppack* pp) {
    }
    static void apply_events(arb_mechanism_ppack*pp) {
    }
    static void post_event(arb_mechanism_ppack*pp) {
    }

#undef PPACK_IFACE_BLOCK
}

extern "C" {
arb_mechanism_interface* make_arb_IOU_catalogue_adex_interface_multicore() {
    static arb_mechanism_interface result;
    result.partition_width = arb::IOU_catalogue::kernel_adex::simd_width_;
    result.backend=arb_backend_kind_cpu;
    result.alignment=1;
    result.init_mechanism  = (arb_mechanism_method)arb::IOU_catalogue::kernel_adex::init;
    result.compute_currents= (arb_mechanism_method)arb::IOU_catalogue::kernel_adex::compute_currents;
    result.apply_events    = (arb_mechanism_method)arb::IOU_catalogue::kernel_adex::apply_events;
    result.advance_state   = (arb_mechanism_method)arb::IOU_catalogue::kernel_adex::advance_state;
    result.write_ions      = (arb_mechanism_method)arb::IOU_catalogue::kernel_adex::write_ions;
    result.post_event      = (arb_mechanism_method)arb::IOU_catalogue::kernel_adex::post_event;
    return &result;
}}
