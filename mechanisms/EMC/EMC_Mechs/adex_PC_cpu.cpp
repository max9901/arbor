//
// Created by lennart/max/elias on 22-07-21.
//

#include <cmath>
#include <arbor/mechanism_abi.h>

namespace arb::EMC_catalogue::kernel_adex_PC {
    static constexpr unsigned simd_width_ = 1;

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
        for (unsigned i = 0; i < _pp_var_width; i++) {
            const auto nidx = _pp_var_node_index[i];
            _pp_var_Vcut[i] = _pp_var_VT[i] + 5. *_pp_var_DeltaT[i];
            _pp_var_vec_i[nidx] = 0;
            //_pp_var_vec_v[nidx] = -50;
        }
    }

    static void compute_currents(arb_mechanism_ppack * pp) {
    }

    static void advance_state(arb_mechanism_ppack* pp) {
        PPACK_IFACE_BLOCK;
        for (unsigned i = 0; i < _pp_var_width; i++) {
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
    static void write_ions(arb_mechanism_ppack* pp) {
    }
    static void apply_events(arb_mechanism_ppack*pp, arb_deliverable_event_stream* stream_ptr) {
        PPACK_IFACE_BLOCK;
        auto ncell = stream_ptr->n_streams;
        for (arb_size_type c = 0; c<ncell; ++c) {
            auto begin  = stream_ptr->events + stream_ptr->begin[c];
            auto end    = stream_ptr->events + stream_ptr->end[c];
            for (auto p = begin; p<end; ++p) {
                auto i_     = p->mech_index;
                auto weight = p->weight;
                if (p->mech_id==_pp_var_mechanism_id) {
                    _pp_var_w[i_] = _pp_var_w[i_] + weight;
                }
            }
        }
    }

    static void post_event(arb_mechanism_ppack*pp) {

    }

#undef PPACK_IFACE_BLOCK
}

extern "C" {
arb_mechanism_interface* make_arb_EMC_catalogue_adex_PC_interface_multicore() {
    static arb_mechanism_interface result;
    result.partition_width = arb::EMC_catalogue::kernel_adex_PC::simd_width_;
    result.backend=arb_backend_kind_cpu;
    result.alignment=1;
    result.init_mechanism  = (arb_mechanism_method)arb::EMC_catalogue::kernel_adex_PC::init;
    result.compute_currents= (arb_mechanism_method)arb::EMC_catalogue::kernel_adex_PC::compute_currents;
    result.apply_events    = (arb_mechanism_method_events)arb::EMC_catalogue::kernel_adex_PC::apply_events;
    result.advance_state   = (arb_mechanism_method)arb::EMC_catalogue::kernel_adex_PC::advance_state;
    result.write_ions      = (arb_mechanism_method)arb::EMC_catalogue::kernel_adex_PC::write_ions;
    result.post_event      = (arb_mechanism_method)arb::EMC_catalogue::kernel_adex_PC::post_event;
    return &result;
}}
