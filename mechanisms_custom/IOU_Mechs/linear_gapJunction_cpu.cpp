//
// Created by max on 09-07-21.
//

#include "iostream" //debugging
#include <arbor/fvm_types.hpp>
#include <arbor/mechanism_abi.h>
#include <cmath>

namespace arb::IOU_catalogue::kernel_linear_gapJunction {
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
[[maybe_unused]] auto* _pp_var_gap_junctions = pp->gap_junctions; \
[[maybe_unused]] auto _pp_var_gap_junction_width = pp->gap_junction_width; \
//End of IFACEBLOCK

// interface methods
    static void init(arb_mechanism_ppack* pp) {
        std::cout << "init linear_gapJunction "<< std::endl;
    }

    static void compute_currents(arb_mechanism_ppack* pp) {
        PPACK_IFACE_BLOCK;
        for (unsigned i = 0; i < _pp_var_gap_junction_width; i++) {
            auto gj = _pp_var_gap_junctions[i];
            auto curr = gj.weight * (_pp_var_vec_v[gj.loc.second] - _pp_var_vec_v[gj.loc.first]); // nA
            _pp_var_vec_i[gj.loc.first] -= curr;      //no need for weiht right now todo is change the way it handles weights put them in the weight vector. fixed in gj.weight
        }
    }

    static void advance_state(arb_mechanism_ppack* pp) {}
    static void write_ions(arb_mechanism_ppack* pp) {}
    static void apply_events(arb_mechanism_ppack*) {}
    static void post_event(arb_mechanism_ppack*) {}

#undef PPACK_IFACE_BLOCK
}

extern "C" {
arb_mechanism_interface* make_arb_IOU_catalogue_linear_gapJunction_interface_multicore() {
    static arb_mechanism_interface result;
    result.partition_width = arb::IOU_catalogue::kernel_linear_gapJunction::simd_width_;
    result.backend=arb_backend_kind_cpu;
    result.alignment=1;
    result.init_mechanism  = (arb_mechanism_method)arb::IOU_catalogue::kernel_linear_gapJunction::init;
    result.compute_currents= (arb_mechanism_method)arb::IOU_catalogue::kernel_linear_gapJunction::compute_currents;
    result.apply_events    = (arb_mechanism_method)arb::IOU_catalogue::kernel_linear_gapJunction::apply_events;
    result.advance_state   = (arb_mechanism_method)arb::IOU_catalogue::kernel_linear_gapJunction::advance_state;
    result.write_ions      = (arb_mechanism_method)arb::IOU_catalogue::kernel_linear_gapJunction::write_ions;
    result.post_event      = (arb_mechanism_method)arb::IOU_catalogue::kernel_linear_gapJunction::post_event;
    return &result;
}}
