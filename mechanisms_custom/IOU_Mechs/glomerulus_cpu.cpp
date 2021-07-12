//
// Created by max on 09-07-21.
//

#include "iostream" //debugging
#include <arbor/mechanism_abi.h>
#include <cmath>

#define S(x) std::cout << #x << "\t\t" << x << std::endl;

namespace arb::IOU_catalogue::kernel_glomerulus {
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
//End of IFACEBLOCK


// interface methods
    static void init(arb_mechanism_ppack* pp) {
        PPACK_IFACE_BLOCK;


        std::cout << "init glumerelus ding\n";

        arb_constraint_partition hoi;


    }

    static void compute_currents(arb_mechanism_ppack* pp) {
        PPACK_IFACE_BLOCK;
    }

    static void advance_state(arb_mechanism_ppack* pp) {}
    static void write_ions(arb_mechanism_ppack* pp) {}
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
