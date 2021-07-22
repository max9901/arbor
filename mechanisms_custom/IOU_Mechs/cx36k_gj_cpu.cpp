//
// Created by lennart on 09-07-21.
//
//
//Test of cx36 differences does it have effect?
// So  just cx36 but 1000 more weight

#include "iostream" //debugging
#include <arbor/fvm_types.hpp>
#include <arbor/mechanism_abi.h>
#include <cmath>

namespace arb::IOU_catalogue::kernel_cx36k_gj {
    static constexpr unsigned simd_width_ = 0;

#define S(x) std::cout << #x << "\t\t" << x << std::endl;

#define PPACK_IFACE_BLOCK \
[[maybe_unused]] auto* _pp_var_vec_v             = pp->vec_v;\
[[maybe_unused]] auto* _pp_var_vec_i             = pp->vec_i;\
[[maybe_unused]] auto* _pp_var_gap_junctions     = pp->gap_junctions; \
[[maybe_unused]] auto _pp_var_gap_junction_width = pp->gap_junction_width;\
//End of IFACEBLOCK

// interface methods
    static void init(arb_mechanism_ppack* pp) {
        std::cout << "init cx36k_gj "<< std::endl;
    }

    static void compute_currents(arb_mechanism_ppack* pp) {
        PPACK_IFACE_BLOCK;
        for (unsigned i = 0; i < _pp_var_gap_junction_width; i++) {
            auto gj = _pp_var_gap_junctions[i];
            const auto vDiff = (_pp_var_vec_v[gj.loc.second] - _pp_var_vec_v[gj.loc.first]);
            const auto fAcc = vDiff  * gj.weight * exp(vDiff * vDiff * (-0.01));
            const auto vAcc = gj.weight * vDiff;
            auto curr = 0.8 * fAcc + 0.2 * vAcc;
            _pp_var_vec_i[gj.loc.first] -= 1000 * curr;
        }
    }

    static void advance_state(arb_mechanism_ppack* pp) {}
    static void write_ions(arb_mechanism_ppack* pp) {}
    static void apply_events(arb_mechanism_ppack*) {}
    static void post_event(arb_mechanism_ppack*) {}

#undef PPACK_IFACE_BLOCK
}

extern "C" {
arb_mechanism_interface* make_arb_IOU_catalogue_cx36k_gj_interface_multicore() {
    static arb_mechanism_interface result;
    result.partition_width = arb::IOU_catalogue::kernel_cx36k_gj::simd_width_;
    result.backend=arb_backend_kind_cpu;
    result.alignment=1;
    result.init_mechanism  = (arb_mechanism_method)arb::IOU_catalogue::kernel_cx36k_gj::init;
    result.compute_currents= (arb_mechanism_method)arb::IOU_catalogue::kernel_cx36k_gj::compute_currents;
    result.apply_events    = (arb_mechanism_method)arb::IOU_catalogue::kernel_cx36k_gj::apply_events;
    result.advance_state   = (arb_mechanism_method)arb::IOU_catalogue::kernel_cx36k_gj::advance_state;
    result.write_ions      = (arb_mechanism_method)arb::IOU_catalogue::kernel_cx36k_gj::write_ions;
    result.post_event      = (arb_mechanism_method)arb::IOU_catalogue::kernel_cx36k_gj::post_event;
    return &result;
}}