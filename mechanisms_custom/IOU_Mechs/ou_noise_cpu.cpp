//
// Created by max on 08-07-21.
//

#include <arbor/mechanism_abi.h>
#include <cmath>
#include <Random123/philox.h>
#include <Random123/boxmuller.hpp>

namespace arb::IOU_catalogue::kernel_ou_noise {
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
                          \
[[maybe_unused]] auto& _pp_var_theta  = pp->globals[0];        \
[[maybe_unused]] auto& _pp_var_sigma  = pp->globals[1];        \
[[maybe_unused]] auto& _pp_var_mu     = pp->globals[2];        \
[[maybe_unused]] auto& _pp_var_alpha  = pp->globals[3];        \
[[maybe_unused]] auto& _pp_var_seed   = pp->globals[4];        \
[[maybe_unused]] auto& _pp_var_cnt    = pp->globals[5];        \
[[maybe_unused]] auto* _pp_var_ouNoise  = pp->state_vars[0];   \
//End of IFACEBLOCK

// interface methods
static void init(arb_mechanism_ppack* pp) {
    PPACK_IFACE_BLOCK;
    _pp_var_cnt = 0.;
    for (arb_size_type i_ = 0; i_ < _pp_var_width; ++i_) {
        _pp_var_ouNoise[i_] = 0;
    }
}

static void compute_currents(arb_mechanism_ppack* pp) {
    PPACK_IFACE_BLOCK;
    philox2x32_key_t k = {{(uint32_t)_pp_var_seed}};
    philox2x32_ctr_t c = {{(uint32_t)_pp_var_cnt}};
    philox2x32_ctr_t cresult = philox2x32(c, k);
    const float rand_global = r123::boxmuller(cresult.v[0],cresult.v[1]).x;
    for (arb_size_type i_ = 0; i_ < _pp_var_width; ++i_) {
        c.v[0] = _pp_var_cnt+i_+1;
        const auto node_indexi_ = _pp_var_node_index[i_];
        const arb_value_type dt = _pp_var_vec_dt[node_indexi_];
        const arb_value_type sqrt_dt = std::sqrt(dt);
        const arb_value_type Iapp_global = _pp_var_sigma * sqrt_dt * rand_global;
        cresult = philox2x32(c, k);
        const float rand_local = r123::boxmuller(cresult.v[0],cresult.v[1]).x;
        _pp_var_ouNoise[i_] +=
                _pp_var_theta * (_pp_var_mu - _pp_var_ouNoise[i_]) * dt +
                (1-_pp_var_alpha) * _pp_var_sigma * sqrt_dt * rand_local
                + _pp_var_alpha * Iapp_global;
        _pp_var_vec_i[node_indexi_] -= _pp_var_ouNoise[i_];
    }
    _pp_var_cnt += _pp_var_width+1;
}

static void advance_state(arb_mechanism_ppack* pp) {}
static void write_ions(arb_mechanism_ppack* pp) {}
static void apply_events(arb_mechanism_ppack*, arb_deliverable_event_stream*) {}
static void post_event(arb_mechanism_ppack*) {}

#undef PPACK_IFACE_BLOCK
}

extern "C" {
  arb_mechanism_interface* make_arb_IOU_catalogue_ou_noise_interface_multicore() {
    static arb_mechanism_interface result;
    result.partition_width = arb::IOU_catalogue::kernel_ou_noise::simd_width_;
    result.backend=arb_backend_kind_cpu;
    result.alignment=1;
    result.init_mechanism  = (arb_mechanism_method)arb::IOU_catalogue::kernel_ou_noise::init;
    result.compute_currents= (arb_mechanism_method)arb::IOU_catalogue::kernel_ou_noise::compute_currents;
    result.apply_events    = (arb_mechanism_method_events)arb::IOU_catalogue::kernel_ou_noise::apply_events;
    result.advance_state   = (arb_mechanism_method)arb::IOU_catalogue::kernel_ou_noise::advance_state;
    result.write_ions      = (arb_mechanism_method)arb::IOU_catalogue::kernel_ou_noise::write_ions;
    result.post_event      = (arb_mechanism_method)arb::IOU_catalogue::kernel_ou_noise::post_event;
    return &result;
  }}
