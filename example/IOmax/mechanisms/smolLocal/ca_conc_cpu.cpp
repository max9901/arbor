#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <arbor/mechanism_abi.h>
#include <arbor/math.hpp>

namespace arb {
namespace smolLocal_catalogue {
namespace kernel_ca_conc {

using ::arb::math::exprelr;
using ::arb::math::safeinv;
using ::std::abs;
using ::std::cos;
using ::std::exp;
using ::std::log;
using ::std::max;
using ::std::min;
using ::std::pow;
using ::std::sin;

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
[[maybe_unused]] auto& _pp_var_index_constraints = pp->index_constraints;\
[[maybe_unused]] auto* _pp_var_extConcentration = pp->state_vars[0];\
[[maybe_unused]] auto* _pp_var_cai = pp->state_vars[1];\
[[maybe_unused]] auto* _pp_var_AREA_SCALE = pp->parameters[0];\
[[maybe_unused]] auto* _pp_var_decayConstant = pp->parameters[1];\
[[maybe_unused]] auto* _pp_var_LENGTH_SCALE = pp->parameters[2];\
[[maybe_unused]] auto* _pp_var_restingConc = pp->parameters[3];\
[[maybe_unused]] auto* _pp_var_rate_concentration = pp->parameters[4];\
[[maybe_unused]] auto* _pp_var_Faraday = pp->parameters[5];\
[[maybe_unused]] auto* _pp_var_shellDepth = pp->parameters[6];\
[[maybe_unused]] auto* _pp_var_initialExtConcentration = pp->parameters[7];\
[[maybe_unused]] auto* _pp_var_surfaceArea = pp->parameters[8];\
[[maybe_unused]] auto* _pp_var_effectiveRadius = pp->parameters[9];\
[[maybe_unused]] auto* _pp_var_eqshellDepth = pp->parameters[10];\
[[maybe_unused]] auto* _pp_var_innerRadius = pp->parameters[11];\
[[maybe_unused]] auto* _pp_var_shellVolume = pp->parameters[12];\
[[maybe_unused]] auto& _pp_var_ion_ca = pp->ion_states[0];\
[[maybe_unused]] auto* _pp_var_ion_ca_index = pp->ion_states[0].index;\
//End of IFACEBLOCK

// procedure prototypes
static void rates(arb_mechanism_ppack* pp, int i_, arb_value_type ica);

// interface methods
static void init(arb_mechanism_ppack* pp) {
    PPACK_IFACE_BLOCK;
    for (arb_size_type i_ = 0; i_ < _pp_var_width; ++i_) {
        auto node_indexi_ = _pp_var_node_index[i_];
        auto ion_ca_indexi_ = _pp_var_ion_ca_index[i_];
        arb_value_type ica = 0.10000000000000001*_pp_var_ion_ca.current_density[ion_ca_indexi_];
        arb_value_type diam = _pp_var_diam_um[node_indexi_];
        arb_value_type cao = _pp_var_ion_ca.external_concentration[ion_ca_indexi_];
        _pp_var_cai[i_] =  3.7151999999999998;
        _pp_var_initialExtConcentration[i_] = cao;
        _pp_var_surfaceArea[i_] =  3.1415899999999999*diam*diam* 0.25;
        rates(pp, i_, ica);
        rates(pp, i_, ica);
        _pp_var_extConcentration[i_] = _pp_var_initialExtConcentration[i_];
    }
    if (!_pp_var_multiplicity) return;
    for (arb_size_type ix = 0; ix < 2; ++ix) {
        for (arb_size_type iy = 0; iy < _pp_var_width; ++iy) {
            pp->state_vars[ix][iy] *= _pp_var_multiplicity[iy];
        }
    }
}

static void advance_state(arb_mechanism_ppack* pp) {
    PPACK_IFACE_BLOCK;
    for (arb_size_type i_ = 0; i_ < _pp_var_width; ++i_) {
        auto ion_ca_indexi_ = _pp_var_ion_ca_index[i_];
        auto node_indexi_ = _pp_var_node_index[i_];
        arb_value_type dt = _pp_var_vec_dt[node_indexi_];
        arb_value_type ica = 0.10000000000000001*_pp_var_ion_ca.current_density[ion_ca_indexi_];
        arb_value_type b_0_;
        rates(pp, i_, ica);
        b_0_ = _pp_var_rate_concentration[i_];
        _pp_var_cai[i_] = _pp_var_cai[i_]+b_0_*dt;
    }
}

static void compute_currents(arb_mechanism_ppack* pp) {
    PPACK_IFACE_BLOCK;
    for (arb_size_type i_ = 0; i_ < _pp_var_width; ++i_) {
        if (_pp_var_cai[i_]< 0.) {
            _pp_var_cai[i_] =  0.;
        }
    }
}

static void write_ions(arb_mechanism_ppack* pp) {
    PPACK_IFACE_BLOCK;
    for (arb_size_type i_ = 0; i_ < _pp_var_width; ++i_) {
        auto ion_ca_indexi_ = _pp_var_ion_ca_index[i_];
        arb_value_type cai_shadowed_ = 0;
        cai_shadowed_ = _pp_var_cai[i_];
        _pp_var_ion_ca.internal_concentration[ion_ca_indexi_] = fma(_pp_var_weight[i_], cai_shadowed_, _pp_var_ion_ca.internal_concentration[ion_ca_indexi_]);
    }
}

static void apply_events(arb_mechanism_ppack*) {}

static void post_event(arb_mechanism_ppack*) {}

// Procedure definitions
static void rates(arb_mechanism_ppack* pp, int i_, arb_value_type ica) {
    PPACK_IFACE_BLOCK;
    arb_value_type iCa;
    iCa =  -0.01*ica*_pp_var_surfaceArea[i_];
    _pp_var_effectiveRadius[i_] = _pp_var_LENGTH_SCALE[i_]*pow(_pp_var_surfaceArea[i_]/(_pp_var_AREA_SCALE[i_]* 12.56636),  0.5);
    _pp_var_eqshellDepth[i_] = _pp_var_shellDepth[i_]-_pp_var_shellDepth[i_]*_pp_var_shellDepth[i_]/_pp_var_effectiveRadius[i_];
    _pp_var_innerRadius[i_] = _pp_var_effectiveRadius[i_]-_pp_var_eqshellDepth[i_];
    _pp_var_shellVolume[i_] =  4.0*(_pp_var_effectiveRadius[i_]*_pp_var_effectiveRadius[i_]*_pp_var_effectiveRadius[i_])* 3.1415899999999999* 0.33333333333333331- 4.0*(_pp_var_innerRadius[i_]*_pp_var_innerRadius[i_]*_pp_var_innerRadius[i_])* 3.1415899999999999* 0.33333333333333331;
    _pp_var_rate_concentration[i_] = iCa/( 2.0*_pp_var_Faraday[i_]*_pp_var_shellVolume[i_])-(_pp_var_cai[i_]-_pp_var_restingConc[i_])/_pp_var_decayConstant[i_];
}
#undef PPACK_IFACE_BLOCK
} // namespace kernel_ca_conc
} // namespace smolLocal_catalogue
} // namespace arb

extern "C" {
  arb_mechanism_interface* make_arb_smolLocal_catalogue_ca_conc_interface_multicore() {
    static arb_mechanism_interface result;
    result.partition_width = arb::smolLocal_catalogue::kernel_ca_conc::simd_width_;
    result.backend=arb_backend_kind_cpu;
    result.alignment=1;
    result.init_mechanism=(arb_mechanism_method)arb::smolLocal_catalogue::kernel_ca_conc::init;
    result.compute_currents=(arb_mechanism_method)arb::smolLocal_catalogue::kernel_ca_conc::compute_currents;
    result.apply_events=(arb_mechanism_method)arb::smolLocal_catalogue::kernel_ca_conc::apply_events;
    result.advance_state=(arb_mechanism_method)arb::smolLocal_catalogue::kernel_ca_conc::advance_state;
    result.write_ions=(arb_mechanism_method)arb::smolLocal_catalogue::kernel_ca_conc::write_ions;
    result.post_event=(arb_mechanism_method)arb::smolLocal_catalogue::kernel_ca_conc::post_event;
    return &result;
  }}

