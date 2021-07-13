#pragma once

#include <cmath>
#include <arbor/mechanism_abi.h>

extern "C" {
  arb_mechanism_type make_arb_IOU_catalogue_linear_gapJunction() {
    static arb_field_info globals[] = {};
    static arb_size_type n_globals = 0;
    static arb_field_info state_vars[] = {};
    static arb_size_type n_state_vars = 0;
    static arb_field_info parameters[] = {};
    static arb_size_type n_parameters = 0;
    static arb_ion_info ions[] = { };
    static arb_size_type n_ions = 0;

    arb_mechanism_type result;
    result.abi_version=ARB_MECH_ABI_VERSION;
    result.fingerprint="<placeholder>";
    result.name="linear_gapJunction";
    result.kind=arb_mechanism_kind_gap_junction;
    result.is_linear=false;
    result.has_post_events=false;
    result.globals=globals;
    result.n_globals=n_globals;
    result.ions=ions;
    result.n_ions=n_ions;
    result.state_vars=state_vars;
    result.n_state_vars=n_state_vars;
    result.parameters=parameters;
    result.n_parameters=n_parameters;
    return result;
  }

  arb_mechanism_interface* make_arb_IOU_catalogue_linear_gapJunction_interface_multicore();
  arb_mechanism_interface* make_arb_IOU_catalogue_linear_gapJunction_interface_gpu();
}
