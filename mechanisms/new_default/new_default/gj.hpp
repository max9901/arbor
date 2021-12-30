#pragma once

#include <cmath>
#include <arbor/mechanism_abi.h>

extern "C" {
  arb_mechanism_type make_arb_new_default_catalogue_gj() {
    // Tables
    static arb_field_info* globals = NULL;
    static arb_size_type n_globals = 0;
    static arb_field_info* state_vars = NULL;
    static arb_size_type n_state_vars = 0;
    static arb_field_info parameters[] = {{ "g", "", 1, -179769313486231570814527423731704356798070567525844996598917476803157260780028538760589558632766878171540458953514382464234321326889464182768467546703537516986049910576551282076245490090389328944075868508455133942304583236903222948165808559332123348274797826204144723168738177180919299881250404026184124858368.000000, 179769313486231570814527423731704356798070567525844996598917476803157260780028538760589558632766878171540458953514382464234321326889464182768467546703537516986049910576551282076245490090389328944075868508455133942304583236903222948165808559332123348274797826204144723168738177180919299881250404026184124858368.000000 } };
    static arb_size_type n_parameters = 1;
    static arb_ion_info* ions = NULL;
    static arb_size_type n_ions = 0;

    arb_mechanism_type result;
    result.abi_version=ARB_MECH_ABI_VERSION;
    result.fingerprint="<placeholder>";
    result.name="gj";
    result.kind=arb_mechanism_kind_gap_junction;
    result.is_linear=true;
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

  arb_mechanism_interface* make_arb_new_default_catalogue_gj_interface_multicore();
  arb_mechanism_interface* make_arb_new_default_catalogue_gj_interface_gpu();
}
