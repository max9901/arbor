#pragma once

#include <cmath>
#include <arbor/mechanism_abi.h>

extern "C" {
  arb_mechanism_type make_arb_smolLocal_catalogue_ca_conc() {
    // Tables
    static arb_field_info globals[] = {  };
    static arb_size_type n_globals = 0;
    static arb_field_info state_vars[] = { { "extConcentration", "mM", NAN, -1.7976931348623157e+308, 1.7976931348623157e+308 },
                                           { "cai", "", NAN, -1.7976931348623157e+308, 1.7976931348623157e+308 } };
    static arb_size_type n_state_vars = 2;
    static arb_field_info parameters[] = { { "AREA_SCALE", "um2", 1000000000000.000000, -1.7976931348623157e+308, 1.7976931348623157e+308 },
                                           { "decayConstant", "ms", 33.333336, -1.7976931348623157e+308, 1.7976931348623157e+308 },
                                           { "LENGTH_SCALE", "um", 1000000.000000, -1.7976931348623157e+308, 1.7976931348623157e+308 },
                                           { "restingConc", "mM", 0.000000, -1.7976931348623157e+308, 1.7976931348623157e+308 },
                                           { "rate_concentration", "", NAN, -1.7976931348623157e+308, 1.7976931348623157e+308 },
                                           { "Faraday", "C / umol", 0.096485, -1.7976931348623157e+308, 1.7976931348623157e+308 },
                                           { "shellDepth", "um", 0.001000, -1.7976931348623157e+308, 1.7976931348623157e+308 },
                                           { "initialExtConcentration", "", NAN, -1.7976931348623157e+308, 1.7976931348623157e+308 },
                                           { "surfaceArea", "", NAN, -1.7976931348623157e+308, 1.7976931348623157e+308 },
                                           { "effectiveRadius", "", NAN, -1.7976931348623157e+308, 1.7976931348623157e+308 },
                                           { "eqshellDepth", "", NAN, -1.7976931348623157e+308, 1.7976931348623157e+308 },
                                           { "innerRadius", "", NAN, -1.7976931348623157e+308, 1.7976931348623157e+308 },
                                           { "shellVolume", "", NAN, -1.7976931348623157e+308, 1.7976931348623157e+308 } };
    static arb_size_type n_parameters = 13;
    static arb_ion_info ions[] = { { "ca", true, false, false, false, false, true, 2} };
    static arb_size_type n_ions = 1;

    arb_mechanism_type result;
    result.abi_version=ARB_MECH_ABI_VERSION;
    result.fingerprint="<placeholder>";
    result.name="ca_conc";
    result.kind=arb_mechanism_kind_density;
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

  arb_mechanism_interface* make_arb_smolLocal_catalogue_ca_conc_interface_multicore();
  arb_mechanism_interface* make_arb_smolLocal_catalogue_ca_conc_interface_gpu();
}
