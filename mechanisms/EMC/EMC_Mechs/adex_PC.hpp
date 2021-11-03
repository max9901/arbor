#pragma once

#include <cmath>
#include <arbor/mechanism_abi.h>

extern "C" {
arb_mechanism_type make_arb_EMC_catalogue_adex_PC() {
    static arb_field_info globals[] = {};
    static arb_size_type n_globals = 0;
    static arb_field_info state_vars[] = {
            { "w"     , "nA", 0                        , -1.7976931348623157e+308, 1.7976931348623157e+308 },
    };
    static arb_size_type n_state_vars = 1;
    static arb_field_info parameters[] = {
            // 0
            { "C"        , "pF", 281e-12               , -1.7976931348623157e+308, 1.7976931348623157e+308 },
            // 1
            { "gL"       , "nS", 30e-9                 , -1.7976931348623157e+308, 1.7976931348623157e+308 },
            // 2
            { "EL"       , "mV", -70.6e-3              , -1.7976931348623157e+308, 1.7976931348623157e+308 },
            // 3
            { "VT"       , "mV", -50.4e-3              , -1.7976931348623157e+308, 1.7976931348623157e+308 },
            // 4
            { "DeltaT"   , "mV", 2e-3                  , -1.7976931348623157e+308, 1.7976931348623157e+308 },
            // 5
            { "Vr"       , "mV", -70.6e-3              , -1.7976931348623157e+308, 1.7976931348623157e+308 },
            // 6
            // Vcut=VT+5*DeltaT
            { "Vcut"     , "pF", NAN                   , -1.7976931348623157e+308, 1.7976931348623157e+308 },
            { "tauw"     , "ms", 144e-3                , -1.7976931348623157e+308, 1.7976931348623157e+308 },
            { "a"        , "nS", 4e-9                  , -1.7976931348623157e+308, 1.7976931348623157e+308 },
            { "b"        , "nA", 0.0805e-9             , -1.7976931348623157e+308, 1.7976931348623157e+308 },
            { "I"        , "nA", 1.5e-9                , -1.7976931348623157e+308, 1.7976931348623157e+308 },
    };
    static arb_size_type n_parameters = 11;
    static arb_ion_info ions[] = { };
    static arb_size_type n_ions = 0;

    arb_mechanism_type result;
    result.abi_version=ARB_MECH_ABI_VERSION;
    result.fingerprint="<placeholder>";
    result.name="adex_PC";
    result.kind=arb_mechanism_kind_point;
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

arb_mechanism_interface* make_arb_EMC_catalogue_adex_PC_interface_multicore();
arb_mechanism_interface* make_arb_EMC_catalogue_adex_PC_interface_gpu();
}
