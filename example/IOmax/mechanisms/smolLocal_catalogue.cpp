#include <arbor/mechcat.hpp>
#include <arbor/mechanism.hpp>
#include <arbor/mechanism_abi.h>

#include "smolLocal/ca_conc.hpp"

namespace arb {

mechanism_catalogue build_smolLocal_catalogue() {
    mechanism_catalogue cat;
    cat.add("ca_conc", make_arb_smolLocal_catalogue_ca_conc());
    cat.register_implementation("ca_conc", std::make_unique<mechanism>(make_arb_smolLocal_catalogue_ca_conc(), *make_arb_smolLocal_catalogue_ca_conc_interface_multicore()));
  return cat;
}

const mechanism_catalogue& global_smolLocal_catalogue() {
    static mechanism_catalogue cat = build_smolLocal_catalogue();
    return cat;
}

} // namespace arb

extern "C" {
    [[gnu::visibility("default")]] const void* get_catalogue() {
        static auto cat = arb::build_smolLocal_catalogue();
        return (void*)&cat;
    }
}


