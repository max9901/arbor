#include <arbor/mechcat.hpp>
#include <arbor/mechanism.hpp>
#include <arbor/mechanism_abi.h>

#include "LocalMechs/ou_noise.hpp"

namespace arb {

mechanism_catalogue build_Local_catalogue() {
    mechanism_catalogue cat;
    cat.add("ou_noise", make_arb_Local_catalogue_ou_noise());
    cat.register_implementation("ou_noise", std::make_unique<mechanism>(make_arb_Local_catalogue_ou_noise(), *make_arb_Local_catalogue_ou_noise_interface_multicore()));
  return cat;
}

const mechanism_catalogue& global_Local_catalogue() {
    static mechanism_catalogue cat = build_Local_catalogue();
    return cat;
}

} // namespace arb

extern "C" {
    [[gnu::visibility("default")]] const void* get_catalogue() {
        static auto cat = arb::build_Local_catalogue();
        return (void*)&cat;
    }
}


