#include <arbor/mechcat.hpp>
#include <arbor/mechanism.hpp>
#include <arbor/mechanism_abi.h>

#include "IOU_Mechs/ou_noise.hpp"
#include "IOU_Mechs/glomerulus.hpp"
#include "IOU_Mechs/linear_gapJunction.hpp"

namespace arb {

mechanism_catalogue build_IOU_catalogue() {
    mechanism_catalogue cat;

    //done !
    cat.add("ou_noise", make_arb_IOU_catalogue_ou_noise());
    cat.register_implementation("ou_noise", std::make_unique<mechanism>(make_arb_IOU_catalogue_ou_noise(), *make_arb_IOU_catalogue_ou_noise_interface_multicore()));
    cat.register_implementation("ou_noise", std::make_unique<mechanism>(make_arb_IOU_catalogue_ou_noise(), *make_arb_IOU_catalogue_ou_noise_interface_gpu()));

    //not tested !
    cat.add("linear_gapJunction", make_arb_IOU_catalogue_linear_gapJunction());
    cat.register_implementation("linear_gapJunction",std::make_unique<mechanism>(make_arb_IOU_catalogue_linear_gapJunction(), *make_arb_IOU_catalogue_linear_gapJunction_interface_multicore()));
    cat.register_implementation("linear_gapJunction",std::make_unique<mechanism>(make_arb_IOU_catalogue_linear_gapJunction(), *make_arb_IOU_catalogue_linear_gapJunction_interface_gpu()));

    //cpu done waiting on validation
    cat.add("glomerulus", make_arb_IOU_catalogue_glomerulus());
    cat.register_implementation("glomerulus", std::make_unique<mechanism>(make_arb_IOU_catalogue_glomerulus(), *make_arb_IOU_catalogue_glomerulus_interface_multicore()));
//    cat.register_implementation("glomerulus", std::make_unique<mechanism>(make_arb_IOU_catalogue_glomerulus(), *make_arb_IOU_catalogue_glomerulus_interface_gpu()));

  return cat;
}

const mechanism_catalogue& global_IOU_catalogue() {
    static mechanism_catalogue cat = build_IOU_catalogue();
    return cat;
}

} // namespace arb
#ifdef STANDALONE
extern "C" {
    [[gnu::visibility("default")]] const void* get_catalogue() {
        static auto cat = arb::build_IOU_catalogue();
        return (void*)&cat;
    }
}
#endif


