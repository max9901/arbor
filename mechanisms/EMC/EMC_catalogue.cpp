#include <arbor/mechcat.hpp>
#include <arbor/mechanism.hpp>

#include "EMC_Mechs/adex.hpp"
#include "EMC_Mechs/ou_noise.hpp"
#include "EMC_Mechs/linear_gj.hpp"
#include "EMC_Mechs/cx36_gj.hpp"

namespace arb {

mechanism_catalogue build_EMC_catalogue() {
    mechanism_catalogue cat;

    //done !
    cat.add("ou_noise", make_arb_EMC_catalogue_ou_noise());
    cat.register_implementation("ou_noise", std::make_unique<mechanism>(make_arb_EMC_catalogue_ou_noise(), *make_arb_EMC_catalogue_ou_noise_interface_multicore()));
#ifdef ARB_HAVE_GPU
    cat.register_implementation("ou_noise", std::make_unique<mechanism>(make_arb_EMC_catalogue_ou_noise(), *make_arb_EMC_catalogue_ou_noise_interface_gpu()));
#endif

    //done !
    cat.add("linear_gj", make_arb_EMC_catalogue_linear_gj());
    cat.register_implementation("linear_gj", std::make_unique<mechanism>(make_arb_EMC_catalogue_linear_gj(), *make_arb_EMC_catalogue_linear_gj_interface_multicore()));
#ifdef ARB_HAVE_GPU
    cat.register_implementation("linear_gj", std::make_unique<mechanism>(make_arb_EMC_catalogue_linear_gj(), *make_arb_EMC_catalogue_linear_gj_interface_gpu()));
#endif

    //done !
    cat.add("cx36_gj", make_arb_EMC_catalogue_cx36_gj());
    cat.register_implementation("cx36_gj", std::make_unique<mechanism>(make_arb_EMC_catalogue_cx36_gj(), *make_arb_EMC_catalogue_cx36_gj_interface_multicore()));
#ifdef ARB_HAVE_GPU
    cat.register_implementation("cx36_gj", std::make_unique<mechanism>(make_arb_EMC_catalogue_cx36_gj(), *make_arb_EMC_catalogue_cx36_gj_interface_gpu()));
#endif

    // elias adex
    cat.add("adex", make_arb_EMC_catalogue_adex());
    cat.register_implementation("adex", std::make_unique<mechanism>(make_arb_EMC_catalogue_adex(), *make_arb_EMC_catalogue_adex_interface_multicore()));
#ifdef ARB_HAVE_GPU
    //todo
//    cat.register_implementation("ou_noise", std::make_unique<mechanism>(make_arb_EMC_catalogue_ou_noise(), *make_arb_EMC_catalogue_ou_noise_interface_gp\u()));
#endif

    return cat;
}

const mechanism_catalogue& global_EMC_catalogue() {
    static mechanism_catalogue cat = build_EMC_catalogue();
    return cat;
}

} // namespace arb

#ifdef STANDALONE
extern "C" {
    [[gnu::visibility("default")]] const void* get_catalogue() {
        static auto cat = arb::build_EMC_catalogue();
        return (void*)&cat;
    }
}
#endif


