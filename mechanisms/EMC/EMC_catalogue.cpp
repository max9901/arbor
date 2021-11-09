#include <arbor/mechcat.hpp>
#include <arbor/mechanism.hpp>

#include "EMC_Mechs/adex.hpp"
#include "EMC_Mechs/adex_DCN.hpp"
#include "EMC_Mechs/adex_PC.hpp"
#include "EMC_Mechs/ou_noise.hpp"
#include "EMC_Mechs/gj_linear.hpp"


namespace arb {

mechanism_catalogue build_EMC_catalogue() {
    mechanism_catalogue cat;

//    done !
    cat.add("ou_noise", make_arb_EMC_catalogue_ou_noise());
    cat.register_implementation("ou_noise", std::make_unique<mechanism>(make_arb_EMC_catalogue_ou_noise(), *make_arb_EMC_catalogue_ou_noise_interface_multicore()));
#ifdef ARB_HAVE_GPU
    cat.register_implementation("ou_noise", std::make_unique<mechanism>(make_arb_EMC_catalogue_ou_noise(), *make_arb_EMC_catalogue_ou_noise_interface_gpu()));
#endif

    //done !
    cat.add("gj_linear", make_arb_EMC_catalogue_gj_linear());
    cat.register_implementation("gj_linear", std::make_unique<mechanism>(make_arb_EMC_catalogue_gj_linear(), *make_arb_EMC_catalogue_gj_linear_interface_multicore()));
#ifdef ARB_HAVE_GPU
    cat.register_implementation("gj_linear", std::make_unique<mechanism>(make_arb_EMC_catalogue_gj_linear(), *make_arb_EMC_catalogue_gj_linear_interface_gpu()));
#endif

    // elias adex
    cat.add("adex", make_arb_EMC_catalogue_adex());
    cat.add("adex_DCN", make_arb_EMC_catalogue_adex_DCN());
    cat.add("adex_PC", make_arb_EMC_catalogue_adex_PC());

    cat.register_implementation("adex", std::make_unique<mechanism>(make_arb_EMC_catalogue_adex(), *make_arb_EMC_catalogue_adex_interface_multicore()));
    cat.register_implementation("adex_DCN", std::make_unique<mechanism>(make_arb_EMC_catalogue_adex_DCN(), *make_arb_EMC_catalogue_adex_DCN_interface_multicore()));
    cat.register_implementation("adex_PC", std::make_unique<mechanism>(make_arb_EMC_catalogue_adex_PC(), *make_arb_EMC_catalogue_adex_PC_interface_multicore()));
#ifdef ARB_HAVE_GPU
    cat.register_implementation("adex", std::make_unique<mechanism>(make_arb_EMC_catalogue_adex(), *make_arb_EMC_catalogue_adex_interface_gpu()));
//    printf("\n\n[WARNING !!!] ADEX not yet implemented on GPU's\n\n");
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


