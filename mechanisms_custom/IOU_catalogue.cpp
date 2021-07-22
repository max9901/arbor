#include <arbor/mechcat.hpp>
#include <arbor/mechanism.hpp>

#include "IOU_Mechs/ou_noise.hpp"
#include "IOU_Mechs/glomerulus_gj.hpp"
#include "IOU_Mechs/linear_gj.hpp"
#include "IOU_Mechs/cx36_gj.hpp"
#include "IOU_Mechs/cx36k_gj.hpp"
#include "IOU_Mechs/smol_dend.hpp"

namespace arb {

mechanism_catalogue build_IOU_catalogue() {
    mechanism_catalogue cat;

    //done !
    cat.add("ou_noise", make_arb_IOU_catalogue_ou_noise());
    cat.register_implementation("ou_noise", std::make_unique<mechanism>(make_arb_IOU_catalogue_ou_noise(), *make_arb_IOU_catalogue_ou_noise_interface_multicore()));
#ifdef ARB_HAVE_GPU
    cat.register_implementation("ou_noise", std::make_unique<mechanism>(make_arb_IOU_catalogue_ou_noise(), *make_arb_IOU_catalogue_ou_noise_interface_gpu()));
#endif

    //not tested !
    cat.add("linear_gj", make_arb_IOU_catalogue_linear_gj());
    cat.register_implementation("linear_gj",std::make_unique<mechanism>(make_arb_IOU_catalogue_linear_gj(), *make_arb_IOU_catalogue_linear_gj_interface_multicore()));
#ifdef ARB_HAVE_GPU
    cat.register_implementation("linear_gj",std::make_unique<mechanism>(make_arb_IOU_catalogue_linear_gj(), *make_arb_IOU_catalogue_linear_gj_interface_gpu()));
#endif

    //not tested !
    cat.add("cx36_gj", make_arb_IOU_catalogue_cx36_gj());
    cat.register_implementation("cx36_gj",std::make_unique<mechanism>(make_arb_IOU_catalogue_cx36_gj(), *make_arb_IOU_catalogue_cx36_gj_interface_multicore()));
#ifdef ARB_HAVE_GPU
    cat.register_implementation("cx36_gj",std::make_unique<mechanism>(make_arb_IOU_catalogue_cx36_gj(), *make_arb_IOU_catalogue_cx36_gj_interface_gpu()));
#endif

    //not tested !
    cat.add("cx36k_gj", make_arb_IOU_catalogue_cx36k_gj());
    cat.register_implementation("cx36k_gj",std::make_unique<mechanism>(make_arb_IOU_catalogue_cx36k_gj(), *make_arb_IOU_catalogue_cx36k_gj_interface_multicore()));

    //cpu/gpu done waiting on validation
    cat.add("glomerulus_gj", make_arb_IOU_catalogue_glomerulus_gj());
    cat.register_implementation("glomerulus_gj", std::make_unique<mechanism>(make_arb_IOU_catalogue_glomerulus_gj(), *make_arb_IOU_catalogue_glomerulus_gj_interface_multicore()));
#ifdef ARB_HAVE_GPU
    cat.register_implementation("glomerulus_gj", std::make_unique<mechanism>(make_arb_IOU_catalogue_glomerulus_gj(), *make_arb_IOU_catalogue_glomerulus_gj_interface_gpu()));
#endif

    //cpu/gpu done waiting on validation
    cat.add("smol_dend", make_arb_IOU_catalogue_smol_dend());
    cat.register_implementation("smol_dend", std::make_unique<mechanism>(make_arb_IOU_catalogue_smol_dend(), *make_arb_IOU_catalogue_smol_dend_interface_multicore()));
#ifdef ARB_HAVE_GPU
    cat.register_implementation("smol_dend", std::make_unique<mechanism>(make_arb_IOU_catalogue_smol_dend(), *make_arb_IOU_catalogue_smol_dend_interface_gpu()));
#endif

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


