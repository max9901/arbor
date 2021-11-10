#include <arbor/mechanism_abi.h>

namespace arb {
    namespace EMC_catalogue {
        void mechanism_adex_DCN_gpu_init_(arb_mechanism_ppack*);
        void mechanism_adex_DCN_gpu_advance_state_(arb_mechanism_ppack*);
        void mechanism_adex_DCN_gpu_compute_currents_(arb_mechanism_ppack*);
        void mechanism_adex_DCN_gpu_write_ions_(arb_mechanism_ppack*);
        void mechanism_adex_DCN_gpu_apply_events_(arb_mechanism_ppack*, arb_deliverable_event_stream*);
        void mechanism_adex_DCN_gpu_post_event_(arb_mechanism_ppack*);

    } // namespace EMC_catalogue
} // namespace arb

extern "C" {
arb_mechanism_interface* make_arb_EMC_catalogue_adex_DCN_interface_gpu() {
    static arb_mechanism_interface result;
    result.backend=arb_backend_kind_gpu;
    result.partition_width=1;
    result.alignment=1;
    result.init_mechanism=       (arb_mechanism_method)arb::EMC_catalogue::mechanism_adex_DCN_gpu_init_;
    result.compute_currents=     (arb_mechanism_method)arb::EMC_catalogue::mechanism_adex_DCN_gpu_compute_currents_;
    result.apply_events=  (arb_mechanism_method_events)arb::EMC_catalogue::mechanism_adex_DCN_gpu_apply_events_;
    result.advance_state=        (arb_mechanism_method)arb::EMC_catalogue::mechanism_adex_DCN_gpu_advance_state_;
    result.write_ions=           (arb_mechanism_method)arb::EMC_catalogue::mechanism_adex_DCN_gpu_write_ions_;
    result.post_event=           (arb_mechanism_method)arb::EMC_catalogue::mechanism_adex_DCN_gpu_post_event_;
    return &result;
}
};

