//
// Created by max on 08-07-21.
//

#include <arbor/mechanism_abi.h>

namespace arb::IOU_catalogue {
        void mechanism_cx36_gj_gpu_init_(arb_mechanism_ppack*);
        void mechanism_cx36_gj_gpu_advance_state_(arb_mechanism_ppack*);
        void mechanism_cx36_gj_gpu_compute_currents_(arb_mechanism_ppack*);
        void mechanism_cx36_gj_gpu_write_ions_(arb_mechanism_ppack*);
        void mechanism_cx36_gj_gpu_apply_events_(arb_mechanism_ppack*);
        void mechanism_cx36_gj_gpu_post_event_(arb_mechanism_ppack*);
    } // namespace arb

extern "C" {
arb_mechanism_interface* make_arb_IOU_catalogue_cx36_gj_interface_gpu() {
    static arb_mechanism_interface result;
    result.backend=arb_backend_kind_gpu;
    result.alignment=1;
    result.init_mechanism  =(arb_mechanism_method)arb::IOU_catalogue::mechanism_cx36_gj_gpu_init_;
    result.compute_currents=(arb_mechanism_method)arb::IOU_catalogue::mechanism_cx36_gj_gpu_compute_currents_;
    result.apply_events    =(arb_mechanism_method)arb::IOU_catalogue::mechanism_cx36_gj_gpu_apply_events_;
    result.advance_state   =(arb_mechanism_method)arb::IOU_catalogue::mechanism_cx36_gj_gpu_advance_state_;
    result.write_ions      =(arb_mechanism_method)arb::IOU_catalogue::mechanism_cx36_gj_gpu_write_ions_;
    result.post_event      =(arb_mechanism_method)arb::IOU_catalogue::mechanism_cx36_gj_gpu_post_event_;
    return &result;
}
};

