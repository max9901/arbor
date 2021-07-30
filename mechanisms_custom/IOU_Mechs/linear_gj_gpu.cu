//
// Created by max on 08-07-21.
//

#include <arbor/gpu/gpu_common.hpp>
#include <arbor/mechanism_abi.h>

namespace arb {
    namespace IOU_catalogue {
        
#define PPACK_IFACE_BLOCK \
auto* _pp_var_vec_i              __attribute__((unused)) = params_.vec_i;\
auto* _pp_var_vec_v              __attribute__((unused)) = params_.vec_v;\
auto& _pp_var_gap_junctions      __attribute__((unused)) = params_.gap_junctions;\
auto& _pp_var_gap_junction_width __attribute__((unused)) = params_.gap_junction_width; \
//End of IFACEBLOCK

namespace {
    __global__ void compute_currents(arb_mechanism_ppack params_) {
        PPACK_IFACE_BLOCK;
        const size_t tid_ = threadIdx.x + blockDim.x*blockIdx.x;
        if(tid_ < _pp_var_gap_junction_width){
            const auto gj   = _pp_var_gap_junctions[tid_];
            const auto curr = gj.weight * (_pp_var_vec_v[gj.loc.second] - _pp_var_vec_v[gj.loc.first]); // nA
            _pp_var_vec_i[gj.loc.first] -= curr;
        }
    }
} // namespace

void mechanism_linear_gj_gpu_init_(arb_mechanism_ppack*) {}
void mechanism_linear_gj_gpu_compute_currents_(arb_mechanism_ppack* p) {
    unsigned block_dim = 128;
    unsigned grid_dim = ::arb::gpu::impl::block_count(p->gap_junction_width, block_dim);
    compute_currents<<<grid_dim, block_dim>>>(*p);
}
void mechanism_linear_gj_gpu_advance_state_(arb_mechanism_ppack*) {}
void mechanism_linear_gj_gpu_write_ions_(arb_mechanism_ppack*) {}
void mechanism_linear_gj_gpu_post_event_(arb_mechanism_ppack*) {}
void mechanism_linear_gj_gpu_apply_events_(arb_mechanism_ppack*, arb_deliverable_event_stream*) {}

} // namespace IOU_catalogue
} // namespace arb
