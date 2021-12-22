#include <memory>
#include <stdexcept>

#include <arbor/arbexcept.hpp>
#include <arbor/common_types.hpp>

#include "backends/multicore/fvm.hpp"
#ifdef ARB_HAVE_GPU
#include "backends/gpu/fvm.hpp"
#endif
#include "fvm_lowered_distributed_cell_impl.hpp"

namespace arb {

fvm_lowered_cell_ptr make_fvm_lowered_distributed_cell(backend_kind p, const execution_context& ctx, const group_description& groupDescription) {
    switch (p) {
    case backend_kind::multicore:
        return fvm_lowered_cell_ptr(new fvm_lowered_distributed_cell_impl<multicore::backend>(ctx, groupDescription));
    case backend_kind::gpu:
#ifdef ARB_HAVE_GPU
        return fvm_lowered_cell_ptr(new fvm_lowered_distributed_cell_impl<gpu::backend>(ctx, groupDescription));
#endif
        ; // fall through
    default:
        throw arbor_internal_error("fvm_lowered_cell: unsupported back-end");
    }
}

} // namespace arb
