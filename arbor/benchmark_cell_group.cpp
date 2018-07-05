#include <chrono>
#include <exception>

#include <arbor/benchmark_cell.hpp>
#include <arbor/recipe.hpp>
#include <arbor/time_sequence.hpp>

#include "cell_group.hpp"
#include "profile/profiler_macro.hpp"
#include "benchmark_cell_group.hpp"

#include "util/span.hpp"

namespace arb {

benchmark_cell_group::benchmark_cell_group(std::vector<cell_gid_type> gids,
                                           const recipe& rec):
    gids_(std::move(gids))
{
    cells_.reserve(gids_.size());
    for (auto gid: gids_) {
        cells_.push_back(util::any_cast<benchmark_cell>(rec.get_cell_description(gid)));
    }

    reset();
}

void benchmark_cell_group::reset() {
    t_ = 0;

    for (auto& c: cells_) {
        c.time_sequence.reset();
    }

    clear_spikes();
}

cell_kind benchmark_cell_group::get_cell_kind() const {
    return cell_kind::benchmark;
}

void benchmark_cell_group::advance(epoch ep,
                                   time_type dt,
                                   const event_lane_subrange& event_lanes)
{
    using std::chrono::high_resolution_clock;
    using duration_type = std::chrono::duration<double, std::micro>;

    PE(advance_bench_cell);
    // Micro-seconds to advance in this epoch.
    auto us = 1e3*(ep.tfinal-t_);
    for (auto i: util::make_span(0, gids_.size())) {
        auto& tseq = cells_[i].time_sequence;
        // Expected time to complete epoch in micro seconds.
        const double duration_us = cells_[i].realtime_ratio*us;
        const auto gid = gids_[i];

        // Start timer.
        auto start = high_resolution_clock::now();

        while (tseq.front()<ep.tfinal) {
            spikes_.push_back({{gid, 0u}, tseq.front()});
            tseq.pop();
        }

        // Wait until the expected time to advance has elapsed. Use a busy-wait
        // so that the resources of this thread are tied up until the interval
        // has elapsed, to emulate a "real" cell.
        while (duration_type(high_resolution_clock::now()-start).count() < duration_us);
    }
    t_ = ep.tfinal;

    PL();
};

const std::vector<spike>& benchmark_cell_group::spikes() const {
    return spikes_;
}

void benchmark_cell_group::clear_spikes() {
    spikes_.clear();
}

void benchmark_cell_group::add_sampler(sampler_association_handle h,
                                   cell_member_predicate probe_ids,
                                   schedule sched,
                                   sampler_function fn,
                                   sampling_policy policy)
{
    std::logic_error("A benchmark_cell group doen't support sampling of internal state!");
}

} // namespace arb