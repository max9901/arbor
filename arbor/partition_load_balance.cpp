#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>


#include <arbor/arbexcept.hpp>
#include <arbor/domain_decomposition.hpp>
#include <arbor/load_balance.hpp>
#include <arbor/recipe.hpp>
#include <arbor/symmetric_recipe.hpp>
#include <arbor/context.hpp>

#include "cell_group_factory.hpp"
#include "execution_context.hpp"
#include "gpu_context.hpp"
#include "util/maputil.hpp"
#include "util/partition.hpp"
#include "util/span.hpp"
#include "util/strprintf.hpp"

#define DEBUGMAX
#ifdef DEBUGMAX
    #include <iostream>
    #include <mpi.h>
#endif
namespace arb {

domain_decomposition partition_load_balance(
    const recipe& rec,
    const context& ctx,
    partition_hint_map hint_map)
{
    const bool gpu_avail = ctx->gpu->has_gpu();

    struct partition_gid_domain {
        partition_gid_domain(const gathered_vector<cell_gid_type>& divs, unsigned domains) {
            auto rank_part = util::partition_view(divs.partition());
            for (auto rank: count_along(rank_part)) {
                for (auto gid: util::subrange_view(divs.values(), rank_part[rank])) {
                    gid_map[gid] = rank;
                }
            }
        }

        int operator()(cell_gid_type gid) const {
            return gid_map.at(gid);
        }

        std::unordered_map<cell_gid_type, int> gid_map;
    };

    struct cell_identifier {
        cell_gid_type id;
        bool is_super_cell;
        bool is_not_so_super_cell;
    };

    using util::make_span;

    unsigned num_domains = ctx->distributed->size();
    unsigned domain_id = ctx->distributed->id();
    auto num_global_cells = rec.num_cells();

    //algoritme too devide the global cellsover the number of domains.
    auto dom_size = [&](unsigned dom) -> cell_gid_type {
        const cell_gid_type B = num_global_cells/num_domains;
        const cell_gid_type R = num_global_cells - num_domains*B;
        return B + (dom<R);
    };


    // Gid_part is all the nodes we will start searching from for this node !

    // Global load balance
    std::vector<cell_gid_type> gid_divisions;
    auto gid_part = make_partition(
        gid_divisions, transform_view(make_span(num_domains), dom_size));

    // Local load balance
    std::vector<std::vector<cell_gid_type>> not_so_super_cells; //cells connected by gj but allowed to be distributed
    std::vector<std::vector<cell_gid_type>> super_cells; //cells connected by gj
    std::vector<cell_gid_type> reg_cells; //independent cells

    // Map to track visited cells (cells that already belong to a group)
    std::unordered_set<cell_gid_type> visited;

    // Connected components algorithm using BFS
    std::queue<cell_gid_type> q;
    for (auto gid: make_span(gid_part[domain_id])) {
        if (!rec.gap_junctions_on(gid).empty()) {
            // If cell hasn't been visited yet, must belong to new super_cell
            // Perform BFS starting from that cell
            if (!visited.count(gid)) {
                visited.insert(gid);
                if(rec.get_cell_kind(gid) == cell_kind::cable) {
                    std::vector<cell_gid_type> cg;
                    q.push(gid);
                    while (!q.empty()) {
                        auto element = q.front();
                        q.pop();
                        cg.push_back(element);
                        // Adjacency list
                        auto conns = rec.gap_junctions_on(element);
                        for (auto c: conns) {
                            if (!visited.count(c.peer.gid)) {
                                visited.insert(c.peer.gid);
                                q.push(c.peer.gid);
                            }
                        }
                    }
                    super_cells.push_back(cg);
                    // if we allow distribution of the cable cell cable_distibuted needs to be selected
                }else if(rec.get_cell_kind(gid) == cell_kind::cable_distributed){
                    std::vector<cell_gid_type> cg;
                    q.push(gid);
                    while (!q.empty()) {
                        auto element = q.front();
                        q.pop();
                        cg.push_back(element);
                        // Adjacency list
                        auto conns = rec.gap_junctions_on(element);
                        for (const auto& c: conns) {
                            if (!visited.count(c.peer.gid)) {
                                visited.insert(c.peer.gid);
                                q.push(c.peer.gid);
                            }
                        }
                    }
                    not_so_super_cells.push_back(cg);
                }
            }
        }
        else {
            // If cell has no gap_junctions, put in separate group of independent cells
            reg_cells.push_back(gid);
        }
    }

    // Sort super_cell groups and only keep those where the first element in the group belongs to domain
    super_cells.erase(std::remove_if(super_cells.begin(), super_cells.end(),
            [gid_part, domain_id](std::vector<cell_gid_type>& cg)
            {
                std::sort(cg.begin(), cg.end());
                return cg.front() < gid_part[domain_id].first;
            }), super_cells.end());

    // Sort not_so_super_cell and keep everything on each and every node
    for(std::vector<cell_gid_type>&nssp : not_so_super_cells) {
        std::sort(nssp.begin(), nssp.end());
    }

    // Collect local gids that belong to this rank, and sort gids into kind lists
    // kind_lists maps a cell_kind to a vector of either:
    // 1. gids of regular cells (in reg_cells)
    // 2. indices of supercells (in super_cells)
    // 3. not_so_super cells only the part that belongs to this rank needs to be added to kind_lists.
    std::vector<cell_gid_type> local_gids;
    std::unordered_map<cell_kind, std::vector<cell_identifier>> kind_lists;
    for (auto gid: reg_cells) {
        local_gids.push_back(gid);
        kind_lists[rec.get_cell_kind(gid)].push_back({gid, false,false});
    }

    for (unsigned i = 0; i < super_cells.size(); i++) {
        auto kind = rec.get_cell_kind(super_cells[i].front());
        for (auto gid: super_cells[i]) {
            if (rec.get_cell_kind(gid) != kind) {
                throw gj_kind_mismatch(gid, super_cells[i].front());
            }
            local_gids.push_back(gid);
        }
        kind_lists[kind].push_back({i, true,false});
    }

    for (unsigned i = 0; i < not_so_super_cells.size(); i++) {
        auto kind = rec.get_cell_kind(not_so_super_cells[i].front());
        for (auto gid: not_so_super_cells[i]) {
            const auto temp = make_span(gid_part[domain_id]);
            if (temp.front() <= gid && gid <= temp.back()) {
                if (rec.get_cell_kind(gid) != kind) {
                    throw gj_kind_mismatch(gid, not_so_super_cells[i].front());
                }
                std::cout << domain_id << ": pushing gid " << gid << " front/back where: " << temp.front() << "/" << temp.back() << std::endl;
                local_gids.push_back(gid);
            }
        }
        kind_lists[kind].push_back({i, false,true});
    }

    // Create a flat vector of the cell kinds present on this node,
    // partitioned such that kinds for which GPU implementation are
    // listed before the others. This is a very primitive attempt at
    // scheduling; the cell_groups that run on the GPU will be executed
    // before other cell_groups, which is likely to be more efficient.
    //
    // TODO: This creates an dependency between the load balancer and
    // the threading internals. We need support for setting the priority
    // of cell group updates according to rules such as the back end on
    // which the cell group is running.

    auto has_gpu_backend = [&ctx](cell_kind c) {
        return cell_kind_supported(c, backend_kind::gpu, *ctx);
    };

    std::vector<cell_kind> kinds;
    for (auto l: kind_lists) {
        kinds.push_back(cell_kind(l.first));
    }
    std::partition(kinds.begin(), kinds.end(), has_gpu_backend);

    //print debug part
    {
#ifdef DEBUGMAX
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        int world_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        for (int worlds = 0; worlds < world_size; worlds++) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (worlds == world_rank) {
                std::cout << "rank: " << world_rank << std::endl;
                std::cout << "domain_id " << domain_id << " : ";
                for (auto gid: make_span(gid_part[domain_id])) {
                    std::cout << " " << gid;
                }
                std::cout << std::endl;

                std::cout << "cells " << " : ";
                for (auto &cell : kind_lists) {
                    std::cout << " " << cell.first;
                }
                std::cout << std::endl;

                size_t count = 0;
                std::cout << "supercells: " << std::endl;
                for (auto &scs : super_cells) {
                    std::cout << count << " : ";
                    for (auto &sc : scs) {
                        std::cout << " " << sc;
                    }
                    std::cout << std::endl;
                    count++;
                }

                count = 0;
                std::cout << "not_so_supercells: " << std::endl;
                for (auto &scs : not_so_super_cells) {
                    std::cout << count << " : ";
                    for (auto &sc : scs) {
                        std::cout << " " << sc;
                    }
                    std::cout << std::endl;
                    count++;
                }

                count = 0;
                std::cout << "local_gids: " << std::endl;
                for (auto &scs : local_gids) {
                    std::cout << scs << " ";
                }
                std::cout << std::endl;
            }
        }
#endif
    }

    //Here we create the groups...
    std::vector<group_description> groups;
    for (auto k: kinds) {
        partition_hint hint;
        if (auto opt_hint = util::value_by_key(hint_map, k)) {
            hint = opt_hint.value();
            if(!hint.cpu_group_size) {
                throw arbor_exception(arb::util::pprintf("unable to perform load balancing because {} has invalid suggested cpu_cell_group size of {}", k, hint.cpu_group_size));
            }
            if(hint.prefer_gpu && !hint.gpu_group_size) {
                throw arbor_exception(arb::util::pprintf("unable to perform load balancing because {} has invalid suggested gpu_cell_group size of {}", k, hint.gpu_group_size));
            }
        }

        backend_kind backend = backend_kind::multicore;
        std::size_t group_size = hint.cpu_group_size;

        if (hint.prefer_gpu && gpu_avail && has_gpu_backend(k)) {
            backend = backend_kind::gpu;
            group_size = hint.gpu_group_size;
        }

        std::vector<cell_gid_type> group_elements;


        // group_elements are sorted such that the gids of all members of a super_cell are consecutive.
            // in the case of distributed gap junctions we need to distribute a super_cell over all members involved.
            // needs to be done before we do this function.
        for (auto cell: kind_lists[k]) {
            if (cell.is_super_cell == true) {

                // als group_element size groter is dan hoe groot de supercell moet zijn voegen we hem toe aan de groeps.
                if (group_elements.size() + super_cells[cell.id].size() > group_size && !group_elements.empty()) {
                    groups.push_back({k, std::move(group_elements), backend});
                    group_elements.clear();
                }

                //toevoegen van de eerste super cell.
                for (auto gid: super_cells[cell.id]) {
                    group_elements.push_back(gid);
                }

            } else if (cell.is_not_so_super_cell == true){

                // als group_element size groter is dan hoe groot de supercell moet zijn voegen we hem toe aan de groeps.
                if (group_elements.size() + not_so_super_cells[cell.id].size() > group_size && !group_elements.empty()) {
                    groups.push_back({k, std::move(group_elements), backend});
                    group_elements.clear();
                }

                //toevoegen van de eerste super cell.
                for (auto gid: not_so_super_cells[cell.id]) {
                    group_elements.push_back(gid);
                }

            } else {
                group_elements.push_back(cell.id);
            }

            //fixen van de hints !
            if (group_elements.size()>=group_size) {
                groups.push_back({k, std::move(group_elements), backend});
                group_elements.clear();
            }
        }
        if (!group_elements.empty()) {
            groups.push_back({k, std::move(group_elements), backend});
        }
    }
    cell_size_type num_local_cells = local_gids.size();

    // Exchange gid list with all other nodes
    // global all-to-all to gather a local copy of the global gid list on each node.
    auto global_gids = ctx->distributed->gather_gids(local_gids);

// MAX TOEVOEGING ::
//    /// for distribution of cell groups exploration work. will be set directly for now
//    size_t my_domain = 0;
//    std::vector<int> domains;
//    std::vector<int> gid_local_offset
    for(auto &group: groups) {
        if(group.kind != arb::cell_kind::cable_distributed){
            //never needed but stil nice to fill it up
            group.my_domain = domain_id;
            group.domains.push_back(domain_id);
            group.gid_local_offsets.push_back(0);
        }else{
            group.my_domain = domain_id;
            for(size_t domain = 0; domain < num_domains ;domain++) {
                for(size_t gid_it = 0; gid_it < group.gids.size() ; gid_it++){
                    auto &i_gid = group.gids[gid_it];
                    if (i_gid > gid_part[domain].second) {
                        break; //break of the search we we are out of range for this domain.
                    } else if (i_gid >= gid_part[domain].first){
                        // first igid that belongs to this domain.
                        // add the domain and this i_gid is the offset,
                        // make that offset relative to the gid vector aka gid_it
                        // and move on to the next domain.
                        group.domains.push_back(domain);
                        group.gid_local_offsets.push_back(gid_it);
                        break;
                    }
                }
            }
        }
    }

    domain_decomposition d;
    d.num_domains = num_domains;
    d.domain_id = domain_id;
    d.num_local_cells = num_local_cells;
    d.num_global_cells = num_global_cells;
    d.groups = std::move(groups);
    d.gid_domain = partition_gid_domain(global_gids, num_domains);
    return d;
}

} // namespace arb

