#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <nlohmann/json.hpp>
#include <arbor/cable_cell.hpp>
#include <arbor/common_types.hpp>
#include <arbor/context.hpp>
#include <arbor/load_balance.hpp>
#include <arbor/profile/meter_manager.hpp>
#include <arbor/simple_sampler.hpp>
#include <arbor/simulation.hpp>

#include <arbor/recipe.hpp>
#include <utility>
#include <arborenv/concurrency.hpp>
#include <arborio/label_parse.hpp>
#include <sup/ioutil.hpp>

//#include <arborenv/gpu_env.hpp>
//#include "mechanisms/Local_catalogue.cpp"
//#include <arbor/assert_macro.hpp>
//#include <arbor/cable_cell.hpp>
//#include <arbor/profile/profiler.hpp>
//#include <arbor/version.hpp>
//#include <arborenv/gpu_env.hpp>
//#include <sup/json_meter.hpp>
//#include <sup/json_params.hpp>


#ifdef ARB_MPI_ENABLED
#include <mpi.h>
#include <arborenv/with_mpi.hpp>
#endif

struct io_params {
    std::string name = "default";
    double simtime   = 10;       //ms
    double timestep  = 0.025;  //ms
    unsigned n_cells = 1;
};

arb::cable_cell IO_cell();
void write_trace_json(const std::vector<arb::trace_vector<double>>& traces, unsigned rank, unsigned cell);
class IO_recipe: public arb::recipe {
public:
    explicit IO_recipe(io_params params,const arb::mechanism_catalogue* catalogue): params_ (std::move(params)),catalogue_(catalogue) {}

    [[nodiscard]] arb::cell_size_type num_cells() const override {
        return params_.n_cells;
    }

    [[nodiscard]] arb::util::unique_any get_cell_description(arb::cell_gid_type gid) const override {
        return IO_cell();
    }

    [[nodiscard]] arb::cell_kind get_cell_kind(arb::cell_gid_type gid) const override {
	    return arb::cell_kind::cable;
    }

    [[nodiscard]] std::any get_global_properties(arb::cell_kind k) const override {
        arb::cable_cell_global_properties a;
        a.default_parameters = arb::neuron_parameter_defaults;
        a.default_parameters.init_membrane_potential = -65;
        a.default_parameters.membrane_capacitance = 0.01;
        a.default_parameters.axial_resistivity = 100;
        a.catalogue = catalogue_;
        return a;
    }

    [[nodiscard]] std::vector<arb::probe_info> get_probes(arb::cell_gid_type gid) const override {
        // Measure membrane voltage at end of soma.
        arb::mlocation loc{4, 1.};
        return {arb::cable_probe_membrane_voltage{loc}};
    }

    [[nodiscard]] std::vector<arb::cell_connection> connections_on(arb::cell_gid_type gid) const override {
        return {};
//        return {arb::cell_connection({gid - 1, 0}, 0, params_.event_weight, params_.event_min_delay)};
    }

//    [[nodiscard]] std::vector<arb::gap_junction_connection> gap_junctions_on(arb::cell_gid_type gid) const override{
//        using policy = arb::lid_selection_policy;
//        std::vector<arb::gap_junction_connection> conns;
//
//        if(gid !=0){
//            conns.push_back(arb::gap_junction_connection({(arb::cell_gid_type)gid -1, "local_0", policy::assert_univalent},
//                                                         {"local_0", policy::assert_univalent},
//                                                         0.015));
//
//            conns.push_back(arb::gap_junction_connection({(arb::cell_gid_type)gid -1, "local_1", policy::assert_univalent},
//                                                         {"local_2", policy::assert_univalent},
//                                                         0.015));
//
//            conns.push_back(arb::gap_junction_connection({(arb::cell_gid_type)gid -1, "local_1", policy::assert_univalent},
//                                                         {"local_2", policy::assert_univalent},
//                                                         0.015));
//
//        }
//        else{
//            conns.push_back(arb::gap_junction_connection({(arb::cell_gid_type)params_.n_cells-1, "local_0", policy::assert_univalent},
//                                                         {"local_0", policy::assert_univalent},
//                                                         0.015));
//
//            conns.push_back(arb::gap_junction_connection({(arb::cell_gid_type)params_.n_cells-1, "local_1", policy::assert_univalent},
//                                                         {"local_2", policy::assert_univalent},
//                                                         0.015));
//
//            conns.push_back(arb::gap_junction_connection({(arb::cell_gid_type)params_.n_cells-1, "local_1", policy::assert_univalent},
//                                                         {"local_2", policy::assert_univalent},
//                                                         0.015));
//
//        }
//
//        if(gid != params_.n_cells-1){
//            conns.push_back(arb::gap_junction_connection({(arb::cell_gid_type)gid +1, "local_0", policy::assert_univalent},
//                                                         {"local_0", policy::assert_univalent},
//                                                         0.015));
//
//            conns.push_back(arb::gap_junction_connection({(arb::cell_gid_type)gid +1, "local_1", policy::assert_univalent},
//                                                         {"local_2", policy::assert_univalent},
//                                                         0.015));
//
//            conns.push_back(arb::gap_junction_connection({(arb::cell_gid_type)gid +1, "local_1", policy::assert_univalent},
//                                                         {"local_2", policy::assert_univalent},
//                                                         0.015));
//        }
//        else{
//            conns.push_back(arb::gap_junction_connection({(arb::cell_gid_type)0, "local_0", policy::assert_univalent},
//                                                         {"local_0", policy::assert_univalent},
//                                                         0.015));
//
//            conns.push_back(arb::gap_junction_connection({(arb::cell_gid_type)0, "local_1", policy::assert_univalent},
//                                                         {"local_2", policy::assert_univalent},
//                                                         0.015));
//
//            conns.push_back(arb::gap_junction_connection({(arb::cell_gid_type)0, "local_1", policy::assert_univalent},
//                                                         {"local_2", policy::assert_univalent},
//                                                         0.015));
//        }
//        return conns;
//    }

//    [[nodiscard]] arb::mechanism_desc gap_junction_mech() const override {
//        return arb::mechanism_desc("glomerulus_gj");
//    }

private:
    io_params params_;
    const arb::mechanism_catalogue* catalogue_;

};

/*
 * MAIN
 ************************************/
int main() {

    arb::proc_allocation resources;
    resources.num_threads = 1;

#ifdef ARB_MPI_ENABLED
    arbenv::with_mpi guard(false);
    auto context = arb::make_context(resources, MPI_COMM_WORLD);
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    }
#else
//     resources.gpu_id = arbenv::default_gpu();
     auto context = arb::make_context(resources);
#endif

#ifdef ARB_PROFILE_ENABLED
    arb::profile::profiler_initialize(context);
#endif

    // Print a banner with information about hardware configuration
    std::cout << "gpu:      " << (has_gpu(context)? "yes": "no") << "\n";
    std::cout << "threads:  " << num_threads(context)            << "\n";
    std::cout << "mpi:      " << (has_mpi(context)? "yes": "no") << "\n";
    std::cout << "ranks:    " << num_ranks(context) << "\n"      << std::endl;

    arb::profile::meter_manager meter;
    meter.start(context);

    io_params params;

    auto cat = (arb::mechanism_catalogue*)&arb::global_smol_catalogue();
    auto IOUCat = (arb::mechanism_catalogue*)&arb::global_IOU_catalogue();
    cat->import(*IOUCat,"");      //actually add the catalogue we are interested in.

    // Create an instance of our recipe.
    IO_recipe recipe(params, cat);
    meter.checkpoint("Create Recipe", context);

    // Greate an decomposition
    auto decomp = arb::partition_load_balance(recipe, context);
    meter.checkpoint("Create Decomp", context);

    // Construct the model.
    arb::simulation sim(recipe, decomp, context);
    meter.checkpoint("Create simmodel", context);

    // Set up the probe that will measure voltage in the cell.
    auto sched = arb::regular_schedule(1);
    std::vector<arb::trace_vector<double>> voltage_traces(decomp.num_local_cells); // This is where the voltage samples will be stored as (time, value) pairs
    // Now attach the sampler to the soma.
    size_t j=0;
    for (const auto& g : decomp.groups) {
        for (auto i : g.gids) {
            sim.add_sampler(arb::one_probe({i, 0}), sched, arb::make_simple_sampler(voltage_traces[j++]));
        }
    }
    meter.checkpoint("set up probes", context);

    // Run the simulation
    sim.run(params.simtime,params.timestep);

    meter.checkpoint("run simulation", context);

    //write the json values to a file
    write_trace_json(voltage_traces, arb::rank(context),0);

    meter.checkpoint("finished", context);
    //create report
    auto report = arb::profile::make_meter_report(meter, context);
    std::cout << report;

//    std::cout << arb::profile::profiler_summary() << "\n";
    return 0;
}


void write_trace_json(const std::vector<arb::trace_vector<double>>& traces, unsigned rank, unsigned cell) {
    std::string path = "./voltages_" + std::to_string(rank) + "_" + std::to_string(cell) + ".json";
    nlohmann::json json;
    json["name"] = "IOmax: cell " + std::to_string(cell);
    json["units"] = "mV";
    json["cell"] = std::to_string(cell);
    json["group"] = std::to_string(rank);
    json["probe"] = "0";
    auto &jt = json["data"]["time"];
    auto &jy = json["data"]["voltage"];

    for (const auto &sample: traces[cell].at(0)) {
        jt.push_back(sample.t);
        jy.push_back(sample.v);
        std::cout <<std::setprecision(64)<< sample.t <<"\t" << sample.v <<"\n";
    }

    std::ofstream file(path);
    file << std::setw(1)  << json << "\n";
}
arb::cable_cell IO_cell(){

    using namespace arborio::literals;

    // Create the sample tree that defines the morphology.
    arb::segment_tree tree;
    double soma_rad = 22.360679775/2.0; // convert diameter to radius in μm
    double dend_rad = 3./2;             // μm
    double axon_rad = 3./2;             // μm

    arb::label_dict labels;
    labels.set("soma", arb::reg::tagged(1));
    labels.set("dend", arb::reg::tagged(2));
    labels.set("axon", arb::reg::tagged(3));

    //append the soma
    tree.append(arb::mnpos, {0,0,0,soma_rad}, {2*soma_rad,0,0,soma_rad}, 1); // soma

    //append the dendriet
    tree.append(0,  {2*soma_rad + 50 ,0  ,0, dend_rad}, 3);  // 1
    tree.append(1,  {2*soma_rad + 100,-10,0, dend_rad}, 3);  // 3
    tree.append(2,  {2*soma_rad + 150,10 ,0, dend_rad}, 3);  // 4
    tree.append(2,  {2*soma_rad + 150,-30,0, dend_rad}, 3);  // 5
    tree.append(3,  {2*soma_rad + 250,20 ,0, dend_rad}, 3);  // 6
    tree.append(4,  {2*soma_rad + 250,0  ,0, dend_rad}, 3);  // 2
    tree.append(4,  {2*soma_rad + 350,60 ,0, dend_rad}, 3);  // 7

    //append the axon
    tree.append(arb::mnpos, {0,0,0, axon_rad}, {-10,0,0, axon_rad}, 2);  // axon

    arb::decor decor;

    //soma mechs
    arb::mechanism_desc na_s("na_s");
    decor.paint("soma"_lab, na_s);

    arb::mechanism_desc kdr("kdr");
    decor.paint("soma"_lab, kdr);

    arb::mechanism_desc ks("k");
    decor.paint("soma"_lab, ks);

    arb::mechanism_desc cal("cal");
    decor.paint("soma"_lab, cal);

    //dend mechs
//    arb::mechanism_desc smol_dend("smol_dend");
//    smol_dend["cah_gmax"]  = 0.010 * 1.7 / 2;
//    smol_dend["kca_gmax"]  = 0.200 * 0.7 * 1.5;
//    smol_dend["h_gmax"]    = 0.025 * 1.7;
//    smol_dend["cacc_gmax"] = 0.007;
//    decor.paint("dend"_lab, smol_dend);

    arb::mechanism_desc cah("cah");
    cah["gmax"] =0.010 * 1.7 / 2;
    decor.paint("dend"_lab, cah);

    arb::mechanism_desc kca("kca");
    kca["gmax"] = 0.200 * 0.7 * 1.5;
    decor.paint("dend"_lab, kca);

    arb::mechanism_desc h("h");
    h["gmax"] = 0.025 * 1.7;
    decor.paint("dend"_lab, h);

    arb::mechanism_desc cacc("cacc");
    cacc["gmax"] = 0.007;
    decor.paint("dend"_lab, cacc);

    //axon mechs
    arb::mechanism_desc na_a("na_a");
    decor.paint("axon"_lab, na_a);

    arb::mechanism_desc k("k");
    decor.paint("axon"_lab, k);

    //global mechs
//    arb::mechanism_desc leak("leak");
//    leak["gmax"] = 1.3e-05;
//    decor.paint("(all)"_reg, leak);

//    arb::mechanism_desc ou_noise("ou_noise/seed=10");
//    decor.paint("(all)"_reg, ou_noise);

    // Add a spike detector to the soma at the beginning
    //    decor.place(arb::mlocation{0,0}, arb::threshold_detector{10});

    // Add gap junction sites at the end of the branch
    decor.place(arb::mlocation{1, 1}, arb::gap_junction_site{},"local_0");
    decor.place(arb::mlocation{3, 1}, arb::gap_junction_site{},"local_1");
    decor.place(arb::mlocation{4, 1}, arb::gap_junction_site{},"local_2");

    // Attach a stimulus
    //    auto stim = arb::i_clamp::box(0, stim_duration, 0.4);
    //    decor.place(arb::mlocation{0, 0.5}, stim);

    // Add a synapse
    //    decor.place(arb::mlocation{0, 0.5}, "expsyn");

    auto out = arb::cable_cell(tree, labels, decor);
    return out;
}
