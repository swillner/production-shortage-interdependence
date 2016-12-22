#include <MRIOTable.h>
#include <nvector.h>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <argp.h>
#include <deque>

typedef unsigned int I; // Index type
typedef float T; // Data type

const char* argp_program_bug_address = "<sven.willner@pik-potsdam.de>";
const char* argp_program_version = "1.0.0";
static char doc[] = "Production Shortage Interdependence measure (PSI).\n\
 Described in:\n\
    N. Glanemann, S.N. Willner, L. Wenz, R. Bierkandt, A. Levermann:\
 Abrupt Events and the Global Supply Network: A Network Measure for Cascading Production Losses.\
 TODO (2016). DOI: TODO\n\n\
Source: https://github.com/swillner/production-shortage-interdependence";
static char args_doc[] = "DATAFILE [INDEXFILE] OUTPUTFILE";
static struct argp_option options[] = {
    { "threshold", 't', "THRESHOLD", 0, "Only read values > THRESHOLD from table (default: 0)" },
    { "gamma", 'g', "GAMMA", 0, "Substitutability 'gamma' (default: 0)" },
    { "max", 'm', "MAX_PSI", 0, "Maximal PSI level to calculate (default: 10)" },
    { "per-region", 'r', "", 0, "Output psi from sector to region" },
    { 0 }
};

struct arguments {
    T threshold;
    T gamma;
    I max_psi;
    bool per_region;
    char* files[3];
};

static error_t parse_opt(int key, char* arg, struct argp_state* state) {
    struct arguments* args = static_cast<arguments*>(state->input);
    switch (key) {
        case 't':
            args->threshold = std::stof(arg);
            break;

        case 'g':
            args->gamma = std::stof(arg);
            break;

        case 'm':
            args->max_psi = std::stod(arg);
            break;

        case 'r':
            args->per_region = true;
            break;

        case ARGP_KEY_ARG:
            if (state->arg_num >= 3) {
                argp_usage(state);
            } else {
                args->files[state->arg_num] = arg;
            }
            break;

        case ARGP_KEY_END:
            if (state->arg_num < 2) {
                argp_usage(state);
            }
            break;

        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static I network_size;
static I sectors_size;
static I regions_size;
static mrio::Table<T, I> table;

static inline I get_index(const I i, const I r) {
    return i + r * sectors_size;
}

static inline I get_sector(const I index) {
    return index % sectors_size;
}

static inline I get_region(const I index) {
    return index / sectors_size;
}

inline const bool ends_with(const std::string& str, const std::string& end) {
    return !(end.size() > str.size()) || std::equal(end.rbegin(), end.rend(), str.rbegin());
}

int main(int argc, char* argv[]) {
#ifndef DEBUG
    try {
#endif
        struct arguments arguments;
        struct argp argp = { options, parse_opt, args_doc, doc };
        arguments.threshold = 1000;
        arguments.gamma = 0;
        arguments.max_psi = 10;
        arguments.files[2] = 0;
        argp_parse(&argp, argc, argv, 0, 0, &arguments);

        if (ends_with(arguments.files[0], ".nc")) {
            table.read_from_netcdf(arguments.files[0], arguments.threshold);
        } else if (ends_with(arguments.files[0], ".mrio")) {
            std::ifstream flows_file(arguments.files[0], std::ios::in | std::ios::binary);
            if (!flows_file) {
                std::cerr << "Could not open flows file '" << arguments.files[0] << "'";
                return -1;
            }
            table.read_from_mrio(flows_file, arguments.threshold);
            flows_file.close();
        } else if (ends_with(arguments.files[0], ".csv")) {
            if (!arguments.files[2]) {
                std::cerr << "Missing index file" << std::endl;
                return -1;
            }
            std::ifstream flows_file(arguments.files[0], std::ios::in);
            if (!flows_file) {
                std::cerr << "Could not open flows file '" << arguments.files[0] << "'";
                return -1;
            }
            std::ifstream index_file(arguments.files[1], std::ios::in);
            if (!index_file) {
                std::cerr << "Could not open index file '" << arguments.files[1] << "'";
                return -1;
            }
            table.read_from_csv(index_file, flows_file, arguments.threshold);
            index_file.close();
            flows_file.close();
        } else {
            std::cerr << "Unknown input file format" << std::endl;
            return -1;
        }

        network_size = table.index_set().size();
        sectors_size = table.index_set().total_sectors_count();
        regions_size = table.index_set().total_regions_count();

        T total = 0;
        nvector<T, 1> in_flow(0, network_size);
        nvector<T, 1> out_flow(0, network_size);
        nvector<T, 2> in_flow_by_sector(0, sectors_size, network_size);
        nvector<T, 1> region_output(0, regions_size);
        nvector<T, 3> psi(0, arguments.max_psi, network_size, network_size);
        nvector<bool, 3, std::deque<bool>> paths_traversed_ir_to_(false, network_size, sectors_size, network_size);
        nvector<bool, 3, std::deque<bool>> next_paths_traversed_ir_to_(false, network_size, sectors_size, network_size);

        for (I ir = 0; ir < network_size; ir++) {
            I i = get_sector(ir);
            I r = get_region(ir);
            for (I js = 0; js < network_size; js++) {
                // if (ir != js) { // TODO self-supply ok?
                in_flow(ir) += table(js, ir);
                out_flow(ir) += table(ir, js);
                in_flow_by_sector(i, js) += table(ir, js);
                //}
            }
            region_output(r) += out_flow(ir);
            total += out_flow(ir);
        }
        for (I ir = 0; ir < network_size; ir++) {
            psi(0, ir, ir) = 1;
            I i = get_sector(ir);
            for (I js = 0; js < network_size; js++) {
                if (ir == js) {
                    for (I level = 1; level < arguments.max_psi; level++) {
                        psi(level, ir, js) = 1;
                    }
                } else if (in_flow_by_sector(i, js) > 0) {
                    psi(1, ir, js) = table(ir, js) / in_flow_by_sector(i, js);
                } else if (table(ir, js) > 0) {
                    return 1; // TODO ?????
                }
            }
        }
        nvector<I, 1> last_max_k(0, network_size);

        for (I ir = 0; ir < network_size; ir++) {
            if (ir > 0) {
                I l;
#pragma omp parallel default(shared) private(l)
                {
#pragma omp for schedule(guided) nowait
                    for (l = 0; l < sectors_size; l++) {
                        for (I js = 0; js < network_size; js++) {
                            for (I mv = 0; mv < network_size; mv++) {
                                paths_traversed_ir_to_(js, l, mv) = false;
                                next_paths_traversed_ir_to_(js, l, mv) = false;
                            }
                        }
                    }
                }
            }
            for (I js = 0; js < network_size; js++) {
                last_max_k(js) = get_sector(js);
            }
            for (I level = 2; level < arguments.max_psi; level++) {
                I js;
#pragma omp parallel default(shared) private(js)
                {
#pragma omp for schedule(guided) nowait
                    for (js = 0; js < network_size; js++) {
                        if (js != ir) {
                            I max_k = last_max_k(js);
                            for (I k = 0; k < sectors_size; k++) {
                                T damage = 0;
                                for (I u = 0; u < regions_size; u++) {
                                    I ku = get_index(k, u);
                                    if (ku != js) {
                                        if (!paths_traversed_ir_to_(ku, k, js) || ku == ir) {
                                            damage += psi(level - 1, ir, ku) * psi(1, ku, js) * (1 - arguments.gamma);
                                        }
                                    }
                                }
                                if (psi(level, ir, js) < damage) {
                                    psi(level, ir, js) = damage;
                                    max_k = k;
                                }
                            }
                            if (level < arguments.max_psi - 1) {
                                next_paths_traversed_ir_to_(js, max_k, js) = true;
                                for (I u = 0; u < regions_size; u++) {
                                    I ku = get_index(max_k, u);
                                    if (psi(level - 1, ir, ku) > 0 && psi(1, ku, js) > 0 && !paths_traversed_ir_to_(ku, max_k, js)) {
                                        for (I l = 0; l < sectors_size; l++) {
                                            for (I mv = 0; mv < network_size; mv++) {
                                                next_paths_traversed_ir_to_(js, l, mv) =
                                                        next_paths_traversed_ir_to_(js, l, mv) || paths_traversed_ir_to_(ku, l, mv);
                                            }
                                        }
                                    }
                                }
                                last_max_k(js) = max_k;
                            }
                        }
                    }
                }
                if (level < arguments.max_psi - 1) {
                    I l2;
#pragma omp parallel default(shared) private(l2)
                    {
#pragma omp for schedule(guided) nowait
                        for (l2 = 0; l2 < sectors_size; l2++) {
                            for (I ir2 = 0; ir2 < network_size; ir2++) {
                                for (I l3 = 0; l3 < network_size; l3++) {
                                    paths_traversed_ir_to_(ir2, l2, l3) = next_paths_traversed_ir_to_(ir2, l2, l3);
                                    next_paths_traversed_ir_to_(ir2, l2, l3) = false;
                                }
                        }
                        }
                    }
                }
            }
        }

        nvector<T, 2> psi_from(0, arguments.max_psi, network_size);
        nvector<T, 2> psi_to(0, arguments.max_psi, network_size);
        nvector<T, 2> psi_from_unweighted(0, arguments.max_psi, network_size);
        for (I level = 0; level < arguments.max_psi; level++) {
            I ir;
#pragma omp parallel default(shared) private(ir)
            {
#pragma omp for schedule(guided) nowait
                for (ir = 0; ir < network_size; ir++) {
                    for (I js = 0; js < network_size; js++) {
                        psi_from(level, ir) += out_flow(js) * psi(level, ir, js);
                        psi_to(level, ir) += out_flow(ir) * psi(level, js, ir);
                        psi_from_unweighted(level, ir) += psi(level, ir, js);
                    }
                    psi_from(level, ir) /= total;
                    psi_to(level, ir) /= total;
                    psi_from_unweighted(level, ir) /= network_size;
                }
            }
        }

        nvector<T, 3> psi_from_sector_to_country(0, arguments.max_psi, network_size, network_size);
        for (I level = 0; level < arguments.max_psi; level++) {
            I ir;
#pragma omp parallel default(shared) private(ir)
            {
#pragma omp for schedule(guided) nowait
                for (ir = 0; ir < network_size; ir++) {
                    for (I s2 = 0; s2 < regions_size; s2++) {
                        for (I j2 = 0; j2 < sectors_size; j2++) {
                            I j2s2 = get_index(j2, s2);
                            psi_from_sector_to_country(level, ir, s2) += out_flow(j2s2) * psi(level, ir, j2s2);
                        }
                        psi_from_sector_to_country(level, ir, s2) /= region_output(s2);
                    }
                }
            }
        }

        char* outfile;
        if (arguments.files[2]) {
            outfile = arguments.files[2];
        } else {
            outfile = arguments.files[1];
        }
        std::ofstream csv_file(outfile);
        csv_file << std::setprecision(5) << std::fixed;
        if (arguments.per_region) {
            for (I ir = 0; ir < network_size; ir++) {
                for (I s = 0; s < regions_size; s++) {
                    csv_file << table.index_set().supersectors()[get_sector(ir)]->name << ",";
                    csv_file << table.index_set().superregions()[get_region(ir)]->name << ",";
                    csv_file << table.index_set().superregions()[s]->name << ",";
                    for (int level = 0; level < arguments.max_psi; level++) {
                        csv_file << psi_from_sector_to_country(level, ir, s);
                        if (level < arguments.max_psi - 1) {
                            csv_file << ",";
                        } else {
                            csv_file << "\n";
                        }
                    }
                }
            }
        } else {
            for (I ir = 0; ir < network_size; ir++) {
                csv_file << table.index_set().supersectors()[get_sector(ir)]->name << ",";
                csv_file << table.index_set().superregions()[get_region(ir)]->name << ",";
                for (int level = 0; level < arguments.max_psi; level++) {
                    csv_file << psi_from(level, ir);
                    if (level < arguments.max_psi - 1) {
                        csv_file << ",";
                    } else {
                        csv_file << "\n";
                    }
                }
            }
        }
        csv_file.close();

        return 0;

#ifndef DEBUG
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        return -1;
    }
#endif
}
