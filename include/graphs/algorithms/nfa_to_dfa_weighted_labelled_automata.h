/*
 * nfa_to_dfa_weighted_labelled_automata.h
 * This file is part of ProbabilisticTraceAlignment
 *
 * Copyright (C) 2020 - Giacomo Bergami
 *
 * ProbabilisticTraceAlignment is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * ProbabilisticTraceAlignment is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ProbabilisticTraceAlignment. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef COMPLETE_MATCHING_NFA_TO_DFA_WEIGHTED_LABELLED_AUTOMATA_H
#define COMPLETE_MATCHING_NFA_TO_DFA_WEIGHTED_LABELLED_AUTOMATA_H

#include <set>
#include <string>
#include <vector>
#include <hashing/hash_combine.h>
#include <graphs/adjacency_graph.h>
#include <basics/trace.h>
#include <ostream>
#include <hashing/pair_hash.h>
#include <hashing/vector_hash.h>

namespace jackbergus {
    namespace fuzzyStringMatching3 {
        namespace graphs {
            namespace algorithms {


/**
 * Defines the global node properties for merging nodes before minimizing the automaton.
 * I.e., a node is defined by its node id, the associated label, and the outgoing edges' labels (from the
 * target node) as well as their associated cost. This will be used to identify when we reached a similar node
 */
                struct determinization_node_eq_safe_heuristic {
                    std::set<size_t> node_id;
                    std::string node_label;
#ifdef DO_STRINGS
                    std::vector<std::pair<std::string, std::string>>
#else
                    std::vector<std::pair<std::string, double>>
#endif
                            outgoing_edges;

                    determinization_node_eq_safe_heuristic(const std::string &nodeLabel);

                    void insert_edge(const std::string &s, double w);

                    void finalize();

                    bool operator==(const determinization_node_eq_safe_heuristic &rhs) const;

                    bool operator!=(const determinization_node_eq_safe_heuristic &rhs) const;

                private:
                    bool sorted = false;

                public:
                    friend std::ostream &operator<<(std::ostream &os, const determinization_node_eq_safe_heuristic &properties) {
                        std::vector<size_t> S{properties.node_id.begin(), properties.node_id.end()};
                        os << "node_id: " << S << " node_label: " << properties.node_label << " outgoing_edges: "
                           << properties.outgoing_edges << " sorted: " << properties.sorted;
                        return os;
                    }
                };


            }
        }
    }
}

namespace std {
    /**
     * Hashing function associated to determinization_node_eq_safe_heuristic
     */
    template <>
    struct hash<jackbergus::fuzzyStringMatching3::graphs::algorithms::determinization_node_eq_safe_heuristic> {
        std::size_t operator()(const jackbergus::fuzzyStringMatching3::graphs::algorithms::determinization_node_eq_safe_heuristic& k) const {
            std::hash<std::string> sh;
#ifndef DO_STRINGS
            std::hash<double> dh;
#endif
            size_t seed = hash_combine(hash_combine(31, k.node_id), k.node_label);

            size_t h = 17;

            for (const auto& cp : k.outgoing_edges) h = hash_combine(h, sh(cp.first)
                                                                        #ifdef DO_STRINGS
                                                                        ^ sh(cp.second));
                                                                        #else
                                                                        ^ dh(cp.second));
#endif

            size_t hash =hash_combine(seed, h);
            return hash;
        }
    };
}


namespace jackbergus{
    namespace fuzzyStringMatching3 {
        namespace graphs {
            namespace algorithms {

/**
 * State information required by the nfa_to_dfa_weighted_labelled_automata algorithm, so to avoid method invocation
 * within a struct/class
 */
                struct determinization_information {
                    ssize_t has_result_well_state = -1;
                    ssize_t is_final_state_inserted_into_result = -1;
                    double  halting_condition = 0.0;
                    size_t  tg_initial, tg_final;

                    determinization_information(double haltingCondition, size_t  tg_initial, size_t  tg_final);

                    determinization_information(const determinization_information& ) = default;
                    determinization_information(determinization_information&& ) = default;
                    determinization_information& operator=(const determinization_information& ) = default;
                    determinization_information& operator=(determinization_information&& ) = default;
                };

                ssize_t nfa_to_dfa_weighted_labelled_automata(weigthed_labelled_automata& graph,
                                                              const probabilisitc_model_trace& MT,
                                                              weigthed_labelled_automata& out,
                                                              determinization_information& info,
                                                              std::unordered_map<determinization_node_eq_safe_heuristic, size_t>& node_compact_info);

            }
        }
    }
}


#endif //COMPLETE_MATCHING_NFA_TO_DFA_WEIGHTED_LABELLED_AUTOMATA_H
