/*
 * nfa_to_dfa_weighted_labelled_automata.cpp
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

#include "graphs/algorithms/nfa_to_dfa_weighted_labelled_automata.h"

jackbergus::fuzzyStringMatching3::graphs::algorithms::determinization_node_eq_safe_heuristic::determinization_node_eq_safe_heuristic(
        const ::std::string &nodeLabel) : node_label(nodeLabel)  {
    sorted = false;
}

void jackbergus::fuzzyStringMatching3::graphs::algorithms::determinization_node_eq_safe_heuristic::insert_edge(const ::std::string &s,
                                                                                                               double w) {
#ifdef DO_STRINGS
    outgoing_edges.emplace_back(s, std::to_string(w));
#else
    outgoing_edges.emplace_back(s, (w));
#endif
    sorted = false;
}

#include <algorithm>

void jackbergus::fuzzyStringMatching3::graphs::algorithms::determinization_node_eq_safe_heuristic::finalize() {
    if (!sorted) {
        ::std::sort(outgoing_edges.begin(), outgoing_edges.end());
        sorted = true;
    }
}

#include <cassert>

bool jackbergus::fuzzyStringMatching3::graphs::algorithms::determinization_node_eq_safe_heuristic::operator==(
        const jackbergus::fuzzyStringMatching3::graphs::algorithms::determinization_node_eq_safe_heuristic &rhs) const {
    assert(sorted);
    assert(rhs.sorted);
    return node_id == rhs.node_id &&
           node_label == rhs.node_label &&
           outgoing_edges == rhs.outgoing_edges;
}

bool jackbergus::fuzzyStringMatching3::graphs::algorithms::determinization_node_eq_safe_heuristic::operator!=(
        const jackbergus::fuzzyStringMatching3::graphs::algorithms::determinization_node_eq_safe_heuristic &rhs) const {
    return !(rhs == *this);
}




jackbergus::fuzzyStringMatching3::graphs::algorithms::determinization_information::determinization_information(
        double haltingCondition, size_t tg_initial, size_t tg_final) : tg_initial(tg_initial), tg_final(tg_final),  halting_condition(haltingCondition), has_result_well_state{-1}, is_final_state_inserted_into_result{-1} {}


ssize_t jackbergus::fuzzyStringMatching3::graphs::algorithms::nfa_to_dfa_weighted_labelled_automata(
        jackbergus::fuzzyStringMatching3::graphs::weigthed_labelled_automata &graph,
        const jackbergus::fuzzyStringMatching3::probabilisitc_model_trace &MT,
        jackbergus::fuzzyStringMatching3::graphs::weigthed_labelled_automata &out,
        jackbergus::fuzzyStringMatching3::graphs::algorithms::determinization_information &info,
        std::unordered_map<determinization_node_eq_safe_heuristic, size_t> &node_compact_info) {

    constexpr double limit = std::numeric_limits<double>::epsilon();
    if ((std::abs(MT.t_with_prob.probability - info.halting_condition) < limit) ||
        (MT.t_with_prob.probability < limit)) return -1;

    //std::unordered_set<size_t> curr;

    std::string curr_node_label = *MT.t_with_prob.t.rbegin();
    determinization_node_eq_safe_heuristic gnp{curr_node_label};
    std::unordered_map<size_t, std::vector<size_t>> node_id_to_sequence_id;
    for (size_t i = 0, N = MT.underlying_sequences.size(); i<N; i++) {
        const auto& eta = MT.underlying_sequences.at(i);
        size_t M = *eta.sequence.rbegin();
        //curr.insert(M);
        node_id_to_sequence_id[M].emplace_back(i);
        gnp.node_id.insert(M);
    }

    std::unordered_map<std::string, probabilisitc_model_trace> MM;
    for (const auto& node_id_to_sequenceId : node_id_to_sequence_id) {
        size_t node_id = node_id_to_sequenceId.first;

        std::unordered_map<std::string, std::vector<std::pair<double, size_t>>> M;
        for (const size_t edge_id : getOutgoingEdgesId((adjacency_graph *) &graph, node_id)) {
            double edge_weight = graph.edge_weight.at(edge_id);
            size_t dst = graph.edge_ids.at(edge_id).second;
            M[node_label(&graph, dst)].emplace_back(edge_weight, dst);
        }

        for (auto& cp : M) {
            probabilisitc_model_trace& trace = MM.emplace(cp.first, 0.0).first->second;
            trace.t_with_prob.t = MT.t_with_prob.t;
            trace.t_with_prob.t.emplace_back(cp.first);
            for (const size_t seqid : node_id_to_sequenceId.second) {
                trace.t_with_prob.probability += MT.underlying_sequences.at(seqid).extend_with_subsequent_step(cp.second, trace.underlying_sequences);
            }
            assert(trace.test());
        }
    }

    // Now, determining the outgoing edges arcs, that are going to be used for
    // determining the node
    auto it = MM.begin();
    while (it != MM.end()) {
        if ((MT.t_with_prob.probability > 0) && ( (it->second.t_with_prob.probability/MT.t_with_prob.probability) > info.halting_condition )) {
            assert(MT.t_with_prob.probability > 0);
            double transition_cost = it->second.t_with_prob.probability / MT.t_with_prob.probability;
            gnp.insert_edge(it->first, transition_cost);
            ///std::cout <<  " -[" << transition_cost << "]-> " << it->first << " : " << it->second.t_with_prob << std::endl << std::endl << std::endl;
            it++;
        } else {
            it = MM.erase(it);
        }
    }
    gnp.finalize();

    // Checking if we arrived in the final state: please remember that, in our model, the final state has no
    // outgoing edges, and so it is not possible to have a state containing both the final state and another state
    bool is_final_state = false;
    if ((gnp.node_id.size() == 1) && (*gnp.node_id.begin() == info.tg_final)) {
        assert(MM.empty());
        is_final_state = true;
    }

    ssize_t src = -1;
    if (is_final_state) {
        // As in this model the final state is just one, when I reach it, then remember that
        // and do not create a novel node, but just use the one that you already have
        if (info.is_final_state_inserted_into_result == -1) {
            info.is_final_state_inserted_into_result = add_node(&out, curr_node_label);
        }
        src = info.is_final_state_inserted_into_result;
    } else {
        ///std::cout << gnp << std::endl;
        // Checking whether we already passed through a node that has the same outgoing edges with the same weight,
        // and whether the node from where we start is indeed the same node. In this case, do not replicate the state
        auto it2 = node_compact_info.emplace(gnp, 0);
        if (!it2.second) {
            src = it2.first->second;
        } else {
            src = it2.first->second = add_node(&out, curr_node_label);
            // Performing the recursive visit of the node towards the adjacent ones
            it = MM.begin();
            while (it != MM.end()) {
                ssize_t dst = nfa_to_dfa_weighted_labelled_automata(graph, it->second, out, info, node_compact_info);
                // See handnotes for the explanation. Roughly speaking,
                // 1. MT.t_with_prob.probability, determines the probability that we accumulated to reach the node,
                //                                after merging the edges pointing towards nodes with the same label
                // 2. it->second.t_with_prob.probability, determines the overall probability that is required to generate
                //                                        the trace by adding a novel character
                // => As (1.) * transition_cost = (2.), then the cost required to perform the actual transition is (2.)/(1.)
                double transition_cost = it->second.t_with_prob.probability / MT.t_with_prob.probability;
                if (dst != -1) {
                    add_edge(&out, src, dst, transition_cost);
                } else {
                    if (info.has_result_well_state == -1) {
                        info.has_result_well_state = add_node(&out, ".");
                        add_edge(&out, info.has_result_well_state, info.has_result_well_state, 1.0);
                    }
                    add_edge(&out, src, info.has_result_well_state, transition_cost);
                }
                it++;
            }
        }
    }



    return src;
}
