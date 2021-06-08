/*
 * adjacency_graph.cpp
 * This file is part of bpm21
 *
 * Copyright (C) 2021 - Giacomo Bergami
 *
 * bpm21 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * bpm21 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with bpm21. If not, see <http://www.gnu.org/licenses/>.
 */


//
// Created by giacomo on 28/02/21.
//

#include "graphs/adjacency_graph.h"
#include <cassert>
using namespace jackbergus::fuzzyStringMatching3::graphs;

adjacency_graph::adjacency_graph() : V_size(0), E_size(0), casusu(ADJACENCY_GRAPH_CASE) {}

/*size_t adjacency_graph::add_node() {
    nodes.emplace_back();
    return V_size++;
}

size_t adjacency_graph::add_edge(size_t src, size_t dst) {
    edge_ids.emplace_back(src, dst);
    ingoing_edges[dst].emplace_back(E_size);
    return nodes.at(src).emplace_back(E_size++);
}*/

size_t jackbergus::fuzzyStringMatching3::graphs::add_node(adjacency_graph *graph, const std::string& label) {
    graph->nodes.emplace_back();
    if (graph->casusu == WEIGHTED_LABELLED_GRAPH_CASE) {
        auto G = ((weigthed_labelled_automata*)graph);
        G->node_label.emplace_back(label);
        G->label_conversion[label].emplace_back(G->count_labels);
        G->inv_label_conversion.insert(std::make_pair(graph->V_size, label));
        auto it = G->final_labels.emplace(label, G->count_labels);
        if (it.second) {
            G->count_labels++;
        }
    }
    return graph->V_size++;
}

size_t jackbergus::fuzzyStringMatching3::graphs::add_edge(adjacency_graph *graph, size_t src, size_t dst, double cost) {
    graph->edge_ids.emplace_back(src, dst);
    graph->ingoing_edges[dst].emplace_back(graph->E_size);
    if (graph->casusu == WEIGHTED_LABELLED_GRAPH_CASE) {
        if (graph->E_size == 0)
            ((weigthed_labelled_automata*)graph)->minimum_edge_weight = cost;
        else
            ((weigthed_labelled_automata*)graph)->minimum_edge_weight = std::min(
                    ((weigthed_labelled_automata*)graph)->minimum_edge_weight, cost);
        ((weigthed_labelled_automata*)graph)->edge_weight.emplace_back(cost);
    }
    return graph->nodes.at(src).emplace_back(graph->E_size++);
}

#include <algorithm>

bool jackbergus::fuzzyStringMatching3::graphs::update_edge_target(adjacency_graph *graph, size_t edge_id, size_t new_dst, double* cost_to_update) {
    if ((!graph) || (graph->E_size <= edge_id) || (graph->V_size <= new_dst))
        return false;
    auto& edge_src_dst = graph->edge_ids.at(edge_id);
    auto& ingoing_dst = graph->ingoing_edges[edge_src_dst.second];
    ingoing_dst.erase(std::remove(ingoing_dst.begin(), ingoing_dst.end(), edge_id), ingoing_dst.end());

    edge_src_dst.second = new_dst;
    graph->ingoing_edges[edge_src_dst.second].emplace_back(edge_id);
    if (cost_to_update && (graph->casusu == WEIGHTED_LABELLED_GRAPH_CASE)) {
        ((weigthed_labelled_automata*)graph)->minimum_edge_weight = std::min(
                    ((weigthed_labelled_automata*)graph)->minimum_edge_weight, *cost_to_update);
        ((weigthed_labelled_automata*)graph)->edge_weight[edge_id] = *cost_to_update;
    }
    return true;
}

bool jackbergus::fuzzyStringMatching3::graphs::remove_edge(adjacency_graph *graph, size_t edge_id) {
    if ((!graph) || (graph->E_size <= edge_id) /*|| (graph->V_size <= new_dst)*/)
        return false;
    auto& edge_src_dst = graph->edge_ids.at(edge_id);

    auto& ingoing_dst = graph->ingoing_edges[edge_src_dst.second];
    ingoing_dst.erase(std::remove(ingoing_dst.begin(), ingoing_dst.end(), edge_id), ingoing_dst.end());
    auto& outgoing_src = graph->nodes[edge_src_dst.first];
    outgoing_src.erase(std::remove(outgoing_src.begin(), outgoing_src.end(), edge_id), outgoing_src.end());
    graph->edge_ids[edge_id] = {-1, -1};

    return true;
}

std::pair<size_t, size_t> jackbergus::fuzzyStringMatching3::graphs::add_undirected_edge(adjacency_graph *graph, size_t src, size_t dst, double cost) {
    return {add_edge(graph, src, dst, cost), add_edge(graph, dst, src, cost)};
}

const std::pair<size_t, size_t> &jackbergus::fuzzyStringMatching3::graphs::edge_from_id(const adjacency_graph *graph, size_t edge_id) {
    return graph->edge_ids.at(edge_id);
}

const std::vector<size_t> &jackbergus::fuzzyStringMatching3::graphs::getOutgoingEdgesId(adjacency_graph *graph, size_t node_id) {
    return graph->nodes.at(node_id);
}

const std::vector<size_t> &jackbergus::fuzzyStringMatching3::graphs::getIngoingEdgesId(const adjacency_graph *graph, size_t node_id) {
    return graph->ingoing_edges.at(node_id);
}

void jackbergus::fuzzyStringMatching3::graphs::DFSUtil(const adjacency_graph *graph, size_t src, std::unordered_set<size_t> &visited) {
    visited.insert(src);

    for (size_t edge_id : graph->nodes.at(src)) {
        size_t dst = graph->edge_ids.at(edge_id).second;
        if (!visited.contains(dst))
            jackbergus::fuzzyStringMatching3::graphs::DFSUtil(graph, dst, visited);
    }
}

void jackbergus::fuzzyStringMatching3::graphs::printAllPathsUtil(const adjacency_graph *graph, size_t u, size_t d, std::unordered_set<size_t> &visited,
                       std::vector<ssize_t> &path, size_t path_index, std::unordered_set<size_t> &visited_src_dst,
                       std::unordered_set<size_t> &global) {
    // Mark the current node and store it in path[]
    global.insert(u);
    visited.insert(u);
    path[path_index] = u;
    path_index++;

    // If current vertex is same as destination, then print
    // current path[]
    if (u == d) {
        for (int i = 0; i < path_index; i++)
            visited_src_dst.insert(path[i]);
    } else {
        for (size_t edge_id : graph->nodes.at(u)) {
            size_t dst = graph->edge_ids.at(edge_id).second;
            if ((!visited.contains(dst)) && (!global.contains(dst)))
                jackbergus::fuzzyStringMatching3::graphs::printAllPathsUtil(graph, dst, d, visited, path, path_index, visited_src_dst, global);
        }
    }

    // Remove current vertex from path[] and mark it as unvisited
    path_index--;
    visited.erase(u);
}

weigthed_labelled_automata::weigthed_labelled_automata() : adjacency_graph{} {
    casusu = WEIGHTED_LABELLED_GRAPH_CASE;
}

double jackbergus::fuzzyStringMatching3::graphs::edge_cost(const adjacency_graph *graph, size_t edge_id) {
    assert(edge_id < graph->E_size);
    if (graph->casusu == ADJACENCY_GRAPH_CASE) {
        return 1.0;
    } else {
        return ((weigthed_labelled_automata*)graph)->edge_weight.at(edge_id);
    }
}

const std::string &jackbergus::fuzzyStringMatching3::graphs::node_label(const adjacency_graph *graph, size_t node_id) {
    assert(node_id < graph->V_size);
    if (graph->casusu == ADJACENCY_GRAPH_CASE) {
        return "";
    } else {
        return ((weigthed_labelled_automata*)graph)->node_label.at(node_id);
    }
}

void from_string(weigthed_labelled_automata &graph, const std::vector<std::string> &trace) {
    assert(!trace.empty());
    size_t n = trace.size();
    size_t curr = 0, prev = 0;
    for (size_t i = 0; i<n; i++) {
        prev = add_node(&graph, trace.at(i));
        std::swap(curr, prev);
        if (i == 0) continue;
        add_edge(&graph, prev, curr);
    }
}

void jackbergus::fuzzyStringMatching3::graphs::dot(adjacency_graph *graph, std::ostream &os) {
    os << "digraph finite_state_machine {" << std::endl;
    os << "    rankdir=LR;" << std::endl;
    os << "    size=\"8,5\"" << std::endl;
    for (int node_id = 0, N = graph->nodes.size(); node_id<N; node_id++) {
        std::string shape = "circle";
        os << "node [shape = circle, label=\"";
        if (graph->casusu == WEIGHTED_LABELLED_GRAPH_CASE) {
            os << ((weigthed_labelled_automata*)graph)->node_label.at(node_id) << " [";
        }
        os << std::to_string(node_id);
        if (graph->casusu == WEIGHTED_LABELLED_GRAPH_CASE) {
            os <<  "] ";
        }
        os << "\", fontsize=10] q" << node_id << ";" << std::endl;
    }
    os << "\n\n";
    for (int node_id = 0, N = graph->nodes.size(); node_id<N; node_id++) {
        std::vector<size_t> outgoing = graph->nodes.at(node_id);
        for (const size_t edge_id : outgoing) {
            size_t dst = graph->edge_ids.at(edge_id).second;
            os << "q" << node_id << " -> q" << dst;
            if (graph->casusu == WEIGHTED_LABELLED_GRAPH_CASE) {
                os << " [ label=\"" << ((weigthed_labelled_automata*)graph)->edge_weight.at(edge_id) << "\" ] ";
            }
            os <<";\n";
        }
    }
    os << "}";
}

void jackbergus::fuzzyStringMatching3::graphs::edge_compacting(weigthed_labelled_automata &graph) {
    for (size_t i = 0; i<graph.V_size; i++) {
        edge_compacting(graph, i);
    }
}

void jackbergus::fuzzyStringMatching3::graphs::edge_compacting(weigthed_labelled_automata &graph, size_t i) {
    std::unordered_map<size_t, std::vector<size_t>> target_node_to_edge_id;
    for (size_t edge_id : graph.nodes.at(i)) {
        auto& ref = graph.edge_ids.at(edge_id);
        if (graph.edge_ids.at(edge_id).second == ((size_t)-1)) continue;
        target_node_to_edge_id[graph.edge_ids.at(edge_id).second].emplace_back(edge_id);
    }

    for (auto& ref : target_node_to_edge_id) {
        if (ref.second.size() == 1) continue; // Do not need to compact if the node has only one outgoing edge
        size_t elected_edge_id = *ref.second.begin();
        double count = 0;
        if (target_node_to_edge_id.size() == 1)
            count = 1;
        else for (size_t edge_id : ref.second)
                count += graph.edge_weight.at(edge_id);
        graph.edge_weight[elected_edge_id] = count;
        for (size_t j = 1, N = ref.second.size(); j < N; j++) {
            assert(remove_edge(&graph, ref.second.at(j)));
        }
    }
}

std::unordered_map<size_t, std::vector<size_t>>::const_iterator jackbergus::fuzzyStringMatching3::graphs::hasIngoingEdges(const adjacency_graph *graph, size_t node_id) {
    return graph->ingoing_edges.find(node_id);
}

void jackbergus::fuzzyStringMatching3::graphs::edge_normalize(weigthed_labelled_automata &graph) {
    for (size_t u = 0; u < graph.V_size; u++) {
        auto& ref = graph.nodes.at(u);
        double total_cost = 0.0;
        for (size_t i = 0, N = ref.size(); i<N; i++) {
            total_cost += graph.edge_weight.at(ref.at(i));
        }
        for (size_t i = 0, N = ref.size(); i<N; i++) {
            graph.edge_weight[ref.at(i)] /= total_cost;
        }
    }
}

