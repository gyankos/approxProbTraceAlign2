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

size_t add_node(adjacency_graph *graph, const std::string& label) {
    graph->nodes.emplace_back();
    if (graph->casusu == WEIGHTED_LABELLED_GRAPH_CASE) {
        auto G = ((weigthed_labelled_automata*)graph);
        G->node_label.emplace_back(label);
        G->inv_label_conversion.insert(std::make_pair(graph->V_size, label)).second;
    }
    return graph->V_size++;
}

size_t add_edge(adjacency_graph *graph, size_t src, size_t dst, double cost) {
    graph->edge_ids.emplace_back(src, dst);
    graph->ingoing_edges[dst].emplace_back(graph->E_size);
    if (graph->casusu == WEIGHTED_LABELLED_GRAPH_CASE) {
        ((weigthed_labelled_automata*)graph)->edge_weight.emplace_back(cost);
    }
    return graph->nodes.at(src).emplace_back(graph->E_size++);
}

std::pair<size_t, size_t> add_undirected_edge(adjacency_graph *graph, size_t src, size_t dst, double cost) {
    return {add_edge(graph, src, dst, cost), add_edge(graph, dst, src, cost)};
}

const std::pair<size_t, size_t> &edge_from_id(const adjacency_graph *graph, size_t edge_id) {
    return graph->edge_ids.at(edge_id);
}

const std::vector<size_t> &getOutgoingEdgesId(adjacency_graph *graph, size_t node_id) {
    return graph->nodes.at(node_id);
}

const std::vector<size_t> &getIngoingEdgesId(const adjacency_graph *graph, size_t node_id) {
    return graph->ingoing_edges.at(node_id);
}

void DFSUtil(const adjacency_graph *graph, size_t src, std::unordered_set<size_t> &visited) {
    visited.insert(src);

    for (size_t edge_id : graph->nodes.at(src)) {
        size_t dst = graph->edge_ids.at(edge_id).second;
        if (!visited.contains(dst))
            DFSUtil(graph, dst, visited);
    }
}

void printAllPathsUtil(const adjacency_graph *graph, size_t u, size_t d, std::unordered_set<size_t> &visited,
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
                printAllPathsUtil(graph, dst, d, visited, path, path_index, visited_src_dst, global);
        }
    }

    // Remove current vertex from path[] and mark it as unvisited
    path_index--;
    visited.erase(u);
}

weigthed_labelled_automata::weigthed_labelled_automata() : adjacency_graph{} {
    casusu = WEIGHTED_LABELLED_GRAPH_CASE;
}

double edge_cost(const adjacency_graph *graph, size_t edge_id) {
    assert(edge_id < graph->E_size);
    if (graph->casusu == ADJACENCY_GRAPH_CASE) {
        return 1.0;
    } else {
        return ((weigthed_labelled_automata*)graph)->edge_weight.at(edge_id);
    }
}

const std::string &node_label(const adjacency_graph *graph, size_t node_id) {
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

void dot(adjacency_graph *graph, std::ostream &os) {
    os << "digraph finite_state_machine {\n"
          "    rankdir=LR;\n"
          "    size=\"8,5\"\n";
    for (int node_id = 0, N = graph->nodes.size(); node_id<N; node_id++) {
        std::string shape = "circle";

        os << "node [shape = circle, label=\"";
        if (graph->casusu == WEIGHTED_LABELLED_GRAPH_CASE) {
            os << ((weigthed_labelled_automata*)graph)->node_label.at(node_id);
        } else {
            os << std::to_string(node_id);
        }
        os << "\", fontsize=10] q" << node_id << ";\n";
    }
    os << "\n\n";
    for (int node_id = 0, N = graph->nodes.size(); node_id<N; node_id++) {
        std::vector<size_t> outgoing = graph->nodes.at(node_id);
        for (const size_t edge_id : outgoing) {
            size_t dst = graph->edge_ids.at(edge_id).second;
            os << "q" << node_id << " -> q" << dst << ";\n";
        }
    }
    os << "}";
}

