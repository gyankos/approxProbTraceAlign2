/*
 * adjacency_graph.h
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

#ifndef CLASSIFIERS_ADJACENCY_GRAPH_H
#define CLASSIFIERS_ADJACENCY_GRAPH_H

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <ostream>


namespace jackbergus {
    namespace fuzzyStringMatching3 {
        namespace graphs {

            enum graph_cases {
                ADJACENCY_GRAPH_CASE,
                WEIGHTED_LABELLED_GRAPH_CASE
            };

            struct adjacency_graph {
                size_t V_size, E_size;
                graph_cases casusu;

                std::vector<std::vector<size_t>> nodes;
                std::vector<std::pair<size_t, size_t>> edge_ids;
                std::unordered_map<size_t, std::vector<size_t>> ingoing_edges;

                adjacency_graph();
                adjacency_graph(const adjacency_graph&) = default;
                adjacency_graph(adjacency_graph&&) = default;
                adjacency_graph& operator=(const adjacency_graph&) = default;
                adjacency_graph& operator=(adjacency_graph&&) = default;
            };

            size_t add_node(adjacency_graph* graph, const std::string& label = "");
            size_t add_edge(adjacency_graph* graph, size_t src, size_t dst, double cost = 1.0);

/**
 *
 * @param graph
 * @param edge_id           Edge to be changed
 * @param new_dst           New destination for the edge
 * @param cost_to_update    New cost for the edge
 * @return
 */
            bool update_edge_target(adjacency_graph *graph, size_t edge_id, size_t new_dst, double* cost_to_update);

            bool remove_edge(adjacency_graph *graph, size_t edge_id);

            std::pair<size_t, size_t> add_undirected_edge(adjacency_graph* graph, size_t src, size_t dst, double cost = 1.0);
            const std::pair<size_t, size_t>& edge_from_id(const adjacency_graph* graph, size_t edge_id);
            const std::vector<size_t> &getOutgoingEdgesId(adjacency_graph* graph, size_t node_id);
            const std::vector<size_t> &getIngoingEdgesId(const adjacency_graph* graph, size_t node_id);
            void DFSUtil(const adjacency_graph* graph, size_t src, std::unordered_set<size_t> &visited);
            void printAllPathsUtil(const adjacency_graph* graph, size_t u, size_t d, std::unordered_set<size_t> &visited, std::vector<ssize_t> &path,
                                   size_t path_index, std::unordered_set<size_t> &visited_src_dst, std::unordered_set<size_t> &global);


            struct weigthed_labelled_automata : public adjacency_graph {
                std::vector<double> edge_weight;
                std::vector<std::string> node_label;
                std::unordered_map<size_t, std::string> inv_label_conversion;          // Maps a node-id to an associated label
                std::unordered_map<std::string, std::vector<size_t>> label_conversion; // Maps a label to the nodes associated to it
                std::unordered_map<std::string, size_t> final_labels;                  // Maps a label to an unique-id associated to it
                size_t count_labels = 0;
                double minimum_edge_weight;

                weigthed_labelled_automata();
                weigthed_labelled_automata(const weigthed_labelled_automata&) = default;
                weigthed_labelled_automata(weigthed_labelled_automata&&) = default;
                weigthed_labelled_automata& operator=(const weigthed_labelled_automata&) = default;
                weigthed_labelled_automata& operator=(weigthed_labelled_automata&&) = default;
            };

            double edge_cost(const adjacency_graph* graph, size_t edge_id);
            const std::string& node_label(const adjacency_graph* graph, size_t node_id);
            void dot(adjacency_graph* graph, std::ostream &os);
            void from_string(weigthed_labelled_automata& graph, const std::vector<std::string>& trace);

            void edge_compacting(weigthed_labelled_automata& graph);
        }
    }
}

#endif //CLASSIFIERS_ADJACENCY_GRAPH_H