/*
 * TransitionGraph.cpp
 * This file is part of ProbabilisticTraceAlignment2
 *
 * Copyright (C) 2021 - Giacomo Bergami
 *
 * ProbabilisticTraceAlignment2 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * ProbabilisticTraceAlignment2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ProbabilisticTraceAlignment. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef COMPLETE_MATCHING_TRANSITIONGRAPH_H
#define COMPLETE_MATCHING_TRANSITIONGRAPH_H


#include <string>
#include <unordered_map>
#include <cstdio>

#include <minauto/Util.hpp>
#include <minauto/Timer.hpp>
#include <minauto/BLAS.hpp>
#include <minauto/Automaton.hpp>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>
#include "adjacency_graph.h"

template <typename T>
struct TransitionGraph {
    typedef boost::numeric::ublas::compressed_matrix<T,boost::numeric::ublas::row_major> Matrix;
    /*struct Triplet {
        size_t src, dst;
        T      value;

        Triplet() : Triplet(0,0,1.0) {}
        Triplet(size_t src, size_t dst, T value) : src(src), dst(dst), value(value) {}
        Triplet(const Triplet& ) = default;
        Triplet(Triplet&& ) = default;
        Triplet& operator=(const Triplet& ) = default;
        Triplet& operator=(Triplet&& ) = default;
    };*/

    //std::unordered_map<size_t, std::string> inv_label_conversion;   // Maps a node-id to an associated label
    //std::vector<std::string>                tmp_labels;
    //std::unordered_map<std::string, size_t> final_labels;           // Maps a label to an unique integer
    //std::vector<Triplet>                    tripletList;
    const weigthed_labelled_automata& graph;

    size_t nodes = 0;
    size_t edges = 0;
    double weight = 1.0;
    size_t source = 0;
    size_t target = 0;
    bool matrices_initialized = false;
    Matrix L, R;

    TransitionGraph(const weigthed_labelled_automata& graph, size_t src, size_t dst, double cost = 1.0, bool im = true) :
        source{src}, target{dst}, weight{cost}, graph{graph} {
        nodes = graph.V_size;
        edges = graph.E_size;

        if (im) {
            initialize_matrices();
            matrices_initialized = true;
        }
    }

    void initialize_matrices() {
        {
            Matrix rectMatrix(graph.final_labels.size()+1, nodes+1);
            for (const auto& cp: graph.inv_label_conversion) {
                rectMatrix(graph.final_labels.at(cp.second), cp.first) = 1;
            }
            std::swap(L, rectMatrix);
        }
        {
            Matrix matrix(nodes+1, nodes+1);
            for (size_t i = 0, N = graph.edge_ids.size(); i<N; i++) {
                const auto& ref = graph.edge_ids.at(i);
                matrix(ref.first, ref.second) = graph.edge_weight.at(i);
            }

            std::swap(R, matrix);
        }
        matrices_initialized = true;
    }
    TransitionGraph(const TransitionGraph& ) = default;
    TransitionGraph(TransitionGraph&& ) = default;
    TransitionGraph& operator=(const TransitionGraph& ) = default;
    TransitionGraph& operator=(TransitionGraph&& ) = default;




#if 0
    /**
     * Initialization of some initial values
     * @param vSize
     * @param eSize
     * @param src
     * @param tgt
     */
    void init(size_t vSize, size_t eSize, size_t src, size_t tgt) {
        nodes = vSize;
        edges = eSize;
        source = src;
        target = tgt;
    }

  /**
  * Adds a id-label association to the graph
  * @param id        node id
  * @param label     label id
  * @return Whether the node was correctly inserted (true) or if the node was already inserted (false)
  */
    bool addNode(size_t id, const std::string& label) {
        nodes = std::max(nodes, id);
        tmp_labels.emplace_back(label);
        return inv_label_conversion.insert(std::make_pair(id, label)).second;
    }

    /**
     * Adds an edge to the matrix (bulk insertion)
     * @param src       node src
     * @param dst       node dst
     * @param weight    edge weight, src->dst
     */
    void addEdge(ssize_t src, ssize_t dst, T weight) {
        assert(0 <= weight && weight <= 1);
        assert(0 <= src && src <= nodes);
        assert(0 <= dst && dst <= nodes);
        if (weight != 0) tripletList.emplace_back(src,dst,weight);
    }

    void finalizeEdgesMatrix(double cost, bool store_matrices = true) {
        weight = cost;
        edges = tripletList.size();
        nodes++;

        // Creating the unique vector tmp_labels
        {
            std::sort(tmp_labels.begin(), tmp_labels.end() );
            tmp_labels.erase(std::unique(tmp_labels.begin(), tmp_labels.end() ), tmp_labels.end() );

            for (size_t i = 0, N = tmp_labels.size(); i < N; i++)
                final_labels[tmp_labels.at(i)] = i;
            tmp_labels.clear();
            if (store_matrices) {
                Matrix rectMatrix(tmp_labels.size()+1, nodes+1);
                for (const auto& cp: inv_label_conversion) {
                    rectMatrix(final_labels.at(cp.second), cp.first) = 1;
                }
                L = rectMatrix;
            }
        }

        // Creating the R matrix
        if (store_matrices) {
            Matrix matrix(nodes+1, nodes+1);
            auto it = tripletList.begin();
            while(it != tripletList.end()) {
                matrix(it->src, it->dst) = it->value;
                it = tripletList.erase(it);
            }

            R = matrix;
        }
    }

#endif

};



#endif //COMPLETE_MATCHING_TRANSITIONGRAPH_H
