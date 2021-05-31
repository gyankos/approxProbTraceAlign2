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
#include <minauto/Automaton.hpp>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>

template <typename T>
struct TransitionGraph {
    typedef boost::numeric::ublas::compressed_matrix<T,boost::numeric::ublas::row_major> Matrix;
    struct Triplet {
        size_t src, dst;
        T      value;

        Triplet() : Triplet(0,0,1.0) {}
        Triplet(size_t src, size_t dst, T value) : src(src), dst(dst), value(value) {}
        Triplet(const Triplet& ) = default;
        Triplet(Triplet&& ) = default;
        Triplet& operator=(const Triplet& ) = default;
        Triplet& operator=(Triplet&& ) = default;
    };

    std::unordered_map<size_t, std::string> inv_label_conversion;
    std::vector<std::string>                tmp_labels;
    std::unordered_map<std::string, size_t> final_labels;
    std::vector<Triplet>                    tripletList;

    size_t nodes = 0;
    size_t edges = 0;
    double weight = 1.0;
    size_t source = 0;
    size_t target = 0;

    Matrix L, R;

    TransitionGraph() {}
    TransitionGraph(const TransitionGraph& ) = default;
    TransitionGraph(TransitionGraph&& ) = default;
    TransitionGraph& operator=(const TransitionGraph& ) = default;
    TransitionGraph& operator=(TransitionGraph&& ) = default;

    TransitionGraph(const std::string &filename) {
        FILE* file = fopen(filename.c_str(), "r");
        bool error = false;
        if (file) {
            int i;
            double w = 1;
            // Reading the number of the nodes
            i = fscanf(file, "nodes: %zd\n", &nodes);
            error = (i == EOF || (i != 1));
            if (error) return;
            i = fscanf(file, "edges: %zd\n", &edges);
            error = (i == EOF || (i != 1));
            if (error) return;
            i = fscanf(file, "source: %zd\n", &source);
            error = (i == EOF || (i != 1));
            if (error) return;
            i = fscanf(file, "target: %zd\n", &target);
            error = (i == EOF || (i != 1));
            if (error) return;
            i = fscanf(file, "weight: %lf\n", &w);
            error = (i == EOF || (i != 1));
            if (error) return;

            for (size_t j = 0; j<nodes; j++) {
                size_t node_no;
                char string[124];
                std::string k;
                i = fscanf(file, "%zd %123s\n", &node_no, string);
                error = (i == EOF || (i != 2));
                if (error) return;
                k = std::string(string);
                if (!addNode(node_no, k)) return;
            }


            for (size_t j = 0; j<edges; j++) {
                size_t src, dst;
                double weight = 0.0;
                i = fscanf(file, "%zd %zd %lf\n", &src, &dst, &weight);
                error = (i == EOF || (i != 3));
                if (error) return;
                assert(std::abs(weight) <= 1.0);
                addEdge(src, dst, weight);
            }
            finalizeEdgesMatrix(w);
            fclose(file);
        }
    }

    TransitionGraph fromString(const std::string &string, double stringWeight) {
        assert(!string.empty());
        TransitionGraph rg;
        size_t n = string.size();
        rg.init(n, n-1, 0, n-1);
        for (size_t i = 0; i<n-1; i++) {
            rg.addNode(i, std::string{string[i]});
            rg.addEdge(i, i+1, 1.0);
        }
        rg.addNode(n-1, std::string{string[n-1]});
        rg.finalizeEdgesMatrix(stringWeight);
        return rg;
    }

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

    void finalizeEdgesMatrix(double cost) {
        weight = cost;
        edges = tripletList.size();
        nodes++;

        // Creating the unique vector tmp_labels
        {
            std::sort(tmp_labels.begin(), tmp_labels.end() );
            tmp_labels.erase(std::unique(tmp_labels.begin(), tmp_labels.end() ), tmp_labels.end() );
            Matrix rectMatrix(tmp_labels.size()+1, nodes+1);

            for (size_t i = 0, N = tmp_labels.size(); i < N; i++)
                final_labels[tmp_labels.at(i)] = i;
            tmp_labels.clear();
            for (const auto& cp: inv_label_conversion) {
                rectMatrix(final_labels.at(cp.second), cp.first) = 1;
            }
            L = rectMatrix;
        }

        // Creating the R matrix
        {
            Matrix matrix(nodes+1, nodes+1);
            auto it = tripletList.begin();
            while(it != tripletList.end()) {
                matrix(it->src, it->dst) = it->value;
                it = tripletList.erase(it);
            }

            R = matrix;
        }
    }
};



#endif //COMPLETE_MATCHING_TRANSITIONGRAPH_H
