#include <iostream>

#include <basics/trace.h>
#include <graphs/TransitionGraph.h>
#include <graphs/adjacency_graph.h>

using namespace jackbergus::fuzzyStringMatching3;


template <typename DBL> probabilisitc_model_trace viterbi_top1(const weigthed_labelled_automata& graph,
                                                                size_t src,
                                                                size_t dst,
                                                                const log_trace& Y) {
    const size_t K = graph.V_size;
    const size_t T = Y.size();
    typename TransitionGraph<DBL>::Matrix T1(graph.nodes, Y.size()),
                                        T2(graph.nodes, Y.size());

    size_t i = src;
    size_t j = 0;
    T1(i, graph.final_labels(Y.at(0))) =
            graph.inv_label_conversion.at(i) == Y.at(0) ? 1.0 : 0.0;

    std::vector<size_t> node_ids;
    for (const auto& cp : graph.inv_label_conversion) {
        node_ids.emplace_back(cp.first);
    }

    std::unordered_map<std::string, std::vector<std::pair<size_t,size_t>>> yVector_labelid_with_offset;
    for (size_t i = 0, NN = Y.size(); i<NN; i++) {
        const std::string& ref = Y.at(i);
        yVector_labelid_with_offset[ref].emplace_back(graph.final_labels.at(ref), i);
    }

    size_t yj;
    for (const auto& cp : graph.inv_label_conversion) {
        i = cp.first;                         // node-id
        j = graph.final_labels.at(cp.second); // label-id

        auto it = yVector_labelid_with_offset.find(cp.second);
        if (it != yVector_labelid_with_offset.end()) {
            for (size_t ii = 0, NN = it->second.size(); ii<NN; ii++) {
                yj = it->second[ii].first;
                j = it->second[ii].second;
                if (j == 0) continue;

                DBL k = -1;
                size_t karg = node_ids.at(0);

                for (const size_t karg_curr : node_ids) {
                    DBL curr = T1(karg_curr, j-1) * graph.R(karg_curr, i);
                    if (curr > k) {
                        k = curr;
                        karg = karg_curr;
                    }
                }

                T1(i, j) = k;
                T2(i, j) = karg;
            }
        } else {

        }
    }


}

int main() {
}
