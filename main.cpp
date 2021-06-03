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

    typename TransitionGraph<DBL>::Matrix T1(graph.V_size, Y.size()),
                                        T2(graph.V_size, Y.size());


    // T2 is a zero matrix

    // Given that T1 is initialized with the initial probability, which is just src,
    // then I just need to check whether the given automata might have a string starting as Y
    {
        const std::string& ref = Y.at(0);
        if (node_label(&graph, src) == ref) {
            T1(src, 0) = 1.0;
        }
    }

    // Mapping each label within the string
    for (size_t j = 1, NN = Y.size(); j < NN; j++) {
        const std::string& alpha = Y.at(j);

        std::cerr << alpha << std::endl;

        auto it = graph.label_conversion.find(alpha);
        if (it != graph.label_conversion.end()) {
            for (const size_t i : it->second) {
                DBL kprev = -1;
                size_t karg = 0;
                for (const size_t edge_k : getIngoingEdgesId(&graph, i)) {
                    size_t k = graph.edge_ids.at(edge_k).first;
                    DBL kval = T1(k, j-1) * graph.edge_weight.at(edge_k);
                    if (kprev < kval) {
                        kprev = kval;
                        karg = k;
                    }
                }
                if (kprev > -1) {
                    T1(i, j) = kprev;
                    T2(i, j) = karg;
                }
            }
        }
    }

    std::cout << T1 << std::endl;
    std::cout << T2 << std::endl;

    DBL    k;
    size_t karg;
    for (size_t karg_curr = 0; karg_curr< graph.V_size; karg_curr++) {
        DBL kcurr = T1(karg_curr, T-1);
        if ((karg_curr == 0) || (kcurr > k)) {
            k = kcurr;
            karg = karg_curr;
        }
    }

    probabilisitc_model_trace pmt;
    auto& ref = pmt.underlying_sequences.emplace_back();
    ref.sequence.emplace_back(karg);
    pmt.t_with_prob.t.emplace_back(graph.inv_label_conversion.at(karg));

    for (size_t j = T-1; j>0; j--) {
        karg = T2(karg, j);
        ref.sequence.emplace_back(karg);
        pmt.t_with_prob.t.emplace_back(graph.inv_label_conversion.at(karg));
    }

    std::reverse(ref.sequence.begin(), ref.sequence.end());
    std::reverse(pmt.t_with_prob.t.begin(), pmt.t_with_prob.t.end());
    pmt.t_with_prob.probability = k;
    return pmt;

#if 0
    size_t i = src;
    size_t j = 0;
    T1(i, graph.final_labels.at(Y.at(0))) =
            graph.inv_label_conversion.at(i) == Y.at(0) ? 1.0 : 0.0;

    std::vector<size_t> node_ids;
    for (const auto& cp : graph.inv_label_conversion) {
        node_ids.emplace_back(cp.first);
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
#endif

}

int main() {

    weigthed_labelled_automata graph;
    from_string(graph, {"ciao", "mamma", "guarda", "come", "mi", "diverto"});
    dot(&graph, std::cout);/*
    std::cout <<
              viterbi_top1<double>(graph, 0, 5, {"ciao", "mamma", "guarda", "come", "mi", "diverto"}) << std::endl;
    std::cout <<
              viterbi_top1<double>(graph, 0, 5, {"ciao", "mamma"}) << std::endl;*/
    std::cout <<
              viterbi_top1<double>(graph, 0, 5, {"ciao", "oddio", "mamma", "guarda", "bene"}) << std::endl;
}
