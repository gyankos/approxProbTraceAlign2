#include <iostream>


#include <hashing/pair_hash.h>
#include <hashing/vector_hash.h>

#include <basics/trace.h>
#include <graphs/TransitionGraph.h>
#include <graphs/adjacency_graph.h>


#include <minauto/RandomValueMap.hpp>
#include <minauto/Worklist.hpp>
#include <minauto/GramSchmidt.hpp>
#include <minauto/WordTree.hpp>
#include <minauto/Minimisation.hpp>
#include <minauto/EquivAlgorithms.hpp>

using namespace jackbergus::fuzzyStringMatching3;
using namespace jackbergus::fuzzyStringMatching3::graphs;

template <typename DBL> probabilisitc_model_trace viterbi_top1(const weigthed_labelled_automata& graph,
                                                                size_t src,
                                                                size_t dst,
                                                                const log_trace& Y,
                                                                bool use_zero_as_minimum_B_cost = true) {
    const size_t K = graph.V_size;
    const size_t T = Y.size();

    typename TransitionGraph<DBL>::Matrix T1(graph.V_size, Y.size()),
                                        T2(graph.V_size, Y.size());
    double minimum_cost = use_zero_as_minimum_B_cost ? 0.0 : graph.minimum_edge_weight / 2.0;

    // T2 is a zero matrix

    // Given that T1 is initialized with the initial probability, which is just src,
    // then I just need to check whether the given automata might have a string starting as Y
    {
        const std::string& ref = Y.at(0);
        if (node_label(&graph, src) == ref) {
            T1(src, 0) = 1.0;
        } else {
            T1(src, 0) = minimum_cost;
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

struct AutomatonDelta {
                    size_t src;
                    size_t target;
                    double weight;
                    std::string label;

    AutomatonDelta() : src{0}, target{0}, weight{1}, label("") {}
    AutomatonDelta(size_t src, size_t target, double weight, const std::string &label) : src(src), target(target),
                                                                                         weight(weight), label(label) {}
    AutomatonDelta(const AutomatonDelta& ) = default;
    AutomatonDelta(AutomatonDelta&& ) = default;
    AutomatonDelta& operator=(const AutomatonDelta& ) = default;
    AutomatonDelta& operator=(AutomatonDelta&& ) = default;
};


template<typename Automaton>
bool check(	MinimisationSettings& msettings,Settings& settings, Automaton& a, Automaton& b, Precision<typename Automaton::Weight>& prec) {

    /*
      typedef ::Automaton<mpfr::mpreal> MPAutomaton ;

        mpfr::mpreal::set_default_prec(256);

        mpfr::mpreal absolute_precision =  (check_absolute_precision == 1e-6) ? std::numeric_limits<mpfr::mpreal>::epsilon() : check_absolute_precision;

        Precision<mpfr::mpreal> prec(absolute_precision,check_relative_precision);
        */
        Automaton c(a);
        Automaton d(b);

        Automaton diff(c,d);

    unsigned k = a.getNrOfStates();

    mpfr::mpreal maximum = kmax<Automaton,true> ( diff, k, settings, prec);

    if(global_setting::verbosity >= 1) {
        std::cout << "Maximum difference " << maximum << " after " << k << " steps " <<std::endl;
    }
    return true; // equiv(c,d,settings,prec);
}

void test_automata_minimization(int n_state, int initial, int final, const std::vector<AutomatonDelta>& delta,
                                int inflation_precision = -1) {

    // Initialization of the automaton
    Automaton<double> automaton, result;

    automaton.setNrOfStates(n_state);
    automaton.setInitial(initial, 1.0);
    automaton.setAccepting(final, 1.0);
    for (const auto& edge: delta) {
        automaton.addTransition(edge.src, edge.target, edge.label, edge.weight);
    }

    // Finalization
    automaton.computeMatrices();


    double absolute_precision = 1E-16,
            relative_precision = 1E-06;
    bool custom_absolute_precision = false;

    MinimisationSettings msettings;
    msettings.lz = false;
    msettings.direction = MinimisationSettings::both;
    msettings.pivot = false;
    msettings.reortho = false;
    msettings.dense = false;

    Settings settings;
    settings.pivot = false;
    settings.reortho = false;
    settings.ce = false;
    settings.onion = false;
    settings.normalise = false;

    if(inflation_precision!=-1) {
        mpfr::mpreal::set_default_prec(inflation_precision);
    }
    absolute_precision =  (custom_absolute_precision == false) ? std::numeric_limits<double>::epsilon() : absolute_precision;
    Precision<double> prec(absolute_precision,relative_precision);

    minimise(automaton, result, msettings, prec);
    std::cout<<"   # states: " << result.getNrOfStates()
             << " # trans: "   << result.getNrOfNonzeros()
             << " # letter " << alphabet.size()<< std::endl;


    {
        result.aiSee(std::cout);
    }

    if(global_setting::check_minimisation) {
        check<Automaton<double>>(msettings,settings,automaton,result,prec);
    }
}

void test_viterbi_1() {
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

#include <cmath>
#include <set>


void generate_my_minimization(double init_prob = 1.0, double stop_probability = 1E-06) {
    weigthed_labelled_automata wla, out;
    std::unordered_map<global_node_properties, size_t> node_compact_info;
#if 1
    size_t n0 = add_node(&wla, "a");
    size_t n1 = add_node(&wla, "a");
    size_t n2 = add_node(&wla, "a");
    size_t n3 = add_node(&wla, "b");
    size_t n4 = add_node(&wla, "c");
    size_t n5 = add_node(&wla, "b");

    add_edge(&wla, n0, n1, 0.5);
    add_edge(&wla, n0, n2, 0.3);
    add_edge(&wla, n0, n3, 0.2);
    add_edge(&wla, n1, n1, 0.5);
    add_edge(&wla, n1, n2, 0.3);
    add_edge(&wla, n1, n3, 0.2);
    add_edge(&wla, n2, n2, 0.3);
    add_edge(&wla, n2, n4, 0.7);
    add_edge(&wla, n3, n5, 1.0);
    add_edge(&wla, n5, n4, 1.0);

    //TransitionGraph<double> tg(wla, n0, n4);
    determinization_information info{stop_probability, n0, n4};
#else
    size_t n0 = add_node(&wla, "a");
    size_t n1 = add_node(&wla, "b");
    size_t n2 = add_node(&wla, "b");
    size_t n3 = add_node(&wla, "c");
    size_t n4 = add_node(&wla, "a");
    size_t n5 = add_node(&wla, "c");
    size_t n6 = add_node(&wla, "a");
    size_t n7 = add_node(&wla, "f");

    bool withLoops = true;

    add_edge(&wla, n0, n1, 0.5);
    add_edge(&wla, n1, n3, 0.5);
    add_edge(&wla, n1, n4, 0.5);
    add_edge(&wla, n4, n7, 1.0);
    add_edge(&wla, n0, n2, 0.5);
    add_edge(&wla, n2, n5, 0.5);
    add_edge(&wla, n2, n6, 0.5);
    add_edge(&wla, n6, n7, 1.0);
    if (withLoops) {
        add_edge(&wla, n3, n1, 0.5);
        add_edge(&wla, n5, n2, 0.5);
    }
    add_edge(&wla, n3, n7, withLoops ? 0.5 : 1.0);
    add_edge(&wla, n5, n7, withLoops ? 0.5 : 1.0);

    //TransitionGraph<double> tg(wla, n0, n4);
    determinization_information info{stop_probability, n0, n7};
#endif

    probabilisitc_model_trace pmt;
    pmt.t_with_prob.probability = init_prob;
    pmt.t_with_prob.t.emplace_back(node_label(&wla, n0));
    pmt.underlying_sequences.emplace_back(init_prob, n0);
    //std::cout << pmt << std::endl;

    // Preliminary minimizing the visit and the complexity by compacting the edge traversing: sum up the edges' probability
    // leading to the same outgoing state
    edge_compacting(wla);

    // Making the NFA into a DFA, by exploiting the prefix strategy
    nfa_to_dfa_weighted_labelled_automata(wla, pmt, out, info, node_compact_info);

    // While creating the wheel state, I might have multiple outgoing edges leading to the same state. Therefore, I
    // could compact those into one single edge by summing up the weights
    edge_compacting(out);



    dot(&out, std::cout);
}

void automata_minimization() {
    test_automata_minimization(4, 0, 3, {{0, 1, 0.5, "a"},
                                         {0, 2, 0.5, "b"},
                                         {1, 1, 0.5, "c"},
                                         {2, 2, 0.5, "c"},
                                         {1, 3, 0.5, "d"},
                                         {2, 3, 0.5, "d"}});
}

int main() {

    generate_my_minimization();
    //automata_minimization();
}
