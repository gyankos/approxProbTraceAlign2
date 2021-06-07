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

/**
 * Defines the global node properties for merging nodes before minimizing the automaton.
 * I.e., a node is defined by its node id, the associated label, and the outgoing edges' labels (from the
 * target node) as well as their associated cost. This will be used to identify when we reached a similar node
 */
struct global_node_properties {
    std::set<size_t>                            node_id;
    std::string                                 node_label;
#ifdef DO_STRINGS
    std::vector<std::pair<std::string, std::string>>
#else
    std::vector<std::pair<std::string, double>>
#endif
    outgoing_edges;

    global_node_properties(const std::string &nodeLabel) : node_label(nodeLabel)  {
        sorted = false;
    }
    void insert_edge(const std::string& s, double w) {
#ifdef DO_STRINGS
        outgoing_edges.emplace_back(s, std::to_string(w));
#else
        outgoing_edges.emplace_back(s, (w));
#endif
        sorted = false;
    }
    void finalize() {
        if (!sorted) {
            std::sort(outgoing_edges.begin(), outgoing_edges.end());
            sorted = true;
        }
    }

    bool operator==(const global_node_properties &rhs) const {
        assert(sorted);
        assert(rhs.sorted);
        return node_id == rhs.node_id &&
               node_label == rhs.node_label &&
               outgoing_edges == rhs.outgoing_edges;
    }
    bool operator!=(const global_node_properties &rhs) const {
        return !(rhs == *this);
    }
              private:
                  bool sorted = false;

public:
    friend std::ostream &operator<<(std::ostream &os, const global_node_properties &properties) {
        std::vector<size_t> S{properties.node_id.begin(), properties.node_id.end()};
        os << "node_id: " << S << " node_label: " << properties.node_label << " outgoing_edges: "
           << properties.outgoing_edges << " sorted: " << properties.sorted;
        return os;
    }
};


namespace std {
    /**
     * Hashing function associated to global_node_properties
     */
    template <>
    struct hash<global_node_properties> {
        std::size_t operator()(const global_node_properties& k) const {
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

struct determinization_information {
    ssize_t has_result_well_state = -1;
    ssize_t is_final_state_inserted_into_result = -1;
    double  halting_condition = 0.0;
    size_t  tg_initial, tg_final;

    determinization_information(double haltingCondition, size_t  tg_initial, size_t  tg_final) : tg_initial(tg_initial), tg_final(tg_final),  halting_condition(haltingCondition), has_result_well_state{-1}, is_final_state_inserted_into_result{-1} {}

    determinization_information(const determinization_information& ) = default;
    determinization_information(determinization_information&& ) = default;
    determinization_information& operator=(const determinization_information& ) = default;
    determinization_information& operator=(determinization_information&& ) = default;
};

ssize_t nfa_to_dfa_weighted_labelled_automata(weigthed_labelled_automata& graph,
                                              const probabilisitc_model_trace& MT,
                                              weigthed_labelled_automata& out,
                                              determinization_information& info,
                                              std::unordered_map<global_node_properties, size_t>& node_compact_info) {

    constexpr double limit = std::numeric_limits<double>::epsilon();
    if ((std::abs(MT.t_with_prob.probability - info.halting_condition) < limit) ||
            (MT.t_with_prob.probability < limit)) return -1;

    //std::unordered_set<size_t> curr;

    std::string curr_node_label = *MT.t_with_prob.t.rbegin();
    global_node_properties gnp{curr_node_label};
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
        std::cout << gnp << std::endl;
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
