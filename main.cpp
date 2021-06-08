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
#include "graphs/algorithms/nfa_to_dfa_weighted_labelled_automata.h"

using namespace jackbergus::fuzzyStringMatching3;
using namespace jackbergus::fuzzyStringMatching3::graphs;
using namespace jackbergus::fuzzyStringMatching3::graphs::algorithms;

#if 0
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
                for (const size_t edge_k : getIngoingEdgesId(&graph, i)) { // TODO: hasIngoingEdges!
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
#endif

#include <cmath>
#include <set>

void extend_set(const std::vector<jackbergus::fuzzyStringMatching3::log_trace> &postfixes,
                const std::string &postfix,
                std::vector<jackbergus::fuzzyStringMatching3::log_trace> &extended) {
    for (auto logtrace : postfixes) {
        logtrace.emplace_back(postfix);
        extended.emplace_back(logtrace);
    }

}

#include <vector>

template<typename T>
void unique_sorted_vector(std::vector<T> &Xyields) {
    if (Xyields.empty()) return;
    std::sort(Xyields.begin(), Xyields.end());
    Xyields.erase(std::unique(Xyields.begin(), Xyields.end()), Xyields.end());
}

#include <vector>

template<typename T, typename F>
void unique_sorted_vector_cmp(std::vector<T> &Xyields, F comp) {
    if (Xyields.empty()) return;
    std::sort(Xyields.begin(), Xyields.end(), comp);
    Xyields.erase(std::unique(Xyields.begin(), Xyields.end()), Xyields.end());
}


#include <stack>

/**
 * Providing an approximated topological sort, where the loops are ignored
 *
 * @param g                 Graph to be approximately visited
 * @param initial_node      The only graph node having no ingoing edges
 * @param vec               output: The vector containing the topologically sorted vertices
 */
void approximated_topo_sort(adjacency_graph *g, size_t initial_node, std::vector<size_t> &vec) {
    std::vector<size_t> order(g->V_size, 0);
    std::vector<bool> visited(g->V_size, false);
    size_t curr_pos = g->V_size;

    std::stack<std::pair<size_t, bool>> visitStack;
    visitStack.emplace(initial_node, false);

    // DFS Visit
    while (!visitStack.empty()) {
        size_t u = 0;
        bool u_visited = false;
        std::tie(u, u_visited) = visitStack.top();
        visitStack.pop();
        if (u_visited) {
            // If it was already visited before continuing the visit, then put it back to the stack
            order[curr_pos--] = u;
        } else {
            if (!visited.at(u)) {
                // Visit the outoing edges
                visited[u] = true;
                const auto &ref = getOutgoingEdgesId(g, u);
                if (!ref.empty()) {
                    // Remember to visit the node back after visiting the adjacent nodes
                    visitStack.emplace(u, true);
                    for (size_t i = 0, N = ref.size(); i < N; i++) {
                        // Visiting the adjacent nodes
                        visitStack.emplace(g->edge_ids.at(ref[i]).second, false);
                    }
                } else {
                    order[curr_pos--] = u;
                }
            }
        }
    }

    std::swap(vec, order);
}

/**
 * Provides the approximately longest path algorithm, by ignoring loops that might diverge the algorithm
 *
 * @param g                 Graph to visit
 * @param initial_node      Initial node having no ingoing edges from which start the visit
 * @return                  Associates to each node id (vector index) its order (value)
 */
std::vector<size_t> approx_longest_path(adjacency_graph *g, size_t initial_node) {
    std::vector<size_t> topo_order;
    std::vector<size_t> neighbour_size(g->V_size, 0);
    approximated_topo_sort(g, initial_node, topo_order);
    for (size_t i = 0; i < (g->V_size); i++) {
        size_t max_size = 0;
        size_t iId = topo_order.at(i);
        ///std::cout << " u = " << iId << std::endl;
        const auto it = hasIngoingEdges(g, iId);
        if (it != g->ingoing_edges.cend()) {
            const auto &ingoing = it->second;
            for (size_t j = 0, N = ingoing.size(); j < N; j++) {
                ///std::cout << "\t" << g->edge_ids.at(ingoing.at(j)) << std::endl;
                size_t nSize = neighbour_size.at(g->edge_ids.at(ingoing.at(j)).first);
                if (nSize > max_size)
                    max_size = nSize;
            }
        }

        max_size++;
        neighbour_size[iId] = max_size;
    }

    return neighbour_size;
}

#include <hashing/pair_hash.h>

void marca_rec(std::pair<size_t, size_t> cp,
               std::unordered_map<std::pair<size_t, size_t>, std::vector<std::pair<size_t, size_t>>> &has_lista,
               boost::numeric::ublas::compressed_matrix<char, boost::numeric::ublas::row_major> &Matrix) {
    if (cp.first == cp.second) return; // imposiburu
    if (cp.first > cp.second)
        std::swap(cp.first, cp.second);

    auto it = has_lista.find(cp);
    if (it != has_lista.end()) {
        auto cpV = it->second;
        has_lista.erase(it);
        Matrix(cp.first, cp.second) = (char) 0;
        for (size_t i = 0, N = cpV.size(); i < N; i++) {
            marca_rec(cpV.at(i), has_lista, Matrix);
        }
    } else {
        Matrix(cp.first, cp.second) = (char) 0;
    }
}


void heuristic_scala(weigthed_labelled_automata &graph, size_t initial_node, size_t final_state) {

    boost::numeric::ublas::compressed_matrix<char, boost::numeric::ublas::row_major> diagonalMatrix(graph.V_size,
                                                                                                    graph.V_size);
    {
        std::unordered_map<std::string, std::unordered_map<jackbergus::fuzzyStringMatching3::graphs::algorithms::determinization_node_eq_safe_heuristic, std::vector<size_t>>> M;
        {
            std::vector<size_t> VertexOrder = approx_longest_path(&graph, initial_node);
            ///std::cerr << "[Init M]" << std::endl;
            for (size_t i = 0, N = graph.V_size; i < N; i++) {
                if (i == final_state) continue;
                const std::string &i_label = graph.node_label.at(i);
                jackbergus::fuzzyStringMatching3::graphs::algorithms::determinization_node_eq_safe_heuristic hm{
                        i_label};
                for (size_t edge_out : getOutgoingEdgesId(&graph, i)) {
                    hm.insert_edge(graph.node_label.at(graph.edge_ids.at(edge_out).second),
                                   graph.edge_weight.at(edge_out));
                }
                hm.finalize();
                M[graph.node_label.at(i)][hm].emplace_back(i);
            }

            ///std::cerr << "[Init diagonalMatrix]" << std::endl;
            auto it = M.begin();
            while (it != M.end()) {
                auto it2 = it->second.begin();
                while (it2 != it->second.end()) {
                    // In this case, I am the only node similar to myself: therefore, I can skip it.
                    if (it2->second.size() <= 1)
                        it2 = it->second.erase(it2);
                    else {
                        size_t N = it2->second.size();
                        std::sort(it2->second.begin(), it2->second.end(),
                                  [&VertexOrder](const size_t u, const size_t v) {
                                      return VertexOrder.at(u) >= VertexOrder.at(v);
                                  });
                        for (size_t i = 0; i < N; i++) {
                            size_t u = it2->second.at(i);
                            for (size_t j = 0; j < i; j++) {
                                size_t v = it2->second.at(j);
                                if (u < v) {
                                    ///std::cout << "[" << u << "," << v << "]=1" << std::endl;
                                    diagonalMatrix(u, v) = (char) 1;
                                } else {
                                    ///std::cout << "[" << v << "," << u << "]=1" << std::endl;
                                    diagonalMatrix(v, u) = (char) 1;
                                }
                            }
                        }
                        it2++;
                    }
                }
                it++;
            }
            VertexOrder.clear();
        }


        // Please observe that the above step is just an heuristic, and I therefore cannot immediately assume that
        // all the nodes that are marked as similar are indeed similar, but we are indeed sure that all the remaining
        // nodes must be different!
        std::unordered_map<std::pair<size_t, size_t>, std::vector<std::pair<size_t, size_t>>> pair_to_remember;
        for (auto cp1 = M.begin(); cp1 != M.end(); cp1 = M.erase(cp1)) {
            for (const auto &cp2 : cp1->second) {
                size_t N = cp2.second.size();
                for (size_t i = 0; i < N; i++) {
                    size_t u = cp2.second.at(i);

                    std::unordered_map<std::string, size_t> edge_map;
                    for (size_t edge_id : getOutgoingEdgesId(&graph, u)) {
                        size_t target = graph.edge_ids.at(edge_id).second;
                        edge_map[graph.node_label.at(target)] = target;
                    }
                    for (size_t j = 0; j < i; j++) {
                        size_t v = cp2.second.at(j);

                        bool are_u_and_v_equivalent = true;
                        std::vector<std::pair<size_t, size_t>> to_emplace_current_elements;
                        for (size_t edge_id : getOutgoingEdgesId(&graph, v)) {
                            size_t target = graph.edge_ids.at(edge_id).second;
                            size_t opponent_target = edge_map.at(graph.node_label.at(target));

                            char val = (char) 1;
                            if (target < opponent_target) {
                                val = diagonalMatrix(target, opponent_target);
                            } else if (opponent_target < target) {
                                val = diagonalMatrix(opponent_target, target);
                            }
                            // If they are not equivalent nodes: please skip the iteration!
                            if (!val) {
                                are_u_and_v_equivalent = false;
                                to_emplace_current_elements.clear();
                                break;
                            } else {
                                if (target < opponent_target) {
                                    to_emplace_current_elements.emplace_back(target, opponent_target);
                                } else if (opponent_target < target) {
                                    to_emplace_current_elements.emplace_back(opponent_target, target);
                                }
                            }
                        }
                        if (!are_u_and_v_equivalent) {
                            std::pair<size_t, size_t> x;
                            if (u < v) {
                                diagonalMatrix(u, v) = (char) 0;
                                x.first = u;
                                x.second = v;
                            } else {
                                diagonalMatrix(v, u) = (char) 0;
                                x.first = v;
                                x.second = u;
                            }
                            marca_rec(x, pair_to_remember, diagonalMatrix);
                        } else {
                            size_t src, dst;
                            if (u < v) {
                                src = u;
                                dst = v;
                            } else {
                                src = v;
                                dst = u;
                            }
                            for (const auto &cp : to_emplace_current_elements) {
                                pair_to_remember[cp].emplace_back(src, dst);
                            }
                        }
                    }
                }
            }
        }
    }


    // Building up the equivalence class
    std::vector<size_t> Rv, Lv, tmpLv;
    std::unordered_map<size_t, size_t> backtrack_to;
    std::vector<std::pair<size_t, size_t>> EqR;
    std::unordered_map<size_t, std::vector<size_t>> eqClassesResolved;
    ///std::cerr << "[Print diagonalMatrix]" << std::endl;
    for (auto it1 = diagonalMatrix.begin1(); it1 != diagonalMatrix.end1(); it1++) {
        for (auto it2 = it1.begin(); it2 != it1.end(); it2++) {
            if (*it2) {
                size_t l = it2.index1(), r = it2.index2();
                Rv.emplace_back(r);
                Lv.emplace_back(l);
                EqR.emplace_back(l, r);
                auto it = backtrack_to.emplace(r, l);
                if (!it.second) {
                    it.first->second = std::min(it.first->second, l);
                }
            }
        }
    }
    unique_sorted_vector(Rv);
    unique_sorted_vector(Lv);
    std::set_difference(Lv.begin(), Lv.end(),
                        Rv.begin(), Rv.end(),
                        std::back_inserter(tmpLv));
    std::vector<bool> W(graph.V_size, false);
    for (size_t u : tmpLv) W[u] = true;
    Rv.clear();
    Lv = std::move(tmpLv); // Minimal-id elements belonging to the equivalence class
    tmpLv.clear();

    std::sort(EqR.begin(), EqR.end());
    ///std::cout << "EqR = " << EqR << std::endl;
    for (auto &cp : backtrack_to) {
        auto it = backtrack_to.find(cp.second);
        while (it != backtrack_to.end()) {
            cp.second = std::min(cp.second, it->second);
            it = backtrack_to.find(cp.second);
        }
    }
    for (auto &cp : EqR) {
        if (W.at(cp.first)) {
            eqClassesResolved[cp.first].emplace_back(cp.second);
            eqClassesResolved[cp.first].emplace_back(cp.first);
        } else {
            assert(backtrack_to.at(cp.first) == backtrack_to.at(cp.second));
            eqClassesResolved[backtrack_to.at(cp.first)].emplace_back(cp.first);
            eqClassesResolved[backtrack_to.at(cp.first)].emplace_back(cp.second);
        }
    }
    EqR.clear();
    W.clear();

    std::unordered_map<size_t, size_t> resolved_classes;
    for (size_t i = 0; i < graph.V_size; i++) {
        auto cp = eqClassesResolved.find(i);
        if (cp != eqClassesResolved.end()) {
            unique_sorted_vector(cp->second);
            size_t new_node = add_node(&graph, graph.node_label.at(i));
            for (size_t ij : cp->second)
                resolved_classes[ij] = new_node;
        }
    }

    for (const auto &classRef_eqClass : eqClassesResolved) {
        std::unordered_set<size_t>  outgoing_met_classes;
        size_t id_new_class = resolved_classes.at(classRef_eqClass.first);

        // Refactoring the outgoing edges
        auto outgoingCopy = graph.nodes.at(classRef_eqClass.first);
        for (size_t edge_id : outgoingCopy) {
            double edge_weight = graph.edge_weight.at(edge_id);
            std::pair<size_t, size_t> edge = graph.edge_ids.at(edge_id);
            if ((edge.first == (size_t) -1) && (edge.second == (size_t) -1)) continue;
            auto it2 = resolved_classes.find(edge.second);
            if ((it2 != resolved_classes.end())) {
                if ((!outgoing_met_classes.contains(it2->second))) {
                    edge.second = it2->second;
                    outgoing_met_classes.insert(it2->second);
                    add_edge(&graph, id_new_class, edge.second, edge_weight);
                }
            } else {
                add_edge(&graph, id_new_class, edge.second, edge_weight);
            }
            remove_edge(&graph, edge_id);
        }
    }

    for (const auto &classRef_eqClass : eqClassesResolved) {
        size_t id_new_class = resolved_classes.at(classRef_eqClass.first);
        for (size_t node_in_class : classRef_eqClass.second) {
            // Refactoring the ingoing edges
            auto it = graph.ingoing_edges.find(node_in_class);
            if (it != graph.ingoing_edges.end()) {
                auto ingoingCopy = it->second;
                for (size_t edge_id : ingoingCopy) {
                    double edge_weight = graph.edge_weight.at(edge_id);
                    std::pair<size_t, size_t> edge = graph.edge_ids.at(edge_id);
                    if ((edge.first == (size_t) -1) && (edge.second == (size_t) -1)) continue;
                    auto it2 = resolved_classes.find(edge.first);
                    if ((it2 == resolved_classes.end())) {
                        add_edge(&graph, edge.first, id_new_class, edge_weight);
                    }
                    remove_edge(&graph, edge_id);
                }
            }
        }
    }

    for (const auto &classRef_eqClass : eqClassesResolved) {
        for (size_t node_in_class : classRef_eqClass.second) {
            remove_node(&graph, node_in_class);
        }
    }

    ///std::cout << "backtrack_to = " << backtrack_to << std::endl;
    ///std::cout << "eqClassesResolved = " << eqClassesResolved << std::endl;
}

#if 0

struct backward_dfa_minimization_safe_heuristic {
    bool sorted;
std::string node_label;
std::vector<std::pair<double, std::set<jackbergus::fuzzyStringMatching3::log_trace>>> outgoing;
std::vector<size_t> loop_nodes;

backward_dfa_minimization_safe_heuristic(const std::string &nodeLabel) : node_label(nodeLabel) {
sorted = false;
}

backward_dfa_minimization_safe_heuristic(const backward_dfa_minimization_safe_heuristic& ) = default;
backward_dfa_minimization_safe_heuristic(backward_dfa_minimization_safe_heuristic&& ) = default;
backward_dfa_minimization_safe_heuristic& operator=(const backward_dfa_minimization_safe_heuristic& ) = default;
backward_dfa_minimization_safe_heuristic& operator=(backward_dfa_minimization_safe_heuristic&& ) = default;
void finalize() {
if (!sorted) {
    std::sort(outgoing.begin(), outgoing.end());
    unique_sorted_vector(loop_nodes);
    sorted = true;
}
}

bool operator==(const backward_dfa_minimization_safe_heuristic &rhs) const {
return node_label == rhs.node_label &&
       outgoing == rhs.outgoing &&
        loop_nodes == rhs.loop_nodes;
}

bool operator!=(const backward_dfa_minimization_safe_heuristic &rhs) const {
return !(rhs == *this);
}
};

namespace std {
template <>
struct hash<backward_dfa_minimization_safe_heuristic> {
std::size_t operator()(const backward_dfa_minimization_safe_heuristic& k) const {
    std::hash<std::string> sh;
    std::hash<double> dh;

    size_t hc = hash_combine(17, sh(k.node_label));

    size_t s = 31;
    for (const auto& ref : k.outgoing) {
        size_t s2 = 13;
        for (const auto& trace : ref.second)
            s2 = hash_combine(s2, trace);
        s = hash_combine(s, dh(ref.first) ^ s2);
    }

    size_t s3 = 19;
    for (const auto& ref : k.loop_nodes)
        s3 = hash_combine(s3, ref);

    return hash_combine(hash_combine(hc, s), s3);
}
};
}

void backward_minimization(weigthed_labelled_automata& graph, const determinization_information& info) {
std::unordered_map<size_t, std::set<jackbergus::fuzzyStringMatching3::log_trace>> postfix_map;
std::vector<bool> visited_nodes(graph.V_size, false);
visited_nodes[info.is_final_state_inserted_into_result] = true;
postfix_map[info.is_final_state_inserted_into_result].insert({graph.node_label.at(info.is_final_state_inserted_into_result)});
std::vector<size_t> VertexOrder = approx_longest_path(&graph, info.tg_initial);
for (size_t i = 0, N = VertexOrder.size(); i<N; i++)
std::cout << "vertex #" << i << " @" << VertexOrder.at(i) << std::endl;
std::cout << std::endl;

std::vector<size_t> frontier;
for (size_t edge_id : graph.ingoing_edges.at(info.is_final_state_inserted_into_result)) {
size_t src = graph.edge_ids.at(edge_id).first;
if (!visited_nodes.at(src)) {
    visited_nodes[src] = true;
    frontier.emplace_back(src);
}
}
unique_sorted_vector_cmp(frontier, [&VertexOrder](const size_t u, const size_t v) {
return VertexOrder.at(u) >= VertexOrder.at(v);
});

while (!frontier.empty()) {
std::vector<size_t> Xyields, uVector;

std::cout << "Frontier F = " << frontier << std::endl;

std::unordered_map<backward_dfa_minimization_safe_heuristic, std::vector<size_t>> Map;
for (size_t node_id : frontier) {
    size_t vertex_order_node_id = VertexOrder.at(node_id);
    std::string postfix{graph.node_label.at(node_id)};
    backward_dfa_minimization_safe_heuristic h{postfix};
    for (size_t outgoing_edges : graph.nodes.at(node_id)) {
        size_t target = graph.edge_ids.at(outgoing_edges).second;
        if (VertexOrder.at(target) >= vertex_order_node_id) {
            assert(postfix_map.contains(target));
            for (auto logtrace : postfix_map.at(target)) {
                logtrace.emplace_back(postfix);
                postfix_map[node_id].insert(logtrace);
            }
            h.outgoing.emplace_back(graph.edge_weight.at(outgoing_edges),
                                    postfix_map.at(target));
        } else {
            // Please observe! Keep track of the nodes actually creating a loop. Making that a limitation,
            // so two equivalent nodes must also have equivalent looping counterparts
            h.loop_nodes.emplace_back(target);
        }
    }
    h.finalize();
    Map[h].emplace_back(node_id);
}
for (auto& ref : Map) {
    std::sort( ref.second.begin(), ref.second.end() );
    ref.second.erase( std::unique( ref.second.begin(), ref.second.end() ), ref.second.end() );
    if (ref.second.size() == 1) continue;
    size_t x_signed = *ref.second.begin();
    Xyields.emplace_back(x_signed);

    for (size_t i = 1, N = ref.second.size(); i<N; i++) {
        std::vector<size_t> IN = graph.ingoing_edges.at(ref.second.at(i));
        for (size_t edge_id : IN) {
            uVector.emplace_back(graph.edge_ids.at(edge_id).first);
            update_edge_target(&graph, edge_id, x_signed, nullptr);
        }
    }
}


unique_sorted_vector(uVector);
for (size_t i : uVector) edge_compacting(graph, i);
uVector.clear();

unique_sorted_vector(Xyields);

frontier.clear();
for (size_t node_id : Xyields) {
    for (size_t edge_id : graph.ingoing_edges.at(node_id)) {
        size_t src = graph.edge_ids.at(edge_id).first;
        if (!visited_nodes.at(src)) {
            visited_nodes[src] = true;
            frontier.emplace_back(src);
        }
    }
}

unique_sorted_vector_cmp(frontier, [&VertexOrder](const size_t u, const size_t v) {
    return VertexOrder.at(u) >= VertexOrder.at(v);
});
}

}
#endif


determinization_information
topology1(weigthed_labelled_automata &wla, bool withLoops = true, double stop_probability = 1E-06) {
    size_t n0 = add_node(&wla, "a");
    size_t n1 = add_node(&wla, "b");
    size_t n2 = add_node(&wla, "b");
    size_t n3 = add_node(&wla, "c");
    size_t n4 = add_node(&wla, "a");
    size_t n5 = add_node(&wla, "c");
    size_t n6 = add_node(&wla, "a");
    size_t n7 = add_node(&wla, "f");

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
    return {stop_probability, n0, n7};
}

determinization_information topology0(weigthed_labelled_automata &wla, double stop_probability = 1E-06) {
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
    return {stop_probability, n0, n4};
}

determinization_information
topology2(weigthed_labelled_automata &wla, bool add_limiting_nodes = true, double stop_probability = 1E-06) {
    size_t n0 = add_node(&wla, ".");
    size_t n1 = add_node(&wla, "a");
    size_t n2 = add_node(&wla, "a");
    size_t n3 = add_node(&wla, "b");
    size_t n4 = add_node(&wla, "c");
    size_t n5 = add_node(&wla, "a");
    size_t n6 = add_node(&wla, "a");
    size_t n7 = add_node(&wla, "b");
    size_t n8 = add_node(&wla, "b");
    size_t n11 = add_node(&wla, "e");
    size_t n13 = add_node(&wla, "a");

    add_edge(&wla, n0, n1, 0.7);
    add_edge(&wla, n1, n3, 1.0);
    add_edge(&wla, n3, n5, 1.0);
    add_edge(&wla, n5, n7, 1.0);
    add_edge(&wla, n7, n13, 0.5);
    add_edge(&wla, n13, n7, 1.0);

    add_edge(&wla, n0, n2, 0.3);
    add_edge(&wla, n2, n4, 1.0);
    add_edge(&wla, n4, n6, 1.0);
    add_edge(&wla, n6, n8, 1.0);
    add_edge(&wla, n8, n6, 0.5);

    if (add_limiting_nodes) {
        size_t n9 = add_node(&wla, "c");
        size_t n10 = add_node(&wla, "d");
        add_edge(&wla, n7, n9, 0.5);
        add_edge(&wla, n9, n11, 1.0);
        add_edge(&wla, n8, n10, 0.5);
        add_edge(&wla, n10, n11, 1.0);
    } else {
        add_edge(&wla, n7, n11, 0.5);
        add_edge(&wla, n8, n11, 0.5);
    }

    return {stop_probability, n0, n11};
}

void generate_my_minimization(double init_prob = 1.0, double stop_probability = 1E-06) {
    weigthed_labelled_automata wla, out;
    std::unordered_map<determinization_node_eq_safe_heuristic, size_t> node_compact_info;

    //determinization_information info = topology0(wla, stop_probability);
    //determinization_information info = topology1(wla, true, stop_probability);
    //determinization_information info = topology1(wla, false, stop_probability);
    //determinization_information info = topology2(wla, true, stop_probability);
    determinization_information info = topology2(wla, false, stop_probability);

    probabilisitc_model_trace pmt;
    pmt.t_with_prob.probability = init_prob;
    pmt.t_with_prob.t.emplace_back(node_label(&wla, info.tg_initial));
    pmt.underlying_sequences.emplace_back(init_prob, info.tg_initial);
    //std::cout << pmt << std::endl;

    // Ensuring that all the edges have 1.0 total sum for their outgoing probability
    // O(|V+E|)
    edge_normalize(wla);

    // Preliminary minimizing the visit and the complexity by compacting the edge traversing: sum up the edges' probability
    // leading to the same outgoing state
    // O(|V+E|)
    edge_compacting(wla);

    // Making the NFA into a DFA, by exploiting the prefix strategy
    // O()
    nfa_to_dfa_weighted_labelled_automata(wla, pmt, out, info, node_compact_info);

    // While creating the wheel state, I might have multiple outgoing edges leading to the same state. Therefore, I
    // could compact those into one single edge by summing up the weights
    // O(|V+E|)
    edge_compacting(out);

    // Providing the backward minimization for DFAs, where edge cost is also considered.
    // This is a sufficient condition, as we just have one single final node.
    heuristic_scala(out, info.tg_initial, info.is_final_state_inserted_into_result);

    dot(&out, std::cout);
}

/*
void automata_minimization() {
    test_automata_minimization(4, 0, 3, {{0, 1, 0.5, "a"},
                                         {0, 2, 0.5, "b"},
                                         {1, 1, 0.5, "c"},
                                         {2, 2, 0.5, "c"},
                                         {1, 3, 0.5, "d"},
                                         {2, 3, 0.5, "d"}});
}
*/

int main() {

    generate_my_minimization();
    //automata_minimization();
}
