/*
 * trace.h
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

#ifndef COMPLETE_MATCHING_TRACE_H
#define COMPLETE_MATCHING_TRACE_H

#include <vector>
#include <string>
#include <ostream>

#include <hashing/vector_hash.h>

namespace jackbergus {
    namespace fuzzyStringMatching3 {

        typedef std::vector<std::string>        log_trace;
        typedef std::vector<unsigned long>      valid_sequence;

        constexpr const double eps = std::numeric_limits<double>::epsilon(); ///< machine epsilon used to ignore some irrelevant values
        constexpr char tau[] = "";                                           ///< empty-skippable event (1)
        constexpr char varepsilon[] = "";                                    ///< empty-skippable event (2)

        /**
         * Definition of a probabilistic log trace, that therefore has no underlying valid sequences.
         * Still, this trace has an associated probability. The default probability value is 1.
         */
        struct probabilistic_log_trace {
            double probability; ///< Trace Probability
            log_trace t;        ///< Actual events outlining the trace

            template<typename... ArgTypes> probabilistic_log_trace(double probability, ArgTypes... Args) :
            probability{probability}, t{std::forward<ArgTypes>(Args)...} {}
            template<typename... ArgTypes> probabilistic_log_trace(ArgTypes... Args) :
            probability{1.0}, t{std::forward<ArgTypes>(Args)...} {}
            probabilistic_log_trace() : probability{1.0}, t{} {}
            probabilistic_log_trace(const probabilistic_log_trace& ) = default;
            probabilistic_log_trace(probabilistic_log_trace&& ) = default;
            probabilistic_log_trace& operator=(const probabilistic_log_trace& ) = default;
            probabilistic_log_trace& operator=(probabilistic_log_trace&& ) = default;

            friend std::ostream &operator<<(std::ostream &os, const probabilistic_log_trace &trace) {
                return os << "ℙ_ℒ(" << trace.t << ") = " << trace.probability;
            }
        };

        /**
         * Definition of a probabilistic valid sequence belonging to a probabilistic model.
         * The model of reference is implicit, and so no explicit conversion from node id and trace label is given for
         * space efficiency reasons.
         * The default probability value is 1.
         */
        struct probabilistic_valid_sequence {
            double probability;       ///< Trace Probability
            valid_sequence sequence;  ///< Actual node succession outlining the valid sequence

            template<typename... ArgTypes> probabilistic_valid_sequence(double probability, ArgTypes... Args) :
                    probability{probability}, sequence{std::forward<ArgTypes>(Args)...} {}
            template<typename... ArgTypes> probabilistic_valid_sequence(ArgTypes... Args) :
                    probability{1.0}, sequence{std::forward<ArgTypes>(Args)...} {}
            probabilistic_valid_sequence() : probability{1.0}, sequence{} {}
            probabilistic_valid_sequence(const probabilistic_valid_sequence& ) = default;
            probabilistic_valid_sequence(probabilistic_valid_sequence&& ) = default;
            probabilistic_valid_sequence& operator=(const probabilistic_valid_sequence& ) = default;
            probabilistic_valid_sequence& operator=(probabilistic_valid_sequence&& ) = default;

            friend std::ostream &operator<<(std::ostream &os, const probabilistic_valid_sequence &trace) {
                return os << "ℙ_N(" << trace.sequence << ") = " << trace.probability;
            }

            double extend_with_subsequent_step(const std::vector<std::pair<double, size_t>>& ref, std::vector<probabilistic_valid_sequence>& V) const {
                size_t N = ref.size();
                double sum = 0;
                for (size_t i = 0; i<N; i++) {
                    auto& refV = V.emplace_back(*this);
                    sum += (refV.probability *= ref.at(i).first);
                    refV.sequence.emplace_back(ref.at(i).second);
                }
                return sum;
            }
        };

        /**
         * Definition of a probabilistic model trace.
         * The model of reference is implicit, and so no explicit conversion from node id and trace label is given for
         * space efficiency reasons. Only the resulting trace after empty-string stripping is given.
         * The default probability value is 1.
         */
        struct probabilisitc_model_trace {
            probabilistic_log_trace  t_with_prob;
            std::vector<struct probabilistic_valid_sequence> underlying_sequences;

            probabilisitc_model_trace(double prob) : t_with_prob(prob) {}
            probabilisitc_model_trace() {}
            probabilisitc_model_trace(const probabilisitc_model_trace& ) = default;
            probabilisitc_model_trace(probabilisitc_model_trace&& ) = default;
            probabilisitc_model_trace& operator=(const probabilisitc_model_trace& ) = default;
            probabilisitc_model_trace& operator=(probabilisitc_model_trace&& ) = default;

            template<typename... ArgTypes> probabilisitc_model_trace(ArgTypes... Args) :
                    t_with_prob{1.0, std::forward<ArgTypes>(Args)...} {}
            template<typename... ArgTypes> probabilisitc_model_trace(double probability, ArgTypes... Args) :
                    t_with_prob{probability, std::forward<ArgTypes>(Args)...} {}
            friend std::ostream &operator<<(std::ostream &os, const probabilisitc_model_trace &trace) {
                os << trace.t_with_prob << std::endl << "having valid sequences: {";
                if (!trace.underlying_sequences.empty()) {
                    os << std::endl;
                    for (const auto& elem : trace.underlying_sequences)
                        os << elem << std::endl;
                }
                return os << "}";
            }

            /**
             * This function tests whether the valid sequences' probability sum up to the trace probability
             * as expected
             * @return Whether the test was passed or not
             */
            bool test();
        };
    }
}

#endif //COMPLETE_MATCHING_TRACE_H
