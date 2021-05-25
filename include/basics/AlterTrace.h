/*
 * AlterString.h
 * This file is part of ProbabilisticTraceAlignment
 *
 * Copyright (C) 2020 - Giacomo Bergami
 *
 * ProbabilisticTraceAlignment is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * ProbabilisticTraceAlignment is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ProbabilisticTraceAlignment. If not, see <http://www.gnu.org/licenses/>.
 */


//
// Created by giacomo on 19/10/20.
//

#ifndef FUZZYSTRINGMATCHING2_ALTERSTRING_H
#define FUZZYSTRINGMATCHING2_ALTERSTRING_H



#include <random>
#include <basics/trace.h>
namespace jackbergus {
    namespace fuzzyStringMatching3 {
        namespace algorithms {

/**
 * Alters a string provided as an input
 */
            struct AlterTrace {
                std::uniform_real_distribution<double>  errorType;
                std::uniform_int_distribution<size_t>  valueType;
                std::mt19937 mersenneError, mersenneValue, mersennePosition;
                std::uniform_int_distribution<size_t> pick;
                jackbergus::fuzzyStringMatching3::log_trace chrs;
                double noiseThreshold;
                size_t seedError = 0;

                AlterTrace(const AlterTrace& ) = default;
                AlterTrace& operator=(const AlterTrace& ) = default;

                /**
                 * Initialization
                 * @param chrs                  Strings to insert
                 * @param noiseThreshold        Threshold under which perform the changes
                 * @param seedError             Seed for generating the error distribution (no change, insertion, deletion, swap)
                 */
                AlterTrace(const jackbergus::fuzzyStringMatching3::log_trace& chrs, double noiseThreshold, size_t seedError = 0);

                /**
                 * Alters a string
                 * @param toAlter   The string to alter, passed by copy
                 * @return          The altered string
                 */
                jackbergus::fuzzyStringMatching3::log_trace alter(const jackbergus::fuzzyStringMatching3::log_trace& toAlter);

            private:
                int doNoise() {
                    double randomNumber =  errorType(mersenneError);
                    return randomNumber <= noiseThreshold ? std::min(3, (int) std::trunc(randomNumber / (noiseThreshold) * 3.0) + 1)
                                                          : 0;
                }

            };


        }
    }
}


#endif //FUZZYSTRINGMATCHING2_ALTERSTRING_H
