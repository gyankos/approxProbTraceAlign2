/*
 * AlterString.cpp
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

#include <string>
#include <cassert>
#include "basics/AlterTrace.h"

using namespace jackbergus::fuzzyStringMatching3::algorithms;

AlterTrace::AlterTrace(const jackbergus::fuzzyStringMatching3::log_trace &chrs, double noiseThreshold, size_t seedError) :
        errorType{0.0,1.0}, pick(0, chrs.size()-2), chrs{chrs},
        noiseThreshold{noiseThreshold}, seedError{seedError} {
    mersenneError.seed(seedError);
}

#include <unordered_map>

jackbergus::fuzzyStringMatching3::log_trace AlterTrace::alter(const jackbergus::fuzzyStringMatching3::log_trace &toAlter) {
    auto it = toAlter.begin();
    size_t originalSize = toAlter.size();
    size_t i = 0;
    jackbergus::fuzzyStringMatching3::log_trace result;
    while ((!toAlter.empty()) && (i< originalSize) &&  it != toAlter.end() ) {
        int casus = doNoise();
        if (toAlter.empty() && (casus != 1)) casus = 1;
        if (casus == 1) {
            // Add character
            std::string c = (std::string)chrs[pick(mersenneValue)];
            result.emplace_back(chrs.at(pick(mersenneValue)));
            result.emplace_back(*it);
            originalSize++;
        } if (casus == 2) {
            // Remove character
            originalSize--;
        } if (casus == 3) {
            // Replace character
            result.emplace_back(*it);
            auto it2 = result.begin()+(result.size()-1);
            size_t d = (size_t)std::distance(result.begin(), it2);
            std::uniform_int_distribution<size_t> randomPosition{0UL, d};
            size_t j = randomPosition(mersennePosition);
            std::iter_swap(it2, result.begin()+j);
            it++;
        } else {
            result.emplace_back(*it);
            it++;
        }
        i++;
    }
    return toAlter;
}
