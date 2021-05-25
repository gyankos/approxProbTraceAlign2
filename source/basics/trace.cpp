/*
 * trace.cpp
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

#include <basics/trace.h>

// This cpp is for the only purpose of implicitly implementing the operator<< for printing some of the structs
bool jackbergus::fuzzyStringMatching3::probabilisitc_model_trace::test() {
    double sum = 0.0;
    for (const auto& x : underlying_sequences)
        sum += x.probability;
    return (sum - t_with_prob.probability) <= jackbergus::fuzzyStringMatching3::eps;
}
