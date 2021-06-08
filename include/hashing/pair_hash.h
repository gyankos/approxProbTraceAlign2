/*
 * pair_hash.h
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

#ifndef FUZZYSTRINGMATCHING2_PAIR_HASH_H
#define FUZZYSTRINGMATCHING2_PAIR_HASH_H

#include <utility>
#include <functional>

namespace std {
    template <class T1, class T2>
    struct hash<std::pair<T1, T2>>
    {
        std::size_t operator()(const std::pair<T1, T2> &pair) const
        {
            return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
        }
    };

}

template < class K, class V > std::ostream& operator << (std::ostream& os, const std::pair<K, V>& v) {
    return os << "⟪" << v.first << ", " << v.second << "⟫";
}

#endif //FUZZYSTRINGMATCHING2_PAIR_HASH_H
