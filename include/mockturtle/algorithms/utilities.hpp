#pragma once
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <tuple>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <array>
#include <fmt/format.h>
#include <thread>
#include "mockturtle/algorithms/nodes.hpp"
#include <algorithm>
#include <locale>

#include <iomanip>
#include <ctime>


#define continue_if(condition) if (condition) continue
#define break_if(condition) if (condition) break
#define shortcut_cost(init_cost) if (init_cost > limit) return std::make_tuple(INF32, false)

#define COMPUTE_IJ(k, N, i, j) \
    uint64_t i = N - 2 - static_cast<uint64_t>(std::sqrt(4*N*(N - 1) - 8*k - 7) / 2.0 - 0.5); \
    uint64_t j = k + i + 1 - N * (N - 1) / 2 + (N - i) * ((N - i) - 1) / 2;

inline uint64_t ceil(uint64_t x, uint64_t y)
{
    return 1 + ((x - 1) / y);
}

/* to be used in the future 
    inline void unique_merge(std::vector<uint64_t>& a, std::vector<uint64_t>& b, std::vector<uint64_t>& c)
    {
        std::vector<uint64_t> c(a.size() + b.size());
        auto it = std::set_union(a.begin(), a.end(), b.begin(), b.end(), c.begin());
        c.resize(it - c.begin());
    }

    inline void unique_insert(std::vector<uint64_t>& v, const uint64_t& element) {
        if (std::find(v.begin(), v.end(), element) == v.end()) {
            v.push_back(element);
        }
    }

    bool is_subset(const std::vector<uint64_t>& v1, const std::vector<uint64_t>& v2) {
        // Check if all elements of v1 are in v2
        return std::all_of(v1.begin(), v1.end(), [&v2](const int& val) {
            return std::find(v2.begin(), v2.end(), val) != v2.end();
        });
    }

    std::vector<uint64_t> get_pi_recursive(uint64_t hash, std::unordered_map<uint64_t, Node>& nodemap, std::unordered_set<uint64_t>& all_hashes)
    {
        auto it = NODE_PIS.find(hash);
        if (it == NODE_PIS.end())
        {
            Node & n = nodemap[hash]; 
            if (n.last_func == fPI) // the node is itself PI, insert hash. No continuation is needed
            {
                std::vector<uint64_t> pis = { hash };
                NODE_PIS[hash] = pis;
                PI_HASHES.emplace(hash);
                all_hashes.erase(hash);
                return pis;
            }
            else if (n.last_func == fNOFUNC)
            {
                // Don't throw, since this may interrupt large scale optimization
                fmt::print("\tWARNING: encountered invalid node.");
                all_hashes.erase(hash);
                return {};
            }
            else // in all other cases, just call the function once again
            {
                std::vector<uint64_t> pis;
                for (uint64_t phash : n.parent_hashes)
                {
                    pis.emplace_back(get_pi_recursive(phash, nodemap, all_hashes));
                }
                NODE_PIS[hash] = pis;
                all_hashes.erase(hash);
                return pis;
            }
        }
        all_hashes.erase(hash);
        return it->second;
    }

    std::unordered_set<uint64_t> get_hashes(const std::unordered_map<uint64_t, Node>& hashmap) 
    {
        std::unordered_set<uint64_t> hashes;
        for (const auto& [hash, value] : hashmap) 
        {
            Node & n = GNM_global[hash];
            if (n.last_func != fNOFUNC && NODE_PIS.find(hash) == NODE_PIS.end())
            {
                hashes.emplace(hash);
            }
        }
        return hashes;
    }

    void update_NODE_PIS()
    {
        std::unordered_set<uint64_t> all_hashes = get_hashes(GNM_global);
        while (!all_hashes.empty())
        {
            get_pi_recursive(*all_hashes.begin(), GNM_global, all_hashes);
        }
    }

    void get_nodemap(std::array<uint8_t, NUM_VARS>& levels, std::unordered_map<uint64_t, Node>& hashmap, std::array<uint64_t, NUM_TT>& local_GEA, std::array<uint64_t, NUM_TT>& local_GEX) 
    // TODO: also UPDATE GEX and GEA
    {
        update_NODE_PIS();
        std::vector<uint64_t> valid_pi_hashes;
        for (const uint64_t pi_hash : PI_HASHES)
        {
            const Node & pi = GNM_global[pi_hash];
            uint8_t pi_idx = INF8;
            switch (pi.func._bits)
            {
                case 0x5555: pi_idx = 0;
                case 0x3333: pi_idx = 1;
                case 0x0F0F: pi_idx = 2;
                case 0x00FF: pi_idx = 3;
            }
            if (levels[pi_idx] == pi.lvl)
            {
                valid_pi_hashes.push_back(pi.hash);
            }
        }
        
        for (auto & [hash, pis] : NODE_PIS)
        {
            if ( std::includes(valid_pi_hashes.begin(), valid_pi_hashes.end(), pis.begin(), pis.end()) )
            {
                hashmap[hash] = GNM_global[hash];
            }
            else 
            {
                hashmap.erase(hash);
            }
        }
    }
*/