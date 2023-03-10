#pragma region imports
#include <iostream>
#include <fstream> //include the filestreamobject as the header files
#include <tuple>
#include <iterator>
#include <array>
#include <algorithm>
#include <vector>
#include <set>
// #include <unordered_map>
#include <unordered_set>
#include <bitset>
// #include <stdio.h>
#include <cstdio>
#pragma endregion imports

#pragma region constants

typedef unsigned short US;
typedef unsigned int   UI;
typedef unsigned long  UL;
typedef const unsigned short CUS;
typedef const unsigned int   CUI;
typedef const unsigned long  CUL;

CUS NVARS = 4u;
CUI BITLEN = 1 << NVARS;
CUL NTT = 1 << BITLEN;
CUI ONES = NTT - 1;
CUI INF = ONES;
const std::array<UI, 2> INF2 = { INF, INF };
const std::array<UI, 3> INF3 = { INF, INF, INF };

typedef std::array<US, NTT> US_NTT;
typedef std::array<UI, NTT> UI_NTT;
typedef std::array<UL, NTT> UL_NTT;
typedef std::array<std::array<UI, 2>, NTT> UI_NTT2;
typedef std::array<std::array<UI, 3>, NTT> UI_NTT3;
typedef std::set<UI> UI_set;
typedef const std::array<US, NTT> CUS_NTT;
typedef const std::array<UI, NTT> CUI_NTT;
typedef const std::array<UL, NTT> CUL_NTT;
typedef const std::array<std::array<UI, 2>, NTT> CUI_NTT2;
typedef const std::array<std::array<UI, 3>, NTT> CUI_NTT3;
typedef const std::set<UI> CUI_set;

CUI kMERGE   = 0;
CUI kMERGE3  = 1;
CUI kXMERGE  = 2;
CUI kNOT     = 3;
CUI kDFF     = 4;
CUI kSPL     = 5;

CUS COSTS[6] = {8, 11, 7, 9, 7, 7};

CUS fDFF   = 0;
CUS fNOT   = 1;
CUS fMERGE = 2;
CUS fOR    = 3;
CUS fAND   = 4;
CUS fXOR   = 5;
CUS fOR3   = 6;
CUS fAND3  = 7;
CUS fMAJ3  = 8;

// CUI PI[4] = {255,  3855, 13107, 21845};
// CUI PI[6] = {0,  255,  3855, 13107, 21845, 65535};
CUI PI[6] = {0, 43690, 52428, 61680, 65280, 65535}; // Reversed - MSB corresponds to 1111, LSB corresponds to 0000
// CUI PI[10] = {0,  255,  3855, 13107, 21845, 43690, 52428, 61680, 65280, 65535};
CUI nPI = sizeof(PI) / sizeof(int);
#pragma endregion constants

#pragma region globals
UI_set frontier(std::begin(PI), std::end(PI));
UI_set xorable_frontier(std::begin(PI), std::end(PI));
UI ci, cj, ck, new_cost, new_depth, OR, AND;
// std::array< std::unordered_map<UI, US> , NTT> spl_map;


US NSPL = 0u;
#pragma endregion globals

// TODO: change costs and depths to signed integers for easier difference calculations

template <typename T>
UI uabs(const T& a, const T& b)
{
    return a > b ? (a - b) : (b - a);
}

UI count_dff(CUI_NTT& depths, CUI tt_a, CUI tt_b = INF, CUI tt_c = INF)
{
    
    if (tt_c == INF)
    {
        if (tt_b == INF)
        {
            return 0u;
        }
        return (depths[tt_b] >= depths[tt_a]) ? (depths[tt_b] - depths[tt_a]) : ((depths[tt_a] - depths[tt_b]));
    }
    else
    {
        if (depths[tt_a] >= depths[tt_b]) // (descending) ABC, ACB, CAB 
        {
            if (depths[tt_a] >= depths[tt_c]) // ABC, ACB - A is the largest
            {
                return 2*depths[tt_a] - depths[tt_b] - depths[tt_c];
            }
            else // CAB - C is the largest
            {
                return 2*depths[tt_c] - depths[tt_b] - depths[tt_a];
            }
        }
        else // BAC, BCA, CBA
        {
            if (depths[tt_b] >= depths[tt_c]) // BAC, BCA - B is the largest
            {
                return 2*depths[tt_b] - depths[tt_c] - depths[tt_c];
            }
            else // CBA - C is the largest
            {
                return 2*depths[tt_c] - depths[tt_b] - depths[tt_a];
            }
        }
    }
}

// // succ - array of unordered maps. For each map, the key is the child TT, the value is the number of times the child is keyed by the parent 
// succ - array of unordered sets. For each parent map, contains the child tt-s preceeded by the parent

bool update(UI_NTT& costs, UI_NTT3& pred, std::array<std::unordered_set<UI>, NTT>& succ, UI_NTT& func, UI_NTT& depths, CUI tt, CUI cost, CUS depth, CUS function, CUI ci, CUI cj = INF, CUI ck = INF, const bool verbose = false)
{
    // frontier.insert(tt);
    if (tt == ci || tt == cj || tt == ck) return false;

    bool status = (depth == depths[tt]) & (cost < costs[tt]);

    if (status)
    {
        if (verbose)
        {
            std::cout << "UPDATE HAPPENS " << cost << " < " << costs[tt] << std::endl;
        }
        auto cost_delta = costs[tt] - cost;
        costs[tt] = cost;
        depths[tt] = depth;
        // DONE: update successors down the line 
        // remove tt from old predecessors
        succ[pred[tt][0]].erase(tt);
        succ[pred[tt][1]].erase(tt);
        succ[pred[tt][2]].erase(tt);
        // update predecessors
        pred[tt][0] = ci;
        pred[tt][1] = cj;
        pred[tt][2] = ck;
        // update successors
        succ[ci].insert(tt);
        succ[cj].insert(tt);
        succ[ck].insert(tt);
        func[tt] = function;

        // Updates the costs of the subsequent nodes
        std::vector<UI> stack( succ[tt].begin(), succ[tt].end() );
        while (!stack.empty())
        {   
            UI child_tt = stack.back();
            costs[child_tt] -= cost_delta;
            stack.insert(stack.end(), pred[child_tt].begin(), pred[child_tt].end());
            stack.pop_back();
        }

        if (verbose)
        {
            std::cout << "\t\tUPDATING: " << tt << ' ' << cost << std::endl;
            std::cout << "\t\tNew pred: " << ci << ' ' << cj << ' ' << std::endl;
            std::cout << "\t\tNew func: " << function << ' ' << std::endl;
        }
    }
    else if (verbose)
    {
        std::cout << "UPDATE DOES NOT HAPPEN " << cost << " >= " << costs[tt] << std::endl;
    }
    return status;
}

void combine2(UI_NTT& costs, UI_NTT3& pred, UI_NTT& func, UI_NTT& depths, UI_NTT& xcosts, UI_NTT3& xpred, UI_NTT& xfunc, UI_NTT& xdepths, CUI ci, CUI cj, const bool verbose = false)
{
    new_depth = std::max(depths[ci], depths[cj]); //depth is incremented at the final stage
    assert(costs[ci] < INF && costs[cj] < INF);
    new_cost = costs[ci] + costs[cj] + COSTS[kMERGE] + COSTS[kDFF] * count_dff(depths, ci, cj);

    //ALL TT-s here are xorable, since they produce a single pulse
    if (verbose) 
    {
        std::cout << std::bitset<16>(ci) << " OR\n" << std::bitset<16>(cj) << "\n----------------\n" << std::bitset<16>(ci | cj) << std::endl;
        std::cout << "\tOriginal cost: " << costs[ci | cj] << std::endl;
        std::cout << "\tNew cost: " << new_cost << std::endl;
    }
    update( costs,  pred,  func,  depths, ci | cj, new_cost, new_depth, fOR, ci, cj, INF, verbose);
    update(xcosts, xpred, xfunc, xdepths, ci | cj, new_cost, new_depth, fOR, ci, cj, INF, verbose);

    if (verbose) 
    {
        std::cout << std::bitset<16>(ci) << " AND\n" << std::bitset<16>(cj) << "\n________________\n" << std::bitset<16>(ci & cj) << std::endl;
        std::cout << "\tOriginal cost: " << costs[ci & cj] << std::endl;
        std::cout << "\tNew cost: " << new_cost << std::endl;
    }
    update( costs,  pred,  func,  depths, ci & cj, new_cost, new_depth, fAND, ci, cj, INF, verbose);
    update(xcosts, xpred, xfunc, xdepths, ci & cj, new_cost, new_depth, fAND, ci, cj, INF, verbose);
}

void combine3(UI_NTT& costs, UI_NTT3& pred, UI_NTT& func, UI_NTT& depths, UI_NTT& xcosts, UI_NTT3& xpred, UI_NTT& xfunc, UI_NTT& xdepths, UI ci, UI cj, UI ck, const bool verbose = false)
{
    // Sort tt-s by depth. ci is now the shallowest, ck is the deepest
    if (depths[ci] > depths[cj]) std::swap(ci, cj);
    if (depths[ci] > depths[ck]) std::swap(ci, ck);
    if (depths[cj] > depths[ck]) std::swap(cj, ck);

    new_depth = depths[ck];//depth is incremented at the final stage
    assert(costs[ci] < INF && costs[cj] < INF && costs[ck] < INF);
    new_cost = costs[ci] + costs[cj] + costs[ck] + COSTS[kMERGE3] + COSTS[kDFF] * count_dff(depths, ci, cj, ck);
    //TODO : NEED TO ADD XORABLES
    if (verbose) 
    {
        std::cout   << std::bitset<16>(ci) << " OR\n"
                    << std::bitset<16>(cj) << " OR\n"
                    << std::bitset<16>(ck) << "\n----------------\n"
                    << std::bitset<16>(ci | cj | ck) << std::endl;
        std::cout << "\tOriginal cost: " << costs[ci | cj | ck] << std::endl;
        std::cout << "\tNew cost: " << new_cost << std::endl;
    }
    update( costs,  pred,  func,  depths, ci | cj | ck, new_cost, new_depth, fOR3, ci, cj, ck, verbose);
    update(xcosts, xpred, xfunc, xdepths, ci | cj | ck, new_cost, new_depth, fOR3, ci, cj, ck, verbose);

    if (verbose) 
    {
        std::cout   << std::bitset<16>(ci) << " AND\n"
                    << std::bitset<16>(cj) << " AND\n"
                    << std::bitset<16>(ck) << "\n----------------\n"
                    << std::bitset<16>(ci & cj & ck) << std::endl;
        std::cout << "\tOriginal cost: " << costs[ci & cj & ck] << std::endl;
        std::cout << "\tNew cost: " << new_cost << std::endl;
    }
    update( costs,  pred,  func,  depths, ci & cj & ck, new_cost, new_depth, fAND3, ci, cj, ck, verbose);
    update(xcosts, xpred, xfunc, xdepths, ci & cj & ck, new_cost, new_depth, fAND3, ci, cj, ck, verbose);

    if (verbose) 
    {
        std::cout   << std::bitset<16>(ci) << " MAJ\n"
                    << std::bitset<16>(cj) << " MAJ\n"
                    << std::bitset<16>(ck) << "\n----------------\n"
                    << std::bitset<16>((ci & cj) | (ci & ck) | (cj & ck)) << std::endl;
        std::cout << "\tOriginal cost: " << costs[(ci & cj) | (ci & ck) | (cj & ck)] << std::endl;
        std::cout << "\tNew cost: " << new_cost << std::endl;
    }
    update( costs,  pred,  func,  depths, (ci & cj) | (ci & ck) | (cj & ck), new_cost, new_depth, fMAJ3, ci, cj, ck, verbose);
    update(xcosts, xpred, xfunc, xdepths, (ci & cj) | (ci & ck) | (cj & ck), new_cost, new_depth, fMAJ3, ci, cj, ck, verbose);
}

UI count_spl(CUI_NTT3& pred, CUI tt_a, CUI tt_b = INF, CUI tt_c = INF)
{

    UI nspl = 0u;
    UI tt;
    std::vector<UI> done(10);
    std::vector<UI> tt_stack(10);
    tt_stack.push_back(tt_a);
    if (tt_b != INF)
    {
        tt_stack.push_back(tt_b);
        if (tt_c != INF)
        {
            tt_stack.push_back(tt_c);
        }
    }
    while (!tt_stack.empty())
    {
        tt = tt_stack.back();
        tt_stack.pop_back();

        if (std::find(done.begin(), done.end(), tt) != done.end()) // if the tt has not been considered yet
        {
            done.push_back(tt);
            for (auto & tt_pred : pred[tt]) //add predecessors to the stack
            {
                if (tt_pred < INF) tt_stack.push_back(tt_pred);
            }
        }
        else
        {
            // std::cout << "Found_splitter at " << tt << " during " << tt_a  << " " << tt_b  << " " << tt_c << std::endl;
            nspl++;
        }
    }
    return nspl;
}

bool merge(UI_NTT& costs, UI_NTT3& pred, UI_NTT& func, UI_NTT& depths, UI_NTT& xcosts, UI_NTT3& xpred, UI_NTT& xfunc, UI_NTT& xdepths, CUI ci, CUI cj, const bool verbose = false)
// TODO: CHECK IF XORABLE IS NEEDED DURING COST CALCULATION
{
    new_depth = std::max(depths[ci], depths[cj]); // depth is incremented at the final stage
    assert(costs[ci] < INF && costs[cj] < INF);
    new_cost = costs[ci] + costs[cj] + COSTS[kMERGE] + COSTS[kDFF] * count_dff(depths, ci, cj); // + COSTS[kSPL] * count_spl(pred, ci, cj)
    UI tt = (ci | cj);
    if (verbose)
    {
        std::cout << std::bitset<16>(ci) << " OR\n" << std::bitset<16>(cj) << "\n----------------\n" << std::bitset<16>(tt) << std::endl;
        std::cout << "\tOriginal cost: " << costs[tt] << std::endl;
        std::cout << "\tNew cost: " << new_cost << std::endl;
    }
    bool status = update(costs, pred, func, depths, tt, new_cost, new_depth, fMERGE, ci, cj, INF, verbose);
    if (status && (verbose))
    {
        std::cout << "Discovered " << std::bitset<16>(ci) << '|' << std::bitset<16>(cj) << '=' << std::bitset<16>(tt) << std::endl;
        std::cout << "Discovered " << ci << '|' << cj << '=' << tt << std::endl;
    }
    if (tt == (ci ^ cj))
    {
        status |= update(xcosts, xpred, xfunc, xdepths, tt, new_cost, new_depth, fMERGE, ci, cj, INF, verbose);
    }
    return status;
}

void print_frontier(CUI_NTT& costs, CUI_NTT3& pred, CUI_NTT& func, CUI_NTT& depths)
{
    UI i = 0;
    for (auto & tt : frontier)
    {
        std::cout << ++i << ": " << std::bitset<16>(tt) << " = " << tt << " (cost= " << costs[tt] << ")" <<  std::endl;
        if (std::find(std::begin(PI), std::end(PI), tt) != std::end(PI))
        {
            std::cout << "\tPrimary input" << std::endl;
        }
        else {
            switch (func[tt]) {
                case fNOT:
                    std::cout << "\tInversion of \n" << std::bitset<16>(pred[tt][0]) << " (cost= " << costs[pred[tt][0]] << ")" << std::endl;
                    break;
                case fMERGE:
                    std::cout << "\tMerging of \n" << std::bitset<16>(pred[tt][0]) << " (cost= " << costs[pred[tt][0]] << ") \n"
                                                   << std::bitset<16>(pred[tt][1]) << " (cost= " << costs[pred[tt][1]] << ")" << std::endl;
                    break;
                case fOR:
                    std::cout << "\tOR-ing of \n" << std::bitset<16>(pred[tt][0]) << " (cost= " << costs[pred[tt][0]] << ") \n"
                                                  << std::bitset<16>(pred[tt][1]) << " (cost= " << costs[pred[tt][1]] << ")" << std::endl;
                    break;
                case fAND:
                    std::cout << "\tAND-ing of \n" << std::bitset<16>(pred[tt][0]) << " (cost= " << costs[pred[tt][0]] << ") \n"
                                                   << std::bitset<16>(pred[tt][1]) << " (cost= " << costs[pred[tt][1]] << ")" << std::endl;
                    break;
                case fXOR:
                    std::cout << "\tXOR-ing of \n" << std::bitset<16>(pred[tt][0]) << " (cost= " << costs[pred[tt][0]] << ") \n"
                                                   << std::bitset<16>(pred[tt][1]) << " (cost= " << costs[pred[tt][1]] << ")" << std::endl;
                    break;
                case fOR3:
                    std::cout << "\tOR-ing of \n"   << std::bitset<16>(pred[tt][0]) << " (cost= " << costs[pred[tt][0]] << ") \n"
                                                    << std::bitset<16>(pred[tt][1]) << " (cost= " << costs[pred[tt][1]] << ") \n"
                                                    << std::bitset<16>(pred[tt][2]) << " (cost= " << costs[pred[tt][2]] << ")" << std::endl;
                    break;
                case fAND3:
                    std::cout << "\tAND-ing of \n"  << std::bitset<16>(pred[tt][0]) << " (cost= " << costs[pred[tt][0]] << ") \n"
                                                    << std::bitset<16>(pred[tt][1]) << " (cost= " << costs[pred[tt][1]] << ") \n"
                                                    << std::bitset<16>(pred[tt][2]) << " (cost= " << costs[pred[tt][2]] << ")" << std::endl;
                    break;
                case fMAJ3:
                    std::cout << "\tMajority of \n" << std::bitset<16>(pred[tt][0]) << " (cost= " << costs[pred[tt][0]] << ") \n"
                                                    << std::bitset<16>(pred[tt][1]) << " (cost= " << costs[pred[tt][1]] << ") \n"
                                                    << std::bitset<16>(pred[tt][2]) << " (cost= " << costs[pred[tt][2]] << ")" << std::endl;
                    break;
            }
        }
    }
}

bool update_frontier(UI_set& fr, CUI_NTT& costs, const bool verbose)
{
    UI ct = 0;
    UI c;
    bool status = false;
    bool temp_status = false;
    for (auto i = 0u; i < NTT; ++i)
    {
        c = costs[i];
        if (c < INF)
        {
            temp_status = fr.insert(i).second;
            status |= temp_status;
            if(temp_status && verbose)
            {
                std::cout << ++ct << ": Found " << i << " with cost " << c << std::endl;
            }
        }
    }
    if (verbose) std::cout << " FRONTIER SIZE IS " << fr.size() << std::endl;
    return status;
}

// no need to count splitters here since only PI-s are considered here
void first_stage(UI_NTT& costs, UI_NTT3& pred, UI_NTT& func, UI_NTT& depths, UI_set& frontier, UI_NTT& xcosts, UI_NTT3& xpred, UI_NTT& xfunc, UI_NTT& xdepths, UI_set& xorable_frontier, const bool three_input = true, const bool verbose = false)
{
    std::cout << "Exploring first functions" << std::endl;
    for (auto ii = frontier.begin(); ii != frontier.end(); ++ii) {
        if (*ii % 1000 == 1) std::cout << "\r" << "ci = " << *ii;
        for (auto jj = std::next(ii); jj != frontier.end(); ++jj) {
            combine2(costs, pred, func, depths, xcosts, xpred, xfunc, xdepths, *ii, *jj, verbose);
            if (three_input)
            {
                for (auto kk = std::next(jj); kk != frontier.end(); ++kk) {
                    combine3(costs, pred, func, depths, xcosts, xpred, xfunc, xdepths, *ii, *jj, *kk, verbose);
                }
            }
        }
    }
    update_frontier(frontier, costs, verbose);
    update_frontier(xorable_frontier, xcosts, verbose);
}

void mergers(UI_NTT& costs, UI_NTT3& pred, UI_NTT& func, UI_NTT& depths, UI_set& frontier, UI_NTT& xcosts, UI_NTT3& xpred, UI_NTT& xfunc, UI_NTT& xdepths, UI_set& xorable_frontier, const bool verbose = false)
{
    US i = 0;
    bool status;
    do {
        status = false;
        std::cout << "Exploring depth " << ++i << std::endl;
        for (auto ii = frontier.begin(); ii != frontier.end(); ++ii) {
            std::cout << "\r" << "ci = " << *ii;
            for (auto jj = std::next(ii); jj != frontier.end(); ++jj) {
                UI tt = *ii | *jj;
                // int old_cost = costs[tt];
                if (verbose)
                {
                    std::cout << "CONSIDERING " << *ii << " | " << *jj << std::endl;
                    std::cout << "CONSIDERING " << std::bitset<16>(*ii) << " | " << std::bitset<16>(*jj) << std::endl;
                    std::cout << "\tBEFORE " << tt << " = " << std::bitset<16>(tt) << " -> " << costs[tt] << std::endl;
                }
                status |= merge(costs, pred, func, depths, xcosts, xpred, xfunc, xdepths, *ii, *jj, verbose);
                if (verbose) 
                {
                    std::cout << "\tAFTER  " << tt << " = " << std::bitset<16>(tt) << " -> " << costs[tt] << std::endl;
                }
                // assert( costs[tt] <= old_cost );
            }
        }
        std::cout << std::endl;
        update_frontier(frontier, costs, verbose);
        // std::cout << " GFRONTIER SIZE IS " << frontier.size() << std::endl;
        update_frontier(xorable_frontier, xcosts, verbose);
        // std::cout << " XFRONTIER SIZE IS " << xorable_frontier.size() << std::endl;
        std::cout << "Found " << frontier.size() << " tt-s at depth " << i << std::endl;
        std::cout << "Found " << xorable_frontier.size() << " XORable tt-s at depth " << i << std::endl;
        if ( (frontier.size() == NTT) && (xorable_frontier.size() == NTT) ) break;
    } while ( status );
}

void final_stage(UI_NTT& costs, UI_NTT3& pred, UI_NTT& func, UI_NTT& depths, UI_NTT& xcosts, UI_NTT3& xpred, UI_NTT& xfunc, UI_NTT& xdepths, const bool verbose = false, const bool debug = false)
{
    UI_NTT not_costs;   not_costs.fill(INF);
    UI_NTT not_depths;  not_depths.fill(INF);
    UI_NTT xor_costs;   xor_costs.fill(INF);
    UI_NTT xor_depths;  xor_depths.fill(INF);
    UI_NTT2 xor_pred;   xor_pred.fill(INF2);
    UI xor_cost, dff_cost, tt;


    for (auto ii = frontier.begin(); ii != frontier.end(); ++ii)
    {   
        tt = *ii ^ ONES;
        not_depths[tt] = depths[*ii];
        assert(costs[*ii] < INF);
        not_costs[tt] = costs[*ii] + COSTS[kNOT];
    }

    std::cout << "Exploring XOR-s" << std::endl;
    for (auto ii = xorable_frontier.begin(); ii != xorable_frontier.end(); ++ii)
    {   
        
        if (verbose)
        {
            std::cout << "Checking:" <<  std::endl;
            std::cout << std::bitset<16>(*ii) << " = " << *ii <<  std::endl;
            std::cout << "FUNC = " << xfunc[*ii] <<  std::endl;
            std::cout << std::bitset<16>(xpred[*ii][0]) << " = " << xpred[*ii][0] <<  std::endl;
            std::cout << std::bitset<16>(xpred[*ii][1]) << " = " << xpred[*ii][1] <<  std::endl;
        }
        assert(xfunc[*ii] != fMERGE || (xpred[*ii][0] | xpred[*ii][1]) == (xpred[*ii][0] ^ xpred[*ii][1]));

        for (auto jj = std::next(ii); jj != xorable_frontier.end(); ++jj) {
            std::cout << "\t\r" << *ii << ' ' << *jj;
            if (verbose)
            {
                std::cout << "Checking:" <<  std::endl;
                std::cout << "FUNC = " << xfunc[*ii] <<  std::endl;
                std::cout << std::bitset<16>(*jj) << " = " << *jj <<  std::endl;
                std::cout << std::bitset<16>(xpred[*jj][0]) << " = " << xpred[*jj][0] <<  std::endl;
                std::cout << std::bitset<16>(xpred[*jj][1]) << " = " << xpred[*jj][1] <<  std::endl;
            }
            assert(xfunc[*jj] != fMERGE || (xpred[*jj][0] | xpred[*jj][1]) == (xpred[*jj][0] ^ xpred[*jj][1]));

            assert(costs[*ii] < INF && costs[*jj] < INF);
            // xor_cost = costs[*ii] + costs[*jj] + COSTS[kXMERGE] + COSTS[kSPL] * count_spl(xpred, *ii, *jj) + COSTS[kDFF] * uabs(depths[*ii], depths[*jj]);
            xor_cost = costs[*ii] + costs[*jj] + COSTS[kXMERGE] + COSTS[kDFF] * count_dff(xdepths, *ii, *jj); //+ COSTS[kSPL] * count_spl(xpred, *ii, *jj) 
            tt = *ii ^ *jj;
            // std::cout << "Discovered " << tt << " with cost " << xor_cost <<
            if (xor_cost < xor_costs[tt])
            {
                xor_depths[tt] = std::max(depths[*ii], depths[*jj]);
                xor_costs[tt] = xor_cost; //need DFFs for path balancing
                xor_pred[tt][0] = *ii;
                xor_pred[tt][1] = *jj;
            }
        }
    }

    for (UI tt = 0u; tt < NTT; ++tt)
    {
        // Prevent overflow
        // if (tt == 0xe800){
        //     true;
        // }
        if (costs[tt] < INF)
        {
            dff_cost = costs[tt];
            if (func[tt] !=  fDFF) dff_cost += COSTS[kDFF];
        }
        else
        {
            dff_cost = INF;
        }
        
        if (costs[tt] == INF && not_costs[tt] == INF && xor_costs[tt] == INF) continue;
        else if (dff_cost <= not_costs[tt] && dff_cost <= xor_costs[tt]) //Adding DFF has the smallest cost
        {
            if (func[tt] !=  fDFF) 
            {
                costs[tt] = dff_cost;
                func[tt] = fDFF;
                depths[tt]++;
                xcosts[tt] = dff_cost; //NOTE: All TT-s become XOR-able at this stage
                xfunc[tt] = fDFF;
                xdepths[tt] = depths[tt];
            }
        }
        else if (not_costs[tt] <= dff_cost && not_costs[tt] <= xor_costs[tt])  //Adding NOT has the smallest cost
        {
            assert(costs[tt] <= xcosts[tt]);
            update( costs,  pred,  func,  depths, tt, not_costs[tt], not_depths[tt] + 1, fNOT, tt ^ ONES);
            update( xcosts,  xpred,  xfunc,  xdepths, tt, not_costs[tt], not_depths[tt] + 1, fNOT, tt ^ ONES);
            // Each merger has to be capped by the 
            // update(xcosts, xpred, xfunc, xdepths, tt, not_costs[tt], not_depths[tt] + 1, fNOT, tt ^ ONES);
            // assert(costs[tt] == xcosts[tt]);
            // assert(pred[tt] == xpred[tt]);
            // assert(func[tt] == xfunc[tt]);
            // assert(depths[tt] == xdepths[tt]);
        }
        else if (xor_costs[tt] <= dff_cost && xor_costs[tt] <= not_costs[tt]) //Adding XOR has the smallest cost
        {
            std::cout << "ADDED XOR: " << std::bitset<16>(tt) << std::endl;
            assert(costs[tt] <= xcosts[tt]);
            update( costs,  pred,  func,  depths, tt, xor_costs[tt], xor_depths[tt] + 1, fXOR, xor_pred[tt][0], xor_pred[tt][1]);
            xcosts[tt] = costs[tt];
            xpred[tt] = pred[tt];
            xfunc[tt] = func[tt];
            xdepths[tt] = depths[tt];
            // update(xcosts, xpred, xfunc, xdepths, tt, xor_costs[tt], xor_depths[tt] + 1, fXOR, xor_pred[tt][0], xor_pred[tt][1]);
            // std::cout << costs[tt] << " <=> " << xcosts[tt] << std::endl;
            // std::cout << pred[tt][0] << " <=> " << xpred[tt][0] << std::endl;
            // std::cout << pred[tt][1] << " <=> " << xpred[tt][1] << std::endl;
            // std::cout << pred[tt][2] << " <=> " << xpred[tt][2] << std::endl;
            // std::cout << func[tt] << " <=> " << xfunc[tt] << std::endl;
            // std::cout << depths[tt] << " <=> " << xdepths[tt] << std::endl;
            // assert(costs[tt] == xcosts[tt]);
            // assert(pred[tt] == xpred[tt]);
            // assert(func[tt] == xfunc[tt]);
            // assert(depths[tt] == xdepths[tt]);
        }
        else
        {
            std::cout << "dff_cost = " << dff_cost;
            std::cout << "not_cost = " << not_costs[tt];
            std::cout << "xor_cost = " << xor_costs[tt];
            throw std::runtime_error("Unsupported case");
        }
    }
    update_frontier(frontier, costs, verbose);
    update_frontier(xorable_frontier, xcosts, verbose);
}

void write_database(const std::string filename, CUI_NTT& costs, CUI_NTT3& pred, CUI_NTT& func, CUI_NTT& depths, CUI_NTT& xcosts, CUI_NTT3& xpred, CUI_NTT& xfunc, CUI_NTT& xdepths, const std::string& header = "", const bool verbose = false)
{
	std::fstream file;
	file.open(filename, std::ios::out);
	if (!file) {
		std::cout << "File not created!";
        return;
	}
    if (verbose) std::cout << "Successfully created " << filename << std::endl;

    file << "/*\n" << header << "\n */\n" ;

    file << "const uint32_t INF = " << INF << ";\n" ;
    file << "const uint32_t NTT = " << NTT << ";\n" ;
    // file << "const uint32_t areas[NTT]  = {[0 ... NTT-1] = INF};\n" ;
    // file << "const uint32_t delays[NTT] = {[0 ... NTT-1] = INF};\n" ;
    file << "// costs, delays, func, pred[0], pred[1], pred[2], xcosts, xdelays, xfunc, xpred[0], xpred[1], xpred[2]\n" ;
    file << "const uint32_t area_delay[NTT][12]  = { \n";

    // auto i = 0u;
    bool pi_or_const;
    for (auto tt = 0u; tt < NTT; ++tt)
    {   
        pi_or_const = (std::find(std::begin(PI), std::end(PI), tt) != std::end(PI)) || (tt == 0 || tt == NTT-1);
        if (pi_or_const)
        {   
            file << "{0, 0, 255, INF, INF, INF, 0, 0, 255, INF, INF, INF}";
            
        }
        else if (costs[tt] < INF)
        {   
            file << "{ " <<  costs[tt] << ", " <<  depths[tt] << ", " <<  func[tt] << ", ";
            file <<  pred[tt][0] << ", " <<  pred[tt][1] << ", " <<  pred[tt][2] << ", ";
            file << xcosts[tt] << ", " << xdepths[tt] << ", " << xfunc[tt] << ", ";
            file << xpred[tt][0] << ", " << xpred[tt][1] << ", " << xpred[tt][2] << "}";
        }
        else
        {
            file << "{INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF}";
        }
        // if (++i % max_entries == 0) file << "\n" ;
        if (tt < NTT - 1) 
        {
            file << ",";
        }
        else
        {
            file << " ";
        }
    file << " // " << std::hex << tt << std::dec << " = " << std::bitset<16>(tt) << " = " << tt << '\n';// fmt::format(", // {0:x} = {0:b} = {0:x}", 
    }
    file << "};\n" ;

    file << "struct lut_sfq_cost\n";
    file << "{\n";
    file << "\t" << "std::pair<uint64_t, uint64_t> operator()( kitty::dynamic_truth_table const& tt )\n";
    file << "\t" << "{\n";
    file << "\t\t" << "uint64_t area = area_delay[tt._bits[0]][0];\n" ;
    file << "\t\t" << "uint64_t delay = area_delay[tt._bits[0]][1];\n" ;
    file << "\t\t" << "return { area , delay };\n" ;
    file << "\t" << "}\n";
    file << "};\n";

    file.close();
}

int main()
{
    
    #pragma region initialization // Initialize costs and depths
    UI_NTT  xorable_costs;  xorable_costs.fill(INF);
    UI_NTT3 xorable_pred;   xorable_pred.fill(INF3);
    UI_NTT  xorable_func;   xorable_func.fill(255);
    UI_NTT  xorable_depths; xorable_depths.fill(INF);

    UI_NTT  global_costs;   global_costs.fill(INF);
    UI_NTT3 global_pred;    global_pred.fill(INF3);
    UI_NTT  global_func;    global_func.fill(255);
    UI_NTT  global_depths;  global_depths.fill(INF);

    for (auto & tt : frontier) {
        global_costs[tt] = 0;
        xorable_costs[tt] = 0;
        global_depths[tt] = 0;
        xorable_depths[tt] = 0;
    }

    UI i = 0u;
    for (auto & d : global_depths)
    {
        if (d < INF) std::cout << "tt : " << i << " --> " << d << std::endl;
        i++;
    }
    #pragma endregion initialization

    US qqq = 0; 
    while (frontier.size() < NTT)
    {
        first_stage(global_costs, global_pred, global_func, global_depths, frontier,
                xorable_costs, xorable_pred, xorable_func, xorable_depths, xorable_frontier, false); //false means 3-input funcs are disabled
        std::cout << "Found " << frontier.size() << " tt-s in first stage" << std::endl;
        print_frontier(global_costs, global_pred, global_func, global_depths);
        
        mergers(global_costs, global_pred, global_func, global_depths, frontier,
                xorable_costs, xorable_pred, xorable_func, xorable_depths, xorable_frontier);
        std::cout << "Found " << frontier.size() << " tt-s in merging stage" << std::endl;
        print_frontier(global_costs, global_pred, global_func, global_depths);

        final_stage(global_costs, global_pred, global_func, global_depths, xorable_costs, xorable_pred, xorable_func, xorable_depths);
        std::cout << "Found " << frontier.size() << " tt-s in final stage" << std::endl;
        print_frontier(global_costs, global_pred, global_func, global_depths);
        if (++qqq == 2) break;
    }

    // const std::string header =  "\t All 1-lvl 4-input TT-s for RSFQ technology \n"
    //                             "\t - first function disabled (e.g. no AND/MAJ)\n"
    //                             "\t - confluence buffer tree (OR-tree\n"
    //                             "\t - completed by DFF-s, INV-s, and XOR-s"
    //                             "\t\t (xorable tracking enabled)";

    const bool write = false;
    if constexpr (write)
    {
        const std::string header =  "\t All 4-input TT-s for RSFQ technology \n"
                                    "\t - first functions:\n"
                                    "\t\t - AND2:\n"
                                    "\t\t - OR2:\n"
                                    // "\t\t - AND3:\n"
                                    // "\t\t - OR3:\n"
                                    // "\t\t - MAJ3:\n"
                                    "\t - confluence buffer tree (OR-tree\n"
                                    "\t - completed by DFF-s, INV-s, and XOR-s"
                                    "\t\t (xorable tracking enabled)";

        write_database("/Users/brainkz/Documents/GitHub/mockturtle/include/mockturtle/utils/bench_files/tt4_first_2_with_xorables.hpp", global_costs, global_pred, global_func, global_depths, xorable_costs, xorable_pred, xorable_func, xorable_depths, header);
    }

}
