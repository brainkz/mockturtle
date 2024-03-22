#pragma once

#include <cstdlib> // Random numbers
#include <chrono>
#include <ctime>
#include <iostream>
#include <fstream> //include the filestreamobject as the header files
#include <tuple>
#include <iterator>
#include <array>
#include <algorithm>
#include <vector>
#include <set>
#include <unordered_map>
#include <bitset>
#include <cstdio>
#include "mockturtle/networks/klut.hpp"

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>

#include "ortools/linear_solver/linear_solver.h"


#include "mockturtle/algorithms/klut_to_graph.hpp" // convert klut to aig
#include "mockturtle/algorithms/simulation.hpp" // simulation engine

#include "BNN/binarize.cpp" //binarization engine

#include <execution>


#pragma region constants
typedef unsigned short US;
typedef unsigned int   UI;
typedef unsigned long  UL;
typedef const unsigned short CUS;
typedef const unsigned int   CUI;
typedef const unsigned long  CUL;
typedef mockturtle::klut_network KLUT;
typedef mockturtle::klut_network::signal SIGNAL;

typedef kitty::dynamic_truth_table DTT;

constexpr CUS MAXVARS = 4u;
CUI RES_SIZE = 2000u;
CUI FINAL_THRES = 2000u;
CUI NUM_PI = 28*28;

#pragma endregion constants



template <typename DTYPE>
std::vector<DTYPE> randomArrayElements(const std::vector<DTYPE>& arr, const std::size_t k) {
    std::vector<DTYPE> result;
    result.reserve(k);

    auto n = arr.size();
    for (std::size_t i = 0; i < k; ++i) {
        std::size_t index = std::rand() % n;
        result.push_back(arr[index]);
    }
    return result;
}

int main ()
{
    // 1. Generate reservoirs
    // 2. Generate connections from initial gates and to final gate
    // 3. Optimize connections based on dataset
    // 4. Measure performance

    KLUT reservoir;

    // Generate PI-s
    std::vector<SIGNAL> inputs;
    inputs.reserve(NUM_PI);

    for (auto i = 0u; i < NUM_PI; ++i)
    {
        SIGNAL s = reservoir.create_pi();
        inputs.push_back(s);
    }
    // std::array<std::string, 5> arr = {"asda","qweq","1241","123412341","sdb45"};
    // std::array<std::string, 3> selected = randomArrayElements<std::string, 5, 3>(arr);
    
    // Generate random TT-s
    std::vector<SIGNAL> res_output_nodes;
    std::vector<UI> res_output_index;
    std::srand(std::time(nullptr));
    for (auto i = 0u; i < RES_SIZE; i++)
    {
        // create random KLUT
        DTT function( MAXVARS );
        UL func = std::rand() & ((1 << MAXVARS) - 1);
        std::array<UL, 1> func_word {func};
        kitty::create_from_words( function, func_word.begin(), func_word.end());

        // select input nodes for KLUT
        std::vector<SIGNAL> selected = randomArrayElements( inputs, MAXVARS );

        // create a KLUT gate
        SIGNAL s = reservoir.create_node( selected, function );

        // assign PO to KLUT output
        res_output_nodes.push_back(s);
        res_output_index.push_back( reservoir.create_po(s) );

        std::cout << "Selected func: " << function._bits[0] << std::endl;;
        std::cout << "Selected elements: ";
        for (auto i : selected) {
            std::cout << i << " ";
        }
        std::cout << std::endl;
    }

    // Convert klut to aig network
    mockturtle::aig_network res_aig = mockturtle::cleanup_dangling ( mockturtle::convert_klut_to_graph<mockturtle::aig_network>( reservoir ) );


    // Import the MNIST dataset 
    std::string  inputImages = "/Users/brainkz/Documents/GitHub/machine-learning-01/data/datasets/mnist/train-images.idx3-ubyte";
    std::string  inputLabels = "/Users/brainkz/Documents/GitHub/machine-learning-01/data/datasets/mnist/train-labels.idx1-ubyte";
    auto[Nimg, Nrow, Ncol, img_vec, lbl_vec] = BNN::binarize_ubyte(inputImages, inputLabels);
    // auto Nimg = img_vec.size();
    auto Npxl = Nrow * Ncol;

    std::cout << "Completed binarization of MNIST" << std::endl;
    
    // Determine bit-vectors corresponding to training dataset
    // Optimize the connections
    UI Nout = res_output_index.size();
    std::vector<std::vector<bool>> outputs;
    outputs.reserve(Nimg);
    for (auto i = 0u; i < img_vec.size(); ++i)
    {
        outputs.push_back( std::vector<bool> {});
        outputs[i].reserve(Npxl);
    }

    std::cout << "Simulating the image vectors" << std::endl;
    auto count = 0u;

    // #pragma omp parallel for schedule(dynamic)
    // #pragma vector aligned
    for (auto & img : img_vec)
    {
        mockturtle::default_simulator<bool> sim( img );
        const auto values = mockturtle::simulate<bool>( res_aig, sim );
        // outputs.push_back(values);
        outputs[count] = values;
        std::cout << ++count << ' ';
        
#ifdef VERBOSE   
        union { int num; char sym; } num_sym;     
        for (auto i = 0u; i < values.size(); i += 4 )
        {
            num_sym.num =   values[i+0] * (1 << 3) | 
                            values[i+1] * (1 << 2) |
                            values[i+2] * (1 << 1) |
                            values[i+3] * (1 << 0) ;
            num_sym.num += ((num_sym.num > 9) ? 55 : 48);
            // std::cout << num_sym.sym;
        }
#endif
        // std::cout << std::endl;
    }
    
    return 0;
}