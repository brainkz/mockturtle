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

#include "mockturtle/io/write_aiger.hpp"
#include "mockturtle/io/aiger_reader.hpp"

#include "BNN/binarize.cpp" //binarization engine
#include "BNN/ga.cpp" //binarization engine

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

template <typename vec>
void print_2d_vec(vec a)  
{
    for (auto & row : a)
    {
        for (auto & item: row)
        {
            std::cout << ' ' << item;
        }
        std::cout << std::endl;
    }
}

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

mockturtle::aig_network generate_reservoir()
{
    KLUT reservoir;

    // Generate PI-s
    std::vector<SIGNAL> inputs;
    inputs.reserve(NUM_PI);

    for (auto i = 0u; i < NUM_PI; ++i)
    {
        SIGNAL s = reservoir.create_pi();
        inputs.push_back(s);
    }
    
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
    return res_aig;
}

template <typename Ntk>
std::vector<std::vector<bool>> sim_reservoir(Ntk res_aig, std::vector<bool> img_vec, UI Nimg, UI Npxl)
{
    UI Nout = res_aig.num_pos();
    std::vector<std::vector<bool>> outputs;
    outputs.reserve(Nimg);
    std::cout << "The image vector size is " << img_vec.size() << std::endl;

    for (auto i = 0u; i < Nimg; ++i)
    {
        outputs.push_back( std::vector<bool> {});
        outputs[i].reserve(Nout);
    }

    std::cout << "Simulating the image vectors" << std::endl;
    auto count = 0u;

    // #pragma omp parallel for schedule(dynamic)
    // #pragma vector aligned
    std::vector<bool> img;
    for (auto i = 0u; i < Nimg; i++)
    {
        img.assign(img_vec.begin() + i*Npxl, img_vec.begin() + i*Npxl + Npxl - 1);
        mockturtle::default_simulator<bool> sim( img );
        const auto values = mockturtle::simulate<bool>( res_aig, sim );
        // outputs.push_back(values);
        outputs[count] = values;
        // std::cout << '(' << ++count << ' ' << outputs.size() << ')';
        std::cout << ++count << '\r';
    }
    return outputs;    
}
/*
std::vector<bool> stringToVector(const std::string& str) {
    std::vector<bool> vec;
    for (char c : str) {
        if (c == '1') {
            vec.push_back(true);
        } else if (c == '0') {
            vec.push_back(false);
        }
    }
    return vec;
*/
// test whether the file is packed correctly

int main()
{
    
    // std::string  Images = "/Users/brainkz/Documents/GitHub/mockturtle/experiments/MNIST.images";
    // std::string  Labels = "/Users/brainkz/Documents/GitHub/mockturtle/experiments/binarized_labels.labels";
    // Import the MNIST dataset 
    std::string  Images = "/Users/brainkz/Documents/GitHub/mockturtle/experiments/binarized_images.data";
    std::string  Labels = "/Users/brainkz/Documents/GitHub/mockturtle/experiments/binarized_labels.data";
    auto[Nimg, Nrow, Ncol, img_vec, lbl_vec] = BNN::read_data(Images, Labels);

    std::string  RES_OUTPUT_BIN = "/Users/brainkz/Documents/GitHub/mockturtle/experiments/RES_OUTPUT.data";
    auto[dims_out, data_out] = BNN::read_ND_vector( RES_OUTPUT_BIN );

    std::cout << "Successfully read binary data" << std::endl;

    assert(Nimg == dims_out[0]);
    auto Nout = dims_out[1];
    std::cout << "Data dimensions: " << std::endl;
    for (const auto & dim : dims_out) { std::cout << dim << std::endl; }

    std::cout << "Reshaping bit data" << std::endl;
    std::vector<std::vector<bool>> outputs(Nimg, std::vector<bool>(Nout)); // create output vector with dimensions 60000 by 2000
    for (int i = 0; i < Nimg; ++i) {
        for (int j = 0; j < Nout; ++j) {
            outputs[i][j] = data_out[i * Nout + j];
        }
    }

    UI Nweights = Nout;
    int wmin = -20;
    int wmax =  20;
    int bmin = -1000;
    int bmax =  1000;
    UI Nspecies = 100;
    UI numdeaths = 60;
    UI num_mutate = 40;
    UI num_xover = 20;
    UI wsigma = 5;
    UI bsigma = 100;
    CUI NUM_ITER = 100;

    std::cout << "Generating the population" << std::endl;
    auto pop = GA::initial_population(Nspecies, Nweights, wmin, wmax, bmin, bmax);
    print_2d_vec(pop);

    std::vector<UI> num_correct(Nspecies, Nimg+1);

    for (auto i = 0u; i < NUM_ITER; i++)
    {
        std::cout << "Iteration " << i << ':' << std::endl;
        GA::ga_iteration(outputs, lbl_vec, pop, num_correct, numdeaths, Nimg, num_mutate, num_xover, wsigma, bsigma, wmin, wmax, bmin, bmax);
    }
}

int main_dummy()
{
    std::string  RES_OUTPUT = "/Users/brainkz/Documents/GitHub/mockturtle/experiments/RES_OUTPUT.txt";
    std::ifstream txt( RES_OUTPUT );
    assert(txt.is_open());

    std::string line;
    std::vector<bool> data_txt;

    std::getline(txt, line);
    size_t dim2 = line.size();
    size_t dim1 = 1;
    for (char c : line) {
        data_txt.push_back( ((c == '1')? true : false) );
    }

    while (std::getline(txt, line))
    {
        if (line.size() == 0) continue;
        for (char c : line) {
            data_txt.push_back( ((c == '1')? true : false) );
        }
#ifdef DEBUG_TXT
        std::cout << line.size() << ' ' << ++dim1 << ' ' << data_txt.size() << std::endl;
#endif
    }
    std::cout << dim1 << " * " << dim2 << " = " << data_txt.size() << std::endl;
    assert(dim1 * dim2 == data_txt.size());

    const std::vector<size_t> dims_txt = { dim1, dim2 };

    const std::string  RES_OUTPUT_BIN = "/Users/brainkz/Documents/GitHub/mockturtle/experiments/RES_OUTPUT.data";
    BNN::write_ND_vector( dims_txt, data_txt, RES_OUTPUT_BIN );

    auto[dims_bin, data_bin] = BNN::read_ND_vector( RES_OUTPUT_BIN );

    std::cout << (dims_bin == dims_txt) << ' ' << (data_bin == data_txt) << std::endl;
    UI WIDTH = 80u;
    for (auto i = 0u; i < data_bin.size(); i+=WIDTH)
    {   
        std::cout << i << " to " << (i + WIDTH - 1) << std::endl;
        for (auto j = 0u; j < WIDTH; j++)
        {
            std::cout << data_bin[WIDTH * i + j];
        }
        std::cout << std::endl;
        for (auto j = 0u; j < WIDTH; j++)
        {
            std::cout << data_txt[WIDTH * i + j];
        }
        std::cout << std::endl;
    }
}

int main_write_test()
{
    // 1. Generate reservoirs
    // 2. Generate connections from initial gates and to final gate
    // 3. Optimize connections based on dataset
    // 4. Measure performance
    const auto res_filename = "/Users/brainkz/Documents/GitHub/mockturtle/experiments/RES_AIG.aiger";
#if false
    mockturtle::aig_network res_aig = generate_reservoir();
    mockturtle::write_aiger(res_aig, res_filename);
#else
    mockturtle::aig_network res_aig;
    lorina::read_aiger( res_filename, mockturtle::aiger_reader( res_aig ) );
#endif
    // Import the MNIST dataset 
    std::string  Images = "/Users/brainkz/Documents/GitHub/mockturtle/experiments/MNIST.images";
    std::string  Labels = "/Users/brainkz/Documents/GitHub/mockturtle/experiments/MNIST.labels";
    auto[Nimg, Nrow, Ncol, img_vec, lbl_vec] = BNN::read_data(Images, Labels);

    std::vector<std::vector<bool>> outputs = sim_reservoir(res_aig, img_vec, Nimg, Nrow * Ncol); // need to import the MNIST
    // create function to write reservoir outputs to file, then work with that file

    std::vector<size_t> dims_in;
    dims_in.push_back(outputs.size());
    dims_in.push_back(outputs[0].size());
    std::vector<bool> data_in;
    auto q = 0ull;
    for (const auto& innerVector : outputs) {
        data_in.insert(data_in.end(), innerVector.begin(), innerVector.end());
#ifdef DEBUG
        for (auto bit : innerVector)  { std::cout << bit; q++; }
        std::cout << ' ' << q << std::endl;
#endif
    }
#ifdef DEBUG
    std::cout << "INPUT DIMENSIONS ARE";
    for (const auto & dim : dims_in)
    {
        std::cout << " " << dim;
    }
    std::cout << std::endl << "TOTAL NUMBER OF BITS IN INPUT IS " << data_in.size() << std::endl;
#endif

    std::string  RES_OUTPUT_BIN = "/Users/brainkz/Documents/GitHub/mockturtle/experiments/RES_OUTPUT.data";
    BNN::write_ND_vector( dims_in, data_in, RES_OUTPUT_BIN );


    auto[dims_out, data_out] = BNN::read_ND_vector( RES_OUTPUT_BIN );
#ifdef DEBUG
    std::cout << "OUTPUT DIMENSIONS ARE";
    for (const auto & dim : dims_out)
    {
        std::cout << " " << dim;
    }
    std::cout << std::endl << "TOTAL NUMBER OF BITS IN OUTPUT IS " << data_out.size() << std::endl;

    std::cout << (dims_out == dims_in) << ' ' << (data_out == data_in) << std::endl;
    UI WIDTH = 200u;
    for (auto i = 0u; i < data_out.size(); i+=WIDTH)
    {   
        std::cout << i << " to " << (i + WIDTH - 1) << std::endl;
        for (auto j = 0u; j < WIDTH; j++)
        {
            std::cout << data_out[i + j];
        }
        std::cout << std::endl;
        for (auto j = 0u; j < WIDTH; j++)
        {
            std::cout << data_in[i + j];
        }
        std::cout << std::endl;
    }
#endif
    return 0;
}   

int main_unused ()
{
    // 1. Generate reservoirs
    // 2. Generate connections from initial gates and to final gate
    // 3. Optimize connections based on dataset
    // 4. Measure performance
    const auto res_filename = "/Users/brainkz/Documents/GitHub/mockturtle/experiments/RES_AIG.aiger";

    mockturtle::aig_network res_aig = generate_reservoir();
    mockturtle::write_aiger(res_aig, res_filename);

    // mockturtle::aig_network res_aig;
    // lorina::read_aiger( res_filename, mockturtle::aiger_reader( res_aig ) );


    // Import the MNIST dataset 
    std::string  Images = "/Users/brainkz/Documents/GitHub/mockturtle/experiments/MNIST.images";
    std::string  Labels = "/Users/brainkz/Documents/GitHub/mockturtle/experiments/MNIST.labels";
    auto[Nimg, Nrow, Ncol, img_vec, lbl_vec] = BNN::read_data(Images, Labels);
    // auto Nimg = img_vec.size();
    auto Npxl = Nrow * Ncol;

    // Determine bit-vectors corresponding to training dataset
    // Optimize the connections
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
    std::vector<bool> img;
    for (auto i = 0u; i < Nimg; i++)
    {
        img.assign(img_vec.begin() + i*Npxl, img_vec.begin() + i*Npxl + Npxl - 1);
        mockturtle::default_simulator<bool> sim( img );
        const auto values = mockturtle::simulate<bool>( res_aig, sim );
        // outputs.push_back(values);
        outputs[count] = values;
        std::cout << ++count << ' ';
    }

    // // OPTIMIZATION PART

    // return 0;
}

// int main()
// {


//     std::string  Images = "/Users/brainkz/Documents/GitHub/mockturtle/experiments/MNIST.images";
//     std::string  Labels = "/Users/brainkz/Documents/GitHub/mockturtle/experiments/MNIST.labels";
//     auto[Nimg, Nrow, Ncol, img_vec, lbl_vec] = BNN::read_data(Images, Labels);
//     // auto Nimg = img_vec.size();
//     auto Npxl = Nrow * Ncol;


// }