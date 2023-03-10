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

// #include "BNN/binarize.cpp" //binarization engine
#include "BNN/ga.cpp" //binarization engine

// #include <armadillo>
#include <ensmallen.hpp>

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

constexpr CUS MAXVARS = 9u;
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

std::vector<uint64_t> generate_random_numbers(size_t N) 
{
    std::vector<uint64_t> numbers(N);
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> dis;

    for (size_t i = 0; i < N; i++) {
        numbers[i] = dis(gen);
    }

    return numbers;
}

template <typename T>
static inline T next_pow_2(T v)
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    return ++v;
}

mockturtle::aig_network random_klut_bin_tree(size_t NUM_PI)
{
    std::random_device rd;
    std::mt19937 g(rd());

    KLUT klut_ntk;
    for (auto i = 0u; i < NUM_PI; ++i)
    {
        SIGNAL s = klut_ntk.create_pi();
    }

    std::vector<SIGNAL> pi_vec(next_pow_2(NUM_PI));
    
    std::iota(pi_vec.begin(), pi_vec.end(), 0);
    std::shuffle(pi_vec.begin(), pi_vec.end(), g);
    
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
        std::cout << "Simulating image " << i << outputs.size() << std::endl;
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

int main()
{
    return 0;
}
/*
int main()
{
    const auto res_filename = "/Users/brainkz/Documents/GitHub/mockturtle/experiments/RES_AIG_new.aiger";
#if true
    mockturtle::aig_network res_aig = generate_reservoir();
    mockturtle::write_aiger(res_aig, res_filename);
#else
    mockturtle::aig_network res_aig;
    lorina::read_aiger( res_filename, mockturtle::aiger_reader( res_aig ) );
#endif
    // Import the MNIST dataset 
    std::string  Images = "/Users/brainkz/Documents/GitHub/mockturtle/experiments/binarized_images.data";
    std::string  Labels = "/Users/brainkz/Documents/GitHub/mockturtle/experiments/binarized_labels.data";
    auto[Nimg, Nrow, Ncol, img_vec, lbl_vec] = BNN::read_data(Images, Labels);

    std::vector<std::vector<bool>> data_out = sim_reservoir(res_aig, img_vec, Nimg, Nrow * Ncol); // need to import the MNIST
    UI Nout = data_out[0].size();
    // std::string  RES_OUTPUT_BIN = "/Users/brainkz/Documents/GitHub/mockturtle/experiments/RES_OUTPUT.data";
    // auto[dims_out, data_out] = BNN::read_ND_vector( RES_OUTPUT_BIN );

    std::cout << "Successfully read binary data" << std::endl;

    // assert(Nimg == dims_out[0]);
    // auto Nout = dims_out[1];
    // std::cout << "Data dimensions: " << std::endl;
    // for (const auto & dim : dims_out) { std::cout << dim << std::endl; }

    std::cout << "Reshaping bit data" << std::endl;
    // std::vector<std::vector<bool>> outputs(Nimg, std::vector<bool>(Nout)); // create output vector with dimensions 60000 by 2000
    arma::imat outputs(Nimg, Nout);
    for (int i = 0; i < Nimg; ++i) {
        for (int j = 0; j < Nout; ++j) {
            // outputs(i,j) = data_out[i * Nout + j];
            outputs(i,j) = data_out[i][j];
        }
    }
    arma::icolvec labels(lbl_vec.size());
    for (int i = 0; i < lbl_vec.size(); ++i) {
        labels(i) = lbl_vec[i];
    }

    std::cout << "Reshaping complete" << std::endl;
    class objFunc
    {
        private:
            // const std::vector<std::vector<bool>> res_output_signals;
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
            arma::imat res_output_signals;
            arma::icolvec labels;
            UI Nout;
            unsigned long long ct = 0;

        public:
            void set_out_sig(arma::imat  sig) {res_output_signals = sig;}
            void set_labels (arma::icolvec  lbl) {labels = lbl;}
            void set_Nout   (arma::uword Nweights) {Nout = Nweights;}

            double Evaluate(const arma::mat& weights)
            {
                auto bias = weights.tail_rows(1)(0);
                auto total = ( (res_output_signals * weights.head_rows(Nout)) > bias ) == labels;
                auto num_correct = arma::sum( total );
                std::cout << ct++ << ' ' << (double) num_correct << std::endl;
                return num_correct;
                // return num_correct;
            };
    };

    std::cout << "Class definition complete" << std::endl;
    // The minimum is at x = [0 0 0].  Our initial point is chosen to be
    // [1.0, -1.0, 1.0].
    arma::mat weights = arma::zeros<arma::mat>( Nout + 1, 1);

    // Create simulated annealing optimizer with default options.
    // The ens::SA<> type can be replaced with any suitable ensmallen optimizer
    // that is able to handle arbitrary functions.

    std::cout << "Creating optimizer object" << std::endl;
    // ens::DE optimizer;
    ens::CNE optimizer;
    objFunc f; // Create function to be optimized.
    std::cout << "Setting internal variables" << std::endl;
    f.set_out_sig(outputs);
    f.set_labels(labels);
    f.set_Nout(Nout);
    std::cout << "Optimizing" << std::endl;
    optimizer.Optimize(f, weights);

    std::cout << "Result " << weights; 
    
}
*/