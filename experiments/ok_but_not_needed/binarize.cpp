#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <vector>

#include "BNN/binarize.cpp" //binarization engine

typedef unsigned int UI;
typedef const unsigned int CUI;

#define VERBOSE

int main()
{
    std::string  inputImages = "/Users/brainkz/Documents/GitHub/machine-learning-01/data/datasets/mnist/train-images.idx3-ubyte";
    std::string outputImages = "/Users/brainkz/Documents/GitHub/mockturtle/experiments/binarized_images.data";
    std::string  inputLabels = "/Users/brainkz/Documents/GitHub/machine-learning-01/data/datasets/mnist/train-labels.idx1-ubyte";
    std::string outputLabels = "/Users/brainkz/Documents/GitHub/mockturtle/experiments/binarized_labels.data";

    auto[Nimg, Nrow, Ncol, img_vec, lbl_vec] = BNN::binarize_ubyte_char(inputImages, inputLabels);
    BNN::write_data(Nimg, Nrow, Ncol, img_vec, lbl_vec, outputImages, outputLabels);
    auto[Nimg_new, Nrow_new, Ncol_new, img_bits, lbl_bits] = BNN::read_data(outputImages, outputLabels);
    assert(Nimg = Nimg_new);
    assert(Nrow = Nrow_new);
    assert(Ncol = Ncol_new);
    // assert(img_vec = img_bits);
    // assert(lbl_vec = lbl_bits);
    std::cout << "#Labels = " << lbl_bits.size() << std::endl;
    for (auto label : lbl_bits)
    {
        std::cout << label;
    }
    std::cout << std::endl;
}