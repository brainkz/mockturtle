#pragma once

#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <array>
#include <bitset>

#include "kitty/partial_truth_table.hpp"

typedef unsigned int UI;
typedef const unsigned int CUI;
typedef unsigned long UL;
typedef const unsigned long CUL;

CUI BITSET_WIDTH = 64;

// #define VERBOSE

// #define cast4byte(charArr) *(reinterpret_cast<int*>(charArr)); // convert 4-byte char array to 32-bit integer
#define print4bytes(charArr) std::cout << charArr[0] << charArr[1] << charArr[2] << charArr[3] << std::endl;
#define char_from8bit(arr, i) (arr[i+7] << 7) | (arr[i+6] << 6) | (arr[i+5] << 5) | (arr[i+4] << 4) | (arr[i+3] << 3) | (arr[i+2] << 2) | (arr[i+1] << 1) | arr[i];
// #define ceiling(x,y) (x + y - 1) / y;
#define ceiling(x,y) x/y + (x % y != 0)
// #define DEBUG_READ

// template<size_t dimcount, typename T>
// struct ND_vector
// {
//     typedef std::vector< typename ND_vector<dimcount-1, T>::type > type;
// };
// template<typename T>
// struct ND_vector<0,T>
// {
//     typedef T type;
// };

namespace BNN
{
inline bool is_big_endian() 
{
    union {
        UI i;
        char c[4];
    } item = {0x01020304};
    return item.c[0] == 1;
}

inline UI cast4byte(char * buf)
{
    if (is_big_endian())
    {
        return (UI)((buf[3] & 0xFF) << 24) | (UI)((buf[2] & 0xFF) << 16) | (UI)((buf[1] & 0xFF) << 8) | (UI)(buf[0] & 0xFF);
    }
    else
    {
        return (UI)((buf[0] & 0xFF) << 24) | (UI)((buf[1] & 0xFF) << 16) | (UI)((buf[2] & 0xFF) << 8) | (UI)(buf[3] & 0xFF);
    }
} 

inline UI get4byte(std::ifstream& input, char * buf)
{   
    UI out;
    input.read(buf, 4);
    if (is_big_endian())
    {
        out = (UI)((buf[3] & 0xFF) << 24) | (UI)((buf[2] & 0xFF) << 16) | (UI)((buf[1] & 0xFF) << 8) | (UI)(buf[0] & 0xFF);
    }
    else
    {
        out = (UI)((buf[0] & 0xFF) << 24) | (UI)((buf[1] & 0xFF) << 16) | (UI)((buf[2] & 0xFF) << 8) | (UI)(buf[3] & 0xFF);
    }
#ifdef VERBOSE
    std::cout << std::hex << out << std::endl;
    print4bytes(buf);
#endif
    return out;
} 

inline void int2char(std::array<unsigned char, 4>& bytes, const UI n)
{
    // note that the function writes in big-endian format for consistency with ubyte format
    bytes[0] = (n >> 24) & 0xFF;
    bytes[1] = (n >> 16) & 0xFF;
    bytes[2] = (n >>  8) & 0xFF;
    bytes[3] =  n        & 0xFF;
}

inline void print_bitvec(std::vector<bool> vec, UI idx, UI len, char false_char = ' ', char true_char = '#')
{   
    for (auto i = 0u; i < len; i++)
    {
        std::cout << (vec[idx + i]) ? true_char : false_char;
    }
    std::cout << std::endl;
}

unsigned int HALF_BYTE = 0x80;

void uintToBytes(unsigned int n, char* bytes) {
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    bytes[3] = static_cast<char>(n & 0xFF);
    bytes[2] = static_cast<char>((n >> 8) & 0xFF);
    bytes[1] = static_cast<char>((n >> 16) & 0xFF);
    bytes[0] = static_cast<char>((n >> 24) & 0xFF);
#else
    bytes[0] = static_cast<char>(n & 0xFF);
    bytes[1] = static_cast<char>((n >> 8) & 0xFF);
    bytes[2] = static_cast<char>((n >> 16) & 0xFF);
    bytes[3] = static_cast<char>((n >> 24) & 0xFF);
#endif
}



// std::pair< std::vector<std::vector<bool>>, std::vector<bool> > 
std::tuple< UI, UI, UI, std::vector<std::vector<bool>>, std::vector<bool>, std::vector<bool>> binarize_ubyte(const std::string& inputImages, const std::string& inputLabels)
{

    std::ifstream img_input (  inputImages, std::ios::binary ); 
    std::ifstream lbl_input (  inputLabels, std::ios::binary ); 

    // Read header information
    char * header = new char [4];
    const UI IMG_Magic = get4byte(img_input, header);
    const UI Nimg      = get4byte(img_input, header);
    const UI Nrow      = get4byte(img_input, header);
    const UI Ncol      = get4byte(img_input, header);
    const UI LBL_Magic = get4byte(lbl_input, header);
    const UI Nlbl      = get4byte(lbl_input, header);

    assert(Nimg == Nlbl);

    std::cout << "Number of images: "  << Nimg << std::endl;
    std::cout << "Number of rows: "    << Nrow << std::endl;
    std::cout << "Number of columns: " << Ncol << std::endl;

    UI BUFFER_SIZE = Nrow * Ncol;
    char * img_buffer = new char [BUFFER_SIZE];
    char * lbl_buffer = new char [1];
    bool bit, label;

    std::vector<bool> all_img_bits;

    std::vector<std::vector<bool>> input_vectors; //outer vector iterates over images, inner vector iterates over image bits
    std::vector<bool> output_labels; //vector iterates over labels 
    input_vectors.reserve(Nimg);
    output_labels.reserve(Nlbl);
    for (auto i = 0u; i < Nimg; ++i)
    {
        //Read image
        std::vector<bool> img_vector;
        img_vector.reserve(BUFFER_SIZE);
        img_input.read(img_buffer, BUFFER_SIZE); // read one image
        for (auto j = 0u; j < BUFFER_SIZE; j++)
        {
            bit = (bool)( *(img_buffer + j) >> 7 ); // get the jth character of the array and determine whether the item is 0 or 1
            all_img_bits.push_back(bit);
            UI newlen = all_img_bits.size();
            if (newlen != 0 && (newlen % 28) == 0) {
                std::cout << '\t';
                print_bitvec(all_img_bits, newlen-28, 28);
            }
            img_vector.push_back(bit);
#ifdef VERBOSE
            char c = bit ? '@' : ' ';
            std::cout << c;
            if ( !(j % Ncol) )  std::cout << std::endl;
#endif
        }
        input_vectors.push_back(img_vector);

        //Read label
        lbl_input.read(lbl_buffer, 1); // read one label
        label = (UI)( *(lbl_buffer) ) >= 5;
#ifdef VERBOSE
        std::cout << ((UI)( *(lbl_buffer) )) << "-->" << label << std::endl;
#endif
        output_labels.push_back(label);
    }
    return std::tuple< UI, UI, UI, std::vector<std::vector<bool>>, std::vector<bool>, std::vector<bool>> (Nimg, Nrow, Ncol, input_vectors, output_labels, all_img_bits);
}


// std::pair< std::vector<std::vector<bool>>, std::vector<bool> > 
std::tuple< UI, UI, UI, std::vector<char>, std::vector<char> > binarize_ubyte_char(const std::string& inputImages, const std::string& inputLabels)
{

    std::ifstream img_input (  inputImages, std::ios::binary ); 
    std::ifstream lbl_input (  inputLabels, std::ios::binary ); 

    // Read header information
    char * header = new char [4];
    CUI IMG_Magic = get4byte(img_input, header); // TODO: just move the cursor by 4 bytes
    CUI Nimg      = get4byte(img_input, header);
    CUI Nrow      = get4byte(img_input, header);
    CUI Ncol      = get4byte(img_input, header);
    CUI LBL_Magic = get4byte(lbl_input, header); // TODO: just move the cursor by 4 bytes
    CUI Nlbl      = get4byte(lbl_input, header);
    CUI Nbits     = Nimg * Nrow * Ncol;

    assert(Nimg == Nlbl);

    std::cout << "Number of images: "  << Nimg << std::endl;
    std::cout << "Number of rows: "    << Nrow << std::endl;
    std::cout << "Number of columns: " << Ncol << std::endl;

    CUI WIDTH = 8;    
    char * buffer = new char [WIDTH];
    char byte;

    // process the image file
    std::vector<char> img_stream;
    UI Nbytes;
    Nbytes = ceiling( Nbits, WIDTH );
    img_stream.reserve( Nbytes );
    auto count = 0u;
    for (auto i = 0u; i < Nbytes * WIDTH; i += WIDTH) // read 8 bytes, pack them into a single 8-bit char
    {   
        byte = 0u; // reset byte to all zeros. 
        img_input.read(buffer, WIDTH); // read 8 bytes
        for (auto j = 0u; j < WIDTH; j++)
        {   
            byte |= ((buffer[j] & 128) >> (7-j)) ; // MSB of byte is placed into a j-th place within a byte
            // std::cout << ((buffer[j] >> 7) ? '@':'.');
            // std::cout << ((byte >> j) ? '@':'.');
            // if (count++ % 28 == 0) std::cout << std::endl;
        }
        img_stream.push_back(byte);
    }
    // process the label file
    std::vector<char> lbl_stream;
    Nbytes = ceiling( Nimg, WIDTH );
    lbl_stream.reserve( Nbytes );
    for (auto i = 0u; i < Nimg; i += WIDTH) // read 8 bytes, pack them into a single 8-bit char
    {   
        // std::cout << std::hex << i;
        byte = 0u; // reset byte to all zeros. 
        lbl_input.read(buffer, WIDTH); // read 8 bytes
        for (auto j = 0u; j < WIDTH; j++)
        {   
            // std::cout << std::hex << (unsigned int) buffer[j];
            byte |= (buffer[j] >= 5) << j ; // false if 01234, true if 56789
        }
        // std::cout << " ---> " << (unsigned int) byte << std::endl;
        lbl_stream.push_back(byte);
    }
    // std::cout << "[BINARIZE] Number of labels : " << std::dec << lbl_stream.size() << std::endl;
    // std::cout << "[BINARIZE] Number of bytes : " << std::dec << Nbytes << std::endl;
    return std::tuple< UI, UI, UI, std::vector<char>, std::vector<char> > (Nimg, Nrow, Ncol, img_stream, lbl_stream);
}


void write_data(UI Nimg, UI Nrow, UI Ncol, std::vector<char> img_stream, std::vector<char> lbl_stream, std::string outputImages, std::string outputLabels)
{
    std::ofstream img_output( outputImages, std::ios::binary );
    assert(img_output.is_open());
    std::ofstream lbl_output( outputLabels, std::ios::binary );
    assert(lbl_output.is_open());


    // Read header information
    char * buffer4byte = new char [4];
    uintToBytes(Nimg, buffer4byte);
    lbl_output.write(buffer4byte, 4);

    lbl_output.write(lbl_stream.data(), lbl_stream.size());
    lbl_output.close(); // Close the file

    img_output.write(buffer4byte, 4);

    uintToBytes(Nrow, buffer4byte);
    img_output.write(buffer4byte, 4);

    uintToBytes(Ncol, buffer4byte);
    img_output.write(buffer4byte, 4);
    img_output.write(img_stream.data(), img_stream.size());
    img_output.close(); // Close the file
}

void write_ND_vector(const std::vector<size_t> & dims, const std::vector<bool> & data, const std::string & filename)
{
    
    std::ofstream out_file( filename, std::ios::binary );
    assert(out_file.is_open());
    // header information
    char * buffer4byte = new char [4];
    uintToBytes(dims.size(), buffer4byte); // number of dimensions
    out_file.write(buffer4byte, 4);
    size_t total_bits = 1u;
    for (const UI & dim : dims)
    {
        total_bits *= dim;
        uintToBytes(dim, buffer4byte); // dimension size
        out_file.write(buffer4byte, 4);
    }

    std::cout << "[WRITE] DIMENSIONS ARE";
    for (const UI & dim : dims)
    {
        std::cout << " " << dim;
    }
    std::cout << " " << std::endl;
    assert(total_bits == data.size());

    // process the image file
    size_t numBytes = (data.size() + 7) / 8;
    std::vector<char> byteData(numBytes); // vector to hold bytes
    for (size_t i = 0; i < data.size(); ++i) {
        if (data[i]) {
            byteData[i / 8] |= 1 << (i % 8);
        }
    }

    out_file.write(byteData.data(), byteData.size());
    out_file.close(); // Close the file
}


// the vector is flat since the number of dimensions should be known at compile time
std::tuple< std::vector<size_t>, std::vector<bool> > read_ND_vector(const std::string & filename)
{
    std::ifstream in_file( filename, std::ios::binary );
    assert(in_file.is_open());
    // Read header information
    char * header = new char [4];
    CUI Ndims = get4byte(in_file, header); 
    std::vector<size_t> dims(Ndims);
    size_t total_bits = 1u;
    std::cout << "[READ] THE DIMENSIONS ARE ";
    for (auto i = 0u; i < Ndims; i++)
    {
        dims[i] = get4byte(in_file, header);
        total_bits *= dims[i];
        std::cout << dims[i] << " ";
    }
    std::cout << std::endl << "[READ] TOTAL NUMBER OF BITS IS " << total_bits << std::endl;
    std::vector<bool> data(total_bits);

    size_t ctr = 0u;
    char byte;
    for (std::size_t i = 0; i < ceiling(total_bits, 8); ++i)
    {
        in_file.read(&byte, 1);
        for (int j = 0; j < 8; ++j)
        {   
            ctr++;
            data[8 * i + j] = (byte >> j) & 1;
            if ((8 * i + j) == total_bits) break;
        }
        // if (ctr % 10000 == 0)
        // {
        //     std::cout << "READ " << ctr << " OUT OF " << total_bits << " BITS SO FAR" << std::endl;
        // }
    }

    // process the image file
    in_file.close(); // Close the file
    return std::make_tuple(dims, data);
}

kitty::partial_truth_table data2ptt(UI Nimg, UI Nrow, UI Ncol, std::vector<bool>)
{
    kitty::partial_truth_table TT;
}

std::tuple<UI, UI, UI, std::vector<bool>, std::vector<bool>> read_data(std::string Images, std::string Labels)
{
    std::ifstream img_file( Images, std::ios::binary );
    assert(img_file.is_open());
    std::ifstream lbl_file( Labels, std::ios::binary );
    assert(lbl_file.is_open());


    // Read header information
    char * header = new char [4];
    CUI Nimg      = get4byte(img_file, header);
    CUI Nrow      = get4byte(img_file, header);
    CUI Ncol      = get4byte(img_file, header);

    std::cout << "Number of images: "  << Nimg << std::endl;
    std::cout << "Number of rows: "    << Nrow << std::endl;
    std::cout << "Number of columns: " << Ncol << std::endl;

    CUI Nlbl      = get4byte(lbl_file, header);
    CUI Nbits     = Nimg * Nrow * Ncol;

    // assert(Nimg == Nlbl);
    std::cout << "Number of labels: "  << Nlbl << std::endl;

    char byte;

    std::vector<bool> img_bits;
    img_bits.reserve(Nbits);
    unsigned long long count = 0u;
    while (img_file.read(&byte, 1)) {
        // std::cout << std::hex << 0 + byte << std::endl;
        for (auto i = 0u; i < 8; i++) {
            img_bits.push_back((bool)((byte >> i) & 1));
#ifdef DEBUG_READ
            std::cout << ((bool)((byte >> i) & 1) ? '#':'.');
            // std::cout << (bool)((byte >> i) & 1);
            if (count++ % Ncol == 0) std::cout << std::endl;
            if (count % (Ncol * Nrow) == 0) std::cout << std::endl;
#endif
        }
        // std::cout << std::endl;
    }
    img_bits.resize(Nbits);

    std::vector<bool> lbl_bits;
    lbl_bits.reserve(Nimg);
    auto ct = 0u;
    auto nbytes = 0u;
    while (lbl_file.read(&byte, 1)) {
        // std::cout << "Number of bytes read from file: " << std::dec << ++nbytes << std::endl;
        // std::cout << std::hex << ct << '\t' << (unsigned int)byte << '\t';
        for (auto i = 0u; i < 8; i++) {
            lbl_bits.push_back((bool)((byte >> i) & 1));
            // std::cout << std::hex << ((byte >> i) & 1);
            ct++;
        }
        // std::cout << std::endl;
    }
    lbl_bits.resize(Nimg);

    // auto ct = 0u;
    // for (bool bit : img_bits)
    // {
    //     if ((ct++ % 28) == 0) std::cout << std::endl;
    //     std::cout << (bit? '#':'.');
    // }


    // // Read header information
    // char * buffer4byte = new char [4];
    // uintToBytes(Nimg, buffer4byte);
    // lbl_output.write(buffer4byte, 4);
    // lbl_output.write(lbl_stream.data(), lbl_stream.size());
    lbl_file.close(); // Close the file

    // img_output.write(buffer4byte, 4);

    // uintToBytes(Nrow, buffer4byte);
    // img_output.write(buffer4byte, 4);

    // uintToBytes(Ncol, buffer4byte);
    // img_output.write(buffer4byte, 4);
    // img_output.write(img_stream.data(), img_stream.size());
    img_file.close(); // Close the file
    return std::make_tuple( Nimg, Nrow, Ncol, img_bits, lbl_bits );
}


/*
template <UI Nimg, UI Nrow, UI Ncol>
std::tuple<std::array<bool, Nimg * Nrow * Ncol>, std::array<bool, Nimg>> read_data_known_dims(const std::string Images, const std::string Labels)
{
    std::ifstream img_file( Images, std::ios::binary );
    assert(img_file.is_open());
    std::ifstream lbl_file( Labels, std::ios::binary );
    assert(lbl_file.is_open());

    // skip headers 
    img_file.seekg(12);
    lbl_file.seekg(4);
    CUI Nbits     = Nimg * Nrow * Ncol;

    char byte;
    std::array<bool, Nbits> img_bits;
    auto bitcount = 0ull;
    auto bytecount = 0ull;
    while (img_file.read(&byte, 1)) {
        for (auto i = 0u; i < 8; i++) {
            img_bits[bytecount * 8 + i] = (byte >> i) & 1;
#ifdef DEBUG_READ
            std::cout << ((bool)((byte >> i) & 1) ? '#':'.');
            // std::cout << (bool)((byte >> i) & 1);
            if (bitcount++ % Ncol == 0) std::cout << std::endl;
            if (bitcount % (Ncol * Nrow) == 0) std::cout << std::endl;
#endif
        }
        bytecount++;
    }

    bytecount = 0ull;
    std::array<bool, Nimg> lbl_bits;
    while (lbl_file.read(&byte, 1)) {
        for (auto i = 0u; i < 8; i++) {
            lbl_bits[bytecount * 8 + i] = (byte >> i) & 1;
        }
        bytecount++;
    }

    // auto ct = 0u;
    // for (bool bit : img_bits)
    // {
    //     if ((ct++ % 28) == 0) std::cout << std::endl;
    //     std::cout << (bit? '#':'.');
    // }


    // // Read header information
    // char * buffer4byte = new char [4];
    // uintToBytes(Nimg, buffer4byte);
    // lbl_output.write(buffer4byte, 4);
    // lbl_output.write(lbl_stream.data(), lbl_stream.size());
    lbl_file.close(); // Close the file

    // img_output.write(buffer4byte, 4);

    // uintToBytes(Nrow, buffer4byte);
    // img_output.write(buffer4byte, 4);

    // uintToBytes(Ncol, buffer4byte);
    // img_output.write(buffer4byte, 4);
    // img_output.write(img_stream.data(), img_stream.size());
    img_file.close(); // Close the file
    return std::make_tuple(img_bits, lbl_bits);
}
*/
}

/*
int main()
{
    std::string  inputImages = "/Users/brainkz/Documents/GitHub/machine-learning-01/data/datasets/mnist/train-images.idx3-ubyte";   // std::string outputImages = "/Users/brainkz/Documents/GitHub/machine-learning-01/data/datasets/mnist/mnist_binarized_image.data";
    std::string  inputLabels = "/Users/brainkz/Documents/GitHub/machine-learning-01/data/datasets/mnist/train-labels.idx1-ubyte";   // std::string outputLabels = "/Users/brainkz/Documents/GitHub/machine-learning-01/data/datasets/mnist/mnist_binarized_label.data";

    auto out = binarize_ubyte(inputImages, inputLabels);
}
*/