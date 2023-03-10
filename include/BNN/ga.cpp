
#include <iostream>
#include <random>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort
#include <vector>

typedef unsigned short  US;
typedef unsigned int    UI;
typedef const unsigned int CUI;

// int wmin = -20;
// int wmax =  20;
// int bmin = -1000;
// int bmax =  1000;
static std::uniform_int_distribution<int> BOOL_DISTR(0, 1);
static std::default_random_engine BOOL_GEN;
namespace GA
{

template <typename T>
std::vector<size_t> argsort(const std::vector<T> &v) 
{
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  std::stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
  return idx;
}
template <typename T>
void remove_at_idx(std::vector<T> & elements, std::vector<T> indices_to_remove) 
{
    // Sort the indices in descending order, so that erasing elements doesn't affect subsequent indices
    std::sort(indices_to_remove.rbegin(), indices_to_remove.rend());

    for (auto idx : indices_to_remove)
    {
        elements.erase(elements.begin() + idx);
    }
}

template<typename T>
void remove_index(std::vector<T>& vector, const std::vector<UI>& to_remove)
{
    auto vector_base = vector.begin();
    // std::vector<T>::size_type 
    size_t down_by = 0;

    for (auto iter = to_remove.cbegin(); 
                iter < to_remove.cend(); 
                iter++, down_by++)
    {
        // std::vector<T>::size_type 
        auto next = (iter + 1 == to_remove.cend() 
                                        ? vector.size() 
                                        : *(iter + 1));

        std::move(vector_base + *iter + 1, 
                vector_base + next, 
                vector_base + *iter - down_by);
    }
    vector.resize(vector.size() - to_remove.size());
}

auto remove_worst(std::vector<std::vector<int>> & pop, std::vector<UI> & num_correct, size_t n) 
{
    // Sort the indices in descending order, so that erasing elements doesn't affect subsequent indices
    
    std::cout << "Initializing " << std::endl;
    std::vector<UI> worst_idx(n);
    std::vector<UI> sorted_num(num_correct); 
    assert(sorted_num == num_correct);

    std::cout << "Partial sort " << std::endl;
    // Use partial_sort to sort the smallest n elements in ascending order.
    // std::partial_sort(sorted_num.begin(), sorted_num.begin() + n, sorted_num.end());
    std::sort(sorted_num.begin(), sorted_num.end());

    // Create a new vector to store the indices of the smallest n elements.

    std::cout << "Partial sort complete" << std::endl;
    // Find the indices of the smallest n elements in the original vector.
    for (int i = 0; i < n; i++) {
        auto it = std::find(num_correct.begin(), num_correct.end(), sorted_num[i]);
        worst_idx[i] = std::distance(num_correct.begin(), it);
    }
    std::cout << "Sorting indices" << std::endl;
    std::sort(worst_idx.begin(), worst_idx.end(), std::greater<int>());

    std::cout << "Erasing "  << std::endl;
    // remove_index(pop, worst_idx);
    std::vector<std::vector<int>> new_pop;
    std::vector<UI> new_num;
    for (auto i = 0u; i < pop.size(); i++)
    {
        if (std::find(worst_idx.begin(), worst_idx.end(), i) == worst_idx.end())
        {
            std::cout << "Erasing num " << i << std::endl;
            new_pop.push_back(pop[i]);
            new_num.push_back(num_correct[i]);
        }
    }
    pop = new_pop;
    num_correct = new_num;
    // return std::make_tuple(new_pop, new_num);

    // std::cout << "Erasing num "  << std::endl;
    // remove_index(num_correct, worst_idx);
    // auto ct = 0u;
    // for (auto idx : worst_idx)
    // {   
    //     pop.erase(pop.begin() + idx);
    //     std::cout << "Erasing num " << ct++ << std::endl;
    //     num_correct.erase(num_correct.begin() + idx);
    // }
}

template<typename T>
std::vector<std::pair<T, T>> select_pairs(const std::vector<T>& vec, const std::vector<UI>& probs, size_t n)
{
    std::vector<std::pair<T, T>> pairs;
    
    // Calculate cumulative sum of probabilities
    std::vector<UI> cum_probs(probs.size());
    std::partial_sum(probs.begin(), probs.end(), cum_probs.begin());
    unsigned int total_prob = cum_probs.back();
    
    // Initialize random number generator
    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<unsigned int> dist(0, total_prob - 1);
    
    // Generate N pairs of elements
    for (size_t i = 0; i < n; i++) {
        // Generate two random indices based on probability distribution
        unsigned int r1 = dist(rng);
        unsigned int r2 = dist(rng);
        auto it1 = std::upper_bound(cum_probs.begin(), cum_probs.end(), r1);
        auto it2 = std::upper_bound(cum_probs.begin(), cum_probs.end(), r2);
        size_t index1 = std::distance(cum_probs.begin(), it1);
        size_t index2 = std::distance(cum_probs.begin(), it2);
        
        // Add pair to result vector
        pairs.emplace_back(vec[index1], vec[index2]);
    }
    
    return pairs;
}

template<typename T>
std::vector<T> select_elements(const std::vector<T>& vec, const std::vector<unsigned int>& probs, size_t n) 
{
    std::vector<T> selected_elems;
    
    // Calculate cumulative sum of probabilities
    std::vector<unsigned int> cum_probs(probs.size());
    std::partial_sum(probs.begin(), probs.end(), cum_probs.begin());
    unsigned int total_prob = cum_probs.back();
    
    // Initialize random number generator
    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<unsigned int> dist(0, total_prob - 1);
    
    // Generate N selected elements
    for (size_t i = 0; i < n; i++) {
        // Generate random index based on probability distribution
        unsigned int r = dist(rng);
        auto it = std::upper_bound(cum_probs.begin(), cum_probs.end(), r);
        size_t index = std::distance(cum_probs.begin(), it);
        
        // Add element to result vector
        selected_elems.push_back(vec[index]);
    }
    
    return selected_elems;
}

int randbetween(int min, int max)
{
    static std::uniform_int_distribution<int> DISTR(min, max);
    static std::default_random_engine GEN;
    return DISTR( GEN );
}
// auto RAND_BOOL = []() { return BOOL_DISTR( BOOL_GEN ); };

auto randgen_funcs(int wmin, int wmax, int bmin, int bmax)
{
    static std::uniform_int_distribution<int> wdistr(wmin, wmax);
    static std::default_random_engine wgen;
    static std::uniform_int_distribution<int> bdistr(bmin, bmax);
    static std::default_random_engine bgen;
    auto WRAND = []() { return wdistr( wgen ); };
    auto BRAND = []() { return bdistr( bgen ); };
    return std::make_tuple(WRAND, BRAND);
}

auto norm_distr_funcs(const UI sigma)
{
    static std::normal_distribution<int> norm_distr(0, sigma);
    static std::default_random_engine norm_gen;
    return []() { return norm_distr( norm_gen ); };
}


std::vector<int> random_individual(UI Nweights, int wmin, int wmax, int bmin, int bmax)
{
    auto[WRAND, BRAND] = randgen_funcs(wmin, wmax, bmin, bmax);
    std::vector<int> weights(Nweights);
    std::generate(weights.begin(), weights.end(), WRAND);
    weights.push_back(BRAND());
    return weights;
}

std::vector<std::vector<int>> initial_population(UI Nspecies, UI Nweights, int wmin, int wmax, int bmin, int bmax) //assumes the bias is added at the end of the vector of weights
{   
    std::vector<std::vector<int>> pop(Nspecies);
    std::generate(pop.begin(), pop.end(), [&]() { return random_individual(Nweights, wmin, wmax, bmin, bmax); });
    return pop;
}

UI eval_thres(const std::vector<std::vector<bool>>& res_output_signals, const std::vector<int>& weights, const std::vector<bool>& labels)
{
    //TODO: try to optimize the bias locally. This way, we got one less variable to worry about
    UI num_correct = 0u;
    const int bias = weights.back();
    const int Nweights = weights.size() - 1;
    for (auto i = 0u; i < res_output_signals.size(); i++)
    {
        int total = 0ul;
        for (auto j = 0u; j < Nweights; j++)
        {
            total += res_output_signals[i][j] * weights[j];
        }
        // std::cout << "\tTotal = " << total << std::endl;
        // std::cout << "\tBias  = " << bias << std::endl;
        // std::cout << "\tLabel = " << labels[i] << std::endl;
        // std::cout << "\tCorrect = " << ((total > bias) == labels[i]) << std::endl;
        num_correct += ((total > bias) == labels[i]);
        std::cout << "IMG #" << i << ": #Correct: " << num_correct << '\r';
    }
    return num_correct;
}
std::vector<int> random_crossover(std::vector<int>& weight_vec1, std::vector<int>& weight_vec2)
{
    // random crossover
    auto Nweights = weight_vec1.size();
    std::vector<int> weights(Nweights);
    for (auto i = 0u; i < Nweights; i++)
    {
        // weights[i] = ( BOOL_DISTR( BOOL_GEN ) ) ? weight_vec1[i] : weight_vec2[i];
        weights[i] = randbetween(weight_vec1[i], weight_vec2[i]);
    }
    return weights;
}
std::vector<int> swap_crossover(std::vector<int>& weight_vec1, std::vector<int>& weight_vec2)
{
    // random crossover
    auto Nweights = weight_vec1.size();
    std::vector<int> weights(Nweights);
    for (auto i = 0u; i < Nweights; i++)
    {
        weights[i] = ( BOOL_DISTR( BOOL_GEN ) ) ? weight_vec1[i] : weight_vec2[i];
        // weights[i] = randbetween(weight_vec1[i], weight_vec2[i]);
    }
    return weights;
}
std::vector<int> mutation(std::vector<int> weights, const UI wsigma, const UI bsigma, int wmin, int wmax, int bmin, int bmax)
{
    // random mutation of single weight vector
    auto WNORM = norm_distr_funcs(wsigma);
    auto BNORM = norm_distr_funcs(bsigma);
    auto Nweights = weights.size()-1;
    for (auto i = 0u; i < Nweights; i++)
    {
        weights[i] = std::clamp(weights[i] + WNORM(), wmin, wmax); 
    }
    weights[Nweights] = std::clamp(weights[Nweights] + BNORM(), bmin, bmax); 
    return weights;
}


void ga_iteration(const std::vector<std::vector<bool>>& res_output_signals, const std::vector<bool>& labels, std::vector<std::vector<int>>& pop, std::vector<UI> & num_correct, CUI Nimg, CUI numdeaths, CUI num_mutate, CUI num_xover, CUI wsigma, CUI bsigma, const int wmin, const int wmax, const int bmin, const int bmax)
{
    UI Nspecies = pop.size();
    // std::vector<UI> num_correct(Nspecies, Nspecies+1);
    std::cout << "Evaluating the fitness" << std::endl;
    for (auto i = 0u; i < Nspecies; i++)
    {
        //TODO: change to using some kind of flag
        if (num_correct[i] > Nimg)
        {
            num_correct[i] = eval_thres(res_output_signals, pop[i], labels );
        }
    }
    std::cout << std::endl;
    for (auto i = 0u; i < Nspecies; i++)
    {
        std::cout << num_correct[i] << ' ';
    }
    std::cout << std::endl;
    
    // Remove weakest species
    std::cout << "Removing the worst" << std::endl;
    remove_worst(pop, num_correct, numdeaths);
    // pop = pair.at(0);
    // num_correct = pair[1];
    // auto[new_pop, new_num]
    // Select <num_xover> pairs for crossover
    std::cout << "Crossover" << std::endl;
    std::vector<std::vector<int>> new_species;
    for (auto[w1, w2] : select_pairs(pop, num_correct, num_xover))
    {
        pop.push_back(random_crossover(w1, w2));
        num_correct.push_back(Nspecies + 1);
    }
    std::cout << "Mutation" << std::endl;
    // Select <num_mutate> species for mutation
    for (auto w : select_elements(pop, num_correct, num_xover))
    {
        pop.push_back(mutation(w, wsigma, bsigma, wmin, wmax, bmin, bmax));
        num_correct.push_back(Nspecies + 1);
    }
    std::cout << "Iteration complete" << std::endl;
}

}



// int main()
// {
//     std::cout << "size of US   = " << sizeof(US)   << std::endl;
//     std::cout << "size of UI   = " << sizeof(UI)   << std::endl;
//     std::cout << "size of UL   = " << sizeof(UL)   << std::endl;
//     std::cout << "size of char = " << sizeof(char) << std::endl;
// }