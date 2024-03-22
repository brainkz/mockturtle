#include <iostream>
#include <vector>
#include <chrono>
#include <thread>
#include <algorithm>
#include <numeric>
#include <random>

constexpr unsigned int VEC_SIZE = 100'000;

void product_calc(const std::vector<unsigned int>& vec, std::vector<unsigned int>& out_vec, int start, int end) {
    for (int i = start; i < end; ++i) {
        for (int j = 0; j < i; ++j) {
            out_vec[i] *= vec[j] * vec[i];
        }
    }
}

int main() {
    std::vector<std::thread> threads;

    for (int num_threads : {1, 5, 10, 50, 100, 500}) {
        std::vector<unsigned int> vec(VEC_SIZE);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, VEC_SIZE);
        for (int i = 0; i < vec.size(); ++i) {
            vec[i] = dis(gen);
        }

        std::vector<unsigned int> out_vec(vec.size(), 1); // Vector to store output
        auto start_time = std::chrono::high_resolution_clock::now();

        // Create threads and divide the work
        int chunk_size = vec.size() / num_threads;
        int remainder = vec.size() % num_threads;
        int start = 0, end = 0;
        for (int i = 0; i < num_threads; ++i) {
            start = end;
            end = start + chunk_size + (i < remainder ? 1 : 0);
            threads.emplace_back(product_calc, std::cref(vec), std::ref(out_vec), start, end);
        }

        // Join threads
        for (auto& t : threads) {
            if (t.joinable())
            {
                t.join();
            }
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

        std::cout << "Time taken for " << num_threads << " threads: " << duration.count() << " ms" << std::endl;
    }

    return 0;
}