#include <algorithm>
// #include <execution>
#include <vector>

int main() {
  std::vector<int> v = {1, 2, 3, 4, 5};

  std::for_each(std::execution::par, v.begin(), v.end(), [](int& x) {
    x *= 2;
  });

  return 0;
}

/*
#include <iostream>
#include <thread>
#include <vector>
#include <algorithm>

const int kNumThreads = 4;

void ProcessChunk(std::vector<int>::iterator begin, std::vector<int>::iterator end) {
  std::sort(begin, end);
}

int main() {
  std::vector<int> data{5, 3, 8, 1, 9, 2, 7, 6, 4, 0};

  // Determine the size of each chunk
  const int chunk_size = data.size() / kNumThreads;

  // Create a vector of threads
  std::vector<std::thread> threads;

  // Start the threads to process each chunk of data
  auto begin = data.begin();
  auto end = begin + chunk_size;
  for (int i = 0; i < kNumThreads; ++i) {
    if (i == kNumThreads - 1) {
      end = data.end();
    }
    threads.emplace_back(ProcessChunk, begin, end);
    begin = end;
    end = begin + chunk_size;
  }

  // Join the threads
  for (auto& t : threads) {
    t.join();
  }

  // Merge the sorted chunks
  std::inplace_merge(data.begin(), data.begin() + chunk_size, data.end());

  // Print the sorted array
  for (auto x : data) {
    std::cout << x << " ";
  }
  std::cout << std::endl;

  return 0;
}
*/

/*
#include <iostream>
#include <thread>
#include <mutex>
#include <unordered_map>

// Global unordered_map
std::unordered_map<int, int> global_map;

// Mutex to protect the global_map
std::mutex mtx;

// Function to increment the value of a key in the global_map
void incrementValue(int key, int increment) {
    // Lock the mutex to protect the global_map
    std::lock_guard<std::mutex> lock(mtx);

    // Increment the value of the key in the global_map
    global_map[key] += increment;
}

int main() {
    // Initialize the global_map with some keys and values
    global_map[1] = 1;
    global_map[2] = 2;
    global_map[3] = 3;

    // Create two threads to increment the value of the key 1 in the global_map
    std::thread t1(incrementValue, 1, 1);
    std::thread t2(incrementValue, 1, 2);

    // Wait for the threads to finish
    t1.join();
    t2.join();

    // Print the final value of the key 1 in the global_map
    std::cout << "Final value of key 1 in global_map: " << global_map[1] << std::endl;

    return 0;
}
*/