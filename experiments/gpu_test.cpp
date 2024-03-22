#include <iostream>
#include <vector>
#include <chrono>
#include <thread>
#include <random>

void outer_bitwise_or(const std::vector<int>& A, const std::vector<int>& B, std::vector<std::vector<int>>& C) {
  const int M = A.size();
  const int N = B.size();

  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      C[i][j] = A[i] | B[j];
    }
  }
}

void threaded_outer_bitwise_or(const std::vector<int>& A, const std::vector<int>& B, std::vector<std::vector<int>>& C, int start_row, int end_row) {
  const int M = A.size();
  const int N = B.size();

  for (int i = start_row; i < end_row; ++i) {
    for (int j = 0; j < N; ++j) {
      C[i][j] = A[i] | B[j];
    }
  }
}

int main() {
  const int M = 10000;
  const int N = 1000;
  std::vector<int> A(M);
  std::vector<int> B(N);
  std::vector<std::vector<int>> C(M, std::vector<int>(N));

  // Initialize input vectors with random numbers
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> dis(1, 100);
  for (int i = 0; i < M; ++i) {
    A[i] = dis(gen);
  }
  for (int j = 0; j < N; ++j) {
    B[j] = dis(gen);
  }

  // Compute outer product without threading
  auto start_time = std::chrono::high_resolution_clock::now();
  outer_bitwise_or(A, B, C);
  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> simple_duration = end_time - start_time;

  // Compute outer product with threading
  const int num_threads = 4;
  std::vector<std::thread> threads(num_threads);
  const int chunk_size = M / num_threads;

  start_time = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < num_threads; ++i) {
    const int start_row = i * chunk_size;
    const int end_row = (i == num_threads - 1) ? M : (i + 1) * chunk_size;
    threads[i] = std::thread(threaded_outer_bitwise_or, std::ref(A), std::ref(B), std::ref(C), start_row, end_row);
  }
  for (int i = 0; i < num_threads; ++i) {
    threads[i].join();
  }
  end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> threaded_duration = end_time - start_time;

  // Print results
  std::cout << "Simple outer product runtime: " << simple_duration.count() << " seconds" << std::endl;
  std::cout << "Threaded outer product runtime: " << threaded_duration.count() << " seconds" << std::endl;

  return 0;
}
