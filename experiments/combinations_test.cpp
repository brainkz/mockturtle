#include <iostream>
#include <vector>

void generate_combinations(int n, int k, int depth, std::vector<int>& current, bool has_zero) {
    if (depth == n) {
        if (has_zero) {
            for (int val : current) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }
        return;
    }

    for (int i = 0; i < k; ++i) {
        if (!current.empty() && i < current.back()) {
            continue; // Skip to maintain sorted order
        }
        current.push_back(i);
        generate_combinations(n, k, depth + 1, current, has_zero || i == 0);
        current.pop_back();
    }
}

int main() {
    int n = 4;
    int k = 5;
    std::vector<int> current;
    generate_combinations(n, k, 0, current, false);

    return 0;
}