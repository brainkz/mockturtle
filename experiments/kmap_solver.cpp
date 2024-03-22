
#include <iostream>
#include <bitset>
#include <vector>
#include <algorithm>
#include <string>

const int NUM_VARS = 4;

struct Minterm {
    int value;
    bool used;
    std::bitset<NUM_VARS> binary;

    Minterm(int v) : value(v), used(false), binary(v) {}
};

bool compareMinterms(const Minterm& a, const Minterm& b) {
    return a.binary.count() < b.binary.count();
}

bool canCombine(const Minterm& a, const Minterm& b, int& diffIndex) {
    int diffCount = 0;

    for (int i = 0; i < NUM_VARS; i++) {
        if (a.binary[i] != b.binary[i]) {
            diffCount++;
            diffIndex = i;
        }
    }

    return diffCount == 1;
}

std::string mintermToExpression(const Minterm& minterm) {
    std::string expr = "";
    char vars[] = {'A', 'B', 'C', 'D'};

    for (int i = 0; i < NUM_VARS; i++) {
        if (minterm.binary[i]) {
            expr += vars[i];
        } else {
            expr += (vars[i] + std::string("'"));
        }
    }

    return expr;
}

int main() {
    int numMinterms;
    std::vector<Minterm> minterms;

    std::cout << "Enter the number of minterms: ";
    std::cin >> numMinterms;

    for (int i = 0; i < numMinterms; i++) {
        int value;
        std::cout << "Enter minterm #" << i + 1 << ": ";
        std::cin >> value;
        minterms.push_back(Minterm(value));
    }

    std::sort(minterms.begin(), minterms.end(), compareMinterms);

    std::vector<Minterm> primeImplicants;
    bool isReduced = false;

    while (!isReduced) {
        isReduced = true;
        std::vector<Minterm> reducedMinterms;

        for (auto it1 = minterms.begin(); it1 != minterms.end(); ++it1) {
            for (auto it2 = it1 + 1; it2 != minterms.end(); ++it2) {
                int diffIndex;

                if (canCombine(*it1, *it2, diffIndex)) {
                    isReduced = false;
                    it1->used = true;
                    it2->used = true;
                    Minterm combinedMinterm(it1->value);
                    combinedMinterm.binary[diffIndex] = 2;
                    reducedMinterms.push_back(combinedMinterm);
                }
            }
        }

        for (const auto& minterm : minterms) {
            if (!minterm.used) {
                primeImplicants.push_back(minterm);
            }
        }

        minterms = reducedMinterms;
    }

    std::cout << "Simplified expression: ";
    for (size_t i = 0; i < primeImplicants.size(); i++) {
        std::cout << mintermToExpression(primeImplicants[i]);
        if (i < primeImplicants.size() - 1) {
            std::cout << " + ";
        }
    }
    std::cout << std::endl;

    return 0;
}
