#include <algorithm>
#include <unordered_set>
#include "sampling.h"

namespace sampling
{

    /*
        Normal sampling.
    */
    std::vector<int> normal_sample(const size_t& n, std::mt19937 gen)
    {
        std::vector<int> result;

        for (int i = 0; i < n; ++i) {
            result.emplace_back(i);
        }
        std::shuffle(result.begin(), result.end(), gen);

        return result;
    }

    /*
        Bob Floyd sampling.
    */
    std::vector<int> bob_floyd_sample(const size_t& k, const size_t& n, std::mt19937 gen)
    {
        std::unordered_set<size_t> elements(k);

        /*
            Tries to add a new random integer to elements.
            If the integer already exists in elements, r is added instead.
        */
        for (size_t r = n - k; r < n; ++r) {
            size_t v = std::uniform_int_distribution<>(1, r)(gen);
            if (!elements.insert(v).second) {
                elements.insert(r);
            }
        }

        std::vector<int> result(elements.begin(), elements.end());

        return result;
    }

    /*
        Returns a vector of k distinct integers, ranging from 0 (inclusive) to n
        (exclusive).
    */
    std::vector<int> sample_range(const size_t& k, const size_t& n, std::mt19937 gen)
    {
        if (k < n) {
            return bob_floyd_sample(k, n, gen);
        }
        return normal_sample(n, gen);
    }

}
