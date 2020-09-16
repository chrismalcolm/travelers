#include <random>
#include <vector>

namespace sampling
{
    std::vector<int> normal_sample(const size_t& n, std::mt19937 gen);
    std::vector<int> bob_floyd_sample(const size_t& k, const size_t& n, std::mt19937 gen);
    std::vector<int> sample_range(const size_t& k, const size_t& n, std::mt19937 gen);
}
