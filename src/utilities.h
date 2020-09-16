#include <vector>

typedef std::vector<float> point;
typedef std::vector<std::vector<float>> matrix;

typedef std::vector<int> permutation;
typedef std::pair<std::vector<int>, float> permutation_cost;

int factorial(int n);
bool sortByCost(const permutation_cost& a, const permutation_cost& b);
float euclideanDistance(
    const point& point_A,
    const point& point_B,
    const unsigned int& n
);
matrix euclideanMatrix(
    const point& source,
    const std::vector<point>& points,
    const point& destination
);
permutation integer_to_permutation(int value, const int& n);
std::vector<int> get_batches(int n, int m);
