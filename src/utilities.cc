#include <algorithm>
#include <random>
#include "utilities.h"

/*
    Factorial function.
*/
int factorial(int n)
{
    return (n == 0 || n == 1) ? 1 : factorial(n - 1) * n;
}

/*
    Sorting function for the permutation costs.
*/
bool sortByCost(const permutation_cost& a, const permutation_cost& b)
{
    return (a.second < b.second);
}

/*
    Returns the distance between two points, in n dimensional euclidean space.
*/
float euclideanDistance(
    const point& point_A,
    const point& point_B,
    const unsigned int& n
)
{
    float squared_sum = 0.0;
    for (int i = 0; i < n; ++i) {
        squared_sum += pow(point_A[i] - point_B[i], 2);
    }
    return sqrt(squared_sum);
}

/*
    Returns a 2d matrix of euclidean distances between points.
    graph[i][j] = the euclidean distance from points[i] to points[j]
*/
matrix euclideanMatrix(const point& source, const std::vector<point>& points, const point& destination)
{
    matrix graph;

    const int count = points.size() + 2;
    const unsigned int dimension = source.size();

    std::vector<float> vector;
    for (int i = 0; i < count; ++i) {
        vector.clear();
        for (int j = 0; j < count; ++j) {
            vector.push_back(euclideanDistance(
                (i != 0 && i != count-1) ? points[i-1] : (i == 0) ? source : destination,
                (j != 0 && j != count-1) ? points[j-1] : (j == 0) ? source : destination,
                dimension
            ));
        }
        graph.push_back(vector);
    }

    return graph;
}

/*
    Returns a permutation of n integers corrosponding to the integer value
    given. The permutation is represented as a vector of integers 1 to n. The
    permutation returned will be unique mod n!.
*/
permutation integer_to_permutation(int value, const int& n)
{
    permutation permutation(n, 0);

    /*
        Sets the permutation for the given value.

        Process is similar to finding the digits of a number in base a, except
        base changes at each step.
    */
    for (int i = n; i > 0; --i) {

        // Get digit 0 <= c < i and mod the value
        int modulo = factorial(i-1);
        int c = value / modulo;
        value %= modulo;

        // Add i into the cth empty slot in the permutation
        for (int j = 0; j < n; ++j) {
            if (permutation[j] != 0) {
                continue;
            }
            if (c == 0) {
                permutation[j] = i;
                break;
            }
            c -= 1;
        }
    }

    return permutation;
}

/*
    Returns the sizes of the batches, where 'n' is the total number of elements
    and 'm' is the maximum batch size. Batches are made to be as similar size
    as possible.
*/
std::vector<int> get_batches(int n, int m)
{
    int a = ((n - 1) / m) + 1;
    int v = n / a;
    int r = n - v * a;

    std::vector<int> result(a, v);
    for (int i = 0; i < r; ++i) {
        result[i] += 1;
    }

    if (r != 0) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::shuffle(result.begin(), result.end(), gen);
    }

    return result;
}
