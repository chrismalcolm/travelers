/*
    Travelers Algorithm

    Authored by Christopher Malcolm (chrismalcolm).

    A set of paths is generated using the n-dimension data points in 'points'.
    They are ordered with a scale with 0 being the path of smallest distance and
    1 being the path of largest distance. On this scale, the path closest to the
    'calibration' is chosen, and the indicies of the points in the path (as they
    appear in 'points') is returned.

    Since it becomes exceedingly difficult to generate all paths for larger sets
    of points, a few shortcuts are made to make computation easier.

    * Instead of applying the algorithm to the whole set of points, the points
    are instead batched into groups no larger than MAX_BATCH and the algorithm
    is applied to each batch instead. Results are returned as if path traverses
    through each batch sequentially.

    * The maximum number of path permutations that can be generated at any stage
    is capped at MAX_PATHS. Any more than this and the set of paths will be
    randomly generated instead.

    Usage:
    Arguments are as follows:
    <count> <dimension> <calibration> [point coordinates]

    where:
    count -> int
    dimension -> int
    calibration -> float
    point coordinates -> floats - list of point coordinates

    Example usage using g++ compiler
    g++ alg.cc && ./a.out 4 2 1.0 0.0 0.0 0.0 0.1 1.0 0.0 1.0 1.0
    1 3 0 2

    Here the algorithm is performed on 4 points in 2 dimensions, calibration 1.0
    with the points (0.0, 0.0), (0.0 0.1), (1.0, 0.0), (1.0, 1.0)
*/

#include <vector>
#include <algorithm>
#include <random>
#include <iostream>
#include <unordered_set>

typedef std::vector<float> point;
typedef std::vector<std::vector<float>> matrix;

typedef std::vector<int> permutation;
typedef std::pair<std::vector<int>, float> permutation_cost;

/*
    Factorial function.
*/
int factorial(int n);

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


/*
    Single unit of the Travelers Algorithm

    A set of paths are generated. These paths start at the 'source' point, go
    through each point in 'points' and end at the 'destination' point.

    The maximum number of paths that can be generated is 'max_paths'. If the
    number of possible permutations for the paths is larger than 'max_paths',
    the paths generated will be chosen at random.

    Once the set of paths has been generated, they are sorted on a scale from
    0 to 1. The 0 indicates the path with the shortest distance and 1 indicates
    the path with the longest. The path chosen is the path which is nearest on
    the scale to 'calibration'.

    This function returns the indicies of the 'points' of the chosen path
*/
std::vector<int> travel(
    const point& source,
    const point& destination,
    const std::vector<point>& points,
    const float& calibration,
    const int& max_paths=40320
)
{
    std::random_device rd;
    std::mt19937 gen(rd());

    /*
        Create a 2d matrix which stores the euclidean distances between each of
        the points.
    */
    matrix graph = euclideanMatrix(source, points, destination);


    int n = points.size();
    int n_factorial = factorial(n);
    int permutation_count;
    std::vector<int> values;

    if (n_factorial <= max_paths) {
        for (int i = 0; i < n_factorial; ++i) {
            values.emplace_back(i);
        }
        std::shuffle(values.begin(), values.end(), gen);
        permutation_count = n_factorial;
    }
    else {
        permutation_count = max_paths;
        values = sample_range(permutation_count, n_factorial, gen);
    }


    permutation permutation;
    std::vector<permutation_cost> permutation_costs;
    for (int value: values) {

        float cost = 0;
        int current_index = 0;

        permutation = integer_to_permutation(value, n);

        for (int next_index: permutation) {
            cost += graph[current_index][next_index];
            current_index = next_index;
        }

        cost += graph[current_index][n+1];

        permutation_costs.push_back(permutation_cost(permutation, cost));
    }

    /*
        Order the permutation_costs from lowest to highest.
        Select the permutation according to the calibration.
    */
    sort(permutation_costs.begin(), permutation_costs.end(), sortByCost);

    return permutation_costs[calibration*(permutation_count-1)].first;
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

/*
     Batching unit of the Travelers Algorithm

     Points are batched the the single unit Travelers Algorithm is applied to
     each batch, returning the indicies in the order as if the path sequentially
     passes through each batch.
*/
std::vector<int> travels(std::vector<point> points, const float& calibration)
{
    int n = points.size();
    std::vector<int> batch_counts = get_batches(n, 13);

    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector<int> indicies = sample_range(n, n, gen);

    std::vector<int> results;

    point source = points[indicies[0]];

    point destination;
    std::vector<point> batch_points;
    int i = 1;
    for (int batch_count: batch_counts) {
        batch_points.clear();
        for (int j = 0; j < batch_count - 1; ++j) {
            batch_points.emplace_back(points[indicies[i++]]);
        }
        destination = points[indicies[i++ % n]];
        std::vector<int> a = travel(source, destination, batch_points, calibration);

        results.emplace_back(indicies[i-batch_count-1]);
        for (int b: a) {
            results.emplace_back(indicies[i-batch_count+b-1]);
        }

        destination = source;
    }

    return results;
}


int main(int argc, char* argv[])
{
    const int count = atoi(argv[1]);
    const int dimension = atoi(argv[2]);
    const float calibration = atof(argv[3]);

    std::vector<point> points;
    point temp;

    for (int i = 0; i < count; ++i) {
        temp.clear();
        for (int j = 0; j < dimension; ++j) {
            temp.emplace_back(atof(argv[i*dimension + j + 4]));
        }
        points.emplace_back(temp);
    }

    std::vector<int> indicies = travels(points, calibration);

    std::cout << indicies[0];
    for (int i = 1; i < indicies.size(); ++i) {
        std::cout << " " << indicies[i];
    }

    return 0;
}
