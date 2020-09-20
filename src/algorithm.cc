#include <algorithm>
#include <iostream>
#include <random>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/thread/thread.hpp>
#include "sampling.h"
#include "utilities.h"

const int MAX_BATCH = 13;
const int MAX_PATHS = 40320;

/*
    Single unit of the Travellers Algorithm

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
    const int& max_paths=MAX_PATHS
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
        values = sampling::sample_range(permutation_count, n_factorial, gen);
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


std::vector<std::pair<int, permutation>> indexed_permutation;
boost::mutex mutex;

void shall(
    int& index,
    const point& source,
    const point& destination,
    const std::vector<point>& points,
    const float& calibration
)
{
    permutation perm = travel(source, destination, points, calibration);

    mutex.lock();
    indexed_permutation.emplace_back(index, perm);
    mutex.unlock();
}

/*
     Batching unit of the Travellers Algorithm

     Points are batched the the single unit Travellers Algorithm is applied to
     each batch, returning the indicies in the order as if the path sequentially
     passes through each batch.
*/
std::vector<int> travels(std::vector<point> points, const float& calibration)
{
    int n = points.size();
    std::vector<int> batch_counts = get_batches(n, MAX_BATCH);

    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector<int> indicies = sampling::sample_range(n, n, gen);
    std::vector<int> results;

    point source = points[indicies[0]];
    point destination;

    std::vector<point> batch_points;
    int i = 1;

    std::vector<boost::thread> threads;

    for (int k = 0; k < batch_counts.size(); ++k) {

        batch_points.clear();
        for (int j = 0; j < batch_counts[k] - 1; ++j) {
            batch_points.emplace_back(points[indicies[i++]]);
        }

        destination = points[indicies[i++ % n]];

        threads.emplace_back(std::move(boost::thread(
            shall,
            i-batch_counts[k]-1,
            source,
            destination,
            batch_points,
            calibration
        )));

        destination = source;
    }

    for (int i = 0; i < threads.size(); ++i) {
        threads[i].join();
    }

    for (std::pair<int, std::vector<int>> pair: indexed_permutation) {
        results.emplace_back(indicies[pair.first]);
        for (int num: pair.second) {
            results.emplace_back(indicies[pair.first+num]);
        }
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
