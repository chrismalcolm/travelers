# Travellers Algorithm
C++ implementation of the Travellers Algorithm.

Authored by Christopher Malcolm (chrismalcolm).

## Description
A set of paths is generated using the n-dimension data points in 'points'.
They are ordered with a scale with 0 being the path of smallest distance and
1 being the path of largest distance. On this scale, the path closest to the
'calibration' is chosen, and the indices of the points in the path (as they
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

## Usage:
Arguments are as follows:

`<count> <dimension> <calibration> [point coordinates]`

where:
* count -> int
* dimension -> int
* calibration -> float
* point coordinates -> floats - list of point coordinates

Example usage using g++ compiler

`g++ -L/usr/include src/sampling.cc src/utilities.cc src/algorithm.cc -o algorithm -lboost_filesystem -lboost_system -lboost_thread && ./algorithm 4 2 1.0 0.0 0.0 0.0 0.1 1.0 0.0 1.0 1.0
1 3 0 2`

Here the algorithm is performed on 4 points in 2 dimensions, calibration 1.0 with the points (0.0, 0.0), (0.0 0.1), (1.0, 0.0), (1.0, 1.0). The directory `/usr/include` should be replaced with the location of the boost library.
