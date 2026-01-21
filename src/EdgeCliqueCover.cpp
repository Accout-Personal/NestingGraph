#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <set>
#include <algorithm>
#include <functional>
#include <chrono>
#include <string>
#include <numeric>
#include <string>
#include<filesystem>

#include <direct.h>
#define GetCurrentDir _getcwd

using namespace std;

std::string getCurrentDirectory() {
    char buffer[FILENAME_MAX];
    GetCurrentDir(buffer, FILENAME_MAX);
    return std::string(buffer);
}

// Graph representation using adjacency matrix
class Graph {
private:
    int numVertices;
    std::vector<std::vector<bool>> adjacencyMatrix;

public:
    // Constructor
    Graph(int n) : numVertices(n) {
        adjacencyMatrix.resize(n, std::vector<bool>(n, false));
    }

    // Add an edge to the graph
    void addEdge(int u, int v) {
        if (u != v) { // No self-loops
            adjacencyMatrix[u][v] = true;
            adjacencyMatrix[v][u] = true;
        }
    }

    // Check if there is an edge between two vertices
    bool hasEdge(int u, int v) const {
        return adjacencyMatrix[u][v];
    }

    // Get number of vertices
    int getNumVertices() const {
        return numVertices;
    }

    // Get neighbors of a vertex
    std::vector<int> getNeighbors(int v) const {
        std::vector<int> neighbors;
        for (int i = 0; i < numVertices; ++i) {
            if (adjacencyMatrix[v][i]) {
                neighbors.push_back(i);
            }
        }
        return neighbors;
    }

    // Get the degree of a vertex
    int getDegree(int v) const {
        int degree = 0;
        for (int i = 0; i < numVertices; ++i) {
            if (adjacencyMatrix[v][i]) {
                degree++;
            }
        }
        return degree;
    }

    // Get edges of the graph
    std::vector<std::pair<int, int>> getEdges() const {
        std::vector<std::pair<int, int>> edges;
        for (int i = 0; i < numVertices; ++i) {
            for (int j = i + 1; j < numVertices; ++j) {
                if (adjacencyMatrix[i][j]) {
                    edges.push_back({ i, j });
                }
            }
        }
        return edges;
    }

    // Get the number of edges in the graph
    int getNumEdges() const {
        int count = 0;
        for (int i = 0; i < numVertices; ++i) {
            for (int j = i + 1; j < numVertices; ++j) {
                if (adjacencyMatrix[i][j]) {
                    count++;
                }
            }
        }
        return count;
    }

    // Get the density of the graph
    double getDensity() const {
        long n = numVertices;
        long m = getNumEdges();
        return (2.0 * m) / (n * (n - 1));
    }
};

// Clique representation
class Clique {
private:
    std::set<int> vertices;

public:
    // Constructor
    Clique() {}

    // Add a vertex to the clique
    void addVertex(int v) {
        vertices.insert(v);
    }

    // Check if vertex is in the clique
    bool containsVertex(int v) const {
        return vertices.find(v) != vertices.end();
    }

    // Get all vertices in the clique
    std::set<int> getVertices() const {
        return vertices;
    }

    // Get size of the clique
    int size() const {
        return vertices.size();
    }

    // Check if the clique is valid (all vertices are connected)
    bool isValid(const Graph& graph) const {
        for (int u : vertices) {
            for (int v : vertices) {
                if (u != v && !graph.hasEdge(u, v)) {
                    return false;
                }
            }
        }
        return true;
    }
};

// Edge Clique Cover representation
class EdgeCliqueCover {
private:
    std::vector<Clique> cliques;

public:
    // Constructor
    EdgeCliqueCover() {}

    // Add a clique to the cover
    void addClique(const Clique& clique) {
        cliques.push_back(clique);
    }

    // Get all cliques in the cover
    std::vector<Clique> getCliques() const {
        return cliques;
    }

    // Get size of the cover (number of cliques)
    int size() const {
        return cliques.size();
    }

    // Check if the cover is valid (all edges are covered)
    bool isValid(const Graph& graph) const {
        std::vector<std::vector<bool>> covered(graph.getNumVertices(),
            std::vector<bool>(graph.getNumVertices(), false));

        // Mark edges covered by cliques
        for (const Clique& clique : cliques) {
            std::set<int> vertices = clique.getVertices();
            for (auto u = vertices.begin(); u != prev(vertices.end()); ++u) {
                for (auto v = next(u); v != vertices.end(); ++v) {
                    if (*u != *v) {
                        covered[*u][*v] = true;
                        covered[*v][*u] = true;
                    }
                }
            }
        }

        // Check if all edges are covered
        for (int i = 0; i < graph.getNumVertices(); ++i) {
            for (int j = i + 1; j < graph.getNumVertices(); ++j) {
                if (graph.hasEdge(i, j) && !covered[i][j]) {
                    return false;
                }
            }
        }

        return true;
    }

    void writeIntoFile(const string& filename)
    {
        //remove the output file if already exists
        try {
            if (filesystem::remove(filename)) {
                cout << "Old clique data has been successfully removed \n" << endl;
            }
        }
        catch (const std::filesystem::filesystem_error& e) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
        std::ofstream outfileClique(filename, std::ios::app);
        for (Clique clique : cliques)
        {
            set<int> vertexes = clique.getVertices();
            for (int vertex : vertexes)
            {
                outfileClique << vertex << " ";
            }
            outfileClique << '\n';
        }

        outfileClique.close();
        cout << "save to file complete \n";
    }

};

// Helper function to intersect two vectors of integers
std::vector<int> intersectVectors(const std::vector<int>& v1, const std::vector<int>& v2) {
    std::vector<int> result;
    std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), std::back_inserter(result));
    return result;
}

// Implementation of Sequence1 heuristic (Algorithm 1 in the paper)
EdgeCliqueCover sequence1Heuristic(const Graph& graph, const std::vector<int>& sequence) {
    EdgeCliqueCover cover;
    int n = graph.getNumVertices();

    // Initialize uncovered edge matrix
    std::vector<std::vector<bool>> uncoveredEdge(n, std::vector<bool>(n, false));
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (graph.hasEdge(i, j)) {
                uncoveredEdge[i][j] = true;
                uncoveredEdge[j][i] = true;
            }
        }
    }

    // Initialize uncovered degree for each vertex
    std::vector<int> uncoveredDegree(n, 0);
    for (int i = 0; i < n; ++i) {
        uncoveredDegree[i] = graph.getDegree(i);
    }

    // Process vertices in the given sequence
    for (int i : sequence) {
        while (uncoveredDegree[i] > 0) {
            // Create a new clique with vertex i
            Clique clique;
            clique.addVertex(i);

            // Initialize iterative neighborhood W
            auto neighbors = graph.getNeighbors(i);
            std::sort(neighbors.begin(), neighbors.end());

            // Initialize local uncovered degree
            std::vector<int> localUncoveredDegree(n, 0);
            for (int j : neighbors) {
                localUncoveredDegree[j] = uncoveredEdge[i][j] ? 1 : 0;
            }

            // Find vertex u that covers maximum uncovered edges
            auto it = std::max_element(neighbors.begin(), neighbors.end(),
                [&localUncoveredDegree](int a, int b) {
                    return localUncoveredDegree[a] < localUncoveredDegree[b];
                });

            while (it != neighbors.end() && localUncoveredDegree[*it] > 0) {
                int u = *it;

                // Cover edges between u and vertices in the clique
                for (int j : clique.getVertices()) {
                    if (uncoveredEdge[u][j]) {
                        uncoveredEdge[u][j] = false;
                        uncoveredEdge[j][u] = false;
                        uncoveredDegree[u]--;
                        uncoveredDegree[j]--;
                    }
                }

                // Add u to the clique
                clique.addVertex(u);

                // Update iterative neighborhood W
                std::vector<int> newNeighbors = graph.getNeighbors(u);
                std::sort(newNeighbors.begin(), newNeighbors.end());
                std::vector<int> intersection;
                std::set_intersection(neighbors.begin(), neighbors.end(),
                    newNeighbors.begin(), newNeighbors.end(),
                    std::back_inserter(intersection));
                neighbors = intersection;

                // Update local uncovered degree
                for (int j : neighbors) {
                    if (uncoveredEdge[u][j]) {
                        localUncoveredDegree[j]++;
                    }
                }

                // Find next vertex to add
                it = std::max_element(neighbors.begin(), neighbors.end(),
                    [&localUncoveredDegree](int a, int b) {
                        return localUncoveredDegree[a] < localUncoveredDegree[b];
                    });

                if (it != neighbors.end() && localUncoveredDegree[*it] == 0) {
                    it = neighbors.end();
                }
            }

            // Add the constructed clique to the cover
            cover.addClique(clique);
        }
    }

    return cover;
}

// Implementation of Sequence2 heuristic (extended Algorithm 1 in the paper)
EdgeCliqueCover sequence2Heuristic(const Graph& graph, const std::vector<int>& sequence) {
    EdgeCliqueCover cover;
    int n = graph.getNumVertices();

    // Initialize uncovered edge matrix
    std::vector<std::vector<bool>> uncoveredEdge(n, std::vector<bool>(n, false));
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (graph.hasEdge(i, j)) {
                uncoveredEdge[i][j] = true;
                uncoveredEdge[j][i] = true;
            }
        }
    }

    // Initialize uncovered degree for each vertex
    std::vector<int> uncoveredDegree(n, 0);
    for (int i = 0; i < n; ++i) {
        uncoveredDegree[i] = graph.getDegree(i);
    }

    // Process vertices in the given sequence
    for (int i : sequence) {
        while (uncoveredDegree[i] > 0) {
            // Create a new clique with vertex i
            Clique clique;
            clique.addVertex(i);

            // Initialize iterative neighborhood W
            auto neighbors = graph.getNeighbors(i);
            std::sort(neighbors.begin(), neighbors.end());

            // Initialize local uncovered degree
            std::vector<int> localUncoveredDegree(n, 0);
            for (int j : neighbors) {
                localUncoveredDegree[j] = uncoveredEdge[i][j] ? 1 : 0;
            }

            // Find vertex u that covers maximum uncovered edges
            auto it = std::max_element(neighbors.begin(), neighbors.end(),
                [&localUncoveredDegree](int a, int b) {
                    return localUncoveredDegree[a] < localUncoveredDegree[b];
                });

            while (it != neighbors.end() && localUncoveredDegree[*it] > 0) {
                int u = *it;

                // Cover edges between u and vertices in the clique
                for (int j : clique.getVertices()) {
                    if (uncoveredEdge[u][j]) {
                        uncoveredEdge[u][j] = false;
                        uncoveredEdge[j][u] = false;
                        uncoveredDegree[u]--;
                        uncoveredDegree[j]--;
                    }
                }

                // Add u to the clique
                clique.addVertex(u);

                // Update iterative neighborhood W
                std::vector<int> newNeighbors = graph.getNeighbors(u);
                std::sort(newNeighbors.begin(), newNeighbors.end());
                std::vector<int> intersection;
                std::set_intersection(neighbors.begin(), neighbors.end(),
                    newNeighbors.begin(), newNeighbors.end(),
                    std::back_inserter(intersection));
                neighbors = intersection;

                // Update local uncovered degree
                for (int j : neighbors) {
                    if (uncoveredEdge[u][j]) {
                        localUncoveredDegree[j]++;
                    }
                }

                // Find next vertex to add
                it = std::max_element(neighbors.begin(), neighbors.end(),
                    [&localUncoveredDegree](int a, int b) {
                        return localUncoveredDegree[a] < localUncoveredDegree[b];
                    });

                if (it != neighbors.end() && localUncoveredDegree[*it] == 0) {
                    // SEQ2 extension: Try to insert the uncovered edge
                    // that maximizes the minimum uncovered degree of its vertices
                    std::pair<int, int> bestEdge = { -1, -1 };
                    int bestMinDegree = -1;

                    for (size_t idx1 = 0; idx1 < neighbors.size(); ++idx1) {
                        for (size_t idx2 = idx1 + 1; idx2 < neighbors.size(); ++idx2) {
                            int u1 = neighbors[idx1];
                            int u2 = neighbors[idx2];
                            if (uncoveredEdge[u1][u2]) {
                                int minDegree = std::min(uncoveredDegree[u1], uncoveredDegree[u2]);
                                if (minDegree > bestMinDegree) {
                                    bestMinDegree = minDegree;
                                    bestEdge = { u1, u2 };
                                }
                            }
                        }
                    }

                    if (bestEdge.first != -1) {
                        // Add the best edge
                        int u = bestEdge.first;
                        int v = bestEdge.second;

                        // Cover the edge
                        uncoveredEdge[u][v] = false;
                        uncoveredEdge[v][u] = false;
                        uncoveredDegree[u]--;
                        uncoveredDegree[v]--;

                        // Add vertices to clique
                        clique.addVertex(u);
                        clique.addVertex(v);

                        // Update iterative neighborhood W
                        std::vector<int> newNeighborsU = graph.getNeighbors(u);
                        std::vector<int> newNeighborsV = graph.getNeighbors(v);
                        std::sort(newNeighborsU.begin(), newNeighborsU.end());
                        std::sort(newNeighborsV.begin(), newNeighborsV.end());

                        std::vector<int> intersection1;
                        std::set_intersection(neighbors.begin(), neighbors.end(),
                            newNeighborsU.begin(), newNeighborsU.end(),
                            std::back_inserter(intersection1));

                        std::vector<int> intersection2;
                        std::set_intersection(intersection1.begin(), intersection1.end(),
                            newNeighborsV.begin(), newNeighborsV.end(),
                            std::back_inserter(intersection2));

                        neighbors = intersection2;

                        // Update local uncovered degree
                        for (int j : neighbors) {
                            if (uncoveredEdge[u][j]) {
                                localUncoveredDegree[j]++;
                            }
                            if (uncoveredEdge[v][j]) {
                                localUncoveredDegree[j]++;
                            }
                        }

                        // Continue with the loop
                        it = std::max_element(neighbors.begin(), neighbors.end(),
                            [&localUncoveredDegree](int a, int b) {
                                return localUncoveredDegree[a] < localUncoveredDegree[b];
                            });

                        if (it != neighbors.end() && localUncoveredDegree[*it] == 0) {
                            it = neighbors.end();
                        }
                    }
                    else {
                        it = neighbors.end();
                    }
                }
            }

            // Add the constructed clique to the cover
            cover.addClique(clique);
        }
    }

    return cover;
}

// Implementation of Maximum1 heuristic (similar to Sequence1 but with vertex ordering by max degree)
EdgeCliqueCover maximum1Heuristic(const Graph& graph) {
    EdgeCliqueCover cover;
    int n = graph.getNumVertices();

    cout << "initialization \n";
    // Initialize uncovered edge matrix
    vector<vector<bool>> uncoveredEdge(n);
    for (int i = 0; i < n; i++) {
        uncoveredEdge[i].resize(i + 1, false); // Initialize with zeros
    }

    // Helper function to get/set value (ensuring triangular access)
    auto getUncoveredEdge = [&uncoveredEdge](int u, int v) -> bool {
        if (u < v) std::swap(u, v); // Ensure u >= v for lower triangular
        return uncoveredEdge[u][v];
        };

    auto setUncoveredEdge = [&uncoveredEdge](int u, int v, bool value) {
        if (u < v) std::swap(u, v); // Ensure u >= v for lower triangular
        uncoveredEdge[u][v] = value;
        };

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (graph.hasEdge(i, j)) {
                setUncoveredEdge(i, j, true);
            }
        }
    }

    // Initialize uncovered degree for each vertex
    std::vector<int> uncoveredDegree(n, 0);
    for (int i = 0; i < n; ++i) {
        uncoveredDegree[i] = graph.getDegree(i);
    }

    // Set of vertices with uncovered edges
    std::set<int> verticesWithUncoveredEdges;
    for (int i = 0; i < n; ++i) {
        if (uncoveredDegree[i] > 0) {
            verticesWithUncoveredEdges.insert(i);
        }
    }

    cout << "main loop \n";
    while (!verticesWithUncoveredEdges.empty()) {
        // Select vertex with maximum uncovered degree
        int i = *std::max_element(verticesWithUncoveredEdges.begin(), verticesWithUncoveredEdges.end(),
            [&uncoveredDegree](int a, int b) {
                return uncoveredDegree[a] < uncoveredDegree[b];
            });

        
        while (uncoveredDegree[i] > 0) {
            // Create a new clique with vertex i
            //cout << uncoveredDegree[i] << "\n";
            Clique clique;
            clique.addVertex(i);

            // Initialize iterative neighborhood W
            auto neighbors = graph.getNeighbors(i);
            std::sort(neighbors.begin(), neighbors.end());

            // Initialize local uncovered degree
            std::vector<int> localUncoveredDegree(n, 0);
            for (int j : neighbors) {
                localUncoveredDegree[j] = getUncoveredEdge(i,j) ? 1 : 0;
            }

            
            // Find vertex u that covers maximum uncovered edges
            auto it = std::max_element(neighbors.begin(), neighbors.end(),
                [&localUncoveredDegree](int a, int b) {
                    return localUncoveredDegree[a] < localUncoveredDegree[b];
                });
            
            while (it != neighbors.end() && localUncoveredDegree[*it] > 0) {
                int u = *it;
                
                // Cover edges between u and vertices in the clique
                for (int j : clique.getVertices()) {
                    if (getUncoveredEdge(u,j)) {
                        setUncoveredEdge(u,j,false);

                        uncoveredDegree[u]--;
                        uncoveredDegree[j]--;
                        
                    }
                }

                // Add u to the clique
                clique.addVertex(u);

                // Update iterative neighborhood W
                std::vector<int> newNeighbors = graph.getNeighbors(u);
                std::sort(newNeighbors.begin(), newNeighbors.end());
                std::vector<int> intersection;
                std::set_intersection(neighbors.begin(), neighbors.end(),
                    newNeighbors.begin(), newNeighbors.end(),
                    std::back_inserter(intersection));
                neighbors = intersection;

                // Update local uncovered degree
                for (int j : neighbors) {
                    if (getUncoveredEdge(u, j)) {
                        localUncoveredDegree[j]++;
                    }
                }

                // Find next vertex to add
                it = std::max_element(neighbors.begin(), neighbors.end(),
                    [&localUncoveredDegree](int a, int b) {
                        return localUncoveredDegree[a] < localUncoveredDegree[b];
                    });

                if (it != neighbors.end() && localUncoveredDegree[*it] == 0) {
                    it = neighbors.end();
                }
            }

            // Add the constructed clique to the cover
            
            cover.addClique(clique);
        }

        // Update the set of vertices with uncovered edges
        verticesWithUncoveredEdges.clear();
        for (int i = 0; i < n; ++i) {
            if (uncoveredDegree[i] > 0) {
                verticesWithUncoveredEdges.insert(i);
            }
        }
    }

    return cover;
}


// Implementation of Minimum1 heuristic (similar to Maximum1 but with vertex ordering by min degree)
EdgeCliqueCover minimum1Heuristic(const Graph& graph) {
    EdgeCliqueCover cover;
    int n = graph.getNumVertices();

    // Initialize uncovered edge matrix
    //std::vector<std::vector<bool>> uncoveredEdge(n, std::vector<bool>(n, false));
    vector<vector<bool>> uncoveredEdge(n);
    for (int i = 0; i < n; i++) {
        uncoveredEdge[i].resize(i + 1, false); // Initialize with zeros
    }

    // Helper function to get/set alpha value (ensuring triangular access)
    auto getUncoveredEdge = [&uncoveredEdge](int u, int v) -> bool {
        if (u < v) std::swap(u, v); // Ensure i >= j for lower triangular
        return uncoveredEdge[u][v];
        };

    auto setUncoveredEdge = [&uncoveredEdge](int u, int v, bool value) {
        if (u < v) std::swap(u, v); // Ensure u >= v for lower triangular
        uncoveredEdge[u][v] = value;
        };

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (graph.hasEdge(i, j)) {
                setUncoveredEdge(i, j,true);
            }
        }
    }

    // Initialize uncovered degree for each vertex
    std::vector<int> uncoveredDegree(n, 0);
    for (int i = 0; i < n; ++i) {
        uncoveredDegree[i] = graph.getDegree(i);
    }

    // Set of vertices with uncovered edges
    std::set<int> verticesWithUncoveredEdges;
    for (int i = 0; i < n; ++i) {
        if (uncoveredDegree[i] > 0) {
            verticesWithUncoveredEdges.insert(i);
        }
    }

    while (!verticesWithUncoveredEdges.empty()) {
        // Select vertex with minimum uncovered degree (but > 0)
        int i = *std::min_element(verticesWithUncoveredEdges.begin(), verticesWithUncoveredEdges.end(),
            [&uncoveredDegree](int a, int b) {
                if (uncoveredDegree[a] == 0) return false;
                if (uncoveredDegree[b] == 0) return true;
                return uncoveredDegree[a] < uncoveredDegree[b];
            });

        while (uncoveredDegree[i] > 0) {
            // Create a new clique with vertex i
            Clique clique;
            clique.addVertex(i);

            // Initialize iterative neighborhood W
            auto neighbors = graph.getNeighbors(i);
            std::sort(neighbors.begin(), neighbors.end());

            // Initialize local uncovered degree
            std::vector<int> localUncoveredDegree(n, 0);
            for (int j : neighbors) {
                localUncoveredDegree[j] = getUncoveredEdge(i,j) ? 1 : 0;
            }

            // Find vertex u that covers maximum uncovered edges
            auto it = std::max_element(neighbors.begin(), neighbors.end(),
                [&localUncoveredDegree](int a, int b) {
                    return localUncoveredDegree[a] < localUncoveredDegree[b];
                });

            while (it != neighbors.end() && localUncoveredDegree[*it] > 0) {
                int u = *it;

                // Cover edges between u and vertices in the clique
                for (int j : clique.getVertices()) {
                    if (getUncoveredEdge(u, j)) {
                        setUncoveredEdge(u, j, false);
                        uncoveredDegree[u]--;
                        uncoveredDegree[j]--;
                    }
                }

                // Add u to the clique
                clique.addVertex(u);

                // Update iterative neighborhood W
                std::vector<int> newNeighbors = graph.getNeighbors(u);
                std::sort(newNeighbors.begin(), newNeighbors.end());
                std::vector<int> intersection;
                std::set_intersection(neighbors.begin(), neighbors.end(),
                    newNeighbors.begin(), newNeighbors.end(),
                    std::back_inserter(intersection));
                neighbors = intersection;

                // Update local uncovered degree
                for (int j : neighbors) {
                    if (getUncoveredEdge(i, j)) {
                        localUncoveredDegree[j]++;
                    }
                }

                // Find next vertex to add
                it = std::max_element(neighbors.begin(), neighbors.end(),
                    [&localUncoveredDegree](int a, int b) {
                        return localUncoveredDegree[a] < localUncoveredDegree[b];
                    });

                if (it != neighbors.end() && localUncoveredDegree[*it] == 0) {
                    it = neighbors.end();
                }
            }

            // Add the constructed clique to the cover
            cover.addClique(clique);
        }

        // Update the set of vertices with uncovered edges
        verticesWithUncoveredEdges.clear();
        for (int i = 0; i < n; ++i) {
            if (uncoveredDegree[i] > 0) {
                verticesWithUncoveredEdges.insert(i);
            }
        }
    }

    return cover;
}

// Implementation of Kou et al. post-processing (Algorithm 3 in the paper)
EdgeCliqueCover kouPostProcessing(const Graph& graph, const EdgeCliqueCover& initialCover) {
    cout << "kou post processing...\n";
    EdgeCliqueCover improvedCover;
    int n = graph.getNumVertices();

    // Initialize edge coefficients (number of cliques covering each edge)
    std::vector<std::vector<int>> alpha(n);
    for (int i = 0; i < n; i++) {
        alpha[i].resize(i + 1, 0); // Initialize with zeros
    }

    // Helper function to get/set alpha value (ensuring triangular access)
    auto getAlpha = [&alpha](int i, int j) -> int& {
        if (i < j) std::swap(i, j); // Ensure i >= j for lower triangular
        return alpha[i][j];
        };

    std::vector<Clique> cliques = initialCover.getCliques();

    // Calculate initial alpha values
    for (const Clique& clique : cliques) {
        std::set<int> vertices = clique.getVertices();
        for (auto u = vertices.begin(); u != prev(vertices.end());++u) {
            for (auto v = next(u); v!=vertices.end();++v) {
                if (*u != *v && graph.hasEdge(*u, *v)) {
                    getAlpha(*u, *v)++;
                    //cout << "alpha value of: " << *u << " " << *v << " " << alpha[*u][*v] << " " << alpha[*v][*u] << "\n";
                }
            }
        }
    }

    // Process each clique
    for (const Clique& clique : cliques) {
        std::set<int> vertices = clique.getVertices();
        bool redundant = true;

        // Check if the clique is redundant
        for (auto u = vertices.begin(); u != prev(vertices.end()); ++u) {
            for (auto v = next(u); v != vertices.end(); ++v) {
                if (*u != *v && graph.hasEdge(*u, *v) && getAlpha(*u, *v) == 1) {
                    redundant = false;
                    break;
                }
            }
            if (!redundant) { 
                break;
            };
        }

        if (!redundant) {
            // Clique is not redundant, add it to the improved cover
            improvedCover.addClique(clique);
        }
        else {
            // Clique is redundant, update alpha values
            for (auto u = vertices.begin(); u != prev(vertices.end()); ++u) {
                for (auto v = next(u); v != vertices.end(); ++v) {
                    if (*u != *v && graph.hasEdge(*u, *v)) {
                        getAlpha(*u, *v)--;
                        //cout << "alpha value of: " << *u << " " << *v << " " << alpha[*u][*v] << " " << alpha[*v][*u] << "\n";
                    }
                }
            }
        }
    }

    return improvedCover;
}

// Implementation of Expansion algorithm (Algorithm 4 in the paper)
Clique expandClique(const Graph& graph, const Clique& initialClique) {

    Clique expandedClique = initialClique;

    // Get the current vertices in the clique
    std::set<int> vertices = expandedClique.getVertices();

    // Compute initial intersection of neighborhoods
    std::vector<int> intersection;
    bool first = true;

    for (int v : vertices) {
        std::vector<int> neighbors = graph.getNeighbors(v);
        std::sort(neighbors.begin(), neighbors.end());

        if (first) {
            intersection = neighbors;
            first = false;
        }
        else {
            std::vector<int> temp;
            std::set_intersection(intersection.begin(), intersection.end(),
                neighbors.begin(), neighbors.end(),
                std::back_inserter(temp));
            intersection = temp;
        }
    }

    // Remove vertices already in the clique
    std::vector<int> candidates;
    for (int v : intersection) {
        if (vertices.find(v) == vertices.end()) {
            candidates.push_back(v);
        }
    }

    // Add more vertices until the clique becomes maximal
    while (!candidates.empty()) {
        int v = candidates.front();
        expandedClique.addVertex(v);
        vertices.insert(v);

        // Update candidates by intersecting with neighborhood of v
        std::vector<int> neighbors = graph.getNeighbors(v);
        std::sort(neighbors.begin(), neighbors.end());

        std::vector<int> newCandidates;
        std::set_intersection(candidates.begin(), candidates.end(),
            neighbors.begin(), neighbors.end(),
            std::back_inserter(newCandidates));

        candidates = newCandidates;
    }

    return expandedClique;
}

// Implementation of Expand-Kou improvement heuristic
EdgeCliqueCover expandKouHeuristic(const Graph& graph, const EdgeCliqueCover& initialCover) {
    
    EdgeCliqueCover expandedCover;

    // Expand each clique in the initial cover
    for (const Clique& clique : initialCover.getCliques()) {
        Clique expandedClique = expandClique(graph, clique);
        expandedCover.addClique(expandedClique);
    }
    
    // Apply Kou post-processing to remove redundant cliques
    return kouPostProcessing(graph, expandedCover);
}

// Helper function to generate sequence of vertices ordered by increasing degree
std::vector<int> generateMinDegreeSequence(const Graph& graph) {
    int n = graph.getNumVertices();
    std::vector<int> sequence(n);
    std::iota(sequence.begin(), sequence.end(), 0); // Fill with 0, 1, ..., n-1

    std::sort(sequence.begin(), sequence.end(),
        [&graph](int a, int b) {
            return graph.getDegree(a) < graph.getDegree(b);
        });

    return sequence;
}

// Helper function to generate sequence of vertices ordered by decreasing degree
std::vector<int> generateMaxDegreeSequence(const Graph& graph) {
    int n = graph.getNumVertices();
    std::vector<int> sequence(n);
    std::iota(sequence.begin(), sequence.end(), 0); // Fill with 0, 1, ..., n-1

    std::sort(sequence.begin(), sequence.end(),
        [&graph](int a, int b) {
            return graph.getDegree(a) > graph.getDegree(b);
        });

    return sequence;
}

tuple<vector<unsigned int>, vector<unsigned int>> readCliquesFromCSV(const string& filename) {

    ifstream file(filename, ios::binary | ios::ate);
    if (!file) throw runtime_error("Cannot open file");

    // Get file size for memory mapping
    streamsize size = file.tellg();
    file.seekg(0, ios::beg);

    // Read file into memory buffer
    vector<char> buffer(size);
    if (!file.read(buffer.data(), size))
        throw runtime_error("Cannot read file");

    // Process buffer
    stringstream ss(string(buffer.data(), size));
    string line;

    // Skip header if exists
    getline(ss, line);

    // Process lines
    vector<unsigned int> clique_start = {};
    vector<unsigned int> clique_end = {};
    while (getline(ss, line)) {
        size_t pos = line.find('\t');
        if (pos != string::npos) {
            clique_start.push_back(stoi(line.substr(0, pos)));
            clique_end.push_back(stoi(line.substr(pos + 1)));
        }
    }

    if (clique_start.size() == 0 || clique_end.size() == 0) { throw runtime_error("size error!"); }

    return { clique_start, clique_end };
}

// Function to read a graph from a file
Graph readGraph(const std::string& edgefile, const std::string& cliquefile ,unsigned int NumberOfNodes) {
    std::ifstream file(edgefile);
    std::string line;
    int n = 0, m = 0;

    cout << "Loading Edges from file: " << edgefile << "\n";
    
    // Create graph
    Graph graph(NumberOfNodes);
    // Read edges from file
    int v1, v2;
    std::vector<int> edge_list;
    vector<tuple<int, int>> edge_list_tuple;
    while (file >> v1 >> v2) {
        graph.addEdge(v1, v2);
    }

    cout << "Generating Cliques \n";
    vector<unsigned int> clique_start;
    vector<unsigned int> clique_end;
    tie(clique_start, clique_end) = readCliquesFromCSV(cliquefile);
    for (unsigned int c = 0; c < clique_start.size(); c++)
    {
        for (int i = clique_start[c]; i < clique_end[c]; i++) {
            for (int j = i + 1; j <= clique_end[c]; j++) {
                graph.addEdge(i, j);
            }
        }
    }

    return graph;
}

tuple<int, int, int> readMetadataFromCSV(const string& filename) {

    ifstream file(filename, ios::binary | ios::ate);
    if (!file) throw runtime_error("Cannot open file");

    // Get file size for memory mapping
    streamsize size = file.tellg();
    file.seekg(0, ios::beg);

    // Read file into memory buffer
    vector<char> buffer(size);
    if (!file.read(buffer.data(), size))
        throw runtime_error("Cannot read file");

    // Process buffer
    stringstream ss(string(buffer.data(), size));
    string line;

    // Skip header if exists
    getline(ss, line);

    // Process lines
    int target_size = 0;
    int edge_size = 0;
    int numberOfNodes = 0;
    while (getline(ss, line)) {
        size_t pos = line.find('\t');
        if (pos != string::npos) {
            string firstPart = line.substr(0, pos);
            if (firstPart == "Total Pieces :") {
                target_size = stoi(line.substr(pos + 1));
            }
            else if (firstPart == "Intra Layer Edges:") {
                edge_size = stoi(line.substr(pos + 1));
            }
            else if (firstPart == "Number of Nodes:") {
                numberOfNodes = stoi(line.substr(pos + 1));
            }

        }
    }

    if (target_size == 0 || edge_size == 0 || numberOfNodes == 0) { throw runtime_error("size error!"); }

    return make_tuple(target_size, edge_size, numberOfNodes);
}



// Main function to test the heuristics
int main(int argc, char* argv[]) {
    vector<string> datasetsThree = { "three", "threep2","threep2w9","threep3","threep3w9" };
    vector<string> datasetsShapes = { "shapes2","shapes4","shapes4-small","shapes5","shapes7","shapes8","shapes9","shapes15" };
    vector<string> datasetsRCO = { "RCO1","RCO2","RCO3","RCO4","RCO5" };
    vector<string> datasetsArtif = { "artif","artif1_2","artif2_4","artif3_6","artif4_8", "artif5_10","artif6_12","artif7_14" };
    vector<string> datasetsBlaz = { "BLAZEWICZ1","BLAZEWICZ2","BLAZEWICZ3","BLAZEWICZ4","BLAZEWICZ5", "blasz2" };
    vector<string> datasetFu = { "fu","fu5","fu6","fu7","fu8","fu9","fu10" };
    vector<string> datasetDagli = { "dagli1" };

    vector<string> datasetsShirts = { "shirts1_2","shirts2_4","shirts3_6","shirts4_8","shirts5_10" };

    vector<string> J1 = { "J1_10_10_0","J1_10_10_1","J1_10_10_2","J1_10_10_3","J1_10_10_4","J1_12_20_0","J1_12_20_1","J1_12_20_2","J1_12_20_3","J1_12_20_4","J1_14_20_0","J1_14_20_1","J1_14_20_2","J1_14_20_3","J1_14_20_4" };
    vector<string> J2 = { "J2_10_35_0","J2_10_35_1","J2_10_35_2","J2_10_35_3","J2_10_35_4","J2_12_35_0","J2_12_35_1","J2_12_35_2","J2_12_35_3","J2_12_35_4","J2_14_35_0","J2_14_35_1","J2_14_35_2","J2_14_35_3","J2_14_35_4" };
    vector<string> poly_jackobs = { "poly1a","poly1b","poly1c","poly1d","poly1e","jakobs1" };
    vector<string> datasetsBlasz2 = {"blasz2" };
    vector<string> datasetsJigSaw = { "JigSaw20" };
    vector<string> datasets;
    /*
    vector<string> datasets = datasetsThree;
    datasets.insert(datasets.end(), datasetsShapes.begin(), datasetsShapes.end());
    datasets.insert(datasets.end(), datasetsRCO.begin(), datasetsRCO.end());
    
    datasets.insert(datasets.end(), datasetsArtif.begin(), datasetsArtif.end());
    datasets.insert(datasets.end(), datasetsBlaz.begin(), datasetsBlaz.end());
    datasets.insert(datasets.end(), datasetFu.begin(), datasetFu.end());
    
    */
    datasets.insert(datasets.end(), datasetDagli.begin(), datasetDagli.end());
    //datasets.insert(datasets.end(), datasetsShirts.begin(), datasetsShirts.end());
    //datasets.insert(datasets.end(), J1.begin(), J1.end());
    //datasets.insert(datasets.end(), J2.begin(), J2.end());
    //datasets.insert(datasets.end(), poly_jackobs.begin(), poly_jackobs.end());

    //datasets = { "three" };
    for (string datasetname : datasetsJigSaw)
    {
        //std::string datasetname = "shapes4";

        std::string CurrDir = getCurrentDirectory();
        std::string masterDirectory = CurrDir + "\\resultsGPU3\\" + datasetname + "\\";
        std::string edgelist = masterDirectory + "graph.csv";
        std::string metadata = masterDirectory + "metadata.csv";
        std::string cliquedata = masterDirectory + "cliques.csv";
        std::string coordsmap = masterDirectory + "pointCoordinate.txt";

        std::string edgeCoverMeta = masterDirectory + "edgeCoverMeta.csv";

        tuple<int, int, int>meta = readMetadataFromCSV(metadata);
        int target_size = get<0>(meta); //this has to be int, otherwise CPLEX would give an error.
        unsigned int NumberofNodes = get<2>(meta);
        // Read graph
        Graph graph = readGraph(edgelist, cliquedata, NumberofNodes);

        std::cout << "Graph loaded: " << graph.getNumVertices() << " vertices, "
            << graph.getNumEdges() << " edges, density: " << graph.getDensity() << std::endl;


        // Run constructive heuristics
        /*
        cout << "running minDegreeSequence...\n";
        auto start = std::chrono::high_resolution_clock::now();
        auto minDegreeSequence = generateMinDegreeSequence(graph);
        EdgeCliqueCover seq1MinCover = sequence1Heuristic(graph, minDegreeSequence);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> seq1Time = end - start;

        cout << "running Sequence2 Heuristic...\n";
        start = std::chrono::high_resolution_clock::now();
        EdgeCliqueCover seq2MinCover = sequence2Heuristic(graph, minDegreeSequence);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> seq2Time = end - start;

        cout << "running minimum1 Heuristic...\n";
        start = std::chrono::high_resolution_clock::now();
        EdgeCliqueCover min1Cover = minimum1Heuristic(graph);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> min1Time = end - start;
        */

        cout << "running maximum1 Heuristic...\n";
        auto start = std::chrono::high_resolution_clock::now();
        EdgeCliqueCover max1Cover = maximum1Heuristic(graph);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> max1Time = end - start;


        /*
        // Run improvement heuristics
        cout << "running Kou Postprocessing for seq1 Min...\n";
        start = std::chrono::high_resolution_clock::now();
        EdgeCliqueCover seq1MinKouCover = kouPostProcessing(graph, seq1MinCover);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> seq1KouTime = end - start;

        cout << "running expand Kou Postprocessing for seq1 Min...\n";
        start = std::chrono::high_resolution_clock::now();
        EdgeCliqueCover seq1MinEKCover = expandKouHeuristic(graph, seq1MinCover);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> seq1EKTime = end - start;

        cout << "running Kou Postprocessing for seq2Min...\n";
        start = std::chrono::high_resolution_clock::now();
        EdgeCliqueCover seq2MinKouCover = kouPostProcessing(graph, seq2MinCover);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> seq2KouTime = end - start;

        cout << "running expand Kou Postprocessing for seq2Min...\n";
        start = std::chrono::high_resolution_clock::now();
        EdgeCliqueCover seq2MinEKCover = expandKouHeuristic(graph, seq2MinCover);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> seq2EKTime = end - start;
        */

        cout << "running expand Kou Postprocessing for MAX1...\n";
        start = std::chrono::high_resolution_clock::now();
        EdgeCliqueCover max1MinEKCover = expandKouHeuristic(graph, max1Cover);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> max1EKTime = end - start;
        /*
        // Print results
        std::cout << "Results:" << std::endl;
        std::cout << "SEQ1-min: " << seq1MinCover.size() << " cliques, time: " << seq1Time.count() << "s" << std::endl;
        std::cout << "SEQ1-min_K: " << seq1MinKouCover.size() << " cliques, time: " << seq1KouTime.count() << "s" << std::endl;
        std::cout << "SEQ1-min_EK: " << seq1MinEKCover.size() << " cliques, time: " << seq1EKTime.count() << "s" << std::endl;
        std::cout << "SEQ2-min: " << seq2MinCover.size() << " cliques, time: " << seq2Time.count() << "s" << std::endl;
        std::cout << "SEQ2-min_K: " << seq2MinKouCover.size() << " cliques, time: " << seq2KouTime.count() << "s" << std::endl;
        std::cout << "SEQ2-min_EK: " << seq2MinEKCover.size() << " cliques, time: " << seq2EKTime.count() << "s" << std::endl;
        std::cout << "MAX1: " << max1Cover.size() << " cliques, time: " << max1Time.count() << "s" << std::endl;
        std::cout << "MIN1: " << min1Cover.size() << " cliques, time: " << min1Time.count() << "s" << std::endl;

        std::cout << "MAX1_EK: " << max1MinEKCover.size() << " cliques, time: " << max1EKTime.count() << "s" << std::endl;

        // Verify solutions
        std::cout << "\nValidation:" << std::endl;
        std::cout << "SEQ1-min: " << (seq1MinCover.isValid(graph) ? "valid" : "invalid") << std::endl;
        std::cout << "SEQ1-min_K: " << (seq1MinKouCover.isValid(graph) ? "valid" : "invalid") << std::endl;
        std::cout << "SEQ1-min_EK: " << (seq1MinEKCover.isValid(graph) ? "valid" : "invalid") << std::endl;
        std::cout << "SEQ2-min: " << (seq2MinCover.isValid(graph) ? "valid" : "invalid") << std::endl;
        std::cout << "SEQ2-min_K: " << (seq2MinKouCover.isValid(graph) ? "valid" : "invalid") << std::endl;
        std::cout << "SEQ2-min_EK: " << (seq2MinEKCover.isValid(graph) ? "valid" : "invalid") << std::endl;
        std::cout << "MAX1: " << (max1Cover.isValid(graph) ? "valid" : "invalid") << std::endl;
        std::cout << "MIN1: " << (min1Cover.isValid(graph) ? "valid" : "invalid") << std::endl;
        std::cout << "MAX1_EK: " << (max1MinEKCover.isValid(graph) ? "valid" : "invalid") << std::endl;
        */

        //Write Run Results
        //remove the output file if already exists
        cout << "save result compilation\n";
        try {
            if (filesystem::remove(edgeCoverMeta)) {
                cout << "Old result has been successfully removed \n" << endl;
            }
        }
        catch (const std::filesystem::filesystem_error& e) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
        std::ofstream outfileMetadata(edgeCoverMeta, std::ios::app);
        
        outfileMetadata << "results:\t" << max1MinEKCover.size() << "\n";
        outfileMetadata << "elapsed time:\t" << (max1Time.count() + max1EKTime.count()) << "\n";
        outfileMetadata << '\n';
        outfileMetadata.close();
        

        max1MinEKCover.writeIntoFile(masterDirectory + "edgecover.txt");
        cout << "save to file complete \n";
    }

    return 0;
}