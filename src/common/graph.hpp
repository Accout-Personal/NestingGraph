#pragma once
#include<vector>
#include <set>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <string>

// Graph representation using adjacency matrix
class Graph {
private:
    uint32_t numVertices;
    uint64_t numberEdges = 0;
    std::vector<std::vector<bool>> adjacencyMatrix;

public:
    // Constructor
    Graph(uint32_t n) : numVertices(n) {
        adjacencyMatrix.resize(n, std::vector<bool>(n, false));
    }

    // Add an edge to the graph
    void addEdge(uint32_t u, uint32_t v) {
        if (u != v) { // No self-loops
            if (!adjacencyMatrix[u][v]) {
                adjacencyMatrix[u][v] = true;
                adjacencyMatrix[v][u] = true;
                numberEdges++;
            }
            
        }
    }

    // Check if there is an edge between two vertices
    bool hasEdge(uint32_t u, uint32_t v) const {
        return adjacencyMatrix[u][v];
    }

    // Get number of vertices
    uint32_t getNumVertices() const {
        return numVertices;
    }

    // Get neighbors of a vertex
    std::vector<uint32_t> getNeighbors(uint32_t v) const {
        std::vector<uint32_t> neighbors;
        for (uint32_t i = 0; i < numVertices; ++i) {
            if (adjacencyMatrix[v][i]) {
                neighbors.push_back(i);
            }
        }
        return neighbors;
    }

    // Get the degree of a vertex
    int getDegree(uint32_t v) const {
        int degree = 0;
        for (uint32_t i = 0; i < numVertices; ++i) {
            if (adjacencyMatrix[v][i]) {
                degree++;
            }
        }
        return degree;
    }

    // Get edges of the graph
    std::vector<std::pair<uint32_t, uint32_t>> getEdges() const {
        std::vector<std::pair<uint32_t, uint32_t>> edges;
        edges.reserve(numberEdges);
        for (uint32_t i = 0; i < numVertices; ++i) {
            for (uint32_t j = i + 1; j < numVertices; ++j) {
                if (adjacencyMatrix[i][j]) {
                    edges.push_back({ i, j });
                }
            }
        }
        return edges;
    }

    void writeEdgesToFile(const std::string& filename) const {
        std::ofstream outFile(filename);
        if (!outFile) {
            throw std::runtime_error("Could not open file for writing: " + filename);
        }

        for (uint32_t i = 0; i < numVertices; ++i) {
            for (uint32_t j = i + 1; j < numVertices; ++j) {
                if (adjacencyMatrix[i][j]) {
                    outFile << i << "\t" << j << "\n";
                }
            }
        }

        outFile.close();
    }


    // Get the number of edges in the graph
    uint32_t getNumEdges() const {
        return numberEdges;
    }

    // Get the density of the graph
    double getDensity() const {
        uint32_t n = numVertices;
        std::cout << "numVertices: " << n << "\n";
        std::cout << "numberEdges: " << numberEdges << "\n";
        uint32_t e = numberEdges;
        return (2.0 * e) / (n * (n - 1));
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
    unsigned int size() const {
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

    void writeIntoFile(const std::string& filename)
    {
        //remove the output file if already exists
        try {
            if (std::filesystem::remove(filename)) {
                std::cout << "Old clique data has been successfully removed \n" << std::endl;
            }
        }
        catch (const std::filesystem::filesystem_error& e) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
        std::ofstream outfileClique(filename, std::ios::app);
        for (Clique clique : cliques)
        {
            std::set<int> vertexes = clique.getVertices();
            for (int vertex : vertexes)
            {
                outfileClique << vertex << " ";
            }
            outfileClique << '\n';
        }

        outfileClique.close();
        std::cout << "save to file complete \n";
    }

};

// Implementation of Maximum1 heuristic (similar to Sequence1 but with vertex ordering by max degree)
EdgeCliqueCover maximum1Heuristic(const Graph& graph) {
    EdgeCliqueCover cover;
    int n = graph.getNumVertices();

    std::cout << "initialization \n";
    // Initialize uncovered edge matrix
    std::vector<std::vector<bool>> uncoveredEdge(n);
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

    std::cout << "main loop \n";
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
                std::vector<uint32_t> newNeighbors = graph.getNeighbors(u);
                std::sort(newNeighbors.begin(), newNeighbors.end());
                std::vector<uint32_t> intersection;
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

// Implementation of Expansion algorithm (Algorithm 4 in the paper)
Clique expandClique(const Graph& graph, const Clique& initialClique) {

    Clique expandedClique = initialClique;

    // Get the current vertices in the clique
    std::set<int> vertices = expandedClique.getVertices();

    // Compute initial intersection of neighborhoods
    std::vector<uint32_t> intersection;
    bool first = true;

    for (int v : vertices) {
        std::vector<uint32_t> neighbors = graph.getNeighbors(v);
        std::sort(neighbors.begin(), neighbors.end());

        if (first) {
            intersection = neighbors;
            first = false;
        }
        else {
            std::vector<uint32_t> temp;
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
        std::vector<uint32_t> neighbors = graph.getNeighbors(v);
        std::sort(neighbors.begin(), neighbors.end());

        std::vector<int> newCandidates;
        std::set_intersection(candidates.begin(), candidates.end(),
            neighbors.begin(), neighbors.end(),
            std::back_inserter(newCandidates));

        candidates = newCandidates;
    }

    return expandedClique;
}

EdgeCliqueCover kouPostProcessing(const Graph& graph, const EdgeCliqueCover& initialCover) {
    std::cout << "kou post processing...\n";
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