#ifndef QUASITHRESHOLDEDITINGLOCALMOVER_H
#define QUASITHRESHOLDEDITINGLOCALMOVER_H

#include <map>
#include <networkit/base/Algorithm.hpp>
#include <networkit/community/QuasiThresholdMover/QuasiThresholdEditingLinear.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {
namespace QuasiThresholdMoving {
class QuasiThresholdEditingLocalMover : public Algorithm {
public:
    enum Initialization { TRIVIAL, EDITING, RANDOM_INSERT, ASC_DEGREE_INSERT, DESC_DEGREE_INSERT, USER_DEFINED_INSERT };

    QuasiThresholdEditingLocalMover(const Graph &G, Initialization initializarion = TRIVIAL,
                                    count maxIterations = 5, bool sortPaths = true,
                                    bool randomness = false, count maxPlateauSize = 4,
                                    bool useBucketQueue = true,
                                    count insertEditCost = 1, count removeEditCost = 1, std::vector<std::vector<int64_t>> editCostMatrix = {});

    void run() override;

    Graph getQuasiThresholdGraph() const;
    count getNumberOfEdits() const;
    count getWeightOfEdits() const;
    count getUsedIterations() const;
    count getPlateauSize() const;
    count getRootEqualBestParents() const;
    std::map<std::string, std::vector<count>> getRunningInfo() const;

    void setInsertionOrder(std::vector<node> order);

private:
    const Graph &G;
    Initialization initialization;
    count maxIterations;
    bool sortPaths;
    bool randomness;
    count maxPlateauSize;
    bool insertRun;
    bool useBucketQueue;
    count insertEditCost;
    count removeEditCost;
    std::vector<std::vector<int64_t>> editCostMatrix;

    std::vector<node> order;

    count usedIterations;
    count numEdits;
    count weightEdits;
    count plateauSize;
    count rootEqualBestParents;

    Graph quasiThresholdGraph;

    std::map<std::string, std::vector<count>> runningInfo;
};
} // namespace QuasiThresholdMoving

} // namespace NetworKit

#endif // QUASITHRESHOLDEDITINGLOCALMOVER_H
