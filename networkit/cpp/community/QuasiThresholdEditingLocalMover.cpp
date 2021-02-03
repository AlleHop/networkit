
#include <networkit/community/QuasiThresholdEditingLocalMover.hpp>
#include <networkit/community/QuasiThresholdMover/EditingRunner.hpp>

namespace NetworKit {
namespace QuasiThresholdMoving {
QuasiThresholdEditingLocalMover::QuasiThresholdEditingLocalMover(
    const Graph &G, Initialization initialization, count maxIterations, bool sortPaths,
    bool randomness, bool moveSubtrees, count maxPlateauSize, bool useBucketQueue, count insertEditCost, count removeEditCost, std::vector<std::vector<int64_t>> editCostMatrix)
    : G(G), initialization(initialization), maxIterations(maxIterations), sortPaths(sortPaths),
      randomness(randomness), moveSubtrees(moveSubtrees), maxPlateauSize(maxPlateauSize), useBucketQueue(useBucketQueue), insertEditCost(insertEditCost), removeEditCost(removeEditCost), editCostMatrix(editCostMatrix),
      usedIterations(0), numEdits(0), weightEdits(0), rootEqualBestParents(0) {}

void QuasiThresholdEditingLocalMover::run() {
    EditingRunner runner(G, initialization, maxIterations, sortPaths, randomness, moveSubtrees, maxPlateauSize,
                         useBucketQueue, order, insertEditCost, removeEditCost, editCostMatrix);
    runner.runLocalMover();
    usedIterations = runner.getUsedIterations();
    numEdits = runner.getNumberOfEdits();
    weightEdits = runner.getWeightOfEdits();
    plateauSize = runner.getPlateauSize();
    rootEqualBestParents = runner.getRootEqualBestParents();
    quasiThresholdGraph = runner.getQuasiThresholdGraph();
    dynamicForestGraph = runner.getDynamicForestGraph();
    runningInfo = runner.getRunningInfo();
    hasRun = true;
}

Graph QuasiThresholdEditingLocalMover::getQuasiThresholdGraph() const {
    assureFinished();
    return quasiThresholdGraph;
}

Graph QuasiThresholdEditingLocalMover::getDynamicForestGraph() const {
    assureFinished();
    return dynamicForestGraph;
}

count QuasiThresholdEditingLocalMover::getNumberOfEdits() const {
    assureFinished();
    return numEdits;
}

count QuasiThresholdEditingLocalMover::getWeightOfEdits() const {
    assureFinished();
    return weightEdits;
}

count QuasiThresholdEditingLocalMover::getUsedIterations() const {
    assureFinished();
    return usedIterations;
}

count QuasiThresholdEditingLocalMover::getPlateauSize() const {
    assureFinished();
    return plateauSize;
}

count QuasiThresholdEditingLocalMover::getRootEqualBestParents() const {
    assureFinished();
    return rootEqualBestParents;
}

std::map<std::string, std::vector<count>> QuasiThresholdEditingLocalMover::getRunningInfo() const {
    assureFinished();
    return runningInfo;
}

void QuasiThresholdEditingLocalMover::setInsertionOrder(std::vector<node> order) {
    this->order = order;
}

} // namespace QuasiThresholdMoving
} // namespace NetworKit
