#ifndef EDITINGRUNNER_H
#define EDITINGRUNNER_H

#include <limits>
#include <map>
#include <networkit/auxiliary/PerfEventCountHardware.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/community/QuasiThresholdEditingLocalMover.hpp>
#include <networkit/community/QuasiThresholdMover/BucketQueue.hpp>
#include <networkit/community/QuasiThresholdMover/DynamicForest.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {
namespace QuasiThresholdMoving {
class EditingRunner {
public:
    EditingRunner(const Graph &G, QuasiThresholdEditingLocalMover::Initialization initialization,
                  count maxIterations, bool sortPaths, bool randomness, bool moveSubtrees, bool subtreeSortPaths, count maxPlateauSize,
                  bool useBucketQueue, std::vector<node> order,
                  count insertEditCost, count removeEditCost, const std::vector<std::vector<int64_t>> &editCostMatrix);

    void runLocalMover();

    count getNumberOfEdits() const { return numEdits; };

    count getCostOfEdits() const { return costEdits; };

    count getUsedIterations() const { return usedIterations; };

    count getPlateauSize() const { return actualMaximumPlateau; };

    count getRootEqualBestParents() const { return rootEqualBestParentsCpy; };

    std::map<std::string, std::vector<count>> getRunningInfo() const { return runningInfo; }

    Graph getQuasiThresholdGraph() const;

    Graph getDynamicForestGraph() const;

private:
    struct TraversalData {
        count generation;
        int64_t scoreMax;
        int64_t childCloseness;
        int64_t scoreMaxCost;
        int64_t childClosenessCost;
        int64_t subtreeEdits;
        int64_t sumPositiveEdits;
        int64_t sumPositiveEditCosts;
        int64_t subtreeEditCosts;
        node bestParentBelow;
        /**
         * Logarithm of the number of equally good choices for this node.
         * The logarithm is used as the number of choices can be exponential in the number of nodes.
         * Therefore, the logarithm is guaranteed to be linear in the number of nodes.
         */
        double logEqualBestChoices;
        count numIndifferentChildren;
        count numCloseChildren;

        TraversalData()
            : generation(none), scoreMax(0), childCloseness(0), scoreMaxCost(0),
              childClosenessCost(0), bestParentBelow(none),
              subtreeEdits(0), subtreeEditCosts(0),
              logEqualBestChoices(-std::numeric_limits<double>::infinity()),
              numIndifferentChildren(0), numCloseChildren(0){};

        void initialize(count currentGeneration) {
            if (currentGeneration != generation) {
                generation = currentGeneration;
                scoreMax = 0;
                childCloseness = 0;
                scoreMaxCost = 0;
                childClosenessCost = 0;
                sumPositiveEdits = 0;
                sumPositiveEditCosts = 0;
                subtreeEdits = 0;
                subtreeEditCosts = 0;
                bestParentBelow = none;
                logEqualBestChoices = -std::numeric_limits<double>::infinity();
                numIndifferentChildren = 0;
                numCloseChildren = 0;
            }
        };

        bool hasChoices() { return logEqualBestChoices > -std::numeric_limits<double>::infinity(); }

        void addEqualChoices(count choices) { addLogChoices(std::log(choices)); }

        void addLogChoices(double logChoices) {
            /* This uses the technique explained on
            https://en.wikipedia.org/w/index.php?title=Log_probability&oldid=954092640#Addition_in_log_space
            to avoid issues when std::exp(logEqualBestChoices) cannot
            be represented as a double anymore. */
            if (logChoices > logEqualBestChoices) {
                logEqualBestChoices =
                    logChoices + std::log1p(std::exp(logEqualBestChoices - logChoices));
            } else {
                logEqualBestChoices += std::log1p(std::exp(logChoices - logEqualBestChoices));
            }
        }

        double calculateOwnWeightForEqualChoices() {
            double ownWeight = numIndifferentChildren * std::log(2);
            // Never adopt only one child, as this is always equivalent to choosing the child as
            // parent and adopting all its children
            if (numCloseChildren != 1
                || numIndifferentChildren > 0) {
                if (numCloseChildren == 0) {
                    if (numIndifferentChildren < 2) {
                        // we need no or at least two indifferent children, so here only 0
                        // remains as option
                        ownWeight = 0;
                    } else {
                        // Subtract the numIndifferentChildren possibilities of choosing exactly
                        // one indifferent child.
                        ownWeight += std::log1p(-std::exp(
                            ownWeight - std::log(numIndifferentChildren)));
                    }
                } else if (numCloseChildren == 1) {
                    // Subtract the one possibility of choosing exactly 0 indifferent children
                    ownWeight += std::log1p(-std::exp(-ownWeight));
                }

                return ownWeight;
            }

            return -std::numeric_limits<double>::infinity();
        }

        std::string toString() {
            std::stringstream ss;
            ss << "\n";
            ss << "scoreMax: " << scoreMax << "\n";
            ss << "childCloseness: " << childCloseness << "\n";
            ss << "scoreMaxCost: " << scoreMaxCost << "\n";
            ss << "childClosenessCost: " << childClosenessCost << "\n";
            ss << "subtreeEdits: " << subtreeEdits << "\n";
            ss << "subtreeEditCosts: " << subtreeEditCosts << "\n";
            ss << "logEqualBestChoices: " << logEqualBestChoices << "\n";
            ss << "bestParentBelow: " << bestParentBelow << "\n";
            return ss.str();
        };
    };

    const Graph &G;
    count maxIterations;
    count usedIterations;
    bool sortPaths;
    bool randomness;
    bool moveSubtrees;
    bool subtreeSortPaths;
    std::vector<node> order;
    count maxPlateauSize;

    bool insertRun;
    bool useBucketQueue;
    bool editMatrixUsed;
    count insertEditCost;
    count removeEditCost;
    const std::vector<std::vector<int64_t>> &editCostMatrix;

    count numEdits;
    count costEdits;

    Aux::SignalHandler handler;
    DynamicForest dynamicForest;
    bool hasMoved;
    std::vector<bool> marker;

    count generation;
    count level;
    std::vector<node> neighborQueue;
    std::vector<node> currentLevel;
    std::vector<node> nextLevel;

    BucketQueue bucketQueue;
    std::vector<node> neighbors;
    count numNeighbors;
    count editCostNeighbors;
    std::vector<int64_t> editCostNodeToMove;
    std::vector<node> touchedNodes;
    std::vector<node> lastVisitedDFSNode;
    std::vector<TraversalData> traversalData;
    std::vector<bool> nodeTouched;

    TraversalData rootData;

    count bestEdits;
    count bestEditCosts;
    count curEdits;
    count curEditCosts;
    int64_t savedEdits;
    count savedEditCosts;
    node curParent;
    std::vector<node> curChildren;
    std::vector<node> bestChildren;

    count numSubtreeMoves;
    count subtreeSize;
    count subtreeOption;
    std::vector<bool> inSubtree;
    std::vector<node> subtreeNeighbors;
    std::vector<count> numNeighborsAll;
    std::vector<node> subtreeNodes;
    std::vector<node> parentCandidates;
    std::vector<node> parentQueue;
    std::vector<int64_t> editCostSubtree;
    count subtreeExtDegree = 0;

    std::vector<bool> existing;

    count rootEqualBestParentsCpy;
    count numNodesMoved;
    
    count editsBefore;
    count editCostsBefore;
    count currentPlateau;
    count actualMaximumPlateau;

    count maxDepth;

    std::mt19937_64 &gen;
    std::uniform_real_distribution<double> realDist;
    std::uniform_int_distribution<count> intDist;

    Aux::Timer timer;
    std::vector<std::pair<std::string, Aux::PerfEventCountHardware>> event_counters;
    std::map<std::string, std::vector<count>> runningInfo;

    void localMove(node nodeToMove);
    void processNode(node u, node nodeToMove);
    void processNodeForSubtree(node u, node nodeToMove);
    void subtreeMove(node nodeToMove);
    Graph getGraphFromEditMatrix();
    std::vector<node> getParentsForTree();
    void compareWithQuadratic(node nodeToMove) const;
    count countNumberOfEdits() const;
    count countCostOfEdits() const;
    count editsIncidentTo(node u) const;

    bool logRandomBool(double logProbability) {
        assert(logProbability <= 0);
        double x = realDist(gen);
        if (x == 0)
            return logProbability > -std::numeric_limits<double>::infinity();
        return std::log(x) < logProbability;
    }

    bool randomBool(count options, count optionsToConsider = 1) {
        assert(options > 0);
        assert(options >= optionsToConsider);
        count x = intDist(gen, std::uniform_int_distribution<count>::param_type(0, options - 1));
        return x < optionsToConsider;
    };
};
} // namespace QuasiThresholdMoving

} // namespace NetworKit

#endif // EDITINGRUNNER_H
