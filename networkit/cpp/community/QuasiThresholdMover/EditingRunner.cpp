/*
 *
 */
#include <tlx/unused.hpp>
#include <networkit/community/QuasiThresholdEditingLocalMover.hpp>
#include <networkit/community/QuasiThresholdMover/EditingRunner.hpp>
#include <networkit/community/QuasiThresholdMover/QuasiThresholdEditingLinear.hpp>
#include <networkit/generators/TreeReachabilityGraphGenerator.hpp>

namespace NetworKit {

namespace QuasiThresholdMoving {
EditingRunner::EditingRunner(const Graph &G,
                             QuasiThresholdEditingLocalMover::Initialization initialization,
                             count maxIterations, bool sortPaths, bool randomness,
                             count maxPlateauSize, bool useBucketQueue, std::vector<node> order,
                             count insertEditCost, count removeEditCost, std::vector<std::vector<count>> editCostMatrix)
    : G(G), maxIterations(maxIterations), usedIterations(0), sortPaths(sortPaths),
      randomness(randomness), maxPlateauSize(maxPlateauSize),
      insertRun(initialization != QuasiThresholdEditingLocalMover::TRIVIAL
                && initialization != QuasiThresholdEditingLocalMover::EDITING),
      useBucketQueue(useBucketQueue), insertEditCost(insertEditCost), removeEditCost(removeEditCost), editCostMatrix(editCostMatrix), handler(), hasMoved(true),
      marker(G.upperNodeIdBound(), false), lastVisitedDFSNode(G.upperNodeIdBound(), none),
      traversalData(G.upperNodeIdBound()), nodeTouched(G.upperNodeIdBound(), false), rootData(),
      existing(G.upperNodeIdBound(), !insertRun), rootEqualBestParentsCpy(0), currentPlateau(0),
      actualMaximumPlateau(0), gen(Aux::Random::getURNG()), realDist(), intDist(0, 1) {

    runningInfo["time"] = std::vector<count>();
    runningInfo["edits"] = std::vector<count>();
    runningInfo["edits_weight"] = std::vector<count>();
    runningInfo["nodes_moved"] = std::vector<count>();

    if (Aux::PerfEventCountHardware::is_available) {
        event_counters.emplace_back("cache_misses", Aux::PerfEventCountHardware::CACHE_MISSES);
        event_counters.emplace_back("cache_references",
                                    Aux::PerfEventCountHardware::CACHE_REFERENCES);
        event_counters.emplace_back("cycles", Aux::PerfEventCountHardware::CPU_CYCLES);
        event_counters.emplace_back("instructions", Aux::PerfEventCountHardware::INSTRUCTIONS);
    }

    for (const auto &ec : event_counters) {
        runningInfo[ec.first] = std::vector<count>();
    }

    for (auto &ec : event_counters) {
        ec.second.enable();
    }

    timer.start();
    switch (initialization) {
    case QuasiThresholdEditingLocalMover::TRIVIAL: {
        if (G.upperNodeIdBound() == G.numberOfNodes()) {
            dynamicForest = DynamicForest(std::vector<node>(G.upperNodeIdBound(), none));
        } else {
            dynamicForest = DynamicForest(G, std::vector<node>(G.upperNodeIdBound(), none));
        }
        break;
    }
    case QuasiThresholdEditingLocalMover::EDITING: {
        QuasiThresholdEditingLinear editing(G);
        editing.run();
        if (G.upperNodeIdBound() == G.numberOfNodes()) {
            dynamicForest = DynamicForest(editing.getParents());
        } else {
            dynamicForest = DynamicForest(G, editing.getParents());
        }
        break;
    }
    case QuasiThresholdEditingLocalMover::RANDOM_INSERT: {
        this->order.reserve(G.numberOfNodes());
        dynamicForest = DynamicForest(std::vector<node>(G.upperNodeIdBound(), none));
        G.forNodesInRandomOrder([&](node u) { this->order.push_back(u); });
        break;
    }
    case QuasiThresholdEditingLocalMover::ASC_DEGREE_INSERT: {
        dynamicForest = DynamicForest(std::vector<node>(G.upperNodeIdBound(), none));

        std::vector<std::vector<node>> buckets(G.numberOfNodes());
        G.forNodes([&](node u) { buckets[G.degree(u)].push_back(u); });
        this->order.reserve(G.numberOfNodes());
        for (const std::vector<node> &bucket : buckets) {
            for (node u : bucket) {
                this->order.push_back(u);
            }
        }
        break;
    }
    case QuasiThresholdEditingLocalMover::USER_DEFINED_INSERT: {
        dynamicForest = DynamicForest(std::vector<node>(G.upperNodeIdBound(), none));
        this->order = std::move(order);
    }
    default:
        break;
    }

    handler.assureRunning();
    if (useBucketQueue) {
        bucketQueue = BucketQueue(G.upperNodeIdBound());
    } else {
        level = 0;
    }

    numEdits = countNumberOfEdits();
    weightEdits = countWeightOfEdits();
    editsBefore = numEdits;

    G.forNodes([&](node u) { lastVisitedDFSNode[u] = u; });

    timer.stop();

    for (auto &ec : event_counters) {
        ec.second.disable();
    }

    for (auto &ec : event_counters) {
        runningInfo[ec.first].push_back(ec.second.readValue());
    }

    runningInfo["time"].push_back(timer.elapsedMicroseconds());
    runningInfo["nodes_moved"].push_back(0);
}

void EditingRunner::runLocalMover() {
    handler.assureRunning();
    if (!insertRun) {
        runningInfo["edits"].push_back(numEdits);
        runningInfo["edits_weight"].push_back(weightEdits);
    }
    for(count i =0; i < editCostMatrix.size(); i++){
        INFO(editCostMatrix[i]);
    }
    count generation = 0;
    //Main loop; check if iteration moved a node and if max iterations were used
    for (count i = insertRun ? 0 : 1; hasMoved && i <= maxIterations; ++i) {
        //check if iteration is on a plateau greater than maxPlateau size
        if (!hasMoved || (randomness && (currentPlateau >= maxPlateauSize)))
            break;
        handler.assureRunning();
        hasMoved = false;
        numNodesMoved = 0;
        timer.start();

        for (auto &ec : event_counters) {
            ec.second.reset();
            ec.second.enable();
        }
        //first run when inserting is enabled. insert all nodes in assigned order
        if (insertRun) {
            for (index j = 0; j < G.numberOfNodes(); j++) {
                node nodeToMove = order[j];
                localMove(nodeToMove, generation++);
                existing[nodeToMove] = 1;
            }
            insertRun = 0;
        } else {
            //local move in random order
            G.forNodesInRandomOrder([&](node nodeToMove) { localMove(nodeToMove, generation++); });
            INFO("Iteration: ", i, " edits: ", numEdits, " moved nodes: ", numNodesMoved);
        }

        timer.stop();
        for (auto &ec : event_counters) {
            ec.second.disable();
        }
        //update running info, time nodes moved, edits
        if (i == 0) {
            runningInfo["nodes_moved"][0] = numNodesMoved;
            runningInfo["time"][0] += timer.elapsedMicroseconds();
            for (auto &ec : event_counters) {
                runningInfo[ec.first][0] += ec.second.readValue();
            }
        } else {
            runningInfo["nodes_moved"].push_back(numNodesMoved);
            runningInfo["time"].push_back(timer.elapsedMicroseconds());
            for (auto &ec : event_counters) {
                runningInfo[ec.first].push_back(ec.second.readValue());
            }
        }
        runningInfo["edits"].push_back(numEdits);
        runningInfo["edits_weight"].push_back(weightEdits);
        usedIterations = i;

        assert(numEdits == countNumberOfEdits());
        assert(weightEdits == countWeightOfEdits());
        //update plateau
        if (numEdits == editsBefore) {
            currentPlateau++;
        } else {
            if (currentPlateau > actualMaximumPlateau) {
                actualMaximumPlateau = currentPlateau;
            }
            currentPlateau = 0;
        }
        editsBefore = numEdits;
    }
}

Graph EditingRunner::getQuasiThresholdGraph() const {
    Graph forest = dynamicForest.toGraph();
    TreeReachabilityGraphGenerator gen(forest);
    gen.run();
    return gen.getGraph();
}

void EditingRunner::localMove(node nodeToMove, count generation) {
    assert(numEdits == countNumberOfEdits());
    assert(weightEdits == countWeightOfEdits());
    INFO("Move node ", nodeToMove);
    handler.assureRunning();
    numNeighbors = 0;
    G.forEdgesOf(nodeToMove, [&](node v) {
        if (!insertRun || existing[v]) {
            ++numNeighbors;
            marker[v] = true;
            neighbors.push_back(v);
        }
    });

    if (!sortPaths) {
        // Do not even attempt to store children when sortPaths is on
        // as due to the sorting they become meaningless.
        curChildren = dynamicForest.children(nodeToMove);
        curParent = dynamicForest.parent(nodeToMove);
    }
    //remove max depth optimization for weighted TODO maybe replace to factor based on editCosts
    maxDepth = none;
    if(insertEditCost == 1 && removeEditCost == 1){
        maxDepth = 2 * numNeighbors;
    }

    //TODO: multiply with weights
    curEditsWeight = removeEditCost * numNeighbors;
    curEdits = numNeighbors;

    if (!insertRun) {
        //check with theory. delete edges for ancestors or children
        // Calculate the old number of edits incident to c to be able to compute the improvement
        // later
        dynamicForest.dfsFrom(
            nodeToMove,
            [&](node c) {
                if (c != nodeToMove) {
                    //for all children in dfs insert one edge; no insertion necessary if child is neighbor of nodeToMove
                    curEdits += 1 - 2 * marker[c];
                    curEditsWeight += 1 * insertEditCost - marker[c] * (insertEditCost + removeEditCost);//TODO weights
                }
            },
            [](node) {});
        dynamicForest.forAncestors(nodeToMove, [&](node p) {
            //for all ancestors insert one edge; no insertion necessary if ancestor is neighbor of nodeToMove
            curEdits += 1 - 2 * marker[p];
            curEditsWeight += 1  * insertEditCost - marker[p] * (insertEditCost + removeEditCost);//TODO weights
         });
    }
    //first step of algorithm; isolate node
    dynamicForest.isolate(nodeToMove);
    //sort path optimization
    if (sortPaths) {
        for (node v : neighbors) {
            dynamicForest.moveUpNeighbor(v, nodeToMove);
        }
    }
    //TODO ignore bucket queue for weighted because maxDepth is not limited
    if (useBucketQueue) {
        bucketQueue.fill(neighbors, dynamicForest);
    } else {
        //put all neighbors in a queue and sort it based on forest depth
        for (node v : neighbors) {
            neighborQueue.emplace_back(v);
        }
        std::stable_sort(neighborQueue.begin(), neighborQueue.end(), [&](node u, node v) {
            return dynamicForest.depth(u) < dynamicForest.depth(v);
        });
    }

    bestChildren.clear();
    rootData.initialize(generation);

    if (useBucketQueue) {
        while (!bucketQueue.empty()) {
            node u = bucketQueue.next();
            processNode(u, nodeToMove, generation);
        }
    } else {
        while (!currentLevel.empty() || !neighborQueue.empty()) {
            //skip empty levels
            if (currentLevel.empty()) {
                level = dynamicForest.depth(neighborQueue.back());
            }
            //process nodes in currentLevel
            for (node u : currentLevel) {
                assert(dynamicForest.depth(u) == level);
                processNode(u, nodeToMove, generation);
            }
            //process neighbors not in current level but in queue
            while (!neighborQueue.empty() && dynamicForest.depth(neighborQueue.back()) == level) {
                node u = neighborQueue.back();
                neighborQueue.pop_back();
                assert(dynamicForest.depth(u) == level);
                if (nodeTouched[u])
                    continue; // if the node was touched in the previous level, it was in
                              // currentLevel and thus has already been processed
                processNode(u, nodeToMove, generation);
            }
            --level;
            currentLevel.clear();
            currentLevel.swap(nextLevel);
        }
    }

    if (!randomness) {
        if (rootData.childClosenessWeight > rootData.scoreMaxWeight) {
            rootData.bestParentBelow = none;
            rootData.scoreMax = rootData.childCloseness;
            rootData.scoreMaxWeight = rootData.childClosenessWeight;
        }
    } else {
        bool coin = false;
        double ownWeight = rootData.numIndifferentChildren * std::log(2);
        if (rootData.childClosenessWeight > rootData.scoreMaxWeight || !rootData.hasChoices()) {
            // INFO("root better");
            rootData.scoreMax = rootData.childCloseness;
            rootData.scoreMaxWeight = rootData.childClosenessWeight;
            rootData.logEqualBestChoices = ownWeight;
            coin = true;
        } else if (rootData.childClosenessWeight == rootData.scoreMaxWeight) {
            ownWeight = rootData.calculateOwnWeightForEqualChoices();
            if (ownWeight > -std::numeric_limits<double>::infinity()) {
                rootData.addLogChoices(ownWeight);
                coin = logRandomBool(ownWeight - rootData.logEqualBestChoices);
            }
            // INFO("root equally good");
        }
        if (coin) {
            rootData.bestParentBelow = none;
        }

        assert(rootData.hasChoices());
    }

    bestEdits = numNeighbors - rootData.scoreMax;
    //TODO fix scoreMax for weighted edits
    bestEditsWeight = (numNeighbors * removeEditCost )- (rootData.scoreMaxWeight);

    // If sortPaths and randomness is on, only adopt children when the chosen parent is the
    // lower end of its path.
    const TraversalData &bestParentData =
        (rootData.bestParentBelow == none) ? rootData : traversalData[rootData.bestParentBelow];
    // Check if there are any close children at all, or if there are no close children, if
    // randomness is enabled and we have at least two indifferent children (one alone would not be
    // adopted).
    if (bestParentData.numCloseChildren > 0
        || (randomness && bestParentData.numIndifferentChildren > 1)) {
        std::vector<node> indifferentChildren;
        for (node u : touchedNodes) {
            if (u != nodeToMove && dynamicForest.parent(u) == rootData.bestParentBelow) {
                if (traversalData[u].childClosenessWeight > 0) {
                    bestChildren.push_back(u);
                } else if (randomness && traversalData[u].childClosenessWeight == 0) {
                    indifferentChildren.push_back(u);
                }
            }
        }

        assert(bestChildren.size() == bestParentData.numCloseChildren);

        if (randomness) { // make sure we adopt either 0 or at least two children
            assert(indifferentChildren.size() == bestParentData.numIndifferentChildren);
            // If we shall adopt one child, there must be another indifferent child that we just
            // sample randomly to get at least two
            if (bestChildren.size() == 1) {
                assert(!indifferentChildren.empty());
                index i = Aux::Random::index(indifferentChildren.size());
                bestChildren.push_back(indifferentChildren[i]);
                indifferentChildren[i] = indifferentChildren.back();
                indifferentChildren.pop_back();
            }

            // If there are no best children, sample either 0 or at least two indifferent children
            if (bestChildren.empty() && !indifferentChildren.empty()) {
                assert(indifferentChildren.size() != 1);
                if (indifferentChildren.size() == 2) {
                    // Sample either 0 or two nodes from indifferentChildren
                    if (randomBool(2)) {
                        for (node u : indifferentChildren) {
                            bestChildren.push_back(u);
                        }
                    }
                } else {
                    // sample either 0 or at least two nodes from indifferentChildren
                    std::vector<node> sample;
                    do {
                        for (node u : indifferentChildren) {
                            sample.clear();
                            if (randomBool(2)) {
                                sample.push_back(u);
                            }
                        }
                    } while (sample.size() == 1);

                    for (node u : sample) {
                        bestChildren.push_back(u);
                    }
                }
            }

            // If there are already two children, just sample randomly from the remaining
            // indifferent children
            if (bestChildren.size() > 1) {
                for (node u : indifferentChildren) {
                    if (randomBool(2)) {
                        bestChildren.push_back(u);
                    }
                }
            }
        }
    }

#ifndef NDEBUG
    compareWithQuadratic(nodeToMove, generation);
#endif

    // calculate the number of saved edits as comparing the absolute number of edits doesn't make
    // sense
    count savedEdits = curEdits - bestEdits;
    count savedEditsWeight = curEditsWeight - bestEditsWeight;

    // cleanup for linear move
    for (node u : touchedNodes) {
        lastVisitedDFSNode[u] = u;
        nodeTouched[u] = false;
    }

    assert(!randomness || rootData.hasChoices());

    if (rootData.logEqualBestChoices < std::log(std::numeric_limits<long long>::max())) {
        rootEqualBestParentsCpy = std::llround(std::exp(rootData.logEqualBestChoices));
    } else {
        rootEqualBestParentsCpy = std::numeric_limits<count>::max();
    }

    for (node v : neighbors) {
        marker[v] = false;
    }
    neighbors.clear();
    touchedNodes.clear();

    if (sortPaths || savedEditsWeight > 0 || randomness) {
        dynamicForest.moveToPosition(nodeToMove, rootData.bestParentBelow, bestChildren);
        hasMoved |= (savedEditsWeight > 0 || (randomness && rootEqualBestParentsCpy > 1));
        numNodesMoved += (savedEditsWeight > 0 || (randomness && rootEqualBestParentsCpy > 1));
        //only consider savedEdits if edits also also save editCosts in weighted case
        if(savedEditsWeight > 0 ){ 
            numEdits -= savedEdits; 
        }
        weightEdits -= savedEditsWeight;
        INFO("Saved Edits: ", savedEdits, ", SavedWeightedEdits: ", savedEditsWeight);
#ifndef NDEBUG
        assert(numEdits == countNumberOfEdits());
        assert(weightEdits == countWeightOfEdits());
#endif
    } else {
        dynamicForest.moveToPosition(nodeToMove, curParent, curChildren);
#ifndef NDEBUG
        assert(numEdits == countNumberOfEdits());
        assert(weightEdits == countWeightOfEdits());
#endif
    }
}

void EditingRunner::processNode(node u, node nodeToMove, count generation) {
    TRACE("Process ", u);
    TRACE("Processing node ", u, " of depth ", dynamicForest.depth(u),
          " (node to move: ", nodeToMove, ")");
    TRACE("Parent: ", dynamicForest.parent(u), ", children: ", dynamicForest.children(u));
    assert(u != nodeToMove);
    tlx::unused(nodeToMove);
    if (useBucketQueue) {
        assert(dynamicForest.depth(u) <= maxDepth);
    }
    if (!nodeTouched[u]) {
        nodeTouched[u] = true;
        touchedNodes.emplace_back(u);
        assert(
            marker[u]); // only marked neighbors may be processed without having been touched before
    }

    traversalData[u].initialize(generation);

    //TODO check if correct with weights
    int64_t sumPositiveEdits = traversalData[u].childCloseness;
    int64_t sumPositiveEditsWeight = traversalData[u].childClosenessWeight;
    assert(traversalData[u].childCloseness >= 0);

    //TODO: multiply with weight
    traversalData[u].childCloseness += marker[u];
    traversalData[u].childCloseness -=
        1 - marker[u]; // if (marker[u]) { ++traversalData[u].childCloseness; } else {
                       // --traversalData[u].childCloseness; }
    traversalData[u].childClosenessWeight += marker[u] * removeEditCost;
    traversalData[u].childClosenessWeight -=
        1 * insertEditCost - marker[u] * insertEditCost;

    TRACE("Edit difference before descending: ", traversalData[u].childCloseness, ", ", traversalData[u].childClosenessWeight);

    assert(!marker[u] || traversalData[u].childClosenessWeight > 0);

    if (traversalData[u].childClosenessWeight >= 0) {
        assert(lastVisitedDFSNode[u] == u);

        node c = dynamicForest.nextDFSNodeOnEnter(u, u);

        while (c != u) {
            //resturcture if 
            if (!nodeTouched[c] || traversalData[c].childClosenessWeight < 0) {

                if (traversalData[u].childCloseness == 0 || dynamicForest.depth(c) > maxDepth) { //TODO depth check ausbauen wenn depth aufwendig
                    traversalData[u].childCloseness = -1;
                } else {
                    --traversalData[u].childCloseness;
                }
                if (nodeTouched[c] && traversalData[c].childClosenessWeight < 0){
                    traversalData[u].childClosenessWeight = traversalData[u].childClosenessWeight + traversalData[c].childClosenessWeight;
                }
                else{
                    traversalData[u].childClosenessWeight -= insertEditCost;
                }
                // advance to the next starting point for the DFS search.
                c = lastVisitedDFSNode[c];

                if (traversalData[u].childClosenessWeight < 0) {
                    lastVisitedDFSNode[u] = c;
                    break;
                }

                c = dynamicForest.nextDFSNodeOnEnter(c, u);
            } else {
                node p = dynamicForest.parent(c);
                c = dynamicForest.nextChild(c, p);

                while (c == p && c != u) {
                    p = dynamicForest.parent(p);
                    c = dynamicForest.nextChild(c, p);
                }
            }
        }
    }

    TRACE("Edit difference after descending: ", traversalData[u].childCloseness, ", ", traversalData[u].childClosenessWeight);

    if (!randomness) {
        if (sumPositiveEditsWeight > traversalData[u].scoreMaxWeight || traversalData[u].scoreMaxWeight == 0) {
            //TODO check for weights
            traversalData[u].scoreMax = sumPositiveEdits;
            traversalData[u].scoreMaxWeight = sumPositiveEditsWeight;
            traversalData[u].bestParentBelow = u;
        }
    } else {
        bool coin = false;
        double ownWeight = traversalData[u].numIndifferentChildren * std::log(2);
        if (sumPositiveEditsWeight > traversalData[u].scoreMaxWeight || !traversalData[u].hasChoices()) {
            // INFO(u, " is better count = 1");
            traversalData[u].scoreMax = sumPositiveEdits;
            traversalData[u].scoreMaxWeight = sumPositiveEditsWeight;
            traversalData[u].logEqualBestChoices = ownWeight;
            // Either we do not adopt children, or we are at the lower end of a path.
            // Otherwise, there must be a node below u that is at least as good.
            assert(!sortPaths || sumPositiveEdits == 0 || dynamicForest.isLowerEnd(u));
            coin = true;
        } else if (sumPositiveEditsWeight == traversalData[u].scoreMaxWeight) {
            ownWeight = traversalData[u].calculateOwnWeightForEqualChoices();
            if (ownWeight > -std::numeric_limits<double>::infinity()) {
                traversalData[u].addLogChoices(ownWeight);
                coin = logRandomBool(ownWeight - traversalData[u].logEqualBestChoices);
            }
            assert(traversalData[u].hasChoices());
            // INFO(u, " equally good count = ", traversalData[u].equalBestParents);
        }
        if (coin) {
            traversalData[u].bestParentBelow = u;
        }
    }

    assert(traversalData[u].scoreMax >= 0);

    traversalData[u].scoreMax += marker[u];
    traversalData[u].scoreMaxWeight += marker[u] * removeEditCost;

    if (traversalData[u].scoreMaxWeight > 0) {
        traversalData[u].scoreMax -= 1 - marker[u];
        traversalData[u].scoreMaxWeight -= 1 *insertEditCost - marker[u] * insertEditCost;
    }
    TRACE("Maximum gain at ", u, ": ", traversalData[u].scoreMax, ", ", traversalData[u].scoreMaxWeight);
    node p = dynamicForest.parent(u);
    TraversalData &parentData = (p == none) ? rootData : traversalData[p];

    parentData.initialize(generation);

    if ((traversalData[u].scoreMaxWeight > 0 || traversalData[u].childClosenessWeight > 0) && p != none) {
        if (useBucketQueue) {
            assert(dynamicForest.depth(p) <= maxDepth);
        }
        if (!nodeTouched[p]) {
            nodeTouched[p] = true;
            touchedNodes.push_back(p);
            if (!useBucketQueue) {
                nextLevel.push_back(p);
            } else if (!marker[p]) { // neighbors already in queue
                bucketQueue.insertParent(p);
            }
        }
    }

    if (traversalData[u].scoreMaxWeight > parentData.scoreMaxWeight) {
        parentData.logEqualBestChoices = traversalData[u].logEqualBestChoices;
        parentData.scoreMax = traversalData[u].scoreMax;
        parentData.scoreMaxWeight = traversalData[u].scoreMaxWeight;
        parentData.bestParentBelow = traversalData[u].bestParentBelow;
        // INFO(u, " better for ", p);
        // INFO("set count to ", parentData.equalBestParents);
    } else if (randomness && traversalData[u].scoreMaxWeight == parentData.scoreMaxWeight) {
        parentData.addLogChoices(traversalData[u].logEqualBestChoices);
        if (logRandomBool(traversalData[u].logEqualBestChoices - parentData.logEqualBestChoices)) {
            parentData.bestParentBelow = traversalData[u].bestParentBelow;
        }
        // INFO(u, " equally good for ", p);
        // INFO("increase count by ", traversalData[u].equalBestParents, " to ",
        // parentData.equalBestParents);
    }

    if (traversalData[u].childClosenessWeight >= 0) {
        assert(traversalData[u].childCloseness <= traversalData[u].scoreMaxWeight);
        parentData.childCloseness += traversalData[u].childCloseness;
        parentData.childClosenessWeight += traversalData[u].childClosenessWeight;
        if (traversalData[u].childClosenessWeight == 0) {
            ++parentData.numIndifferentChildren;
        } else {
            ++parentData.numCloseChildren;
        }
    }

    assert(!dynamicForest.children(u).empty() || traversalData[u].childCloseness == 1);
}

void EditingRunner::compareWithQuadratic(node nodeToMove, count generation) const {
    std::vector<int64_t> missingBelow, missingAbove, existingBelow, existingAbove;
    missingBelow.resize(G.upperNodeIdBound(), 0);
    missingAbove.resize(G.upperNodeIdBound(), 0);
    existingBelow.resize(G.upperNodeIdBound(), 0);
    existingAbove.resize(G.upperNodeIdBound(), 0);
    std::vector<bool> usingDeepNeighbors(G.upperNodeIdBound(), false);
    dynamicForest.forChildrenOf(none, [&](node r) {
        if (existing[r]) {
            dynamicForest.dfsFrom(
                r,
                [&](node u) {
                    if (dynamicForest.depth(u) > maxDepth)
                        usingDeepNeighbors[u] = true;
                    if (u != nodeToMove) {
                        missingBelow[u] = missingAbove[u] = 1 - marker[u];
                        existingBelow[u] = existingAbove[u] = marker[u];
                    }
                    node p = dynamicForest.parent(u);
                    if (p != none) {
                        missingAbove[u] += missingAbove[p];
                        existingAbove[u] += existingAbove[p];
                    }
                },
                [&](node u) {
                    node p = dynamicForest.parent(u);
                    if (p != none) {
                        missingBelow[p] += missingBelow[u];
                        existingBelow[p] += existingBelow[u];
                        if (usingDeepNeighbors[u])
                            usingDeepNeighbors[p] = true;
                    }
                });
        }
    });

    assert(missingBelow[nodeToMove] == 0);
    assert(existingBelow[nodeToMove] == 0);

    if (!sortPaths) {
        bool exactValue = true;
        for (node c : curChildren) {
            missingBelow[nodeToMove] += missingBelow[c];
            existingBelow[nodeToMove] += existingBelow[c];
            if (usingDeepNeighbors[c])
                exactValue = false;
        }

        if (curParent != none) {
            missingAbove[nodeToMove] = missingAbove[curParent];
            existingAbove[nodeToMove] = existingAbove[curParent];
            if (usingDeepNeighbors[curParent])
                exactValue = false;
        }
        if (exactValue) {
            assert(curEdits
                   == numNeighbors - existingAbove[nodeToMove] - existingBelow[nodeToMove]
                          + missingAbove[nodeToMove] + missingBelow[nodeToMove]);
        }
    }

    count minEdits = std::numeric_limits<count>::max();
    std::vector<node> minChildren;
    node minParent = curParent;
    G.forNodes([&](node u) {
        if (u == nodeToMove || usingDeepNeighbors[u] || !existing[u])
            return;
        if (existingBelow[u] >= missingBelow[u]
            || (traversalData[u].generation == generation && traversalData[u].childCloseness > 0)) {
            assert(traversalData[u].childCloseness == existingBelow[u] - missingBelow[u]);
        } else if (nodeTouched[u]) {
            assert(traversalData[u].childCloseness < 0);
        }
    });

    G.forNodes([&](node u) {
        if (dynamicForest.children(u).empty() && marker[u] && !usingDeepNeighbors[u]
            && existing[u]) {
            assert(traversalData[u].childCloseness == 1);
        }
    });

    auto tryEditBelow = [&](node p) {
        if (p == nodeToMove)
            return;

        count edits = numNeighbors;
        if (p != none) {
            edits += missingAbove[p];
            edits -= existingAbove[p];
        }

        std::vector<node> children;
        dynamicForest.forChildrenOf(p, [&](node c) {
            if (c == nodeToMove || usingDeepNeighbors[c] || !existing[c])
                return;
            if (existingBelow[c] > missingBelow[c]) { // TODO try >= (more children...)
                if (dynamicForest.children(c).empty() && marker[c]) {
                    assert(traversalData[c].childCloseness == 1);
                }
                assert(traversalData[c].childCloseness == existingBelow[c] - missingBelow[c]);

                children.emplace_back(c);
                edits -= existingBelow[c] - missingBelow[c];
            }
        });

        if (edits < minEdits) {
            minEdits = edits;
            minChildren = std::move(children);
            minParent = p;
        }
    };

    dynamicForest.dfsFrom(
        none, [](node) {}, tryEditBelow);
    tryEditBelow(none);

    assert(minEdits >= bestEdits);

    count childClosenessControl = numNeighbors;
    if (rootData.bestParentBelow != none) {
        childClosenessControl -=
            (existingAbove[rootData.bestParentBelow] - missingAbove[rootData.bestParentBelow]);
    }
    for (node u : bestChildren) {
        childClosenessControl -= (existingBelow[u] - missingBelow[u]);
    }
    if (!sortPaths) {
        TRACE("Current edits: ", curEdits, " (with parent ", curParent, " and current children ",
              curChildren, "), minimum edits: ", minEdits);
    } else {
        TRACE("Current edits: ", curEdits, ", minimum edits: ", minEdits);
    }
    TRACE("Quadratic algorithm wants to have new parent ", minParent, " and new children ",
          minChildren);
    TRACE("Linear algorithm wants to have new parent ", rootData.bestParentBelow,
          " and new children ", bestChildren, " edits: ", childClosenessControl);
    assert(minEdits >= childClosenessControl);

    G.forNodes([&](node u) {
        tlx::unused(u);
        TRACE("Node", marker[u] ? "[x]" : "", " ", u, ", parent ", dynamicForest.parent(u));
    });
}

count EditingRunner::countNumberOfEdits() const {
    // count number of edits that are needed with the initial given forest
    count numExistingEdges = 0;
    count numMissingEdges = 0;
    std::vector<bool> marker(G.upperNodeIdBound());

    dynamicForest.forChildrenOf(none, [&](node r) {
        count depth = 0;
        dynamicForest.dfsFrom(
            r,
            [&](node u) { // on enter
                count upperNeighbors = 0;

                G.forNeighborsOf(u, [&](node v) {
                    if (marker[v])
                        ++upperNeighbors;
                });

                numExistingEdges += upperNeighbors;
                numMissingEdges += depth - upperNeighbors;
                marker[u] = true;
                depth += 1;
            },
            [&](node u) { // on exit
                marker[u] = false;
                depth -= 1;
            });
    });

    return numMissingEdges + (G.numberOfEdges() - numExistingEdges);
}

count EditingRunner::countWeightOfEdits() const {
    // count weight of edits that are needed with the initial given forest
    count numExistingEdges = 0;
    count numMissingEdges = 0;
    std::vector<bool> marker(G.upperNodeIdBound());

    dynamicForest.forChildrenOf(none, [&](node r) {
        count depth = 0;
        dynamicForest.dfsFrom(
            r,
            [&](node u) { // on enter
                count upperNeighbors = 0;

                G.forNeighborsOf(u, [&](node v) {
                    if (marker[v])
                        ++upperNeighbors;
                });

                numExistingEdges += upperNeighbors;
                numMissingEdges += depth - upperNeighbors;
                marker[u] = true;
                depth += 1;
            },
            [&](node u) { // on exit
                marker[u] = false;
                depth -= 1;
            });
    });
    TRACE("Missing Edges:  ", numMissingEdges,", NumEdges-ExistingEdges: " ,(G.numberOfEdges() - numExistingEdges));
    return ((numMissingEdges * insertEditCost) + ((G.numberOfEdges() - numExistingEdges) * removeEditCost));
}

count EditingRunner::editsIncidentTo(node u) const {
    std::vector<bool> visited(G.upperNodeIdBound());
    count edits = 0;
    node parent = dynamicForest.parent(u);
    while (parent != none) {
        visited[parent] = 1;
        if (!G.hasEdge(parent, u))
            edits++;
        parent = dynamicForest.parent(parent);
    }
    dynamicForest.forChildrenOf(u, [&](node c) {
        dynamicForest.dfsFrom(
            c,
            [&](node w) {
                visited[w] = 1;
                if (!G.hasEdge(w, u))
                    edits++;
            },
            [&](node) {});
    });

    G.forNeighborsOf(u, [&](node w) {
        if (!visited[w])
            edits++;
    });

    return edits;
}
} // namespace QuasiThresholdMoving

} // namespace NetworKit
