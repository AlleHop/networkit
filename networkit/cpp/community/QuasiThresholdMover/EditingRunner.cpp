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
                             count maxIterations, bool sortPaths, bool randomness, bool moveSubtrees,
                             count maxPlateauSize, bool useBucketQueue, std::vector<node> order,
                             count insertEditCost, count removeEditCost, std::vector<std::vector<int64_t>> editCostMatrix)
    : G(G), maxIterations(maxIterations), usedIterations(0), sortPaths(sortPaths),
      randomness(randomness), moveSubtrees(moveSubtrees), maxPlateauSize(maxPlateauSize),
      insertRun(initialization != QuasiThresholdEditingLocalMover::TRIVIAL
                && initialization != QuasiThresholdEditingLocalMover::EDITING),
      useBucketQueue(useBucketQueue
      && !(editCostMatrix.size() != 0 && editCostMatrix.size() == G.upperNodeIdBound()) 
      && (insertEditCost == 1 && removeEditCost ==1)), 
      insertEditCost(insertEditCost), removeEditCost(removeEditCost), editCostMatrix(editCostMatrix), handler(), hasMoved(true),
      marker(G.upperNodeIdBound(), false), inSubtree(G.upperNodeIdBound(), false),  numNeighborsAll(G.upperNodeIdBound(), 0),lastVisitedDFSNode(G.upperNodeIdBound(), none),
      traversalData(G.upperNodeIdBound()), nodeTouched(G.upperNodeIdBound(), false), nodeTouchedSubtree(G.upperNodeIdBound(), false), rootData(),
      existing(G.upperNodeIdBound(), !insertRun), rootEqualBestParentsCpy(0), currentPlateau(0),
      actualMaximumPlateau(0), gen(Aux::Random::getURNG()), realDist(), intDist(0, 1) {

    runningInfo["time"] = std::vector<count>();
    runningInfo["edits"] = std::vector<count>();
    runningInfo["edits_weight"] = std::vector<count>();
    runningInfo["nodes_moved"] = std::vector<count>();
    runningInfo["subtrees_moved"] = std::vector<count>();

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

    //check if edit Matrix has correct size
    editMatrixUsed = (editCostMatrix.size() != 0 && editCostMatrix.size() == G.upperNodeIdBound() );
    if(editMatrixUsed){
        //G = getGraphFromEditMatrix();
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
    case QuasiThresholdEditingLocalMover::DESC_DEGREE_INSERT: {
        dynamicForest = DynamicForest(std::vector<node>(G.upperNodeIdBound(), none));

        std::vector<std::vector<node>> buckets(G.numberOfNodes());
        G.forNodes([&](node u) { buckets[G.degree(u)].push_back(u); });
        this->order.reserve(G.numberOfNodes());
        for (auto it = buckets.rbegin(); it != buckets.rend(); ++it) {
            for (node u : *it) {
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
    editsWeightBefore = weightEdits;

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
    runningInfo["subtrees_moved"].push_back(0);
}

void EditingRunner::runLocalMover() {
    handler.assureRunning();
    if (!insertRun) {
        runningInfo["edits"].push_back(numEdits);
        runningInfo["edits_weight"].push_back(weightEdits);
    }

    generation = 0;
    //Main loop; check if iteration moved a node and if max iterations were used
    for (count i = insertRun ? 0 : 1; hasMoved && i <= maxIterations; ++i) {
        //check if iteration is on a plateau greater than maxPlateau size
        if (!hasMoved || (randomness && (currentPlateau >= maxPlateauSize)))
            break;
        handler.assureRunning();
        hasMoved = false;
        numNodesMoved = 0;
        numSubtreeMoves = 0;
        timer.start();

        for (auto &ec : event_counters) {
            ec.second.reset();
            ec.second.enable();
        }
        //first run when inserting is enabled. insert all nodes in assigned order
        if (insertRun) {
            for (index j = 0; j < G.numberOfNodes(); j++) {
                node nodeToMove = order[j];
                localMove(nodeToMove);
                existing[nodeToMove] = 1;
            }
            insertRun = 0;
        } else {
            //local move in random order
            G.forNodesInRandomOrder([&](node nodeToMove) { localMove(nodeToMove); });
            INFO("Iteration: ", i, " edits: ", numEdits, " moved nodes: ", numNodesMoved, " subtrees moved: ", numSubtreeMoves);
        }

        timer.stop();
        for (auto &ec : event_counters) {
            ec.second.disable();
        }
        //update running info, time nodes moved, edits
        if (i == 0) {
            runningInfo["nodes_moved"][0] = numNodesMoved;
            runningInfo["subtrees_moved"][0] = numSubtreeMoves;
            runningInfo["time"][0] += timer.elapsedMicroseconds();
            for (auto &ec : event_counters) {
                runningInfo[ec.first][0] += ec.second.readValue();
            }
        } else {
            runningInfo["nodes_moved"].push_back(numNodesMoved);
            runningInfo["subtrees_moved"].push_back(numSubtreeMoves);
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
        if (weightEdits == editsWeightBefore) {
            currentPlateau++;
        } else {
            if (currentPlateau > actualMaximumPlateau) {
                actualMaximumPlateau = currentPlateau;
            }
            currentPlateau = 0;
        }
        editsWeightBefore = weightEdits;
    }
}

Graph EditingRunner::getQuasiThresholdGraph() const {
    Graph forest = dynamicForest.toGraph();
    TreeReachabilityGraphGenerator gen(forest);
    gen.run();
    return gen.getGraph();
}

Graph EditingRunner::getDynamicForestGraph() const {
    Graph forest = dynamicForest.toGraph();
    return forest;
}

void EditingRunner::localMove(node nodeToMove) {
    generation ++;
    assert(numEdits == countNumberOfEdits());
    assert(weightEdits == countWeightOfEdits());
    INFO("Move node ", nodeToMove);
    handler.assureRunning();
    if(editMatrixUsed){
        editCostNodeToMove = editCostMatrix[nodeToMove];
    }
    numNeighbors = 0;
    editCostNeighbors = 0;
    G.forEdgesOf(nodeToMove, [&](node v) {
        if (!insertRun || existing[v]) {
            ++numNeighbors;
            if(editMatrixUsed){
                editCostNeighbors += editCostNodeToMove[v];
            }
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
    //TODO maybe replace to factor based on editCosts
    maxDepth = none;
    if(insertEditCost == 1 && removeEditCost == 1 && editMatrixUsed == false){
        maxDepth = 2 * numNeighbors;
    }   

    if(editMatrixUsed){
        curEditsWeight = editCostNeighbors;
    }
    else{
        curEditsWeight = removeEditCost * numNeighbors;
    }
    curEdits = numNeighbors;
    if (!insertRun) {
        // Calculate the old number of edits incident to c to be able to compute the improvement
        // later
        dynamicForest.dfsFrom(
            nodeToMove,
            [&](node c) {
                if (c != nodeToMove) {
                    //for all children in dfs insert one edge; no insertion necessary if child is neighbor of nodeToMove
                    curEdits += 1 - 2 * marker[c];
                    if(editMatrixUsed){
                        curEditsWeight -= editCostNodeToMove[c];
                    }
                    else{
                        curEditsWeight += 1 * insertEditCost - marker[c] * (insertEditCost + removeEditCost);
                    }
                }
            },
            [](node) {});
        dynamicForest.forAncestors(nodeToMove, [&](node p) {
            //for all ancestors insert one edge; no insertion necessary if ancestor is neighbor of nodeToMove
            curEdits += 1 - 2 * marker[p];
            if(editMatrixUsed){
                curEditsWeight -= editCostNodeToMove[p];
            }
            else{
                curEditsWeight += 1  * insertEditCost - marker[p] * (insertEditCost + removeEditCost);  
            }
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
            processNode(u, nodeToMove);
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
                processNode(u, nodeToMove);
            }
            //process neighbors not in current level but in queue
            while (!neighborQueue.empty() && dynamicForest.depth(neighborQueue.back()) == level) {
                node u = neighborQueue.back();
                neighborQueue.pop_back();
                assert(dynamicForest.depth(u) == level);
                if (nodeTouched[u])
                    continue; // if the node was touched in the previous level, it was in
                              // currentLevel and thus has already been processed
                processNode(u, nodeToMove);
            }
            --level;
            currentLevel.clear();
            currentLevel.swap(nextLevel);
        }
    }

    //save sum positive editWeights for root
    rootData.sumPositiveEdits = rootData.childCloseness;
    rootData.sumPositiveEditsWeight = rootData.childClosenessWeight;
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
    if(editMatrixUsed){
        bestEditsWeight = editCostNeighbors - (rootData.scoreMaxWeight);
    }
    else{
        bestEditsWeight = (numNeighbors * removeEditCost ) - (rootData.scoreMaxWeight);
    }


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
                        sample.clear();
                        for (node u : indifferentChildren) {
                            if (randomBool(2)) {
                                sample.push_back(u);
                            }
                        }
                    } while (sample.size() == 1);

                    for (node u : sample) {
                        bestChildren.push_back(u);
                    }
                }
            } else if (bestChildren.size() > 1) {
                // If there are already two children, just sample randomly from the remaining
                // indifferent children
                for (node u : indifferentChildren) {
                    if (randomBool(2)) {
                        bestChildren.push_back(u);
                    }
                }
            }
        }
    }

#ifndef NDEBUG
    compareWithQuadratic(nodeToMove);
#endif

    // calculate the number of saved edits as comparing the absolute number of edits doesn't make
    // sense
    int64_t savedEdits = curEdits - bestEdits;
    count savedEditsWeight = curEditsWeight - bestEditsWeight;

    assert(!randomness || rootData.hasChoices());

    if (rootData.logEqualBestChoices < std::log(std::numeric_limits<long long>::max())) {
        rootEqualBestParentsCpy = std::llround(std::exp(rootData.logEqualBestChoices));
    } else {
        rootEqualBestParentsCpy = std::numeric_limits<count>::max();
    }

    if (sortPaths || savedEditsWeight > 0 || randomness) {
        dynamicForest.moveToPosition(nodeToMove, rootData.bestParentBelow, bestChildren);
        hasMoved |= (savedEditsWeight > 0 || (randomness && rootEqualBestParentsCpy > 1));
        numNodesMoved += (savedEditsWeight > 0 || (randomness && rootEqualBestParentsCpy > 1));
        numEdits -= savedEdits;
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

    if(moveSubtrees){
        //increase generation by one to reuse childCloseness but reset everything else
        generation++;
        subtreeSize = 0;
        subtreeExtDegree = 0;
        editCostNeighbors = 0;
        curChildren = dynamicForest.children(nodeToMove);
        curParent = dynamicForest.parent(nodeToMove);

        //calculate subtree size and nodes in subtree
        dynamicForest.dfsFrom(nodeToMove, [&](node c) {
            subtreeSize++;
            inSubtree[c] = true;
            subtreeNodes.push_back(c);
		}, [](node) {});
        TRACE("Subtreesize:", subtreeSize);
        TRACE("Children: ", curChildren);
        TRACE("Parent: ", curParent);

    }

    //subtreeMove when subtree is more than one node
    if(moveSubtrees && subtreeSize > 1 && !insertRun){
        std::vector<int64_t> editCostNodeU = {};
        editCostNeighbors = 0;
        TRACE(dynamicForest.printPaths());
        //calculate editCost for Subtree
        dynamicForest.dfsFrom(nodeToMove, [&](node u) {
            if(editMatrixUsed){
                editCostNodeU = editCostMatrix[u];
            }
		    G.forNeighborsOf(u, [&](node v) {
				if (!inSubtree[v]) {
					++numNeighborsAll[v];
                    ++subtreeExtDegree;
                    if(editMatrixUsed){
                        editCostNeighbors += editCostNodeU[v];
                    }
				}
			});
		}, [](node) {});

        //add up editCost in matrix for nodes in subtree
        if(editMatrixUsed){
            editCostSubtree = editCostMatrix[nodeToMove];
            for ( node v : subtreeNodes){
                if(v != nodeToMove){
                    std::transform(editCostSubtree.begin(), editCostSubtree.end(), editCostMatrix[v].begin(),
                        editCostSubtree.begin(), std::plus<int64_t>());
                }
            }
        }

        if(editMatrixUsed){
            curEditsWeight = editCostNeighbors;
        }
        else{
            curEditsWeight = removeEditCost * subtreeExtDegree;
        }
        curEdits = subtreeExtDegree;
        if (!insertRun) {
            // Calculate the old number of edits incident to NodeToMove to be able to compute the improvement
            // later
            dynamicForest.forAncestors(nodeToMove, [&](node p) {
                //for all ancestors insert edge to subtree; no insertions necessary if ancestor is neighbor of node in subtree
                curEdits += subtreeSize - 2 * numNeighborsAll[p];
                if(editMatrixUsed){
                   curEditsWeight -= editCostSubtree[p];
                }
                else{
                    curEditsWeight += subtreeSize * insertEditCost - numNeighborsAll[p] * (insertEditCost + removeEditCost); 
                }
            });
            //not necessary for children because they are moved and edits dont change
        }

        //correct child closenesee above subtree
        int64_t sumChildCloseness = 0;
        int64_t sumChildClosenessWeight = 0;
        for( node c : curChildren){
            sumChildCloseness += traversalData[c].childCloseness;
            sumChildClosenessWeight += traversalData[c].childClosenessWeight;
            assert(dynamicForest.parent(c)== nodeToMove);
        }
        //sum childCloseness has to be positive because only then we adopt children
        assert(sumChildClosenessWeight >= 0);

        //correct rootdata
        //rootData.childCloseness -= sumChildCloseness;
        //rootData.childClosenessWeight -= sumChildClosenessWeight;
        if( sumChildClosenessWeight > rootData.sumPositiveEditsWeight){
            rootData.sumPositiveEdits = 0;
            rootData.sumPositiveEditsWeight = 0;
        } else {
            rootData.sumPositiveEdits -= sumChildCloseness;
            rootData.sumPositiveEditsWeight -= sumChildClosenessWeight;
        }
    
        dynamicForest.forAncestors(nodeToMove, [&](node p) {
            traversalData[p].childCloseness -= sumChildCloseness;
            traversalData[p].childClosenessWeight -= sumChildClosenessWeight;
            if( sumChildClosenessWeight > traversalData[p].sumPositiveEditsWeight){
                traversalData[p].sumPositiveEdits = 0;
                traversalData[p].sumPositiveEditsWeight = 0;
            } else {
                traversalData[p].sumPositiveEdits -= sumChildCloseness;
                traversalData[p].sumPositiveEditsWeight -= sumChildClosenessWeight;
            }
        });

        //add parent candidates to queue
        G.forNodes([&](node v) {
            if(!inSubtree[v]) {
                int64_t editCost = ((numNeighborsAll[v] * removeEditCost) + ((numNeighborsAll[v] - subtreeSize) * insertEditCost));
                //fill queue with candidates with positive editCost or positive childcloseness
                if ((traversalData[v].childClosenessWeight > 0 && nodeTouched[v] ) ||
                    (!editMatrixUsed && (editCost >= 0)) ||
                    (editMatrixUsed && editCostSubtree[v] >= 0)) {
                    parentCandidates.push_back(v);
                }
            }
        });

        if (useBucketQueue) {
            bucketQueue.fill(parentCandidates, dynamicForest);
        } else {
            //put all candidates in a queue and sort it based on forest depth
            for (node v : parentCandidates) {
                parentQueue.emplace_back(v);
            }
            std::stable_sort(parentQueue.begin(), parentQueue.end(), [&](node u, node v) {
                return dynamicForest.depth(u) < dynamicForest.depth(v);
            });
        }

        bestChildren.clear();
        rootData.initializeForSubtree(generation);

        if (useBucketQueue) {
            while (!bucketQueue.empty()) {
                node u = bucketQueue.next();
                processNodeForSubtree(u, nodeToMove);
            }
        } else {
            while (!currentLevel.empty() || !parentQueue.empty()) {
                //skip empty levels
                if (currentLevel.empty()) {
                   level = dynamicForest.depth(parentQueue.back());
                }
                //process nodes in currentLevel
                for (node u : currentLevel) {
                    assert(dynamicForest.depth(u) == level);
                   processNodeForSubtree(u, nodeToMove);
                }
                //process neighbors not in current level but in queue
                while (!parentQueue.empty() && dynamicForest.depth(parentQueue.back()) == level) {
                    node u = parentQueue.back();
                    parentQueue.pop_back();
                    assert(dynamicForest.depth(u) == level);
                    if (nodeTouchedSubtree[u])
                        continue; // if the node was touched in the previous level, it was in
                              // currentLevel and thus has already been processed
                    processNodeForSubtree(u, nodeToMove);
                }
                --level;
                currentLevel.clear();
                currentLevel.swap(nextLevel);
            }
        }

        if (!randomness) {
            if (rootData.sumPositiveEditsWeight > rootData.scoreMaxWeight) {
                rootData.bestParentBelow = none;
                rootData.scoreMax = rootData.sumPositiveEdits;
                rootData.scoreMaxWeight = rootData.sumPositiveEditsWeight;
            }
        } else {
            bool coin = false;
            double ownWeight = rootData.numIndifferentChildren * std::log(2);
            if (rootData.sumPositiveEditsWeight > rootData.scoreMaxWeight || !rootData.hasChoices()) {
                // INFO("root better");
                rootData.scoreMax = rootData.sumPositiveEdits;
                rootData.scoreMaxWeight = rootData.sumPositiveEditsWeight;
                rootData.logEqualBestChoices = ownWeight;
                coin = true;
            } else if (rootData.sumPositiveEditsWeight == rootData.scoreMaxWeight) {
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

        bestEdits = subtreeExtDegree - rootData.scoreMax;
        if(editMatrixUsed){
            bestEditsWeight = editCostNeighbors - (rootData.scoreMaxWeight);
        }
        else{
            bestEditsWeight = (removeEditCost * subtreeExtDegree) - (rootData.scoreMaxWeight);
        }

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
            assert(bestChildren.size() == bestParentData.numCloseChildren || curParent == bestParentData.bestParentBelow);

            if (randomness) { // make sure we adopt either 0 or at least two children
                assert(indifferentChildren.size() == bestParentData.numIndifferentChildren);
                // If we shall adopt one child, there must be another indifferent child that we just
                // sample randomly to get at least two
                if (bestChildren.size() == 1 && !indifferentChildren.empty()) {
                    //for subtreemove adopting one child should be ok
                    //assert(!indifferentChildren.empty());
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
                            sample.clear();
                            for (node u : indifferentChildren) {
                                if (randomBool(2)) {
                                    sample.push_back(u);
                                }
                            }
                        } while (sample.size() == 1);

                        for (node u : sample) {
                            bestChildren.push_back(u);
                        }
                    }
                } else if (bestChildren.size() > 1) {
                    // If there are already two children, just sample randomly from the remaining
                    // indifferent children
                    for (node u : indifferentChildren) {
                        if (randomBool(2)) {
                            bestChildren.push_back(u);
                        }
                    }
                }
            }
        }

        // calculate the number of saved edits as comparing the absolute number of edits doesn't make
        // sense
        int64_t savedEdits = curEdits - bestEdits;
        count savedEditsWeight = curEditsWeight - bestEditsWeight;

        assert(!randomness || rootData.hasChoices());

        if (rootData.logEqualBestChoices < std::log(std::numeric_limits<long long>::max())) {
            rootEqualBestParentsCpy = std::llround(std::exp(rootData.logEqualBestChoices));
        } else {
            rootEqualBestParentsCpy = std::numeric_limits<count>::max();
        }

        TRACE("Best Parent for Subtree ", rootData.bestParentBelow);
        TRACE("Best adopted Children for Subtree ", bestChildren);
        if (savedEditsWeight > 0 || randomness) {
            std::vector<node> vect{ nodeToMove };
            TRACE("SubtreeMove: current Parent: ", curParent, " new Parent: ", rootData.bestParentBelow);
            dynamicForest.moveToAnyPosition(rootData.bestParentBelow, vect);
            dynamicForest.moveToAnyPosition(nodeToMove, bestChildren);
            hasMoved |= (savedEditsWeight > 0 || (randomness && rootEqualBestParentsCpy > 1));
            numNodesMoved += (savedEditsWeight > 0 || (randomness && rootEqualBestParentsCpy > 1));
            numSubtreeMoves += (savedEditsWeight > 0 || (randomness && rootEqualBestParentsCpy > 1));
            numEdits -= savedEdits;
            weightEdits -= savedEditsWeight;
            INFO("Saved Edits: ", savedEdits, ", SavedWeightedEdits: ", savedEditsWeight);
#ifndef NDEBUG
            assert(numEdits == countNumberOfEdits());
            assert(weightEdits == countWeightOfEdits());
#endif
        }
    }
    // cleanup for linear move
    for (node u : touchedNodes) {
        lastVisitedDFSNode[u] = u;
        nodeTouched[u] = false;
        nodeTouchedSubtree[u] = false;
    }
    for (node v : neighbors) {
        marker[v] = false;
    }
    for (node v : subtreeNodes) {
        inSubtree[v] = false;
    }
    neighbors.clear();
    parentCandidates.clear();
    touchedNodes.clear();
    subtreeNodes.clear();
    //TODO Optimize: only reset necessary values
    if(moveSubtrees){
        std::fill(numNeighborsAll.begin(), numNeighborsAll.end(), 0);
    }

}

void EditingRunner::processNode(node u, node nodeToMove) {
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
        // only marked neighbors may be processed without having been touched before
        assert(marker[u]); 
    }

    traversalData[u].initialize(generation);

    //save sum of positive edits for subtreemove
    int64_t sumPositiveEdits = traversalData[u].childCloseness;
    int64_t sumPositiveEditsWeight = traversalData[u].childClosenessWeight;
    traversalData[u].sumPositiveEdits = traversalData[u].childCloseness;
    traversalData[u].sumPositiveEditsWeight = traversalData[u].childClosenessWeight;
    assert(traversalData[u].childClosenessWeight >= 0);

    if (marker[u]) { 
        ++traversalData[u].childCloseness;
    } else {
        --traversalData[u].childCloseness; 
    }
    if(editMatrixUsed){
        //increase childCloseness if neighbor, decrease if not neighbor
        traversalData[u].childClosenessWeight +=  editCostNodeToMove[u];
    }
    else{
        traversalData[u].childClosenessWeight += marker[u] * removeEditCost;
        traversalData[u].childClosenessWeight -=
        1 * insertEditCost - marker[u] * insertEditCost;
    }


    TRACE("Edit difference before descending: ", traversalData[u].childCloseness, ", ", traversalData[u].childClosenessWeight);

    assert(!marker[u] || traversalData[u].childClosenessWeight > 0);

    if (traversalData[u].childClosenessWeight >= 0) {
        assert(lastVisitedDFSNode[u] == u);

        node c = dynamicForest.nextDFSNodeOnEnter(u, u);

        while (c != u) {
            //resturcture if 
            if (!nodeTouched[c] || traversalData[c].childClosenessWeight < 0) {
                if (nodeTouched[c] ){
                    assert(traversalData[c].childClosenessWeight < 0);
                    traversalData[u].childClosenessWeight += traversalData[c].childClosenessWeight;
                    traversalData[u].childCloseness+= traversalData[c].childCloseness;
                    // advance to the next starting point for the DFS search.
                    c = lastVisitedDFSNode[c];
                }
                else{
                    if(editMatrixUsed){
                        //decrease childCloseness because editCost is here negative
                        traversalData[u].childClosenessWeight += editCostNodeToMove[c];
                    }
                    else{
                        traversalData[u].childClosenessWeight -= insertEditCost;
                    }
                    --traversalData[u].childCloseness;
                }

                if (traversalData[u].childClosenessWeight < 0 || dynamicForest.depth(c) > maxDepth) {
                    lastVisitedDFSNode[u] = c;
                    break;
                }

                c = dynamicForest.nextDFSNodeOnEnter(c, u);
            } else {
                if(traversalData[c].childClosenessWeight == 0){
                    traversalData[u].childCloseness += traversalData[c].childCloseness;
                }
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
            traversalData[u].scoreMax = sumPositiveEdits;
            traversalData[u].scoreMaxWeight = sumPositiveEditsWeight;
            traversalData[u].bestParentBelow = u;
        }
    } else {
        bool coin = false;
        double ownWeight = traversalData[u].numIndifferentChildren * std::log(2);
        if (sumPositiveEditsWeight > traversalData[u].scoreMaxWeight || !traversalData[u].hasChoices()) {
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
        }
        if (coin) {
            traversalData[u].bestParentBelow = u;
        }
    }

    assert(traversalData[u].scoreMaxWeight >= 0);

    if(marker[u]){
        ++traversalData[u].scoreMax;
    } else {
        --traversalData[u].scoreMax;
    }

    if(editMatrixUsed){
        traversalData[u].scoreMaxWeight += editCostNodeToMove[u];     
    }
    else{
        traversalData[u].scoreMaxWeight += marker[u] * removeEditCost;
        traversalData[u].scoreMaxWeight -= 1 * insertEditCost - marker[u] * insertEditCost;
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
    } else if (randomness && traversalData[u].scoreMaxWeight == parentData.scoreMaxWeight) {
        parentData.addLogChoices(traversalData[u].logEqualBestChoices);
        if (logRandomBool(traversalData[u].logEqualBestChoices - parentData.logEqualBestChoices)) {
            parentData.bestParentBelow = traversalData[u].bestParentBelow;
        }
    }
    
    //propagate positive childCloseness
    if (traversalData[u].childClosenessWeight >= 0) {
        assert(traversalData[u].childClosenessWeight <= traversalData[u].scoreMaxWeight);
        if (traversalData[u].childClosenessWeight == 0) {
            ++parentData.numIndifferentChildren;
        } else {
            ++parentData.numCloseChildren;
            parentData.childCloseness += traversalData[u].childCloseness;
            parentData.childClosenessWeight += traversalData[u].childClosenessWeight;
        }
    }

    assert(!dynamicForest.children(u).empty() || traversalData[u].childCloseness == 1);
}

void EditingRunner::processNodeForSubtree(node u, node nodeToMove) {
    TRACE("Process for Subtree ", u);
    TRACE("Processing node for Subtree ", u, " of depth ", dynamicForest.depth(u),
          " (node to move: ", nodeToMove, ")");
    TRACE("Parent: ", dynamicForest.parent(u), ", children: ", dynamicForest.children(u));
    assert(u != nodeToMove);
    tlx::unused(nodeToMove);
    if (useBucketQueue) {
        assert(dynamicForest.depth(u) <= maxDepth);
    }
    if (!nodeTouchedSubtree[u]) {
        nodeTouchedSubtree[u] = true;
        if(!nodeTouched[u]){
            touchedNodes.emplace_back(u);
        }
        assert( (traversalData[u].childClosenessWeight > 0 && nodeTouched[u]) || 
        (editMatrixUsed && editCostSubtree[u] >= 0) || 
        (!editMatrixUsed && ((2 * numNeighborsAll[u] - subtreeSize) * insertEditCost) >= 0)); 
        // only node with editCosts >= 0 are not touched and are processed
    }

    traversalData[u].initializeForSubtree(generation);

    //required edits to this node for subtree
    traversalData[u].subtreeEdits = (2 * numNeighborsAll[u] - subtreeSize);
    if( editMatrixUsed){
        traversalData[u].subtreeEditCosts = editCostSubtree[u];
    } else {
        traversalData[u].subtreeEditCosts = (numNeighborsAll[u] * removeEditCost) + ((numNeighborsAll[u] - subtreeSize) * insertEditCost);
    }

    //sum of Positive childCloseness below node
    int64_t sumPositiveEdits = traversalData[u].sumPositiveEdits;
    int64_t sumPositiveEditsWeight = traversalData[u].sumPositiveEditsWeight;
    assert(sumPositiveEditsWeight >=0);

    if (!randomness) {
        if (sumPositiveEditsWeight > traversalData[u].scoreMaxWeight || traversalData[u].scoreMaxWeight == 0) {
            traversalData[u].scoreMax = sumPositiveEdits;
            traversalData[u].scoreMaxWeight = sumPositiveEditsWeight;
            traversalData[u].bestParentBelow = u;
        }
    } else {
        bool coin = false;
        double ownWeight = traversalData[u].numIndifferentChildren * std::log(2);
        if (sumPositiveEditsWeight > traversalData[u].scoreMaxWeight || !traversalData[u].hasChoices()) {
            traversalData[u].scoreMax = sumPositiveEdits;
            traversalData[u].scoreMaxWeight = sumPositiveEditsWeight;
            traversalData[u].logEqualBestChoices = ownWeight;
            // Either we do not adopt children, or we are at the lower end of a path.
            // Otherwise, there must be a node below u that is at least as good.
            //TODO check is that correct because we dont sort path before subtreeMove
            //assert(!sortPaths || sumPositiveEdits == 0 || dynamicForest.isLowerEnd(u));
            coin = true;
        } else if (sumPositiveEditsWeight == traversalData[u].scoreMaxWeight) {
            ownWeight = traversalData[u].calculateOwnWeightForEqualChoices();
            if (ownWeight > -std::numeric_limits<double>::infinity()) {
                traversalData[u].addLogChoices(ownWeight);
                coin = logRandomBool(ownWeight - traversalData[u].logEqualBestChoices);
            }
            assert(traversalData[u].hasChoices());
        }
        if (coin) {
            traversalData[u].bestParentBelow = u;
        }
    }

    //not the case if we have a node with childCloseness < 0 and subtreeEditCost  = 0
    //assert(traversalData[u].scoreMaxWeight >= 0);

    traversalData[u].scoreMax += (2 * numNeighborsAll[u] - subtreeSize);
    if(editMatrixUsed){
        traversalData[u].scoreMaxWeight += editCostSubtree[u];      
    }
    else{
        traversalData[u].scoreMaxWeight += (numNeighborsAll[u] * removeEditCost) + ((numNeighborsAll[u] - subtreeSize) * insertEditCost);
    }

    TRACE("Maximum score at ", u, ": ", traversalData[u].scoreMax, ", ", traversalData[u].scoreMaxWeight);
    node p = dynamicForest.parent(u);
    TraversalData &parentData = (p == none) ? rootData : traversalData[p];

    parentData.initializeForSubtree(generation);
    int64_t editCostParent = ((numNeighborsAll[p] * removeEditCost) + ((numNeighborsAll[p] - subtreeSize) * insertEditCost));

    if ((traversalData[u].scoreMaxWeight > 0 || traversalData[u].subtreeEditCosts > 0) && p != none) {
        if (useBucketQueue) {
            assert(dynamicForest.depth(p) <= maxDepth);
        }
        if (!nodeTouchedSubtree[p]) {
            nodeTouchedSubtree[p] = true;
            if(!nodeTouched[p]){
                touchedNodes.push_back(p);
            }
            if (!useBucketQueue) {
                nextLevel.push_back(p);
            } else if (!(traversalData[p].childClosenessWeight > 0 && nodeTouched[p]) && 
                !(editMatrixUsed && editCostSubtree[u] >= 0) && 
                !(!editMatrixUsed && editCostParent >= 0)) { 
                // positive editcost already in queue
                bucketQueue.insertParent(p);
            }
        }
    }

    if (traversalData[u].scoreMaxWeight > parentData.scoreMaxWeight) {
        parentData.logEqualBestChoices = traversalData[u].logEqualBestChoices;
        parentData.scoreMax = traversalData[u].scoreMax;
        parentData.scoreMaxWeight = traversalData[u].scoreMaxWeight;
        parentData.bestParentBelow = traversalData[u].bestParentBelow;
    } else if (randomness && traversalData[u].scoreMaxWeight == parentData.scoreMaxWeight) {
        parentData.addLogChoices(traversalData[u].logEqualBestChoices);
        if (logRandomBool(traversalData[u].logEqualBestChoices - parentData.logEqualBestChoices)) {
            parentData.bestParentBelow = traversalData[u].bestParentBelow;
        }
    }

    if (traversalData[u].childClosenessWeight >= 0) {
        //not true in this case: positive childCloseness when node neighbor of nodeToMove but  not too subtreenodes
        //assert(traversalData[u].childClosenessWeight <= traversalData[u].scoreMaxWeight);
        if (traversalData[u].childClosenessWeight == 0) {
            ++parentData.numIndifferentChildren;
        } else {
            ++parentData.numCloseChildren;
            //dont propagate childcloseness, already happend in local move
            //parentData.childCloseness += traversalData[u].childCloseness;
            //parentData.childClosenessWeight += traversalData[u].childClosenessWeight;
        }
    }
}

void EditingRunner::compareWithQuadratic(node nodeToMove) const {
    std::vector<int64_t> missingBelow, missingAbove, existingBelow, existingAbove;
    std::vector<int64_t> missingBelowWeighted, missingAboveWeighted, existingBelowWeighted, existingAboveWeighted;
    missingBelow.resize(G.upperNodeIdBound(), 0);
    missingAbove.resize(G.upperNodeIdBound(), 0);
    existingBelow.resize(G.upperNodeIdBound(), 0);
    existingAbove.resize(G.upperNodeIdBound(), 0);
    missingBelowWeighted.resize(G.upperNodeIdBound(), 0);
    missingAboveWeighted.resize(G.upperNodeIdBound(), 0);
    existingBelowWeighted.resize(G.upperNodeIdBound(), 0);
    existingAboveWeighted.resize(G.upperNodeIdBound(), 0);
    std::vector<bool> usingDeepNeighbors(G.upperNodeIdBound(), false);
    if (editMatrixUsed)
    {
        std::vector<int64_t> editCostNodeToMove = editCostMatrix[nodeToMove];
    }
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
                        if(editMatrixUsed){
                            missingBelowWeighted[u] = missingAboveWeighted[u] = 1 * (-1) * editCostNodeToMove[u] - marker[u] * (-1) * editCostNodeToMove[u];
                            existingBelowWeighted[u] = existingAboveWeighted[u] = marker[u] * editCostNodeToMove[u];
                        }
                        else{
                            missingBelowWeighted[u] = missingAboveWeighted[u] = 1 * insertEditCost - marker[u] * insertEditCost;
                            existingBelowWeighted[u] = existingAboveWeighted[u] = marker[u] * removeEditCost;
                        }
                    }
                    node p = dynamicForest.parent(u);
                    if (p != none) {
                        missingAbove[u] += missingAbove[p];
                        existingAbove[u] += existingAbove[p];
                        missingAboveWeighted[u] += missingAboveWeighted[p];
                        existingAboveWeighted[u] += existingAboveWeighted[p];
                    }
                },
                [&](node u) {
                    node p = dynamicForest.parent(u);
                    if (p != none) {
                        missingBelow[p] += missingBelow[u];
                        existingBelow[p] += existingBelow[u];
                        missingBelowWeighted[p] += missingBelowWeighted[u];
                        existingBelowWeighted[p] += existingBelowWeighted[u];
                        if (usingDeepNeighbors[u])
                            usingDeepNeighbors[p] = true;
                    }
                });
        }
    });

    assert(missingBelow[nodeToMove] == 0);
    assert(existingBelow[nodeToMove] == 0);
    assert(missingBelowWeighted[nodeToMove] == 0);
    assert(existingBelowWeighted[nodeToMove] == 0);

    if (!sortPaths) {
        bool exactValue = true;
        for (node c : curChildren) {
            missingBelow[nodeToMove] += missingBelow[c];
            existingBelow[nodeToMove] += existingBelow[c];
            missingBelowWeighted[nodeToMove] += missingBelowWeighted[c];
            existingBelowWeighted[nodeToMove] += existingBelowWeighted[c];
            if (usingDeepNeighbors[c])
                exactValue = false;
        }

        if (curParent != none) {
            missingAbove[nodeToMove] = missingAbove[curParent];
            existingAbove[nodeToMove] = existingAbove[curParent];
            missingAboveWeighted[nodeToMove] = missingAboveWeighted[curParent];
            existingAboveWeighted[nodeToMove] = existingAboveWeighted[curParent];
            if (usingDeepNeighbors[curParent])
                exactValue = false;
        }
        if (exactValue && !editMatrixUsed) {
            assert(curEdits
                   == numNeighbors - existingAbove[nodeToMove] - existingBelow[nodeToMove]
                          + missingAbove[nodeToMove] + missingBelow[nodeToMove]);
            assert(curEditsWeight
                   == numNeighbors * removeEditCost - existingAboveWeighted[nodeToMove] - existingBelowWeighted[nodeToMove]
                          + missingAboveWeighted[nodeToMove] + missingBelowWeighted[nodeToMove]);
        }
    }

    count minEdits = std::numeric_limits<count>::max();
    count minEditsWeight = std::numeric_limits<count>::max();
    count editsWeightOffset = 0;
    std::vector<node> minChildren;
    node minParent = curParent;
    G.forNodes([&](node u) {
        if (u == nodeToMove || usingDeepNeighbors[u] || !existing[u])
            return;
        if (existingBelowWeighted[u] >= missingBelowWeighted[u]
            || (traversalData[u].generation == generation && traversalData[u].childClosenessWeight > 0)) {
            assert(traversalData[u].childClosenessWeight == existingBelowWeighted[u] - missingBelowWeighted[u]);
            //
            assert(traversalData[u].childCloseness == existingBelow[u] - missingBelow[u]);
        } else if (nodeTouched[u]) {
            assert(traversalData[u].childClosenessWeight < 0);
        }
    });

    G.forNodes([&](node u) {
        if (dynamicForest.children(u).empty() && marker[u] && !usingDeepNeighbors[u]
            && existing[u]) {
            assert(traversalData[u].childCloseness == 1);
        }
        if(marker[u] && editMatrixUsed)
        {       
            assert(editCostNodeToMove[u]>=0);
            editsWeightOffset += editCostNodeToMove[u];
        }
    });

    auto tryEditBelow = [&](node p) {
        if (p == nodeToMove)
            return;
        count editsWeight = editsWeightOffset;
        count edits = numNeighbors;
        if(!editMatrixUsed)
        {       
            editsWeight = numNeighbors * removeEditCost;
        }
        if (p != none) {
            edits += missingAbove[p];
            edits -= existingAbove[p];
            editsWeight += missingAboveWeighted[p];
            editsWeight -= existingAboveWeighted[p];
        }

        std::vector<node> children;
        dynamicForest.forChildrenOf(p, [&](node c) {
            if (c == nodeToMove || usingDeepNeighbors[c] || !existing[c])
                return;
            if (existingBelowWeighted[c] > missingBelowWeighted[c]) { // TODO try >= (more children...)
                if (dynamicForest.children(c).empty() && marker[c]) {
                    assert(traversalData[c].childCloseness == 1);
                    if(editMatrixUsed){
                        assert(traversalData[c].childClosenessWeight == 1 * editCostNodeToMove[c]);
                    }
                    else{
                        assert(traversalData[c].childClosenessWeight == 1 * removeEditCost);
                    }

                }
                assert(traversalData[c].childCloseness == existingBelow[c] - missingBelow[c]);//
                assert(traversalData[c].childClosenessWeight == existingBelowWeighted[c] - missingBelowWeighted[c]);

                children.emplace_back(c);
                edits -= existingBelow[c] - missingBelow[c];
                editsWeight -= existingBelowWeighted[c] - missingBelowWeighted[c];
            }
        });

        if (editsWeight < minEditsWeight) {
            minEdits = edits;
            minEditsWeight = editsWeight;
            minChildren = std::move(children);
            minParent = p;
        }
    };

    dynamicForest.dfsFrom(
        none, [](node) {}, tryEditBelow);
    tryEditBelow(none);
    //correct assertion? assert(minEdits >= bestEdits);
    //assert(minEdits <= bestEdits);
    assert(minEditsWeight == bestEditsWeight);

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
    //correct assertion? assert(minEdits >= childClosenessControl);
    assert(minEdits <= childClosenessControl);

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
    count weightOfEdits = 0;
    count weightOfRemoveEdits = 0;
    count weightOfInsertEdits = 0;
    count numExistingEdges = 0;
    count numMissingEdges = 0;
    std::vector<int64_t> editCostNodeU = {};
    std::vector<int64_t> editCostNodeP = {};
    std::vector<bool> marker(G.upperNodeIdBound());

    dynamicForest.forChildrenOf(none, [&](node r) {
        count depth = 0;
        dynamicForest.dfsFrom(
            r,
            [&](node u) { // on enter
                count upperNeighbors = 0;

                if(editMatrixUsed){
                    editCostNodeU = editCostMatrix[u];
                }
                G.forNeighborsOf(u, [&](node v) {
                    if (marker[v]){
                        ++upperNeighbors;
                        //-remove
                        if(editMatrixUsed){
                            weightOfRemoveEdits -= editCostNodeU[v];
                        }
                    }
                    else if(editMatrixUsed) {
                        //+remove
                        weightOfRemoveEdits += editCostNodeU[v];
                    }
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
    if(editMatrixUsed){
        G.forNodes([&](node u) {
            editCostNodeU = editCostMatrix[u];
            dynamicForest.forAncestors(u, [&](node p) {
                if(editCostNodeU[p] < 0 && u != p){
                    weightOfInsertEdits += (-1) * editCostNodeU[p];
                }
            });
        });
    }
    TRACE("Missing Edges:  ", numMissingEdges,", NumEdges-ExistingEdges: " ,(G.numberOfEdges() - numExistingEdges));
    if(editMatrixUsed){
        weightOfEdits = (weightOfInsertEdits) + (weightOfRemoveEdits / 2);
    }
    else{
        weightOfEdits = ((numMissingEdges * insertEditCost) + ((G.numberOfEdges() - numExistingEdges) * removeEditCost));
    }
    return weightOfEdits;
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
Graph EditingRunner::getGraphFromEditMatrix() {
    Graph editGraph = Graph(editCostMatrix.size(),false,false);
    std::vector<int64_t> editCostNodeU = {};
    NetworKit::Unsafe unsafe;
    for( count u; u < editCostMatrix.size(); u++){
        editCostNodeU = editCostMatrix[u];
        for( count v; v < editCostNodeU.size(); v++){
            if(editCostNodeU[v] > 0){
                editGraph.addPartialEdge(unsafe, u,v);
            }
        }
    }
    return editGraph;
}
} // namespace QuasiThresholdMoving

} // namespace NetworKit
