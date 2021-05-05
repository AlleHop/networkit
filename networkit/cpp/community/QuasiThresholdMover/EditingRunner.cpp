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
                             count maxIterations, bool sortPaths, bool randomness, bool moveSubtrees, bool subtreeSortPaths,
                             count maxPlateauSize, bool useBucketQueue, std::vector<node> order,
                             count insertEditCost, count removeEditCost, const std::vector<std::vector<int64_t>> &editCostMatrix)
    : G(G), maxIterations(maxIterations), usedIterations(0), sortPaths(sortPaths),
      randomness(randomness), moveSubtrees(moveSubtrees), subtreeSortPaths(subtreeSortPaths), maxPlateauSize(maxPlateauSize),
      insertRun(initialization != QuasiThresholdEditingLocalMover::TRIVIAL
                && initialization != QuasiThresholdEditingLocalMover::EDITING
                && initialization != QuasiThresholdEditingLocalMover::RANDOM_TREE),
      useBucketQueue(useBucketQueue
      && !(editCostMatrix.size() != 0 && editCostMatrix.size() == G.upperNodeIdBound()) 
      && (insertEditCost == 1 && removeEditCost ==1)), 
      insertEditCost(insertEditCost), removeEditCost(removeEditCost), editCostMatrix(editCostMatrix), handler(), hasMoved(true),
      marker(G.upperNodeIdBound(), false), inSubtree(G.upperNodeIdBound(), false),  numNeighborsAll(G.upperNodeIdBound(), 0),lastVisitedDFSNode(G.upperNodeIdBound(), none),
      traversalData(G.upperNodeIdBound()), nodeTouched(G.upperNodeIdBound(), false), rootData(),
      existing(G.upperNodeIdBound(), !insertRun), rootEqualBestParentsCpy(0), currentPlateau(0),
      actualMaximumPlateau(0), gen(Aux::Random::getURNG()), realDist(), intDist(0, 1) {

    runningInfo["time"] = std::vector<count>();
    runningInfo["edits"] = std::vector<count>();
    runningInfo["edit_costs"] = std::vector<count>();
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
    /*if(editMatrixUsed){
        G = getGraphFromEditMatrix();
    }*/

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
    case QuasiThresholdEditingLocalMover::RANDOM_TREE: {
        std::vector<node> parents;
        parents = getParentsForTree();
        if (G.upperNodeIdBound() == G.numberOfNodes()) {
            dynamicForest = DynamicForest(parents);
        } else {
            dynamicForest = DynamicForest(G, parents);
        }
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

    ratioParentCandidates = 0;
    ratioSubtreeNeighbors = 0;
    countSubtreeMoves = 0;

    numEdits = countNumberOfEdits();
    costEdits = countCostOfEdits();
    editsBefore = numEdits;
    editCostsBefore = costEdits;

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
        runningInfo["edit_costs"].push_back(costEdits);
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
            INFO("RatioSubtreeNeighbors: ", ratioSubtreeNeighbors, " RatioParentCandidates: ", ratioParentCandidates);
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
        runningInfo["edit_costs"].push_back(costEdits);
        usedIterations = i;

        assert(numEdits == countNumberOfEdits());
        assert(costEdits == countCostOfEdits());
        //update plateau
        if (costEdits == editCostsBefore) {
            currentPlateau++;
        } else {
            if (currentPlateau > actualMaximumPlateau) {
                actualMaximumPlateau = currentPlateau;
            }
            currentPlateau = 0;
        }
        editCostsBefore = costEdits;
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
    assert(costEdits == countCostOfEdits());
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
        curEditCosts = editCostNeighbors;
    }
    else{
        curEditCosts = removeEditCost * numNeighbors;
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
                        curEditCosts -= editCostNodeToMove[c];
                    }
                    else{
                        curEditCosts += 1 * insertEditCost - marker[c] * (insertEditCost + removeEditCost);
                    }
                }
                return DynamicForest::ReturnState::CONTINUE;
            },
            [](node) {return DynamicForest::ReturnState::CONTINUE;});
        dynamicForest.forAncestors(nodeToMove, [&](node p) {
            //for all ancestors insert one edge; no insertion necessary if ancestor is neighbor of nodeToMove
            curEdits += 1 - 2 * marker[p];
            if(editMatrixUsed){
                curEditCosts -= editCostNodeToMove[p];
            }
            else{
                curEditCosts += 1  * insertEditCost - marker[p] * (insertEditCost + removeEditCost);  
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

    //save sum positive editCosts for root
    rootData.sumPositiveEdits = rootData.childCloseness;
    rootData.sumPositiveEditCosts = rootData.childClosenessCost;
    if (!randomness) {
        if (rootData.sumPositiveEditCosts > rootData.scoreMaxCost) {
            rootData.bestParentBelow = none;
            rootData.scoreMax = rootData.sumPositiveEdits;
            rootData.scoreMaxCost = rootData.sumPositiveEditCosts;
        }
    } else {
        bool coin = false;
        double ownWeight = rootData.numIndifferentChildren * std::log(2);
        if (rootData.sumPositiveEditCosts > rootData.scoreMaxCost || !rootData.hasChoices()) {
            // INFO("root better");
            rootData.scoreMax = rootData.sumPositiveEdits;
            rootData.scoreMaxCost = rootData.sumPositiveEditCosts;
            rootData.logEqualBestChoices = ownWeight;
            coin = true;
        } else if (rootData.sumPositiveEditCosts == rootData.scoreMaxCost) {
            ownWeight = rootData.calculateOwnWeightForEqualChoices();
            if (ownWeight > -std::numeric_limits<double>::infinity()) {
                rootData.addLogChoices(ownWeight);
                coin = logRandomBool(ownWeight - rootData.logEqualBestChoices);
            }
            // INFO("root equally good");
        }
        if (coin) {
            //correct scoremax incase root is equally good for weighted case
            rootData.scoreMax = rootData.sumPositiveEdits;
            rootData.bestParentBelow = none;
        }

        assert(rootData.hasChoices());
    }

    bestEdits = numNeighbors - rootData.scoreMax;
    if(editMatrixUsed){
        bestEditCosts = editCostNeighbors - (rootData.scoreMaxCost);
    }
    else{
        bestEditCosts = (numNeighbors * removeEditCost ) - (rootData.scoreMaxCost);
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
                if (traversalData[u].childClosenessCost > 0) {
                    bestChildren.push_back(u);
                } else if (randomness && traversalData[u].childClosenessCost == 0) {
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
                //child closeness is not necessarily 0 if weighted childcloseness is 0
                bestEdits -= traversalData[indifferentChildren[i]].childCloseness;
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
                            //child closeness is not necessarily 0 if weighted childcloseness is 0
                            bestEdits -= traversalData[u].childCloseness;
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
                        //child closeness is not necessarily 0 if weighted childcloseness is 0
                        bestEdits -= traversalData[u].childCloseness;
                    }
                }
            } else if (bestChildren.size() > 1) {
                // If there are already two children, just sample randomly from the remaining
                // indifferent children
                for (node u : indifferentChildren) {
                    if (randomBool(2)) {
                        bestChildren.push_back(u);
                        //child closeness is not necessarily 0 if weighted childcloseness is 0
                        bestEdits -= traversalData[u].childCloseness;
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
    savedEdits = curEdits - bestEdits;
    savedEditCosts = curEditCosts - bestEditCosts;

    assert(!randomness || rootData.hasChoices());

    if (rootData.logEqualBestChoices < std::log(std::numeric_limits<long long>::max())) {
        rootEqualBestParentsCpy = std::llround(std::exp(rootData.logEqualBestChoices));
    } else {
        rootEqualBestParentsCpy = std::numeric_limits<count>::max();
    }

    if (sortPaths || savedEditCosts > 0 || randomness) {
        dynamicForest.moveToPosition(nodeToMove, rootData.bestParentBelow, bestChildren);
        hasMoved |= (savedEditCosts > 0 || (randomness && rootEqualBestParentsCpy > 1));
        numNodesMoved += (savedEditCosts > 0 || (randomness && rootEqualBestParentsCpy > 1));
        numEdits -= savedEdits;
        costEdits -= savedEditCosts;
        INFO("Saved Edits: ", savedEdits, ", Saved EditCosts: ", savedEditCosts);
#ifndef NDEBUG
        assert(numEdits == countNumberOfEdits());
        assert(costEdits == countCostOfEdits());
#endif
    } else {
        dynamicForest.moveToPosition(nodeToMove, curParent, curChildren);
#ifndef NDEBUG
        assert(numEdits == countNumberOfEdits());
        assert(costEdits == countCostOfEdits());
#endif
    }
    // cleanup for linear move
    for (node u : touchedNodes) {
        lastVisitedDFSNode[u] = u;
        nodeTouched[u] = false;
    }
    touchedNodes.clear();

    if(moveSubtrees && !insertRun){
        subtreeMove(nodeToMove);
    }
    // cleanup for subtreemove move
    for (node u : touchedNodes) {
        lastVisitedDFSNode[u] = u;
        nodeTouched[u] = false;
    }
    for (node v : neighbors) {
        marker[v] = false;
    }
    neighbors.clear();
    touchedNodes.clear();
}

void EditingRunner::subtreeMove(node nodeToMove){
    savedEdits =  0;
    savedEditCosts = 0;
    for( subtreeOption = 0; subtreeOption < 4 && savedEditCosts == 0 ; subtreeOption ++){
        subtreeSize = 0;
        bestEdits = 0;
        bestEditCosts = 0;
        //subtreeOption = 1;
        if(subtreeOption == 2 || subtreeOption == 3){
            dynamicForest.moveNodeToLowerEnd(nodeToMove, marker);
        }
        curChildren = dynamicForest.children(nodeToMove);
        curParent = dynamicForest.parent(nodeToMove);
        if((subtreeOption == 1 || subtreeOption == 3)&& !curChildren.empty()){
            std::shuffle(curChildren.begin(), curChildren.end(), Aux::Random::getURNG());
            curChildren.pop_back();
            dynamicForest.moveToAnyPosition(curParent, curChildren);
            if (subtreeOption == 1 || subtreeOption == 3) {
                for(node v : curChildren){
                bestEdits += traversalData[v].childCloseness;
                bestEditCosts += traversalData[v].childClosenessCost;
                }
            }
        }
        //calculate subtree size and nodes in subtree
        dynamicForest.dfsFrom(nodeToMove, [&](node c) {
            subtreeSize++;
            inSubtree[c] = true;
            subtreeNodes.push_back(c);
            if(subtreeSize > std::max(50.0,sqrt(G.numberOfNodes()))){
                return DynamicForest::ReturnState::BREAK;
            } else {
                return DynamicForest::ReturnState::CONTINUE;
            }
		}, [](node) {return DynamicForest::ReturnState::CONTINUE;});
        TRACE("Subtreesize:", subtreeSize);
        if(subtreeSize > std::max(50.0,sqrt(G.numberOfNodes()))){
            dynamicForest.moveToAnyPosition(nodeToMove, curChildren);
        }

    //subtreeMove when subtree is more than one node
    if(subtreeSize > 1  && subtreeSize <= std::max(50.0,sqrt(G.numberOfNodes()))){
        generation++;
        subtreeExtDegree = 0;
        editCostNeighbors = 0;
        TRACE("Children: ", curChildren);
        TRACE("Parent: ", curParent);
        editCostNeighbors = 0;
        //calculate editCost for Subtree
        dynamicForest.dfsFrom(nodeToMove, [&](node u) {
		    G.forEdgesOf(u, [&](node v) {
				if (!inSubtree[v]) {
					++numNeighborsAll[v];
                    ++subtreeExtDegree;
                    if(numNeighborsAll[v] == 1){
                        subtreeNeighbors.push_back(v);
                    }
                    if(editMatrixUsed){
                        editCostNeighbors += editCostMatrix[u][v];
                    }
				}
			});
            return DynamicForest::ReturnState::CONTINUE;
		}, [](node) {return DynamicForest::ReturnState::CONTINUE;});

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
            curEditCosts = editCostNeighbors;
        }
        else{
            curEditCosts = removeEditCost * subtreeExtDegree;
        }
        curEdits = subtreeExtDegree;
        if (!insertRun) {
            // Calculate the old number of edits incident to NodeToMove to be able to compute the improvement
            // later
            dynamicForest.forAncestors(nodeToMove, [&](node p) {
                //for all ancestors insert edge to subtree; no insertions necessary if ancestor is neighbor of node in subtree
                curEdits += subtreeSize - 2 * numNeighborsAll[p];
                if(editMatrixUsed){
                   curEditCosts -= editCostSubtree[p];
                }
                else{
                    curEditCosts += subtreeSize * insertEditCost - numNeighborsAll[p] * (insertEditCost + removeEditCost); 
                }
            });
            //not necessary for children because they are moved and edits dont change
        }

        //subtree sortPath
        if (subtreeSortPaths) {
            node nonSubtreeNeighbor;
            for (node v : subtreeNeighbors) {
                nonSubtreeNeighbor = dynamicForest.moveUpSubtreeNeighbor(v, nodeToMove);
            }
        }
        //add parent candidates to queue
        //TODO get rid of forNodes?
        for (node v : subtreeNeighbors)  {
                int64_t testEditCost = ((numNeighborsAll[v] * removeEditCost) + ((numNeighborsAll[v] - subtreeSize) * insertEditCost));
                //fill queue with candidates with positive editCost or neighbors of nodetomove
                if ( marker[v] || 
                    (!editMatrixUsed && (testEditCost >= 0)) ||
                    (editMatrixUsed && editCostSubtree[v] >= 0)) {
                    parentCandidates.push_back(v);
                }
        };

        double curRatioSubtreeNeighbors = subtreeNeighbors.size() * 1.0 / G.numberOfNodes();
        double curRatioParentCandidates = parentCandidates.size() * 1.0 / G.numberOfNodes();
        countSubtreeMoves ++;
        ratioSubtreeNeighbors = ratioSubtreeNeighbors + (curRatioSubtreeNeighbors - ratioSubtreeNeighbors) / countSubtreeMoves;
        ratioParentCandidates = ratioParentCandidates + (curRatioParentCandidates - ratioParentCandidates) / countSubtreeMoves;

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
        rootData.initialize(generation);

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
                    if (nodeTouched[u])
                        continue; // if the node was touched in the previous level, it was in
                              // currentLevel and thus has already been processed
                    processNodeForSubtree(u, nodeToMove);
                }
                --level;
                currentLevel.clear();
                currentLevel.swap(nextLevel);
            }
        }

        //save sum positive editCosts for root
        rootData.sumPositiveEdits = rootData.childCloseness;
        rootData.sumPositiveEditCosts = rootData.childClosenessCost;
        if (!randomness) {
            if (rootData.sumPositiveEditCosts > rootData.scoreMaxCost) {
                rootData.bestParentBelow = none;
                rootData.scoreMax = rootData.sumPositiveEdits;
                rootData.scoreMaxCost = rootData.sumPositiveEditCosts;
            }
        } else {
            bool coin = false;
            double ownWeight = rootData.numIndifferentChildren * std::log(2);
            if (rootData.sumPositiveEditCosts > rootData.scoreMaxCost || !rootData.hasChoices()) {
                // INFO("root better");
                rootData.scoreMax = rootData.sumPositiveEdits;
                rootData.scoreMaxCost = rootData.sumPositiveEditCosts;
                rootData.logEqualBestChoices = ownWeight;
                coin = true;
            } else if (rootData.sumPositiveEditCosts == rootData.scoreMaxCost) {
                rootData.addLogChoices(ownWeight);
                coin = logRandomBool(ownWeight - rootData.logEqualBestChoices);
                // INFO("root equally good");
            }
            if (coin) {
                //correct scoremax incase root is equally good for weighted case
                rootData.scoreMax = rootData.sumPositiveEdits;
                rootData.bestParentBelow = none;
            }
            assert(rootData.hasChoices());
        }

        bestEdits += subtreeExtDegree - rootData.scoreMax;
        if(editMatrixUsed){
            bestEditCosts += editCostNeighbors - (rootData.scoreMaxCost);
        }
        else{
            bestEditCosts += (removeEditCost * subtreeExtDegree) - (rootData.scoreMaxCost);
        }
        // If sortPaths and randomness is on, only adopt children when the chosen parent is the
        // lower end of its path.
        const TraversalData &bestParentData =
            (rootData.bestParentBelow == none) ? rootData : traversalData[rootData.bestParentBelow];
        // For subtree move adopting one child is okay, because we have at least two children because of subtree
        if (bestParentData.numCloseChildren > 0
            || (randomness && bestParentData.numIndifferentChildren > 0)) {
            for (node u : touchedNodes) {
                if (u != nodeToMove && dynamicForest.parent(u) == rootData.bestParentBelow && nodeTouched[u]) {
                    if (traversalData[u].childClosenessCost > 0) {
                        bestChildren.push_back(u);
                    } else if (randomness && traversalData[u].childClosenessCost == 0) {
                        //direclty insert into bestChildren based on coin toss
                        if (randomBool(2)) {
                            //child closeness is not necessarily 0 if weighted childcloseness is 0
                            bestEdits -= traversalData[u].childCloseness;
                            bestChildren.push_back(u);
                        }              
                    }
                }
            }
            assert(randomness || bestChildren.size() == bestParentData.numCloseChildren);
        }

        // calculate the number of saved edits as comparing the absolute number of edits doesn't make
        // sense
        savedEdits = curEdits - bestEdits;
        savedEditCosts = curEditCosts - bestEditCosts;

        assert(!randomness || rootData.hasChoices());

        if (rootData.logEqualBestChoices < std::log(std::numeric_limits<long long>::max())) {
            rootEqualBestParentsCpy = std::llround(std::exp(rootData.logEqualBestChoices));
        } else {
            rootEqualBestParentsCpy = std::numeric_limits<count>::max();
        }

        TRACE("Best Parent for Subtree ", rootData.bestParentBelow);
        TRACE("Best adopted Children for Subtree ", bestChildren);
        if (savedEditCosts > 0 || randomness || subtreeOption == 1 || subtreeOption == 3) {
            assert(savedEditCosts >= 0);
            std::vector<node> nodeToMoveVec { nodeToMove };
            TRACE("SubtreeMove: current Parent: ", curParent, " new Parent: ", rootData.bestParentBelow);
            dynamicForest.moveToAnyPosition(rootData.bestParentBelow, nodeToMoveVec);
            dynamicForest.moveToAnyPosition(nodeToMove, bestChildren);
            hasMoved |= (savedEditCosts > 0 || (randomness && rootEqualBestParentsCpy > 1));
            numNodesMoved += (savedEditCosts > 0 || (randomness && rootEqualBestParentsCpy > 1));
            numSubtreeMoves += (savedEditCosts > 0 || (randomness && rootEqualBestParentsCpy > 1));
            numEdits -= savedEdits;
            costEdits -= savedEditCosts;
            INFO("Saved Edits: ", savedEdits, ", Saved EditCosts: ", savedEditCosts);
#ifndef NDEBUG
            assert(numEdits == countNumberOfEdits());
            assert(costEdits == countCostOfEdits());
#endif
        }

        for (node u : touchedNodes) {
            lastVisitedDFSNode[u] = u;
            nodeTouched[u] = false;
        }
        for (node u : subtreeNeighbors) {
            numNeighborsAll[u] = 0;
        }

        parentCandidates.clear();
        touchedNodes.clear();
        subtreeNeighbors.clear();
    }
    for (node v : subtreeNodes) {
        inSubtree[v] = false;
    }
    subtreeNodes.clear();
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
    traversalData[u].sumPositiveEdits = traversalData[u].childCloseness;
    traversalData[u].sumPositiveEditCosts = traversalData[u].childClosenessCost;
    assert(traversalData[u].childClosenessCost >= 0);

    if (marker[u]) { 
        ++traversalData[u].childCloseness;
    } else {
        --traversalData[u].childCloseness; 
    }
    if(editMatrixUsed){
        //increase childCloseness if neighbor, decrease if not neighbor
        traversalData[u].childClosenessCost +=  editCostNodeToMove[u];
    }
    else{
        traversalData[u].childClosenessCost += marker[u] * removeEditCost;
        traversalData[u].childClosenessCost -=
        1 * insertEditCost - marker[u] * insertEditCost;
    }


    TRACE("Edit difference before descending: ", traversalData[u].childCloseness, ", ", traversalData[u].childClosenessCost);

    assert(!marker[u] || traversalData[u].childClosenessCost > 0 || (editMatrixUsed && editCostNodeToMove[u]== 0));

    if (traversalData[u].childClosenessCost >= 0) {
        assert(lastVisitedDFSNode[u] == u);

        node c = dynamicForest.nextDFSNodeOnEnter(u, u);

        while (c != u) {
            //resturcture if 
            if (!nodeTouched[c] || traversalData[c].childClosenessCost < 0) {
                if (nodeTouched[c] ){
                    assert(traversalData[c].childClosenessCost < 0);
                    traversalData[u].childClosenessCost += traversalData[c].childClosenessCost;
                    traversalData[u].childCloseness+= traversalData[c].childCloseness;
                    // advance to the next starting point for the DFS search.
                    c = lastVisitedDFSNode[c];
                }
                else{
                    if(editMatrixUsed){
                        //decrease childCloseness because editCost is here negative
                        traversalData[u].childClosenessCost += editCostNodeToMove[c];
                    }
                    else{
                        traversalData[u].childClosenessCost -= insertEditCost;
                    }
                    --traversalData[u].childCloseness;
                }

                if (traversalData[u].childClosenessCost < 0 || dynamicForest.depth(c) > maxDepth) {
                    lastVisitedDFSNode[u] = c;
                    break;
                }

                c = dynamicForest.nextDFSNodeOnEnter(c, u);
            } else {
                if(traversalData[c].childClosenessCost == 0){
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

    TRACE("Edit difference after descending: ", traversalData[u].childCloseness, ", ", traversalData[u].childClosenessCost);

    if (!randomness) {
        if (traversalData[u].sumPositiveEditCosts > traversalData[u].scoreMaxCost || traversalData[u].scoreMaxCost == 0) {
            traversalData[u].scoreMax = traversalData[u].sumPositiveEdits;
            traversalData[u].scoreMaxCost = traversalData[u].sumPositiveEditCosts;
            traversalData[u].bestParentBelow = u;
        }
    } else {
        bool coin = false;
        double ownWeight = traversalData[u].numIndifferentChildren * std::log(2);
        if (traversalData[u].sumPositiveEditCosts > traversalData[u].scoreMaxCost || !traversalData[u].hasChoices()) {
            traversalData[u].scoreMax = traversalData[u].sumPositiveEdits;
            traversalData[u].scoreMaxCost = traversalData[u].sumPositiveEditCosts;
            traversalData[u].logEqualBestChoices = ownWeight;
            // Either we do not adopt children, or we are at the lower end of a path.
            // Otherwise, there must be a node below u that is at least as good.
            assert(!sortPaths || traversalData[u].sumPositiveEdits == 0 || dynamicForest.isLowerEnd(u));
            coin = true;
        } else if (traversalData[u].sumPositiveEditCosts == traversalData[u].scoreMaxCost) {
            ownWeight = traversalData[u].calculateOwnWeightForEqualChoices();
            if (ownWeight > -std::numeric_limits<double>::infinity()) {
                traversalData[u].addLogChoices(ownWeight);
                coin = logRandomBool(ownWeight - traversalData[u].logEqualBestChoices);
            }
            assert(traversalData[u].hasChoices());
        }
        if (coin) {
            //update scoremax if scoremax weighted is identical to sumPositiveEditCosts, scoremax can be not identical to traversalData[u].sumPositiveEdits
            traversalData[u].scoreMax = traversalData[u].sumPositiveEdits;
            traversalData[u].bestParentBelow = u;
        }
    }

    assert(traversalData[u].scoreMaxCost >= 0);

    if(marker[u]){
        ++traversalData[u].scoreMax;
    } else {
        --traversalData[u].scoreMax;
    }

    if(editMatrixUsed){
        traversalData[u].scoreMaxCost += editCostNodeToMove[u];     
    }
    else{
        traversalData[u].scoreMaxCost += marker[u] * removeEditCost;
        traversalData[u].scoreMaxCost -= 1 * insertEditCost - marker[u] * insertEditCost;
    }

    TRACE("Maximum gain at ", u, ": ", traversalData[u].scoreMax, ", ", traversalData[u].scoreMaxCost);

    node p = dynamicForest.parent(u);
    TraversalData &parentData = (p == none) ? rootData : traversalData[p];
    parentData.initialize(generation);

    if ((traversalData[u].scoreMaxCost > 0 || traversalData[u].childClosenessCost > 0) && p != none) {
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

    if (traversalData[u].scoreMaxCost > parentData.scoreMaxCost) {
        parentData.logEqualBestChoices = traversalData[u].logEqualBestChoices;
        parentData.scoreMax = traversalData[u].scoreMax;
        parentData.scoreMaxCost = traversalData[u].scoreMaxCost;
        parentData.bestParentBelow = traversalData[u].bestParentBelow;
    } else if (randomness && traversalData[u].scoreMaxCost == parentData.scoreMaxCost) {
        parentData.addLogChoices(traversalData[u].logEqualBestChoices);
        if (logRandomBool(traversalData[u].logEqualBestChoices - parentData.logEqualBestChoices)) {
            //update scoremax of parent if scoremaxCost is equal and current node offers more choices
            parentData.scoreMax = traversalData[u].scoreMax;
            parentData.bestParentBelow = traversalData[u].bestParentBelow;
        }
    }
    
    //propagate positive childCloseness
    if (traversalData[u].childClosenessCost >= 0) {
        assert(traversalData[u].childClosenessCost <= traversalData[u].scoreMaxCost);
        if (traversalData[u].childClosenessCost == 0) {
            ++parentData.numIndifferentChildren;
        } else {
            ++parentData.numCloseChildren;
            parentData.childCloseness += traversalData[u].childCloseness;
            parentData.childClosenessCost += traversalData[u].childClosenessCost;
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
    if (!nodeTouched[u]) {
        nodeTouched[u] = true;
        touchedNodes.emplace_back(u);
        assert(marker[u] || 
        (editMatrixUsed && editCostSubtree[u] >= 0) || 
        (!editMatrixUsed && (((numNeighborsAll[u] * removeEditCost) + ((numNeighborsAll[u] - subtreeSize) * insertEditCost)) >= 0))); 
        // only node with editCosts >= 0 are not touched and are processed
    }

    traversalData[u].initialize(generation);

    //required edits to this node for subtree
    traversalData[u].subtreeEdits = (2 * numNeighborsAll[u] - subtreeSize);
    if( editMatrixUsed){
        traversalData[u].subtreeEditCosts = editCostSubtree[u];
    } else {
        traversalData[u].subtreeEditCosts = (numNeighborsAll[u] * removeEditCost) + ((numNeighborsAll[u] - subtreeSize) * insertEditCost);
    }

    //sum of Positive childCloseness below node
    traversalData[u].sumPositiveEdits = traversalData[u].childCloseness;
    traversalData[u].sumPositiveEditCosts = traversalData[u].childClosenessCost;
    assert(traversalData[u].sumPositiveEditCosts >=0);

    if (marker[u]) { 
        ++traversalData[u].childCloseness;
    } else {
        --traversalData[u].childCloseness; 
    }
    if(editMatrixUsed){
        //increase childCloseness if neighbor, decrease if not neighbor
        traversalData[u].childClosenessCost +=  editCostNodeToMove[u];
    }
    else{
        traversalData[u].childClosenessCost += marker[u] * removeEditCost;
        traversalData[u].childClosenessCost -=
        1 * insertEditCost - marker[u] * insertEditCost;
    }


    TRACE("Edit difference before descending: ", traversalData[u].childCloseness, ", ", traversalData[u].childClosenessCost);

    assert(!marker[u] || traversalData[u].childClosenessCost > 0 || (editMatrixUsed && editCostNodeToMove[u]== 0));

    if (traversalData[u].childClosenessCost >= 0) {
        assert(lastVisitedDFSNode[u] == u);

        node c = dynamicForest.nextDFSNodeOnEnter(u, u);

        while (c != u) {
            //resturcture if 
            if ((!nodeTouched[c] || traversalData[c].childClosenessCost < 0) && c != nodeToMove) {
                if (nodeTouched[c] ){
                    assert(traversalData[c].childClosenessCost < 0);
                    traversalData[u].childClosenessCost += traversalData[c].childClosenessCost;
                    traversalData[u].childCloseness+= traversalData[c].childCloseness;
                    // advance to the next starting point for the DFS search.
                    c = lastVisitedDFSNode[c];
                }
                else{
                    if(editMatrixUsed){
                        //decrease childCloseness because editCost is here negative
                        traversalData[u].childClosenessCost += editCostNodeToMove[c];
                    }
                    else{
                        traversalData[u].childClosenessCost -= insertEditCost;
                    }
                    --traversalData[u].childCloseness;
                }

                if (traversalData[u].childClosenessCost < 0 || dynamicForest.depth(c) > maxDepth) {
                    lastVisitedDFSNode[u] = c;
                    break;
                }

                c = dynamicForest.nextDFSNodeOnEnter(c, u);
            } else {
                if(traversalData[c].childClosenessCost == 0){
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

    TRACE("Edit difference after descending: ", traversalData[u].childCloseness, ", ", traversalData[u].childClosenessCost);

    if (!randomness) {
        if (traversalData[u].sumPositiveEditCosts > traversalData[u].scoreMaxCost || traversalData[u].scoreMaxCost == 0) {
            traversalData[u].scoreMax = traversalData[u].sumPositiveEdits;
            traversalData[u].scoreMaxCost = traversalData[u].sumPositiveEditCosts;
            traversalData[u].bestParentBelow = u;
        }
    } else {
        bool coin = false;
        double ownWeight = traversalData[u].numIndifferentChildren * std::log(2);
        if (traversalData[u].sumPositiveEditCosts > traversalData[u].scoreMaxCost || !traversalData[u].hasChoices()) {
            traversalData[u].scoreMax = traversalData[u].sumPositiveEdits;
            traversalData[u].scoreMaxCost = traversalData[u].sumPositiveEditCosts;
            traversalData[u].logEqualBestChoices = ownWeight;
            // Either we do not adopt children, or we are at the lower end of a path.
            // Otherwise, there must be a node below u that is at least as good.
            //TODO check is that correct because we dont sort path before subtreeMove
            //assert(!sortPaths || sumPositiveEdits == 0 || dynamicForest.isLowerEnd(u));
            coin = true;
        } else if (traversalData[u].sumPositiveEditCosts == traversalData[u].scoreMaxCost) {
            //ownWeight = traversalData[u].calculateOwnWeightForEqualChoices();
            traversalData[u].addLogChoices(ownWeight);
            coin = logRandomBool(ownWeight - traversalData[u].logEqualBestChoices);
            assert(traversalData[u].hasChoices());
        }
        if (coin) {
            //update scoremax if scoremax weighted is identical to sumPositiveEditCosts, scoremax can be not identical to sumPositiveEdits
            traversalData[u].scoreMax = traversalData[u].sumPositiveEdits;
            traversalData[u].bestParentBelow = u;
        }
    }

    assert(traversalData[u].scoreMaxCost >= 0);

    traversalData[u].scoreMax += traversalData[u].subtreeEdits;
    traversalData[u].scoreMaxCost += traversalData[u].subtreeEditCosts;      

    TRACE("Maximum score at ", u, ": ", traversalData[u].scoreMax, ", ", traversalData[u].scoreMaxCost);
    node p = dynamicForest.parent(u);
    TraversalData &parentData = (p == none) ? rootData : traversalData[p];

    parentData.initialize(generation);

    int64_t parentEditCost = 0;
    if ((traversalData[u].scoreMaxCost > 0 || traversalData[u].subtreeEditCosts > 0 || traversalData[u].childClosenessCost > 0) && p != none) {
        if (useBucketQueue) {
            assert(dynamicForest.depth(p) <= maxDepth);
        }
        if (!nodeTouched[p]) {
            nodeTouched[p] = true;
            if(!editMatrixUsed){
                parentEditCost = ((numNeighborsAll[p] * removeEditCost) + ((numNeighborsAll[p] - subtreeSize) * insertEditCost));
            }
            touchedNodes.push_back(p);
            if (!useBucketQueue) {
                nextLevel.push_back(p);
            } else if (!(marker[p]) && 
                !(editMatrixUsed && editCostSubtree[p] >= 0) && 
                !(!editMatrixUsed && parentEditCost >= 0)) { 
                // positive editcost already in queue
                bucketQueue.insertParent(p);
            }
        }
    }

    if (traversalData[u].scoreMaxCost > parentData.scoreMaxCost) {
        parentData.logEqualBestChoices = traversalData[u].logEqualBestChoices;
        parentData.scoreMax = traversalData[u].scoreMax;
        parentData.scoreMaxCost = traversalData[u].scoreMaxCost;
        parentData.bestParentBelow = traversalData[u].bestParentBelow;
    } else if (randomness && traversalData[u].scoreMaxCost == parentData.scoreMaxCost) {
        parentData.addLogChoices(traversalData[u].logEqualBestChoices);
        if (logRandomBool(traversalData[u].logEqualBestChoices - parentData.logEqualBestChoices)) {
            //update scoremax of parent if scoremaxCost is equal and current node offers more choices
            parentData.scoreMax = traversalData[u].scoreMax;
            parentData.bestParentBelow = traversalData[u].bestParentBelow;
        }
    }

    if (traversalData[u].childClosenessCost >= 0) {
        //not true in this case: positive childCloseness when node neighbor of nodeToMove but  not too subtreenodes
        //assert(traversalData[u].childClosenessCost <= traversalData[u].scoreMaxCost);
        if (traversalData[u].childClosenessCost == 0) {
            ++parentData.numIndifferentChildren;
        } else {
            ++parentData.numCloseChildren;
            //dont propagate childcloseness, already happend in local move
            parentData.childCloseness += traversalData[u].childCloseness;
            parentData.childClosenessCost += traversalData[u].childClosenessCost;
        }
    }
}

void EditingRunner::compareWithQuadratic(node nodeToMove) const {
    std::vector<int64_t> missingBelow, missingAbove, existingBelow, existingAbove;
    std::vector<int64_t> missingBelowCosted, missingAboveCosted, existingBelowCosted, existingAboveCosted;
    missingBelow.resize(G.upperNodeIdBound(), 0);
    missingAbove.resize(G.upperNodeIdBound(), 0);
    existingBelow.resize(G.upperNodeIdBound(), 0);
    existingAbove.resize(G.upperNodeIdBound(), 0);
    missingBelowCosted.resize(G.upperNodeIdBound(), 0);
    missingAboveCosted.resize(G.upperNodeIdBound(), 0);
    existingBelowCosted.resize(G.upperNodeIdBound(), 0);
    existingAboveCosted.resize(G.upperNodeIdBound(), 0);
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
                            missingBelowCosted[u] = missingAboveCosted[u] = 1 * (-1) * editCostNodeToMove[u] - marker[u] * (-1) * editCostNodeToMove[u];
                            existingBelowCosted[u] = existingAboveCosted[u] = marker[u] * editCostNodeToMove[u];
                        }
                        else{
                            missingBelowCosted[u] = missingAboveCosted[u] = 1 * insertEditCost - marker[u] * insertEditCost;
                            existingBelowCosted[u] = existingAboveCosted[u] = marker[u] * removeEditCost;
                        }
                    }
                    node p = dynamicForest.parent(u);
                    if (p != none) {
                        missingAbove[u] += missingAbove[p];
                        existingAbove[u] += existingAbove[p];
                        missingAboveCosted[u] += missingAboveCosted[p];
                        existingAboveCosted[u] += existingAboveCosted[p];
                    }
                    return DynamicForest::ReturnState::CONTINUE;
                },
                [&](node u) {
                    node p = dynamicForest.parent(u);
                    if (p != none) {
                        missingBelow[p] += missingBelow[u];
                        existingBelow[p] += existingBelow[u];
                        missingBelowCosted[p] += missingBelowCosted[u];
                        existingBelowCosted[p] += existingBelowCosted[u];
                        if (usingDeepNeighbors[u])
                            usingDeepNeighbors[p] = true;
                    }
                    return DynamicForest::ReturnState::CONTINUE;
                });
        }
    });

    assert(missingBelow[nodeToMove] == 0);
    assert(existingBelow[nodeToMove] == 0);
    assert(missingBelowCosted[nodeToMove] == 0);
    assert(existingBelowCosted[nodeToMove] == 0);

    if (!sortPaths) {
        bool exactValue = true;
        for (node c : curChildren) {
            missingBelow[nodeToMove] += missingBelow[c];
            existingBelow[nodeToMove] += existingBelow[c];
            missingBelowCosted[nodeToMove] += missingBelowCosted[c];
            existingBelowCosted[nodeToMove] += existingBelowCosted[c];
            if (usingDeepNeighbors[c])
                exactValue = false;
        }

        if (curParent != none) {
            missingAbove[nodeToMove] = missingAbove[curParent];
            existingAbove[nodeToMove] = existingAbove[curParent];
            missingAboveCosted[nodeToMove] = missingAboveCosted[curParent];
            existingAboveCosted[nodeToMove] = existingAboveCosted[curParent];
            if (usingDeepNeighbors[curParent])
                exactValue = false;
        }
        if (exactValue && !editMatrixUsed) {
            assert(curEdits
                   == numNeighbors - existingAbove[nodeToMove] - existingBelow[nodeToMove]
                          + missingAbove[nodeToMove] + missingBelow[nodeToMove]);
            assert(curEditCosts
                   == numNeighbors * removeEditCost - existingAboveCosted[nodeToMove] - existingBelowCosted[nodeToMove]
                          + missingAboveCosted[nodeToMove] + missingBelowCosted[nodeToMove]);
        }
    }

    count minEdits = std::numeric_limits<count>::max();
    count minEditCosts = std::numeric_limits<count>::max();
    count editCostsOffset = 0;
    std::vector<node> minChildren;
    node minParent = curParent;
    G.forNodes([&](node u) {
        if (u == nodeToMove || usingDeepNeighbors[u] || !existing[u])
            return;
        if (existingBelowCosted[u] >= missingBelowCosted[u]
            || (traversalData[u].generation == generation && traversalData[u].childClosenessCost > 0)) {
            assert(traversalData[u].childClosenessCost == existingBelowCosted[u] - missingBelowCosted[u]);
            //
            assert(traversalData[u].childCloseness == existingBelow[u] - missingBelow[u]);
        } else if (nodeTouched[u]) {
            assert(traversalData[u].childClosenessCost < 0);
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
            editCostsOffset += editCostNodeToMove[u];
        }
    });

    auto tryEditBelow = [&](node p) {
        if (p == nodeToMove)
            return DynamicForest::ReturnState::CONTINUE;;
        count editCosts = editCostsOffset;
        count edits = numNeighbors;
        if(!editMatrixUsed)
        {       
            editCosts = numNeighbors * removeEditCost;
        }
        if (p != none) {
            edits += missingAbove[p];
            edits -= existingAbove[p];
            editCosts += missingAboveCosted[p];
            editCosts -= existingAboveCosted[p];
        }

        std::vector<node> children;
        dynamicForest.forChildrenOf(p, [&](node c) {
            if (c == nodeToMove || usingDeepNeighbors[c] || !existing[c])
                return;
            if (existingBelowCosted[c] > missingBelowCosted[c]) { // TODO try >= (more children...)
                if (dynamicForest.children(c).empty() && marker[c]) {
                    assert(traversalData[c].childCloseness == 1);
                    if(editMatrixUsed){
                        assert(traversalData[c].childClosenessCost == 1 * editCostNodeToMove[c]);
                    }
                    else{
                        assert(traversalData[c].childClosenessCost == 1 * removeEditCost);
                    }

                }
                assert(traversalData[c].childCloseness == existingBelow[c] - missingBelow[c]);//
                assert(traversalData[c].childClosenessCost == existingBelowCosted[c] - missingBelowCosted[c]);

                children.emplace_back(c);
                edits -= existingBelow[c] - missingBelow[c];
                editCosts -= existingBelowCosted[c] - missingBelowCosted[c];
            }
        });

        if (editCosts < minEditCosts) {
            minEdits = edits;
            minEditCosts = editCosts;
            minChildren = std::move(children);
            minParent = p;
        }
        return DynamicForest::ReturnState::CONTINUE;
    };

    dynamicForest.dfsFrom(
        none, [](node) {return DynamicForest::ReturnState::CONTINUE;}, tryEditBelow);
    tryEditBelow(none);
    //correct assertion? assert(minEdits >= bestEdits);
    //assert(minEdits <= bestEdits);
    assert(minEditCosts == bestEditCosts);

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
    //assertion only correct for weighted childCLoseness
    //assert(minEdits <= childClosenessControl);

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

                G.forEdgesOf(u, [&](node v) {
                    if (marker[v])
                        ++upperNeighbors;
                });

                numExistingEdges += upperNeighbors;
                numMissingEdges += depth - upperNeighbors;
                marker[u] = true;
                depth += 1;
                return DynamicForest::ReturnState::CONTINUE;
            },
            [&](node u) { // on exit
                marker[u] = false;
                depth -= 1;
                return DynamicForest::ReturnState::CONTINUE;
            });
    });

    return numMissingEdges + (G.numberOfEdges() - numExistingEdges);
}

count EditingRunner::countCostOfEdits() const {
    // count weight of edits that are needed with the initial given forest
    count costOfEdits = 0;
    count costOfRemoveEdits = 0;
    count costOfInsertEdits = 0;
    count numExistingEdges = 0;
    count numMissingEdges = 0;
    std::vector<bool> marker(G.upperNodeIdBound());

    dynamicForest.forChildrenOf(none, [&](node r) {
        count depth = 0;
        dynamicForest.dfsFrom(
            r,
            [&](node u) { // on enter
                count upperNeighbors = 0;

                G.forEdgesOf(u, [&](node v) {
                    if (marker[v]){
                        ++upperNeighbors;
                        //-remove
                        if(editMatrixUsed){
                            costOfRemoveEdits -= editCostMatrix[u][v];
                        }
                    }
                    else if(editMatrixUsed) {
                        //+remove
                        costOfRemoveEdits += editCostMatrix[u][v];
                    }
                });
                numExistingEdges += upperNeighbors;
                numMissingEdges += depth - upperNeighbors;
                marker[u] = true;  
                depth += 1;
                return DynamicForest::ReturnState::CONTINUE;
            },
            [&](node u) { // on exit
                marker[u] = false;
                depth -= 1;
                return DynamicForest::ReturnState::CONTINUE;
            });
    });
    if(editMatrixUsed){
        G.forNodes([&](node u) {
            dynamicForest.forAncestors(u, [&](node p) {
                if(editCostMatrix[u][p] < 0 && u != p){
                    costOfInsertEdits += (-1) * editCostMatrix[u][p];
                }
            });
        });
    }
    TRACE("Missing Edges:  ", numMissingEdges,", NumEdges-ExistingEdges: " ,(G.numberOfEdges() - numExistingEdges));
    if(editMatrixUsed){
        costOfEdits = (costOfInsertEdits) + (costOfRemoveEdits / 2);
    }
    else{
        costOfEdits = ((numMissingEdges * insertEditCost) + ((G.numberOfEdges() - numExistingEdges) * removeEditCost));
    }
    return costOfEdits;
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
                return DynamicForest::ReturnState::CONTINUE;
            },
            [&](node) {return DynamicForest::ReturnState::CONTINUE;});
    });

    G.forEdgesOf(u, [&](node w) {
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

std::vector<node> EditingRunner::getParentsForTree() {
    std::vector<node> parents = std::vector<node>(G.upperNodeIdBound(), none);
    std::vector<bool> nodeHandled = std::vector<bool>(G.upperNodeIdBound(), 0);
    G.forNodesInRandomOrder([&](node u) { order.push_back(u); });
        for (index j = 0; j < G.numberOfNodes(); j++) {
                node u = order[j];
                G.forEdgesOf(u, [&](node v) {
                    if(nodeHandled[v]){
                        parents[u] = v;
                    }
                });        
                nodeHandled[u] = 1;
        }
    return parents;
}
} // namespace QuasiThresholdMoving

} // namespace NetworKit
