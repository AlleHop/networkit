/*
 * CentralityGTest.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include "DynTopClosenessGTest.h"
#include "../../generators/DorogovtsevMendesGenerator.h"
#include "../DynTopHarmonicCloseness.h"
#include "../HarmonicCloseness.h"
#include "../TopHarmonicCloseness.h"

namespace NetworKit {

TEST_F(DynTopClosenessGTest, testDynTopHarmonicCloseness) {

  Aux::Random::setSeed(42, false);
  Graph G = DorogovtsevMendesGenerator(500).generate();

  count k = 10;

  DynTopHarmonicCloseness centrality(G, k, false);
  centrality.run();

  HarmonicCloseness reference(G, false);
  reference.run();

  auto scores = centrality.ranking();
  auto refScores = reference.ranking();

  for (count j = 0; j < k; ++j) {
    EXPECT_TRUE(scores[j].first == refScores[j].first);
    EXPECT_FLOAT_EQ(scores[j].second, refScores[j].second);
  }

  count numInsertions = 1;

  std::vector<GraphEvent> deletions;
  std::vector<GraphEvent> insertions;

  for (count i = 0; i < numInsertions; i++) {

    node u = G.upperNodeIdBound();
    node v = G.upperNodeIdBound();

    do {
      u = G.randomNode();
      v = G.randomNode();

    } while (G.hasEdge(u, v));

    GraphEvent edgeAddition(GraphEvent::EDGE_ADDITION, u, v);
    insertions.insert(insertions.begin(), edgeAddition);

    GraphEvent edgeDeletion(GraphEvent::EDGE_REMOVAL, u, v);
    deletions.push_back(edgeDeletion);

    G.addEdge(u, v);
  }

  for (auto e : insertions) {
    G.removeEdge(e.u, e.v);
  }

  for (GraphEvent edgeAddition : insertions) {

    node u = edgeAddition.u;
    node v = edgeAddition.v;

    G.addEdge(u, v);

    centrality.update(edgeAddition);
    reference.run();

    scores = centrality.ranking();
    refScores = reference.ranking();

    for (count j = 0; j < k; ++j) {
      EXPECT_TRUE(scores[j].first == refScores[j].first);
      EXPECT_FLOAT_EQ(scores[j].second, refScores[j].second);
    }
  }

  for (GraphEvent edgeDeletion : deletions) {

    node u = edgeDeletion.u;
    node v = edgeDeletion.v;

    G.removeEdge(u, v);

    centrality.update(edgeDeletion);
    reference.run();

    scores = centrality.ranking();
    refScores = reference.ranking();

    for (count j = 0; j < k; ++j) {
      EXPECT_TRUE(scores[j].first == refScores[j].first);
      EXPECT_FLOAT_EQ(scores[j].second, refScores[j].second);
    }
  }
}
} // namespace NetworKit
