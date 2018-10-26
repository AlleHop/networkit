#ifndef CKBDYNAMIC_COMMUNITY_H_
#define CKBDYNAMIC_COMMUNITY_H_

#include <tlx/counting_ptr.hpp>
#include <tlx/math/integer_log2.hpp>
#include "../../Globals.h"
#include "../../auxiliary/SamplingSet.h"
#include "NodePairHash.h"

namespace NetworKit {
	namespace CKBDynamicImpl {
		class CKBDynamicImpl;

		class Community : public tlx::ReferenceCounter {
		public:
			/**
			 * Create an exact copy of a given community @a o
			 *
			 * @param o The community to copy.
			 */
			Community(const Community& o);

			Community(double edgeProbability, CKBDynamicImpl& generator);

			/**
			 * Remove the given node @a u from the community.
			 *
			 * This removes all edges incident to @a u from the community.
			 */
			void removeNode(node u);

			/**
			 * Remove a random node from the community.
			 */
			void removeRandomNode();

			/**
			 * Adds a node to the community and generates
			 * edges from @a u to the community with the
			 * current edge probability of the community.
			 *
			 * @param u The node to add.
			 */
			void addNode(node u);

			double getEdgeProbability() const { return edgeProbability; }

			/**
			 * Change the probability of edges to @a prob.
			 *
			 * This adds or removes random edges to achieve the desired edge probability.
			 *
			 * @param u The node to remove.
			 */
			void changeEdgeProbability(double prob);

			/**
			 * Combine the edges of this community with the edges of another community.
			 * This assumes that the nodes sets are identical.
			 * Overlapping edges are removed and the other community will be empty afterwards.
			 *
			 * @param other The community to combine the nodes with.
			 */
			void combineWith(Community& other);

			/**
			 * Perturb edges with a certain probability. 0
			 * changes no edges, 1 removes all edges and
			 * completely re-generates the community.
			 *
			 * FIXME possibly we want to call this with an number of edges to perturb.
			 *
			 * @param p The probability/percentage with which edges shall be removed or added.
			 */
			void perturbEdges(double prob);


			count getNumberOfNodes() const { return nodes.size(); };

			count getNumberOfEdges() const { return edges.size(); };

			count getMaximumNumberOfEdges() const { return (nodes.size() * (nodes.size() - 1) / 2); }

			const Aux::SamplingSet<node>& getNodes() const { return nodes; };

			bool hasNode(node u) const { return nodes.contains(u); };

			index getId() const { return id; };

			bool isAvailable() const { return available; }

			void setAvailable(bool available);
		private:
			void removeEdge(node u, node v);
			void addEdge(node u, node v);

			void removeRandomEdges(count k);
			void addRandomEdges(count k);

			count drawDesiredNumberOfEdges(double prob) const;

			void verifyInvariants() const;

			std::pair<node, node> edgeFromIndex(index i) const;

			index id;
			Aux::SamplingSet<std::pair<node, node>, NodePairHash> edges;
			// only used if edgeProbability > 0.5.
			Aux::SamplingSet<std::pair<node, node>, NodePairHash> nonEdges;
			Aux::SamplingSet<node> nodes;
			std::unordered_map<node, Aux::SamplingSet<node>> neighbors;
			double edgeProbability;
			bool available;
			bool storeNonEdges;
			CKBDynamicImpl& generator;
		};

		using CommunityPtr = tlx::CountingPtr<Community>;
	}
}

namespace std {
	// inspired by https://stackoverflow.com/questions/20953390/what-is-the-fastest-hash-function-for-pointers
	template <>
	struct hash<NetworKit::CKBDynamicImpl::CommunityPtr> {
		std::size_t operator()(const NetworKit::CKBDynamicImpl::CommunityPtr& k) const {
			using std::size_t;
			static const size_t shift = tlx::integer_log2_floor(1 + sizeof(NetworKit::CKBDynamicImpl::Community));
			return reinterpret_cast<size_t>(k.get()) >> shift;
		};
	};
}

#endif