#include "CKBDynamicImpl.h"
#include "CommunityDeathEvent.h"
#include "CommunityBirthEvent.h"
#include "CommunitySplitEvent.h"
#include "CommunityMergeEvent.h"
#include "PowerlawCommunitySizeDistribution.h"
#include "../../auxiliary/SignalHandling.h"
#include <tlx/unused.hpp>

namespace NetworKit {
	namespace CKBDynamicImpl {
		namespace {
			/**
			 * Returns number of steps you need to wait until the next success (edge) occurs.
			 */
			count get_next_edge_distance(const double log_cp) {
				return static_cast<count>(1 + floor(log(1.0 - Aux::Random::probability()) / log_cp));
			}

		}

		void CKBDynamicImpl::addEdge(node u, node v) {
			auto e = Community::canonicalEdge(u, v);
			auto it = edgesAlive.find(e);

			if (it == edgesAlive.end()) {
				edgesAlive.insert(std::make_pair(e, 1ul));
				auto evIt = currentEdgeEvents.find(e);
				if (evIt != currentEdgeEvents.end()) {
					assert(evIt->second.type == GraphEvent::EDGE_REMOVAL);
					currentEdgeEvents.erase(evIt);
				} else {
					currentEdgeEvents.emplace(e, GraphEvent(GraphEvent::EDGE_ADDITION, e.first, e.second));
				}
			} else {
				it->second += 1;
			}
		}

		void CKBDynamicImpl::removeEdge(node u, node v) {
			auto e = Community::canonicalEdge(u, v);
			auto it = edgesAlive.find(e);
			if (it == edgesAlive.end()) {
				throw std::runtime_error("Error, removing edge that does not exist");
			}

			if (it->second > 1) {
				--it->second;
			} else {
				edgesAlive.erase(it);
				auto evIt = currentEdgeEvents.find(e);
				if (evIt != currentEdgeEvents.end()) {
					assert(evIt->second.type == GraphEvent::EDGE_ADDITION);
					currentEdgeEvents.erase(evIt);
				} else {
					currentEdgeEvents.emplace(e, GraphEvent(GraphEvent::EDGE_REMOVAL, e.first, e.second));
				}
			}
		}

		void CKBDynamicImpl::addNodeToCommunity(node u, CommunityPtr com) {
			if (com != globalCommunity) {
				nodeCommunities[u].insert(com);
				auto comNode = std::make_pair(com->getId(), u);
				auto it = currentCommunityEvents.find(comNode);
				if (it != currentCommunityEvents.end()) {
					assert(it->second.type == CommunityEvent::NODE_LEAVES_COMMUNITY);
					currentCommunityEvents.erase(it);
				} else {
					currentCommunityEvents.emplace(comNode, CommunityEvent(CommunityEvent::NODE_JOINS_COMMUNITY, u, com->getId()));
				}

				++currentCommunityMemberships;
			}
		}

		void CKBDynamicImpl::removeNodeFromCommunity(node u, CommunityPtr com) {
			if (com != globalCommunity) {
				nodeCommunities[u].erase(com);
				auto comNode = std::make_pair(com->getId(), u);
				auto it = currentCommunityEvents.find(comNode);
				if (it != currentCommunityEvents.end()) {
					assert(it->second.type == CommunityEvent::NODE_JOINS_COMMUNITY);
					currentCommunityEvents.erase(it);
				} else {
					currentCommunityEvents.emplace(comNode, CommunityEvent(CommunityEvent::NODE_LEAVES_COMMUNITY, u, com->getId()));
				}
				--currentCommunityMemberships;
			}
		}

		void CKBDynamicImpl::addCommunity(CommunityPtr com) {
			if (com->isAvailable()) {
				// If community is too small, remove community again!!
				if (com->getNumberOfNodes() < communitySizeSampler->getMinSize()) {
					INFO("community has only ", com->getNumberOfNodes(), " nodes, destroying.");
					currentEvents.emplace_back(new CommunityDeathEvent(com, 0, 1, *this));
				} else {
					availableCommunities.insert(com);
				}
			} else {
				availableCommunities.erase(com);
			}
			communities.insert(com);
		}

		void CKBDynamicImpl::removeCommunity(CommunityPtr com) {
			availableCommunities.erase(com);
			communities.erase(com);
		}

		index CKBDynamicImpl::nextCommunityId() {
			index result = maxCommunityId;
			++maxCommunityId;
			return result;
		}

		CKBDynamicImpl::CKBDynamicImpl(const CKBDynamic::param_type &params) :
			communitySizeSampler(new PowerlawCommunitySizeDistribution(params.minCommunitySize, params.maxCommunitySize, params.communitySizeExponent, params.intraCommunityEdgeProbability, params.intraCommunityEdgeExponent, params.minSplitRatio)),
			membershipDistribution(params.minCommunityMembership, params.maxCommunityMembership, params.communityMembershipExponent),
			maxCommunityId(0),
			sumOfDesiredMemberships(0),
			n(params.n),
			communityEventProbability(params.communityEventProbability),
			nodeEventProbability(params.nodeEventProbability),
			perturbationProbability(params.perturbationProbability),
			epsilon(params.epsilon),
			numTimesteps(params.numTimesteps),
			currentCommunityMemberships(0) {
			membershipDistribution.run();
		}

		std::vector<GraphEvent> CKBDynamicImpl::getGraphEvents() const {
			this->assureFinished();
			return std::move(graphEvents);
		}

		std::vector<CommunityEvent> CKBDynamicImpl::getCommunityEvents() const {
			this->assureFinished();
			return std::move(communityEvents);
		}

		void CKBDynamicImpl::generateNode() {
			node u = desiredMemberships.size();
			desiredMemberships.push_back(membershipDistribution.getDegree());
			sumOfDesiredMemberships += desiredMemberships.back();
			nodesAlive.insert(u);
			nodeCommunities.emplace_back();
			globalCommunity->addNode(u);
			graphEvents.emplace_back(GraphEvent::NODE_ADDITION, u);
		}

		void CKBDynamicImpl::eraseNode() {
			node u = nodesAlive.at(Aux::Random::index(nodesAlive.size()));
			sumOfDesiredMemberships -= desiredMemberships[u];
			desiredMemberships[u] = 0;

			while (nodeCommunities[u].size() > 0) {
				CommunityPtr com = *nodeCommunities[u].begin();
				com->removeNode(u);
				// if a community becomes too small, erase it
				if (com->isAvailable() && com->getNumberOfNodes() < communitySizeSampler->getMinSize()) {
					INFO("Available community has only ", com->getNumberOfNodes(), " nodes, destroying.");
					currentEvents.emplace_back(new CommunityDeathEvent(com, 0, 1, *this));
				}
			}

			assert(nodesAlive.contains(u));
			nodesAlive.erase(u);
			globalCommunity->removeNode(u);
			currentErasedNodes.push_back(u);
		}

		void CKBDynamicImpl::finishTimeStep() {
			for (const auto &it : currentEdgeEvents) {
				graphEvents.push_back(it.second);
			}
			currentEdgeEvents.clear();

			for (const auto &it : currentCommunityEvents) {
				communityEvents.push_back(it.second);
			}
			currentCommunityEvents.clear();

			// Ensure that node removals are after edge events.
			for (node u : currentErasedNodes) {
				graphEvents.emplace_back(GraphEvent::NODE_REMOVAL, u);
			}
			currentErasedNodes.clear();

			communityEvents.emplace_back(CommunityEvent::TIME_STEP);
			graphEvents.emplace_back(GraphEvent::TIME_STEP);
		}

		count CKBDynamicImpl::sampleNumSteps() const {
			return Aux::Random::integer(5, 15);
		}


		void CKBDynamicImpl::run() {
			if (hasRun) throw std::runtime_error("Error, run has already been called");

			Aux::SignalHandler handler;

			// initialization
			globalCommunity = CommunityPtr(new Community(epsilon, *this));
			communities.erase(globalCommunity);
			assert(!globalCommunity->isAvailable());

			for (node u = 0; u < n; ++u) {
				generateNode();
			}

			const count initialNumberOfNodes = nodesAlive.size();

			while (currentCommunityMemberships < sumOfDesiredMemberships) {
				handler.assureRunning();
				count communitySize;
				double edgeProbability;
				std::tie(communitySize, edgeProbability) = communitySizeSampler->drawCommunity();

				CommunityPtr com(new Community(edgeProbability, *this));
				com->setDesiredNumberOfNodes(communitySize);
				com->setAvailable(true);
			}

			assignNodesToCommunities();

			// Finish initial graph generation.
			finishTimeStep();

			std::binomial_distribution<count> numEventDistribution;
			double deathProbability = 0.25, birthProbability = 0.25, splitProbability = 0.25, mergeProbability = 0.25;
			tlx::unused(mergeProbability);

			for (count timestep = 0; timestep < numTimesteps; ++timestep) {
				handler.assureRunning();
				numEventDistribution.param(std::binomial_distribution<count>::param_type(communities.size(), communityEventProbability));
				const count numCommunityEvents = numEventDistribution(Aux::Random::getURNG());

				numEventDistribution.param(std::binomial_distribution<count>::param_type(nodesAlive.size(), nodeEventProbability));
				const count numNodeEvents = numEventDistribution(Aux::Random::getURNG());

				INFO("Timestep ", timestep, " generating ", numCommunityEvents, " community events and ", numNodeEvents, " node events");

				for (count i = 0; i < numCommunityEvents; ++i) {
					handler.assureRunning();
					count numSteps = sampleNumSteps();
					double r = Aux::Random::real();
					if (r < birthProbability) {
						// generate new community
						currentEvents.emplace_back(new CommunityBirthEvent(numSteps, *this));
					} else if (r < birthProbability + deathProbability) {
						// let a community die
						if (availableCommunities.size() > 0) {
							CommunityPtr com = availableCommunities.at(Aux::Random::index(availableCommunities.size()));
							count coreSize = std::max<count>(0.1 * com->getNumberOfNodes(), communitySizeSampler->getMinSize());
							currentEvents.emplace_back(new CommunityDeathEvent(com, coreSize, numSteps, *this));
							assert(!com->isAvailable());
						} else {
							WARN("No community available for death event.");
						}
					} else if (r < birthProbability + deathProbability + splitProbability) {
						// Split a community
						if (availableCommunities.size() > 0) {
							CommunityPtr com = availableCommunities.at(Aux::Random::index(availableCommunities.size()));
							auto comSizeProbA = communitySizeSampler->drawCommunity();
							auto comSizeProbB = communitySizeSampler->drawCommunity();
							currentEvents.emplace_back(new CommunitySplitEvent(com, comSizeProbA.first, comSizeProbA.second, comSizeProbB.first, comSizeProbB.second, numSteps, *this));
							assert(!com->isAvailable());
						} else {
							WARN("No community available for splitting.");
						}
					} else {
						// merge two communities
						if (availableCommunities.size() > 1) {
							index ia = Aux::Random::integer(availableCommunities.size() - 1);
							index ib = Aux::Random::integer(1, availableCommunities.size() - 1);
							if (ia == ib) {
								ib = 0;
							}

							CommunityPtr comA = availableCommunities.at(ia);
							CommunityPtr comB = availableCommunities.at(ib);

							count targetSize;
							double targetEdgeProbability;
							std::tie(targetSize, targetEdgeProbability) = communitySizeSampler->drawCommunity();
							currentEvents.emplace_back(new CommunityMergeEvent(comA, comB, targetSize, targetEdgeProbability, numSteps, *this));
							assert(!comA->isAvailable());
							assert(!comB->isAvailable());
						} else {
							WARN("No two communities available for merge.");
						}
					}
				} // generated all new community events

				// generate node events
				const double wantedNodeFraction = initialNumberOfNodes * 1.0 / nodesAlive.size();
				const double nodeBirthProbability = wantedNodeFraction / (1 + wantedNodeFraction);

				// First generate all death events, then all birth events.
				// This ensures that no node that is born in this time step dies again in this time step.
				numEventDistribution.param(std::binomial_distribution<count>::param_type(numNodeEvents, nodeBirthProbability));
				const count nodesBorn = numEventDistribution(Aux::Random::getURNG());

				for (count j = 0; j < (numNodeEvents - nodesBorn); ++j) {
					eraseNode();
				}
				for (count j = 0; j < nodesBorn; ++j) {
					generateNode();
				}

				// Trigger all current events
				for (size_t e = 0; e < currentEvents.size();) {
					handler.assureRunning();
					currentEvents[e]->nextStep();

					if (!currentEvents[e]->isActive()) {
						std::swap(currentEvents[e], currentEvents.back());
						currentEvents.pop_back();
					} else {
						++e;
					}
				}

				// edge perturbations
				if (perturbationProbability > 0) {
					globalCommunity->perturbEdges(perturbationProbability);

					const double sqrtPerturbationProbability = std::sqrt(perturbationProbability);

					const double log_perturb = std::log(1.0 - sqrtPerturbationProbability);

					for (count ci = get_next_edge_distance(log_perturb) - 1; ci < communities.size(); ci += get_next_edge_distance(log_perturb)) {
						handler.assureRunning();
						communities.at(ci)->perturbEdges(sqrtPerturbationProbability);
					}
				}

				assignNodesToCommunities();

				// adjust event probabilities
				{
					const double x = sumOfDesiredMemberships * 1.0 / currentCommunityMemberships;
					birthProbability = 0.5 * x / (1 + x);
					deathProbability = 0.5 - birthProbability;
					INFO("Current memberships: ", currentCommunityMemberships, " desired: ", sumOfDesiredMemberships, " number of communities: ", communities.size(), " available: ", availableCommunities.size(), " active events ", currentEvents.size());
					INFO("Current nodes ", nodesAlive.size(), " current edges: ", edgesAlive.size(), " total graph events ", graphEvents.size(), " total community events ", communityEvents.size());
				}

				finishTimeStep();

			}

			hasRun = true;
		}

		void CKBDynamicImpl::assignNodesToCommunities() {
			count totalMissingMembers = 0;

			// TODO: Do something about imbalances between totalMissingMemberships and totalMissingMembers.
			// If totalMissingMemberships > totalMissingMembers, a) find nodes that have too many memberships and remove those nodes from some of them and b) create empty communities of size 1 that are removed afterwards.
			// If totalMissingMembers > totalMissingMemberships, find nodes to which we can assign additional memberships, i.e., nodes that are not much over their allocated memberships.
			std::vector<std::pair<CommunityPtr, count>> communitiesByDesiredMembers;

			{
				std::vector<count> numCommunitiesWithDesired;

				auto getMissing = [](const CommunityPtr &com) {
					const count desired = com->getDesiredNumberOfNodes();
					const count actual = com->getNumberOfNodes();
					assert(actual <= desired);
					return desired - actual;
				};

				for (CommunityPtr com : communities) {
					if (getMissing(com) == 0) continue;
					const count desired = com->getDesiredNumberOfNodes();

					if (desired >= numCommunitiesWithDesired.size()) {
						numCommunitiesWithDesired.resize(desired + 1);
					}

					++numCommunitiesWithDesired[desired];
				}

				if (numCommunitiesWithDesired.size() == 0) {
					return;
				}

				count sum = 0;
				for (auto it = numCommunitiesWithDesired.rbegin(); it != numCommunitiesWithDesired.rend(); ++it) {
					const count tmp = *it;
					*it = sum;
					sum += tmp;
				}

				communitiesByDesiredMembers.resize(sum);

				for (CommunityPtr com : communities) {
					const count missing = getMissing(com);
					const count desired = com->getDesiredNumberOfNodes();
					if (missing == 0) continue;

					totalMissingMembers += missing;
					communitiesByDesiredMembers[numCommunitiesWithDesired[desired]] = {com, missing};
					++numCommunitiesWithDesired[desired];
				}
			}

			count totalMissingMemberships = 0;

			// This map records for every node for which the desired number of memberships can impossibly be satisfied how many can be satisfied.
			std::unordered_map<node, count> satisfiableMemberships;
			for (node u : nodesAlive) {
				const count desired = desiredMemberships[u];
				const count actual = nodeCommunities[u].size();

				if (desired > actual) {
					count missing = desired - actual;
					count satisfiable = 0;
					for (auto it = communitiesByDesiredMembers.begin(); satisfiable < missing && it != communitiesByDesiredMembers.end(); ++it) {
						if (!it->first->hasNode(u)) {
							++satisfiable;
						}
					}

					if (satisfiable < missing) {
						satisfiableMemberships[u] = satisfiable;
					}

					totalMissingMemberships += satisfiable;
				}
			}

			if (totalMissingMemberships < totalMissingMembers) {
				// TODO: find additional nodes to assign
			}

			std::vector<node> nodesByDesiredMemberships;

			{
				std::vector<count> nodesPerDesired;

				for (node u : nodesAlive) {
					const count desired = desiredMemberships[u];
					const count actual = nodeCommunities[u].size();

					if (desired > actual) {
						if (nodesPerDesired.size() < desired) {
							nodesPerDesired.resize(desired + 1);
						}

						++nodesPerDesired[desired];
					}
				}
				count sum = 0;
				// Reverse prefix sum so the actual order is reversed
				for (auto it = nodesPerDesired.rbegin(); it != nodesPerDesired.rend(); ++it) {
					const count temp = *it;
					*it = sum;
					sum += temp;
				}

				nodesByDesiredMemberships.resize(sum);

				for (node u : nodesAlive) {
					const count desired = desiredMemberships[u];
					const count actual = nodeCommunities[u].size();

					if (desired > actual) {
						nodesByDesiredMemberships[nodesPerDesired[desired]] = u;
						++nodesPerDesired[desired];
					}
				}
			}


			std::vector<std::vector<CommunityPtr>> freshAssignments;

			for (node lu = 0; lu < nodesByDesiredMemberships.size(); ++lu) {
				const node u = nodesByDesiredMemberships[lu];
				count communitiesToFind = desiredMemberships[u] - nodeCommunities[u].size();
				for (auto cit = communitiesByDesiredMembers.begin(); cit != communitiesByDesiredMembers.end() && communitiesToFind > 0; ++cit) {
					// TODO: actually remove communities with missing = 0 to speed this up.
					const CommunityPtr &com = cit->first;
					count &missing = cit->second;
					if (!com->hasNode(u) && missing > 0) {
						--missing;
						freshAssignments[lu].push_back(com);
						--communitiesToFind;
						if (communitiesToFind == 0) {
							break;
						}
					}
				}
			}

			// TODO: collect all nodes that have not found enough communities and add them to remaining communities, creating duplicates.

			// TODO: Randomize freshAssignments using curveball or edge switches and thereby try to eliminate duplicates.

			// TODO: Actually assign nodes to communities, thereby eliminating any remaining duplicates.
		}
	}
}
