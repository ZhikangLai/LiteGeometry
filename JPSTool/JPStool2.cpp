#include "JPStool2.h"
//*********************************************************************************************//
/**
 * @brief Performs Jump Point Search (JPS) to find the shortest path in a grid from a start point to an end point.
 *
 * This function implements the JPS algorithm to efficiently search for a path in a grid. It starts by initializing
 * the search with the provided start and end points, then iteratively explores the grid nodes while applying the
 * JPS pruning technique to skip unnecessary nodes. The algorithm uses an open list, closed list, and a node pool
 * to track the nodes being processed. When the end point is reached, the path is reconstructed from the
 * `parent` pointers of the nodes.
 *
 * @param[in] _startPoint The starting point of the search, provided as a 3D vector.
 * @param[in] _endPoint The target end point of the search, provided as a 3D vector.
 * @return A vector of 3D grid points representing the best path from the start point to the end point.
 *         If no path is found, an empty vector is returned.
 */
std::vector<Eigen::RowVector3i> JPSPathFinder2::JPSGraphSearch(
    const Eigen::RowVector3d& _startPoint,
    const Eigen::RowVector3d& _endPoint
) {
    endPoint = _endPoint.array().round().cast<int>();
    openList.clear();
    closedList.clear();
    openSet.clear();
    nodePool.clear();

    GridNode startNode;
    startNode.dirs = dirDecsMap[Eigen::RowVector3i::Zero()];
    startNode.position = _startPoint.array().round().cast<int>();
    addJumpNode(nullptr, startNode);

    int iterCount = 0;
    std::vector<Eigen::RowVector3i> bestPaths;
    while (!openList.empty()) {
        GridNode* currentNode = openList.top();
        openList.pop();
        openSet.erase(currentNode->position);
        if (closedList.find(currentNode->position) != closedList.end()) {
            continue;
        }
        closedList.emplace(currentNode->position);

        if (currentNode->position == endPoint) {
            const GridNode* node = currentNode;
            while (node != nullptr) {
                bestPaths.emplace_back(node->position);
                node = node->parent;
            }
            std::reverse(bestPaths.begin(), bestPaths.end());
            break;
        }

        for (const auto& dir : currentNode->dirs) {
            findJumpNodes(currentNode, dir);
        }

        if (iterCount > maxiter) {
            std::cout << "Reached maximum iterations!" << std::endl;
            return {};
        }
        iterCount++;
    }

    return bestPaths;
}
//*********************************************************************************************//

//*********************************************************************************************//
bool JPSPathFinder2::_searchAlongLine(
    const Eigen::RowVector3i& point,
    const Eigen::RowVector3i& dir,
    GridNode& jumpNode
) {
    const auto& ortDirsSet = ortDirsMap.at(dir);
    Eigen::RowVector3i nextPoint1d = point + dir;
    std::vector<Eigen::RowVector3i> forcedNeighbors;
    forcedNeighbors.reserve(8);

    while (this->isMovable(nextPoint1d)) {
        if (nextPoint1d == endPoint) {
            jumpNode.position = nextPoint1d;
            return true;
        }
        for (const auto& ortDir : ortDirsSet) {
            Eigen::RowVector3i blockPoint = nextPoint1d + ortDir;
            //if (this->isObstacle(blockPoint) && this->isMovable(blockPoint + dir)) {
            if (!this->isMovable(blockPoint) && this->isMovable(blockPoint + dir)) {
                forcedNeighbors.emplace_back(ortDir + dir);
            }
        }
        if (!forcedNeighbors.empty()) {
            forcedNeighbors.emplace_back(dir);
            jumpNode.position = nextPoint1d;
            jumpNode.dirs = forcedNeighbors;
            return true;
        }
        nextPoint1d += dir;
    }
    return false;
}
//*********************************************************************************************//

//*********************************************************************************************//
void JPSPathFinder2::searchAlongLine(
    GridNode* currentNode,
    const Eigen::RowVector3i& dir
) {
    GridNode jumpNode;
    if (_searchAlongLine(currentNode->position, dir, jumpNode)) {
        addJumpNode(currentNode, jumpNode);
    }
}
//*********************************************************************************************//

//*********************************************************************************************//
bool JPSPathFinder2::_searchAlongDiagonal2d(
    const Eigen::RowVector3i& point,
    const Eigen::RowVector3i& dir,
    std::vector<GridNode>& jumpNodes
) {
    const auto& dirDecs = dirDecsMap[dir];
    Eigen::RowVector3i dirDec1 = dirDecs[0];
    Eigen::RowVector3i dirDec2 = dirDecs[1];
    Eigen::RowVector3i nextPoint2d = point;
    std::vector<Eigen::RowVector3i> searchDirs;
    searchDirs.reserve(8);

    GridNode currentJumpNode;
    auto processStep = [&](
        const Eigen::RowVector3i& B,
        const Eigen::RowVector3i& dirDec,
        bool needForcedExplore) -> bool {

            if (nextPoint2d == endPoint) {
                currentJumpNode.position = nextPoint2d;
                jumpNodes.emplace_back(currentJumpNode);
                return true;
            }

            if (needForcedExplore) {
                Eigen::RowVector3i N = B + dirDec;
                if (this->isMovable(N)) {
                    searchDirs.emplace_back(N - nextPoint2d);
                }
            }

            GridNode nextJumpNode1, nextJumpNode2;
            bool isExitJumpNode1 = _searchAlongLine(nextPoint2d, dirDec1, nextJumpNode1);
            bool isExitJumpNode2 = _searchAlongLine(nextPoint2d, dirDec2, nextJumpNode2);
            if (isExitJumpNode1 || isExitJumpNode2 || !searchDirs.empty()) {
                searchDirs.emplace_back(dir);
                currentJumpNode.position = nextPoint2d;
                currentJumpNode.dirs = searchDirs;
                jumpNodes.emplace_back(currentJumpNode);

                if (isExitJumpNode1) {
                    jumpNodes.emplace_back(nextJumpNode1);
                }

                if (isExitJumpNode2) {
                    jumpNodes.emplace_back(nextJumpNode2);
                }
                return true;
            }
            return false;
        };

    while (true) {
        Eigen::RowVector3i B1 = nextPoint2d + dirDec1;
        Eigen::RowVector3i B2 = nextPoint2d + dirDec2;
        nextPoint2d += dir;
        if (!this->isMovable(nextPoint2d)) return false;

        bool isB1Movable = this->isMovable(B1);
        bool isB2Movable = this->isMovable(B2);

        if (isB1Movable) { // B1 is movable
            if (isB2Movable) { // Case 1: B1 and B2 are movable
                if (processStep(B1, dirDec1, false)) {
                    return true;
                }
            }
            else { // Case 2: B1 is movable, but B2 is not
                if (processStep(B2, dirDec2, true)) {
                    return true;
                }
            }
        }
        else { // B1 is not movable
            if (isB2Movable) {
                // Case 3: B2 is movable, but B1 is not
                if (processStep(B1, dirDec1, true)) {
                    return true;
                }
            }
            else { // Case 4: B1 and B2 are not movable
                // Neither is movable, nothing happens, and the diagonal movement in this direction is terminated
                return false;
            }
        }
    }
}
//*********************************************************************************************//

//*********************************************************************************************//
void JPSPathFinder2::searchAlongDiagonal2d(
    GridNode* currentNode,
    const Eigen::RowVector3i& dir
) {
    std::vector<GridNode> jumpNodes;
    jumpNodes.reserve(3);
    if (_searchAlongDiagonal2d(currentNode->position, dir, jumpNodes)) {

        if (!addJumpNode(currentNode, jumpNodes[0])) return;

        auto& handle = openSet.at(jumpNodes[0].position);
        GridNode* nextJumpNode = *handle;
        for (size_t i = 1; i < jumpNodes.size(); ++i) {
            addJumpNode(nextJumpNode, jumpNodes[i]);
        }
    }
}
//*********************************************************************************************//

//*********************************************************************************************//

bool JPSPathFinder2::_searchAlongDiagonal3d(
    const Eigen::RowVector3i& point,
    const Eigen::RowVector3i& dir,
    std::pair<GridNode, std::vector<std::vector<GridNode>>>& jumpNodesPair
) {
    Eigen::RowVector3i dirDec1(dir.x(), 0, 0);
    Eigen::RowVector3i dirDec2(0, dir.y(), 0);
    Eigen::RowVector3i dirDec3(0, 0, dir.z());

    Eigen::RowVector3i dirDiag1(dir.x(), dir.y(), 0);
    Eigen::RowVector3i dirDiag2(dir.x(), 0, dir.z());
    Eigen::RowVector3i dirDiag3(0, dir.y(), dir.z());
    Eigen::RowVector3i nextPoint3d = point;

    std::vector<Eigen::RowVector3i> searchDirs;
    searchDirs.reserve(8);
    GridNode currentJumpNode;
    auto processStep = [&](
        const Eigen::RowVector3i& blockPoint,
        const Eigen::RowVector3i& dirDiag,
        bool needForcedExplore) -> bool {
            if (nextPoint3d == endPoint) {
                currentJumpNode.position = nextPoint3d;
                jumpNodesPair.first = currentJumpNode;
                return true;
            }
            if (needForcedExplore) {
                Eigen::RowVector3i N = blockPoint + dirDiag;
                if (this->isMovable(N)) {
                    searchDirs.emplace_back(N - nextPoint3d);
                }
            }
            std::vector<GridNode> nextJumpNodes1, nextJumpNodes2, nextJumpNodes3;
            bool isExitJumpNode1 = _searchAlongDiagonal2d(nextPoint3d, dirDiag1, nextJumpNodes1);
            bool isExitJumpNode2 = _searchAlongDiagonal2d(nextPoint3d, dirDiag2, nextJumpNodes2);
            bool isExitJumpNode3 = _searchAlongDiagonal2d(nextPoint3d, dirDiag3, nextJumpNodes3);
            if (isExitJumpNode1 || isExitJumpNode2 || isExitJumpNode3 || !searchDirs.empty()) {
                searchDirs.emplace_back(dir);
                currentJumpNode.position = nextPoint3d;
                currentJumpNode.dirs = searchDirs;
                jumpNodesPair.first = currentJumpNode;
                std::vector<std::vector<GridNode>> jumpNodesVec;
                jumpNodesVec.reserve(3);
                if (isExitJumpNode1) {
                    jumpNodesVec.emplace_back(nextJumpNodes1);
                }
                if (isExitJumpNode2) {
                    jumpNodesVec.emplace_back(nextJumpNodes2);
                }
                if (isExitJumpNode3) {
                    jumpNodesVec.emplace_back(nextJumpNodes3);
                }
                jumpNodesPair.second = jumpNodesVec;
                return true;
            }
            return false;
        };

    while (true) {
        Eigen::RowVector3i B1 = nextPoint3d + dirDec1;
        Eigen::RowVector3i B2 = nextPoint3d + dirDec2;
        Eigen::RowVector3i B3 = nextPoint3d + dirDec3;
        nextPoint3d += dir;
        if (!this->isMovable(nextPoint3d)) return false;

        bool isB1Movable = this->isMovable(B1);
        bool isB2Movable = this->isMovable(B2);
        bool isB3Movable = this->isMovable(B3);

        bool isDiag1Movable = isB1Movable || isB2Movable;
        bool isDiag2Movable = isB1Movable || isB3Movable;
        bool isDiag3Movable = isB2Movable || isB3Movable;

        if (isDiag1Movable) {
            if (isDiag2Movable) {
                if (isDiag3Movable) {
                    // Case 1: B1, B2 and B3 are movable
                    if (processStep(B1, dirDiag1, false)) {
                        return true;
                    }
                }
                else {
                    // Case 2: B1 and B2 are movable, but B3 is not
                    if (processStep(B3, dirDiag3, true)) {
                        return true;
                    }
                }
            }
            else {
                if (isDiag3Movable) {
                    // Case 3: B1 and B3 are movable, but B2 is not
                    if (processStep(B2, dirDiag2, true)) {
                        return true;
                    }
                }
                else {
                    // Case 4: B1 is movable, but B2 and B3 are not movable
                    if (processStep(B2, dirDec2, true)) {
                        return true;
                    }
                    if (processStep(B3, dirDec3, true)) {
                        return true;
                    }
                }
            }
        }
        else {
            if (isDiag2Movable) {
                if (isDiag3Movable) {
                    // Case 5: B1 is not movable, but B2 and B3 are movable
                    if (processStep(B1, dirDec1, true)) {
                        return true;
                    }
                }
                else {
                    // Case 6: B1 and B3 are not movable, but B2 is movable
                    if (processStep(B1, dirDec1, true)) {
                        return true;
                    }
                    if (processStep(B3, dirDec3, true)) {
                        return true;
                    }
                }
            }
            else {
                if (isDiag3Movable) {
                    // Case 7: B1 and B2 are not movable, but B3 is movable
                    if (processStep(B1, dirDec1, true)) {
                        return true;
                    }
                    if (processStep(B2, dirDec2, true)) {
                        return true;
                    }
                }
                else {
                    // Case 8: B1, B2 and B3 are not movable
                    return false;
                }
            }
        }
    }
}
//*********************************************************************************************//

//*********************************************************************************************//
void JPSPathFinder2::searchAlongDiagonal3d(
    GridNode* currentNode,
    const Eigen::RowVector3i& dir
) {
    std::pair<GridNode, std::vector<std::vector<GridNode>>> jumpNodesPair;
    if (_searchAlongDiagonal3d(currentNode->position, dir, jumpNodesPair)) {

        if (!addJumpNode(currentNode, jumpNodesPair.first)) return;
        auto& jumpNodesVec = jumpNodesPair.second;

        if (!jumpNodesVec.empty()) {

            auto& handle = openSet.at(jumpNodesPair.first.position);
            GridNode* nextJumpNode = *handle;

            for (const auto& jumpNodes : jumpNodesVec) {
                if (jumpNodes.empty()) continue;

                if (!addJumpNode(nextJumpNode, jumpNodes[0])) return;

                auto& child_handle = openSet.at(jumpNodes[0].position);
                GridNode* nextNextJumpNode = *child_handle;

                for (size_t i = 1; i < jumpNodes.size(); ++i) {
                    addJumpNode(nextNextJumpNode, jumpNodes[i]);
                }
            }
        }

    }
}
//*********************************************************************************************//

//*********************************************************************************************//
void JPSPathFinder2::findJumpNodes(
    GridNode* currentNode,
    const Eigen::RowVector3i& dir
) {
    int normDir = dir.array().abs().sum();

    switch (normDir) {
    case 1: // 1D-line search
        //std::cout << "\n1D-line search" << std::endl;
        searchAlongLine(currentNode, dir);
        //std::cout << "1D-line search\n" << std::endl;
        return;
    case 2: // 2D-diagonal search 
        //std::cout << "\n2D-line search" << std::endl;
        searchAlongDiagonal2d(currentNode, dir);
        //std::cout << "2D-line search\n" << std::endl;
        return;
    case 3: // 3D-diagonal search
        //std::cout << "\n3D-line search" << std::endl;
        searchAlongDiagonal3d(currentNode, dir);
        //std::cout << "3D-line search\n" << std::endl;
        return;
    }
}
//*********************************************************************************************//

//*********************************************************************************************//
inline bool JPSPathFinder2::addJumpNode(
    GridNode* parent,
    const GridNode& jumpNode
) {
    const auto& point = jumpNode.position;
    if (closedList.find(point) != closedList.end()) {
        return false;
    }
    double tentative_gCost = (parent == nullptr)
        ? 0.0
        : parent->gCost + static_cast<double>((point - parent->position).squaredNorm());

    auto it = openSet.find(point);
    if (it == openSet.end()) { // create new node
        nodePool.emplace_back(std::make_unique<GridNode>(jumpNode));
        GridNode* newNode = nodePool.back().get();
        newNode->gCost = tentative_gCost;
        newNode->fCost = tentative_gCost + getHeu(point);
        newNode->parent = parent;
        newNode->heapkey = openList.push(newNode);
        openSet.emplace(point, newNode->heapkey);

    }
    else { // update exit node
        GridNode* existNode = *(it->second);
        if (existNode->gCost > tentative_gCost) {
            existNode->gCost = tentative_gCost;
            existNode->fCost = tentative_gCost + getHeu(point);
            existNode->parent = parent;
            existNode->dirs = jumpNode.dirs;
            openList.update(existNode->heapkey);
        }
    }
    return true;
}
//*********************************************************************************************//

//*********************************************************************************************//
inline double JPSPathFinder2::getHeu(const Eigen::RowVector3i& point) const {
    double squared_distance = static_cast<double>((point - endPoint).squaredNorm());
    double w = (squared_distance > DISTANCE_THRESHOLD_SQUARED) ? WEIGHT_HIGH : WEIGHT_LOW;
    return w * squared_distance;
}
//*********************************************************************************************//