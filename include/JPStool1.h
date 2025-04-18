#ifndef JPSTOOL1_H_
#define JPSTOOL1_H_
#include "JPSUtils.h"
/**​
* @class JPSPathFinder1
* @brief A path finding class implementing Jump Point Search(JPS) algorithm in 3D space.
*
*Inherits from JPSenv1 environment class and provides graph search capabilities
* using memory - efficient JPS algorithm with obstacle checking and path optimization.
*/
class JPSPathFinder1 : public JPSenv1
{
public:
	JPSPathFinder1(
		const std::unordered_map<size_t, SCData>& SCSet,
		const std::unordered_map<size_t, EFData>& EFSet)
		:JPSenv1(SCSet, EFSet),
		endPoint(Eigen::RowVector3i::Zero()),
		openList(),
		closedList(),
		openSet(),
		dirDecsMap(getDirDecsMap()),
		ortDirsMap(getOrtDirsMap()),
		nodePool()
	{
	}
	/**​
	* @brief Main graph search function to find path between two points
	* @param _startPoint 3D starting position in continuous coordinates
	* @param _endPoint 3D target position in continuous coordinates
	* @return Vector of 3D grid points representing the found path
	*/
	std::vector<Eigen::RowVector3i> JPSGraphSearch(
		const Eigen::RowVector3d& _startPoint,
		const Eigen::RowVector3d& _endPoint
	);

private:
	Eigen::RowVector3i endPoint;
	priorityQueue openList; // Priority queue for open nodes
	absl::flat_hash_set<Eigen::RowVector3i, EigenRowVectoriHash> closedList; // closed list
	absl::flat_hash_map<Eigen::RowVector3i, priorityQueue::handle_type, EigenRowVectoriHash> openSet; // Fast lookup for open nodes

	DirDecsMap dirDecsMap; // Precomputed directional decomposition map
	OrtDirsMap ortDirsMap; // Orthogonal directions mapping for JPS pruning
	std::vector<std::unique_ptr<GridNode>> nodePool; // Memory pool for node allocation

	static constexpr double DISTANCE_THRESHOLD_SQUARED = 25.0; // 5.0^2
	static constexpr double WEIGHT_HIGH = 8.0;
	static constexpr double WEIGHT_LOW = 1.0;
	static constexpr int maxiter = 200000;

	//*********************************************************************************************//
	/**​
		* @brief Adds a jump node to the open list
		* @param parent Pointer to parent grid node
		* @param jumpNode Found jump node candidate
		* @return True if node was successfully added, False otherwise
	*/
	bool addJumpNode(GridNode* parent, const GridNode& jumpNode);

	//*********************************************************************************************//
	/**​
		* @brief Calculates heuristic cost between given point and target
		* @param point Current grid position
		* @return Heuristic cost estimate to target
	*/
	double getHeu(const Eigen::RowVector3i& point) const;
	//*********************************************************************************************//

	//*********************************************************************************************//

	bool _searchAlongLine(
		const Eigen::RowVector3i& point,
		const Eigen::RowVector3i& dir,
		GridNode& jumpNode
	);

	/**​
		* @brief Wrapper for straight - line search that handles node expansion
		* @param currentNode Pointer to current node being processed
		* @param dir Search direction vector
	*/
	void searchAlongLine(
		GridNode* currentNode,
		const Eigen::RowVector3i& dir
	);

	//*********************************************************************************************//

	bool _searchAlongDiagonal2d(
		const Eigen::RowVector3i& point,
		const Eigen::RowVector3i& dir,
		std::vector<GridNode>& jumpNodes
	);

	/**​
		* @brief Wrapper for 2D diagonal search with node expansion
		* @param currentNode Pointer to current node being processed
		* @param dir Diagonal direction vector
	*/
	void searchAlongDiagonal2d(
		GridNode* currentNode,
		const Eigen::RowVector3i& dir
	);

	//*********************************************************************************************//

	bool _searchAlongDiagonal3d(
		const Eigen::RowVector3i& point,
		const Eigen::RowVector3i& dir,
		std::pair<GridNode, std::vector<std::vector<GridNode>>>& jumpNodesPair
	);

	/**​
		* @brief Wrapper for 3D diagonal search with node expansion
		* @param currentNode Pointer to current node being processed
		* @param dir 3D diagonal direction vector
	*/
	void searchAlongDiagonal3d(
		GridNode* currentNode,
		const Eigen::RowVector3i& dir
	);

	//*********************************************************************************************//

	/**​
		* @brief Main jump point discovery function
		* @param currentNode: Pointer to current node being processed
		* @param dir: Current search direction vector
	*/
	void findJumpNodes(
		GridNode* currentNode,
		const Eigen::RowVector3i& dir
	);

};

#endif
