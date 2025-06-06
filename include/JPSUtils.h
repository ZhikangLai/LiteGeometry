#ifndef JPSUTILS_H_
#define JPSUTILS_H_
#include "Utility.h"
#include "Prepare.h"
#include "RasterizationTool.h"
#include "LightGeoBase.h"
#include "AbslEigenRowVectoriHash.h"
#include <xxhash.hpp>
#include <absl/container/flat_hash_map.h>
#include <boost/heap/d_ary_heap.hpp> 

using DirDecsMap = std::unordered_map<Eigen::RowVector3i, std::vector<Eigen::RowVector3i>, EigenRowVectoriHash>;
/**
 * @brief Constructs a mapping for decomposing diagonal 3D directions.
 *
 * This function creates and returns a `DirDecsMap` that associates each primary 3D integer vector 
 * (`Eigen::RowVector3i`) with its valid decomposition directions. 
 * 
 * Moreover, when the key is (0, 0, 0), the returned decomposition vectors include all 26 possible directions
 *
 * @return A DirDecsMap mapping each diagonal 3D direction vector (`Eigen::RowVector3i`) to its corresponding
 *         decomposition vectors (`std::vector<Eigen::RowVector3i>`).
 */
DirDecsMap getDirDecsMap();

using OrtDirsMap = std::unordered_map<Eigen::RowVector3i, std::vector<Eigen::RowVector3i>, EigenRowVectoriHash>;
/**
 * @brief Generates a mapping from main directions to their orthogonal directions.
 *
 * This function creates and returns an `OrtDirsMap` that links each of the six principal 3D 
 * direction vectors (e.g., (1,0,0), (-1,0,0), etc.) with a collection of all 3D integer vectors 
 * orthogonal to their corresponding main direction.
 *
 * @return An OrtDirsMap mapping each main 3D direction vector (`Eigen::RowVector3i`) to a vector of 3D integer vectors
 *         (`std::vector<Eigen::RowVector3i>`) that are orthogonal to it.
 */
OrtDirsMap getOrtDirsMap();

//************************************************************************************************************************//
struct GridNode;
template <typename T>
struct GridNodePtrComparator {
	bool operator()(const T& n1, const T& n2) const {
		return n1->fCost > n2->fCost;
	}
};

using priorityQueue = boost::heap::d_ary_heap<
	GridNode*,
	boost::heap::mutable_<true>,
	boost::heap::arity<2>,
	boost::heap::compare<GridNodePtrComparator<GridNode*>>
>;

struct GridNode {
	Eigen::RowVector3i position = Eigen::RowVector3i::Zero();
	std::vector<Eigen::RowVector3i> dirs;
	double gCost = inf;
	double fCost = inf;
	GridNode* parent = nullptr;
	priorityQueue::handle_type heapkey;
};
//************************************************************************************************************************//

struct DetectionObjects {
	Mesh scUnion;
	std::unique_ptr<Tree> scTree;
	std::unique_ptr<Point_inside> insideSCDetecter;

	DetectionObjects(Mesh&& sc)
		: scUnion(std::move(sc))
	{
		if (!scUnion.is_empty()) {
			scTree = std::make_unique<Tree>(faces(scUnion).first, faces(scUnion).second, scUnion);
			insideSCDetecter = std::make_unique<Point_inside>(*scTree);
		}
	}
	~DetectionObjects() = default;
};

class JPSAABBEnv
{
public:
	JPSAABBEnv(
		const std::unordered_map<size_t, SCData>& SCSet,
		const std::unordered_map<size_t, EFData>& EFSet)
		:
		detector(nullptr)
	{
		initGridMap(SCSet, EFSet);
	}

protected:
	inline bool isMovable(const Eigen::RowVector3i& point) const {
		return (*detector->insideSCDetecter)(Point_3(point.x(), point.y(), point.z())) == CGAL::ON_BOUNDED_SIDE;
	}

private:
	std::unique_ptr<DetectionObjects> detector;
	/**
	 * @brief Initializes the navigability map using Boolean mesh operations.
	 *
	 * Constructs a collision detection structure by combining navigable regions (SCData)
	 * and subtracting obstacle regions (EFData) through CGAL mesh processing.
	 * 
	 * Process flow:
	 *   1. For each entry in SCSet:
	 *      - Generate convex hull from OBB vertices (sc.obbBoxVertices).
	 *      - Merge into cumulative SC union mesh via Boolean union.
	 *
	 *   2. If EFSet is non-empty:
	 *      - For each entry in EFSet:
	 *        - Build obstacle mesh from expanded OBB vertices (ef.obbBoxExpand1Vertices).
	 *        - Merge into cumulative EF union mesh via Boolean union.
	 *		- Subtract EF union from SC union via Boolean difference.
	 *
	 *   3. Finalize Collision Detection: 
	 *      - Construct AABB tree acceleration structure from final mesh
	 *		- Initialize detector with spatial query optimization
	 *
	 * @param[in] SCSet: `unordered_map<size_t, SCData>` where each value contains
	 *             original OBB vertices for constructing navigable convex hulls.
	 * @param[in] EFSet: `unordered_map<size_t, EFData>` where each value provides
	 *             expanded OBB vertices defining obstacle regions to exclude.
	 */
	void initGridMap(
		const std::unordered_map<size_t, SCData>& SCSet,
		const std::unordered_map<size_t, EFData>& EFSet
	);
};
//************************************************************************************************************************//



//************************************************************************************************************************//
// Hash functor for Eigen::RowVector4i.
// Computes the hash using only the first three components (3D coordinate),
// ignoring the fourth component (label).
struct EigenRowVector4iHash {
	size_t operator()(const Eigen::RowVector4i& v) const {
		return static_cast<size_t>(xxh::xxhash<64>(v.data(), sizeof(int) * 3, kSeed));
	}
};

// Equality functor for Eigen::RowVector4i.
// Compares only the first three components, so vectors with the same 3D coordinates
// are considered equal regardless of their label.
struct EigenRowVector4iEqual {
	bool operator()(const Eigen::RowVector4i& v1, const Eigen::RowVector4i& v2) const {
		return v1.head<3>() == v2.head<3>(); 
	}
};

using RowVector4iSet = absl::flat_hash_set<Eigen::RowVector4i, EigenRowVector4iHash, EigenRowVector4iEqual>;

class JPSVoxelEnv
{
public:
	JPSVoxelEnv(
		const std::unordered_map<size_t, SCData>& SCSet,
		const std::unordered_map<size_t, EFData>& EFSet)
		:
		GridPointsSet()
	{
		initGridMap(SCSet, EFSet);
	}

protected:

	//inline bool isObstacle(const Eigen::RowVector3i& point) const {
	//	auto it = GridPointsSet.find(Eigen::RowVector4i(point.x(), point.y(), point.z(), 1));
	//	if (it != GridPointsSet.end()) return static_cast<bool>((*it)(3));
	//	return false;
	//}

	inline bool isMovable(const Eigen::RowVector3i& point) const {
		auto it = GridPointsSet.find(Eigen::RowVector4i(point.x(), point.y(), point.z(), 0));
		if (it != GridPointsSet.end()) return false;
		return true;
	}
private:
	RowVector4iSet GridPointsSet;
	/**
	 * @brief Initializes a discrete surface voxel map for UAV collision detection using surface rasterization.
	 * 
	 * Constructs prohibited surface markers through convex hull operations and mesh rasterization. 
	 * Optimized for surface-only storage in large 3D environments.
	 * 
	 * Process flow:
	 *   1. Flight Zone polyhedrons Surface Construction:
	 *      - For each entry in SCSet:
	 *		  * Generate convex hull from OBB vertices (sc.obbBoxVertices).
	 *        * Merge into cumulative SC union mesh via Boolean union.
	 *
	 *   2. If EFSet is non-empty:
	 *      - For each entry in EFSet:
	 *        * Build obstacle mesh from expanded OBB vertices (ef.obbBoxExpand1Vertices).
	 *        * Merge into cumulative EF union mesh via Boolean union.
	 *		- Compute SC กษ EF intersection of polyhedrons.
	 *		- Rasterize intersection surfaces to voxels labeled 1 (prohibited obstacle surfaces).
	 * 
	 *   3. Final Surface Voxel Map:
	 *		- Rasterize SC union polyhedron surface to voxels labeled 0 (prohibited flight zone boundaries).
	 *      - Combine all surface voxels into GridPointsSet with XYZW format:
	 *        * W=0: Flight zone boundary surface
	 *        * W=1: Obstacle intrusion surface
	 *
	 * Key Constraints:
	 *    - Only polyhedron surface voxels are stored - polyhedron interiors are implicit free space.
	 *	  - All surface voxels (W=0/1) represent non-traversable space.
	 *    - Valid paths must stay strictly within polyhedron interiors without intersecting any surface voxels.
	 * 
	 * @param[in] SCSet: `unordered_map<size_t, SCData>` where each value contains
	 *             original OBB vertices for constructing flight zone polyhedrons.
	 * @param[in] EFSet: `unordered_map<size_t, EFData>` where each value provides
	 *             expanded OBB vertices defining obstacle polyhedrons to exclude.
	 */
	void initGridMap(
		const std::unordered_map<size_t, SCData>& SCSet,
		const std::unordered_map<size_t, EFData>& EFSet
	);
};

//************************************************************************************************************************//
#endif