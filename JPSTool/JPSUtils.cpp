#include "JPSUtils.h"
//************************************************************************************************************************//
DirDecsMap getDirDecsMap() {
	DirDecsMap dirDecsMap;
	dirDecsMap.reserve(20);
	//////////////////////////////// ID0: dir(-1, -1, -1) ///////////////////////////////////
	{
		Eigen::RowVector3i dir(-1, -1, -1);
		std::vector<Eigen::RowVector3i> dirDecs = {
			Eigen::RowVector3i(-1,  0,  0),
			Eigen::RowVector3i(0, -1,  0),
			Eigen::RowVector3i(0,  0, -1),
			Eigen::RowVector3i(-1, -1,  0),
			Eigen::RowVector3i(-1,  0, -1),
			Eigen::RowVector3i(0, -1, -1)
		};
		dirDecsMap[dir] = dirDecs;
	}

	//////////////////////////////// ID1: dir(0, -1, -1) ////////////////////////////////////
	{
		Eigen::RowVector3i dir(0, -1, -1);
		std::vector<Eigen::RowVector3i> dirDecs = {
			Eigen::RowVector3i(0,  0,  -1),
			Eigen::RowVector3i(0, -1,  0)
		};
		dirDecsMap[dir] = dirDecs;
	}

	//////////////////////////////// ID2: dir(1, -1, -1) ////////////////////////////////////
	{
		Eigen::RowVector3i dir(1, -1, -1);
		std::vector<Eigen::RowVector3i> dirDecs = {
			Eigen::RowVector3i(1,  0,  0),
			Eigen::RowVector3i(0, -1,  0),
			Eigen::RowVector3i(0,  0, -1),
			Eigen::RowVector3i(1, -1,  0),
			Eigen::RowVector3i(1,  0, -1),
			Eigen::RowVector3i(0, -1, -1)
		};
		dirDecsMap[dir] = dirDecs;
	}

	//////////////////////////////// ID3: dir(-1, 0, -1) ////////////////////////////////////
	{
		Eigen::RowVector3i dir(-1, 0, -1);
		std::vector<Eigen::RowVector3i> dirDecs = {
			Eigen::RowVector3i(0,  0, -1),
			Eigen::RowVector3i(-1,  0,  0)
		};
		dirDecsMap[dir] = dirDecs;
	}

	//////////////////////////////// ID4: dir(1, 0, -1) ////////////////////////////////////
	{
		Eigen::RowVector3i dir(1, 0, -1);
		std::vector<Eigen::RowVector3i> dirDecs = {
			Eigen::RowVector3i(0,  0, -1),
			Eigen::RowVector3i(1,  0,  0)
		};
		dirDecsMap[dir] = dirDecs;
	}

	//////////////////////////////// ID5: dir(-1, 1, -1) ////////////////////////////////////
	{
		Eigen::RowVector3i dir(-1, 1, -1);
		std::vector<Eigen::RowVector3i> dirDecs = {
			Eigen::RowVector3i(-1,  0,  0),
			Eigen::RowVector3i(0,  1,  0),
			Eigen::RowVector3i(0,  0, -1),
			Eigen::RowVector3i(-1,  1,  0),
			Eigen::RowVector3i(-1,  0, -1),
			Eigen::RowVector3i(0,  1, -1)
		};
		dirDecsMap[dir] = dirDecs;
	}

	//////////////////////////////// ID6: dir(0, 1, -1) ////////////////////////////////////
	{
		Eigen::RowVector3i dir(0, 1, -1);
		std::vector<Eigen::RowVector3i> dirDecs = {
			Eigen::RowVector3i(0,  0, -1),
			Eigen::RowVector3i(0,  1,  0)
		};
		dirDecsMap[dir] = dirDecs;
	}

	//////////////////////////////// ID7: dir(1, 1, -1) ////////////////////////////////////
	{
		Eigen::RowVector3i dir(1, 1, -1);
		std::vector<Eigen::RowVector3i> dirDecs = {
			Eigen::RowVector3i(1,  0,  0),
			Eigen::RowVector3i(0,  1,  0),
			Eigen::RowVector3i(0,  0, -1),
			Eigen::RowVector3i(1,  1,  0),
			Eigen::RowVector3i(1,  0, -1),
			Eigen::RowVector3i(0,  1, -1)
		};
		dirDecsMap[dir] = dirDecs;
	}

	//////////////////////////////// ID8: dir(-1, -1, 0) ////////////////////////////////////
	{
		Eigen::RowVector3i dir(-1, -1, 0);
		std::vector<Eigen::RowVector3i> dirDecs = {
			Eigen::RowVector3i(0, -1,  0),
			Eigen::RowVector3i(-1,  0,  0)
		};
		dirDecsMap[dir] = dirDecs;
	}


	//////////////////////////////// ID9: dir(1, -1, 0) ////////////////////////////////////
	{
		Eigen::RowVector3i dir(1, -1, 0);
		std::vector<Eigen::RowVector3i> dirDecs = {
			Eigen::RowVector3i(0, -1,  0),
			Eigen::RowVector3i(1,  0,  0)
		};
		dirDecsMap[dir] = dirDecs;
	}

	//////////////////////////////// ID10: dir(0, 0, 0) ////////////////////////////////////
	{
		Eigen::RowVector3i dir(0, 0, 0);
		std::vector<Eigen::RowVector3i> dirDecs = {
			Eigen::RowVector3i(1,  0,  0),
			Eigen::RowVector3i(-1,  0,  0),
			Eigen::RowVector3i(0,  1,  0),
			Eigen::RowVector3i(1,  1,  0),
			Eigen::RowVector3i(-1,  1,  0),
			Eigen::RowVector3i(0, -1,  0),
			Eigen::RowVector3i(1, -1,  0),
			Eigen::RowVector3i(-1, -1,  0),
			Eigen::RowVector3i(0,  0,  1),
			Eigen::RowVector3i(1,  0,  1),
			Eigen::RowVector3i(-1,  0,  1),
			Eigen::RowVector3i(0,  1,  1),
			Eigen::RowVector3i(1,  1,  1),
			Eigen::RowVector3i(-1,  1,  1),
			Eigen::RowVector3i(0, -1,  1),
			Eigen::RowVector3i(1, -1,  1),
			Eigen::RowVector3i(-1, -1,  1),
			Eigen::RowVector3i(0,  0, -1),
			Eigen::RowVector3i(1,  0, -1),
			Eigen::RowVector3i(-1,  0, -1),
			Eigen::RowVector3i(0,  1, -1),
			Eigen::RowVector3i(1,  1, -1),
			Eigen::RowVector3i(-1,  1, -1),
			Eigen::RowVector3i(0, -1, -1),
			Eigen::RowVector3i(1, -1, -1),
			Eigen::RowVector3i(-1, -1, -1)
		};
		dirDecsMap[dir] = dirDecs;
	}


	//////////////////////////////// ID11: dir(-1, 1, 0) ////////////////////////////////////
	{
		Eigen::RowVector3i dir(-1, 1, 0);
		std::vector<Eigen::RowVector3i> dirDecs = {
			Eigen::RowVector3i(0,  1,  0),
			Eigen::RowVector3i(-1,  0,  0)
		};
		dirDecsMap[dir] = dirDecs;
	}

	//////////////////////////////// ID12: dir(1, 1, 0) ////////////////////////////////////
	{
		Eigen::RowVector3i dir(1, 1, 0);
		std::vector<Eigen::RowVector3i> dirDecs = {
			Eigen::RowVector3i(0, 1, 0),
			Eigen::RowVector3i(1, 0, 0)
		};
		dirDecsMap[dir] = dirDecs;
	}

	//////////////////////////////// ID13: dir(-1, -1, 1) ////////////////////////////////////
	{
		Eigen::RowVector3i dir(-1, -1, 1);
		std::vector<Eigen::RowVector3i> dirDecs = {
			Eigen::RowVector3i(-1,  0,  0),
			Eigen::RowVector3i(0, -1,  0),
			Eigen::RowVector3i(0,  0,  1),
			Eigen::RowVector3i(-1, -1,  0),
			Eigen::RowVector3i(-1,  0,  1),
			Eigen::RowVector3i(0, -1,  1)
		};
		dirDecsMap[dir] = dirDecs;
	}

	//////////////////////////////// ID14: dir(0, -1, 1) ////////////////////////////////////
	{
		Eigen::RowVector3i dir(0, -1, 1);
		std::vector<Eigen::RowVector3i> dirDecs = {
			Eigen::RowVector3i(0,  0,  1),
			Eigen::RowVector3i(0, -1,  0)
		};
		dirDecsMap[dir] = dirDecs;
	}

	//////////////////////////////// ID15: dir(1, -1, 1) ////////////////////////////////////
	{
		Eigen::RowVector3i dir(1, -1, 1);
		std::vector<Eigen::RowVector3i> dirDecs = {
			Eigen::RowVector3i(1,  0,  0),
			Eigen::RowVector3i(0, -1,  0),
			Eigen::RowVector3i(0,  0,  1),
			Eigen::RowVector3i(1, -1,  0),
			Eigen::RowVector3i(1,  0,  1),
			Eigen::RowVector3i(0, -1,  1)
		};
		dirDecsMap[dir] = dirDecs;
	}

	//////////////////////////////// ID16: dir(-1, 0, 1) ////////////////////////////////////
	{
		Eigen::RowVector3i dir(-1, 0, 1);
		std::vector<Eigen::RowVector3i> dirDecs = {
			Eigen::RowVector3i(0,  0,  1),
			Eigen::RowVector3i(-1,  0,  0)
		};
		dirDecsMap[dir] = dirDecs;
	}

	//////////////////////////////// ID17: dir(1, 0, 1) ////////////////////////////////////
	{
		Eigen::RowVector3i dir(1, 0, 1);
		std::vector<Eigen::RowVector3i> dirDecs = {
			Eigen::RowVector3i(0, 0, 1),
			Eigen::RowVector3i(1, 0, 0)
		};
		dirDecsMap[dir] = dirDecs;
	}

	//////////////////////////////// ID18: dir(-1, 1, 1) ////////////////////////////////////
	{
		Eigen::RowVector3i dir(-1, 1, 1);
		std::vector<Eigen::RowVector3i> dirDecs = {
			Eigen::RowVector3i(-1,  0,  0),
			Eigen::RowVector3i(0,  1,  0),
			Eigen::RowVector3i(0,  0,  1),
			Eigen::RowVector3i(-1,  1,  0),
			Eigen::RowVector3i(-1,  0,  1),
			Eigen::RowVector3i(0,  1,  1)
		};
		dirDecsMap[dir] = dirDecs;
	}

	//////////////////////////////// ID19: dir(0, 1, 1) ////////////////////////////////////
	{
		Eigen::RowVector3i dir(0, 1, 1);
		std::vector<Eigen::RowVector3i> dirDecs = {
			Eigen::RowVector3i(0, 0, 1),
			Eigen::RowVector3i(0, 1, 0)
		};
		dirDecsMap[dir] = dirDecs;
	}

	//////////////////////////////// ID20: dir(1, 1, 1) ////////////////////////////////////
	{
		Eigen::RowVector3i dir(1, 1, 1);
		std::vector<Eigen::RowVector3i> dirDecs = {
			Eigen::RowVector3i(1, 0, 0),
			Eigen::RowVector3i(0, 1, 0),
			Eigen::RowVector3i(0, 0, 1),
			Eigen::RowVector3i(1, 1, 0),
			Eigen::RowVector3i(1, 0, 1),
			Eigen::RowVector3i(0, 1, 1)
		};
		dirDecsMap[dir] = dirDecs;
	}
	return dirDecsMap;
}
//************************************************************************************************************************//

//************************************************************************************************************************//
OrtDirsMap getOrtDirsMap() {
	OrtDirsMap _OrtDirsMap;
	_OrtDirsMap.reserve(6);

	std::vector<Eigen::RowVector3i> mainDirs = {
		{1, 0, 0}, {-1, 0, 0},
		{0, 1, 0}, {0, -1, 0},
		{0, 0, 1}, {0, 0, -1}
	};

	auto generateOrthogonal = [](const Eigen::RowVector3i& dir) -> std::vector<Eigen::RowVector3i> {
		std::vector<Eigen::RowVector3i> ortVecs;
		ortVecs.reserve(8);
		for (int x = -1; x <= 1; ++x) {
			for (int y = -1; y <= 1; ++y) {
				for (int z = -1; z <= 1; ++z) {
					if (x == 0 && y == 0 && z == 0) continue;

					int dot = x * dir.x() + y * dir.y() + z * dir.z();

					if (dot == 0) ortVecs.emplace_back(x, y, z);
				}
			}
		}
		return ortVecs;
		};

	for (const auto& dir : mainDirs) {
		_OrtDirsMap.emplace(dir, generateOrthogonal(dir));
	}

	return _OrtDirsMap;
}
//************************************************************************************************************************//

//************************************************************************************************************************//
void JPSAABBEnv::initGridMap(
	const std::unordered_map<size_t, SCData>& SCMap,
	const std::unordered_map<size_t, EFData>& EFMap
) {
	Mesh scUnionHull, efUnionHull;
	for (const auto& [_, sc] : SCMap) {
		Mesh currentHull = buildConvexHullMesh(sc.obbBoxVertices);
		Mesh temp;
		if (PMP::corefine_and_compute_union(
			scUnionHull, currentHull, temp,
			PMP::parameters::allow_self_intersections(true),
			PMP::parameters::throw_on_self_intersection(false)
		)) {
			scUnionHull = std::move(temp);
		}
	}

	if (!EFMap.empty()) {
		for (const auto& [_, ef] : EFMap) {
			Mesh currentHull = buildConvexHullMesh(ef.obbBoxExpand1Vertices);
			Mesh temp;
			if (PMP::corefine_and_compute_union(
				efUnionHull, currentHull, temp,
				PMP::parameters::allow_self_intersections(true),
				PMP::parameters::throw_on_self_intersection(false)
			)) {
				efUnionHull = std::move(temp);
			}
		}

		Mesh tempDiff;
		if (PMP::corefine_and_compute_difference(
			scUnionHull, efUnionHull, tempDiff,
			PMP::parameters::allow_self_intersections(true),
			PMP::parameters::throw_on_self_intersection(false)
		)) {
			scUnionHull = std::move(tempDiff);
		}
	}
	detector = std::make_unique<DetectionObjects>(std::move(scUnionHull));
}
//************************************************************************************************************************//

//************************************************************************************************************************//
void JPSVoxelEnv::initGridMap(
	const std::unordered_map<size_t, SCData>& SCMap,
	const std::unordered_map<size_t, EFData>& EFMap
) {
	RowVector4iSet efGridPointsSet, scGridPointsSet;
	Mesh scCombinedHull, efCombinedHull, combinedHullsIntersection;
	for (const auto& [_, sc] : SCMap) {
		Mesh currentHull = buildConvexHullMesh(sc.obbBoxVertices);
		Mesh temp;
		if (PMP::corefine_and_compute_union(
			scCombinedHull, currentHull, temp,
			PMP::parameters::allow_self_intersections(true),
			PMP::parameters::throw_on_self_intersection(false)
		)) {
			scCombinedHull = std::move(temp);
		}
	}

	if (!EFMap.empty()) {
		for (const auto& [_, ef] : EFMap) {
			Mesh currentHull = buildConvexHullMesh(ef.obbBoxExpand1Vertices);
			Mesh temp;
			if (PMP::corefine_and_compute_union(
				efCombinedHull, currentHull, temp,
				PMP::parameters::allow_self_intersections(true),
				PMP::parameters::throw_on_self_intersection(false)
			)) {
				efCombinedHull = std::move(temp);
			}
		}

		Mesh tempIntersect;
		if (PMP::corefine_and_compute_intersection(
			scCombinedHull, efCombinedHull, tempIntersect,
			PMP::parameters::allow_self_intersections(true),
			PMP::parameters::throw_on_self_intersection(false)
		)) {
			combinedHullsIntersection = std::move(tempIntersect);
		}

		RowVector3iSet efGridPoints = rasterizePolyhedron(combinedHullsIntersection);
		for (const auto& point : efGridPoints) {
			GridPointsSet.emplace(point(0), point(1), point(2), 1);
		}
	}

	RowVector3iSet scGridPoints = rasterizePolyhedron(scCombinedHull);
	for (const auto& point : scGridPoints) {
		GridPointsSet.emplace(point(0), point(1), point(2), 0);
	}

	//std::ofstream file("GridPoints.csv", std::ios::binary);
	//if (file.is_open()) {
	//	for (const auto& point : GridPointsSet) {
	//		file << point[0] << "," << point[1] << "," << point[2] << "," << point[3] << "\n";
	//	}
	//	file.close();
	//}
}
//************************************************************************************************************************//