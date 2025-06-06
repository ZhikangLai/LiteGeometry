#include "Prepare.h"
#include "JPSAABBEnv.h"
#include "JPSVoxelEnv.h"
#include <chrono>  
int main() {
    const auto& [SCMap,EFMap] = loadVolumeMaps("./testData/testdata_simple.json");

	Eigen::RowVector3d startPoint(748876.1520, 2564890.4544, 65.5);
	Eigen::RowVector3d endPoint1(748678.2699, 2564651.6236, 50.0761);
	Eigen::RowVector3d endPoint2(748378.7634, 2564889.0147, 45.0761);

	JPSAABBPathFinder _JPSAABBPathFinder(SCMap, EFMap);
	auto start1 = std::chrono::high_resolution_clock::now();
	const auto& bestPaths1 = _JPSAABBPathFinder.JPSGraphSearch(startPoint, endPoint1);
	auto end1 = std::chrono::high_resolution_clock::now();

	auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1);
	std::cout << "JPSAABBPathFinder execution time: " << duration1.count() << " ms" << std::endl;

	JPSVoxelPathFinder _JPSVoxelPathFinder(SCMap, EFMap);
	auto start2 = std::chrono::high_resolution_clock::now();
	const auto& bestPaths2 = _JPSVoxelPathFinder.JPSGraphSearch(startPoint, endPoint1);
	auto end2 = std::chrono::high_resolution_clock::now();

	auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2);
	std::cout << "JPSVoxelPathFinder execution time: " << duration2.count() << " ms" << std::endl;

	if (!bestPaths1.empty()) {
		std::ofstream file1("./testResult/test_jps/JPSAABBPath.csv", std::ios::binary);
		if (file1.is_open()) {
			for (const auto& point : bestPaths1) {
				file1 << point[0] << "," << point[1] << "," << point[2] << "\n";
			}
			file1.close();
		}
	}

	if (!bestPaths2.empty()) {
		std::ofstream file2("./testResult/test_jps/JPSVoxelPath.csv", std::ios::binary);
		if (file2.is_open()) {
			for (const auto& point : bestPaths2) {
				file2 << point[0] << "," << point[1] << "," << point[2] << "\n";
			}
			file2.close();
		}
	}

	return 0;
}

