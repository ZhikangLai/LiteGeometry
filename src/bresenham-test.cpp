#include "Prepare.h"
#include "bresenham.h"

int main() {
	const auto& [SCSet, EFSet] = getDataSet("./testData/testdata_simple.json");
	std::string savePath = "./matlab_show/bresenham";
	//const auto& testPolyhedron = SCSet.at(0).obbBoxVertices;
	const auto& testPolyhedron = EFSet.at(3).obbBoxExpand1Vertices;
	//************************************** test 3D line ****************************************//
	{
		RowVector3iSet linePointsSet = rasterizeBresenhamLine3D(testPolyhedron.row(0), testPolyhedron.row(6));
		std::ofstream file(savePath + "/linePoints.csv", std::ios::binary);
		if (file.is_open()) {
			for (const auto& point : linePointsSet) {
				file << point[0] << "," << point[1] << "," << point[2] << "\n";
			}
			file.close();
		}
	}
	//************************************** test 3D line ****************************************//

	//************************************** test 3D polygon ****************************************//
	{
		RowVector3iSet facePointsSet = rasterizePolygon3D(testPolyhedron.block<4, 3>(2, 0));
		std::ofstream file(savePath + "/facePoints.csv", std::ios::binary);
		if (file.is_open()) {
			for (const auto& point : facePointsSet) {
				file << point[0] << "," << point[1] << "," << point[2] << "\n";
			}
			file.close();
		}
	}
	//************************************** test 3D polygon ****************************************//

	//************************************** test polyhedron ****************************************//
	{
		RowVector3iSet polyhedronPointsSet = rasterizePolyhedron(testPolyhedron);
		std::ofstream file(savePath + "/polyhedronPoints.csv", std::ios::binary);
		if (file.is_open()) {
			for (const auto& point : polyhedronPointsSet) {
				file << point[0] << "," << point[1] << "," << point[2] << "\n";
			}
			file.close();
		}
	}
	//************************************** test polyhedron ****************************************//
	return 0;
}
