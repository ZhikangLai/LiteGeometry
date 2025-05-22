#include "Prepare.h"
#include "lightGeo.h"

int main() {
	const auto& [SCMap, EFMap] = loadVolumeMaps("./testData/testdata_simple.json");
	std::string savePath = "./matlab_show/lightGeo";
	const auto& polyhedronVertices = EFMap.at(3).obbBoxExpand1Vertices;

	Eigen::Matrix<double,4,3> unorderedVertices;
	unorderedVertices << 748694.4250704022, 2564734.3476669602, 49.5,
		748674.4539194419, 2564739.5306861168, 81.5,
		748674.4154279609, 2564734.5007915981, 49.5,
		748694.4635618832, 2564739.3775614789, 81.5;

	Eigen::MatrixX3d polygon = generateClosedPolygon(unorderedVertices);
	Eigen::RowVector3d testPoint(748674.435, 2564737.016, 65.5);

	{// test generateClosedPolygon

		std::cout << "\n========== Test: generateClosedPolygon ==========\n";
		Eigen::MatrixX3d polygon = generateClosedPolygon(unorderedVertices);
		std::cout << "[Input] Original unordered vertices (n x 3):\n";
		std::cout << unorderedVertices.format(Eigen::FullPrecision) << "\n\n";

		std::cout << "[Output] Generated closed polygon (n x 3):\n";
		std::cout << polygon.format(Eigen::FullPrecision) << "\n";

		std::cout << "========== End of Test: generateClosedPolygon ==========\n";

	}

	{// test computePlaneNormal and compute3Dto2DTransformMatrix

		std::cout << "\n========== Test: computePlaneNormal & compute3Dto2DTransformMatrix ==========\n";

		const auto& [v1, v2] = computePlaneNormal(unorderedVertices);
		std::cout << "[computePlaneNormal] Result:\n";
		std::cout << "  Normal Vector 1 (v1): " << v1.format(Eigen::FullPrecision) << "\n";
		std::cout << "  Normal Vector 2 (v2): " << v2.format(Eigen::FullPrecision) << "\n\n";

		Eigen::Matrix<double, 3, 2> T = compute3Dto2DTransformMatrix(v1);
		std::cout << "[compute3Dto2DTransformMatrix] Projection matrix from normal v1 (3x2):\n";
		std::cout << T.format(Eigen::FullPrecision) << "\n";

		std::cout << "========== End of Test: computePlaneNormal & compute3Dto2DTransformMatrix ==========\n";

	}

	{// test isCoplanar

		std::cout << "\n========== Test: isCoplanar ==========\n";
		bool _isCoplanar = isCoplanar(unorderedVertices);
		if (_isCoplanar) {
			std::cout << "The provided points are coplanar" << std::endl;
		}
		else {
			std::cout << "The provided points are not coplanar" << std::endl;
		}
		std::cout << "========== End of Test: isCoplanar ==========\n";

	}

	{// test isPointOnLine3D

		std::cout << "\n========== Test: isPointOnLine3D ==========\n";
	    Eigen::RowVector3d A(0, 0, 0);
		Eigen::RowVector3d B(1, 1, 1);
		Eigen::RowVector3d C(0.5, 0.5, 0.5);
		Eigen::RowVector3d D(2, 2, 2);
		Eigen::RowVector3d E(-1, -1, -1);
		Segment3D segAB = Segment3D{ A,B };
		Ray3D rayAB = Ray3D{ A,B };
		Line3D lineAB = Line3D{ A,B };


		std::cout << "--- On Segment AB ---\n";
		std::cout << "  C on Segment AB: " << isPointOnLine3D(segAB, C) << "\n";
		std::cout << "  D on Segment AB: " << isPointOnLine3D(segAB, D) << "\n";
		std::cout << "  E on Segment AB: " << isPointOnLine3D(segAB, E) << "\n\n";

		std::cout << "--- On Ray AB ---\n";
		std::cout << "  C on Ray AB: " << isPointOnLine3D(rayAB, C) << "\n";
		std::cout << "  D on Ray AB: " << isPointOnLine3D(rayAB, D) << "\n";
		std::cout << "  E on Ray AB: " << isPointOnLine3D(rayAB, E) << "\n\n";

		std::cout << "--- On Line AB ---\n";
		std::cout << "  C on Line AB: " << isPointOnLine3D(lineAB, C) << "\n";
		std::cout << "  D on Line AB: " << isPointOnLine3D(lineAB, D) << "\n";
		std::cout << "  E on Line AB: " << isPointOnLine3D(lineAB, E) << "\n";

		std::cout << "========== End of Test: isPointOnLine3D ==========\n";

	}

	{// test isPointInPolygon3D

		std::cout << "\n========== Test: isPointInPolygon3D ==========\n";
		bool isInsideWithEdge = isPointInPolygon3D(unorderedVertices, testPoint);
		bool isInsideStrict = isPointInPolygon3D(unorderedVertices, testPoint, false);
		std::cout << "Including boundary (point on face/edge is considered inside): "
			<<  isInsideWithEdge << "\n\n";
		std::cout << "Strict interior only (point on boundary is considered outside): "
			<<  isInsideStrict << "\n";
		std::cout << "========== End of Test: isPointInPolygon3D ==========\n";

	}

	{// test isPointInPolyhedron

		std::cout << "\n========== Test: isPointInPolyhedron ==========\n";

		Eigen::RowVector3d bottleFaceCenter(748684.4399944221, 2564736.4396765385, 49);
		bool isInsideWithEdge = isPointInPolyhedron(polyhedronVertices, bottleFaceCenter);
		bool isInsideStrict = isPointInPolyhedron(polyhedronVertices, bottleFaceCenter, false);
		std::cout << "Including boundary (point on face/edge is considered inside): "
			      << isInsideWithEdge << "\n\n";
		std::cout << "Strict interior only (point on boundary is considered outside): "
				  << isInsideStrict << "\n";
		std::cout << "========== End of Test: isPointInPolyhedron ==========\n";

	}

	{// test computeLinesDistance
		std::cout << "\n========== Test: computeLinesDistance ==========\n";
		{
			Eigen::RowVector3d A(1, 1, 0);
			Eigen::RowVector3d B(0, -1, 0.4);
			Eigen::RowVector3d C(1, -1, 0);
			Eigen::RowVector3d D(0, 0, 0.5);
			Segment3D AB{ A ,B };
			Segment3D CD{ C ,D };
			const auto& [closestP1, closestP2, dist] = computeLinesDistance(AB, CD);
			std::cout << "Closest Point on Segment AB: " << closestP1.format(Eigen::FullPrecision) << "\n";
			std::cout << "Closest Point on Segment CD: " << closestP2.format(Eigen::FullPrecision) << "\n";
			std::cout << "Distance: " << dist << "\n";
		}

		{
			Eigen::RowVector3d A(1, 1, 0);
			Eigen::RowVector3d B(0, 0.5, 0.55);
			Eigen::RowVector3d C(1, -1, 0);
			Eigen::RowVector3d D(0, 0, 0.45);
			Ray3D AB{ A ,B };
			Ray3D CD{ C ,D };
			std::tuple<Eigen::RowVector3d, Eigen::RowVector3d, double> Intersection;
			const auto& [closestP1, closestP2, dist] = computeLinesDistance(AB, CD);
			std::cout << "Closest Point on Ray AB: " << closestP1.format(Eigen::FullPrecision) << "\n";
			std::cout << "Closest Point on Ray CD: " << closestP2.format(Eigen::FullPrecision) << "\n";
			std::cout << "Distance: " << dist << "\n";
		}
		std::cout << "========== End of Test: computeLinesDistance ==========\n";
	}

	{// test isLinesIntersection3D

		std::cout << "\n========== Test: isLinesIntersection3D ==========\n";
		{
			Eigen::RowVector3d A(1, 1, 0);
			Eigen::RowVector3d B(0, -1, 0.4);
			Eigen::RowVector3d C(1, -1, 0);
			Eigen::RowVector3d D(0, 0, 0.5);
			Segment3D AB{ A ,B };
			Segment3D CD{ C ,D };
			std::tuple<Eigen::RowVector3d, Eigen::RowVector3d, double> Intersection;
			bool isLinesIntersection = isLinesIntersection3D(AB, CD, Intersection, 0.2);
			std::cout << "\n--- Test 1: Segment Intersection ---\n";
			if (isLinesIntersection) {
				std::cout << "Closest Point on Segment AB: " << std::get<0>(Intersection).format(Eigen::FullPrecision) << "\n";
				std::cout << "Closest Point on Segment CD: " << std::get<1>(Intersection).format(Eigen::FullPrecision) << "\n";
				std::cout << "Distance: " << std::get<2>(Intersection) << "\n";
			}
		}

		{
			Eigen::RowVector3d A(1, 1, 0);
			Eigen::RowVector3d B(0, 0.5, 0.55);
			Eigen::RowVector3d C(1, -1, 0);
			Eigen::RowVector3d D(0, 0, 0.45);
			Ray3D AB{ A ,B };
			Ray3D CD{ C ,D };
			std::tuple<Eigen::RowVector3d, Eigen::RowVector3d, double> Intersection;
			bool isLinesIntersection = isLinesIntersection3D(AB, CD, Intersection, 0.2);
			std::cout << "\n--- Test 2: Ray Intersection ---\n";
			if (isLinesIntersection) {
				std::cout << "Closest Point on Ray AB: " << std::get<0>(Intersection).format(Eigen::FullPrecision) << "\n";
				std::cout << "Closest Point on Ray CD: " << std::get<1>(Intersection).format(Eigen::FullPrecision) << "\n";
				std::cout << "Distance: " << std::get<2>(Intersection) << "\n";
			}
		}
		std::cout << "========== End of Test: isLinesIntersection3D ==========\n";

	}

	{// test isLinePolygonIntersection3D

		std::cout << "\n========== Test: isLinePolygonIntersection3D ==========\n";
		Eigen::RowVector3d A(748674.6211, 2564712.4947, 64.5593);
		Eigen::RowVector3d B(748691.5686, 2564750.4571, 61.6250);
		Segment3D AB{ A ,B };
		std::vector<Eigen::RowVector3d> intersections1;
		bool isIntersection1 = isLinePolygonIntersection3D(polygon, AB, intersections1, true, false);
		std::cout << "\n-- Segment AB Intersection --\n";
		if (isIntersection1) {
			std::cout << "Intersection Points:\n"
					  << intersections1[0].format(Eigen::FullPrecision) << "\n";
		}
		
		std::cout << "\n--------------------------------------------------------\n";

		Eigen::RowVector3d C(748705.94909344427, 2564736.774573423, 65.5);
		Eigen::RowVector3d D(748697.44422830443, 2564736.8396572643, 65.5);
		Ray3D CD{ C ,D };
		std::vector<Eigen::RowVector3d> intersections2;
		bool isIntersection2 = isLinePolygonIntersection3D(polygon, CD, intersections2, true, false);
		std::cout << "\n-- Ray CD Intersection --\n";
		if (isIntersection2) {
			std::cout << "Intersection Points:\n";
			for (const auto& point : intersections2) {
				std::cout << point.format(Eigen::FullPrecision) << "\n";
			}
		}
		std::cout << "========== End of Test: isLinePolygonIntersection3D ==========\n";
	}

	{// test isLinePolyhedronIntersection

		std::cout << "\n========== Test: isLinePolyhedronIntersection ==========\n";
		{
			Eigen::RowVector3d A(748674.6211, 2564712.4947, 64.5593);
			Eigen::RowVector3d B(748708.8839, 2564753.8769, 59.6075);
			Segment3D AB{ A ,B };
			std::vector<Eigen::RowVector3d> intersections;
			bool isIntersection = isLinePolyhedronIntersection(polyhedronVertices, AB, intersections);
			std::cout << "\n-- Segment AB Intersection with Polyhedron --\n";
			if (isIntersection) {
				std::cout << "Intersection Points:\n";
				for (const auto& point : intersections) {
					std::cout << point.format(Eigen::FullPrecision) << "\n";
				}
			}
		}
		std::cout << "\n--------------------------------\n";
		{
			Eigen::RowVector3d A(748674.6211, 2564712.4947, 64.5593);
			Eigen::RowVector3d B(748680.9715683657, 2564720.1647197301, 63.641504336681471);
			Ray3D AB{ A ,B };
			std::vector<Eigen::RowVector3d> intersections;
			std::cout << "\n-- Ray AB Intersection with Polyhedron --\n";
			bool isIntersection = isLinePolyhedronIntersection(polyhedronVertices, AB, intersections);
			if (isIntersection) {
				std::cout << "Intersection Points:\n";
				for (const auto& point : intersections) {
					std::cout << point.format(Eigen::FullPrecision) << "\n";
				}
			}
		}
		std::cout << "========== End of Test: isLinePolyhedronIntersection ==========\n";

	}

	{// test filterPointsByPolyhedron

		std::cout << "\n========== Test: filterPointsByPolyhedron & filterPointsByPolygon==========\n";
		std::vector<Eigen::RowVector3d> testPoints3D;
		std::vector<Eigen::RowVector2d> testPoints2D;
		testPoints3D.reserve(200);
		testPoints2D.reserve(200);
		std::ifstream file("./testData/randomPoints.csv");
		std::string line;
		while (getline(file, line)) {
			std::istringstream lineStream(line);
			std::string c1, c2, c3;
			if (std::getline(lineStream, c1, ',') &&
				std::getline(lineStream, c2, ',') &&
				std::getline(lineStream, c3)) {

				Eigen::RowVector3d vec3d;
				vec3d << stod(c1), stod(c2), stod(c3);
				testPoints3D.emplace_back(vec3d);

				Eigen::RowVector2d vec2d;
				vec2d << stod(c1), stod(c2);
				testPoints2D.emplace_back(vec2d);

			}
		}
		std::cout << "Loaded " << testPoints3D.size() << " test points.\n";

		Eigen::Matrix<double, 8, 3> testPolyhedronVertices;
		testPolyhedronVertices << 20.785733063246152, 76.974749045767595, 20.082613534316067,
			79.748658377982537, 77.104254429716875, 20.082613534316067,
			81.118930088650117, 18.361256017841775, 20.082613534316067,
			22.15600477391374, 18.231750633892496, 20.082613534316067,
			20.785733063246152, 76.974749045767595, 80.113871877439578,
			79.748658377982537, 77.104254429716875, 80.113871877439578,
			81.118930088650117, 18.361256017841775, 80.113871877439578,
			22.15600477391374, 18.231750633892496, 80.113871877439578;

		std::cout << "\n[Polyhedron Filtering]";
		std::vector<Eigen::RowVector3d> insidePoints3D = filterPointsByPolyhedron(testPolyhedronVertices, testPoints3D, false, false);
		std::cout << "\nNumber of points inside the polyhedron: " << insidePoints3D.size() << "\n";

		std::vector<Eigen::RowVector3d> outsidePoints3D = filterPointsByPolyhedron(testPolyhedronVertices, testPoints3D,false, true);
		std::cout << "\nNumber of points outside the polyhedron: " << outsidePoints3D.size() << "\n";

		std::cout << "\n[Polygon Filtering]";
		std::vector<Eigen::RowVector2d> insidePoints2D = filterPointsByPolygon(testPolyhedronVertices.block<4, 2>(0, 0), testPoints2D, true, false, true);
		std::cout << "\nNumber of points outside the polygon: " << insidePoints2D.size() << "\n";

		std::vector<Eigen::RowVector2d> outsidePoints2D = filterPointsByPolygon(testPolyhedronVertices.block<4,2>(0,0), testPoints2D, true, true, true);
		std::cout << "\nNumber of points inside the polygon: " << outsidePoints2D.size() << "\n";


		std::cout << "========== End of Test: filterPointsByPolyhedron & filterPointsByPolygon ==========\n";

	}


	{// test computePCABox and computeMinBoundBox

		std::cout << "\n========== Test: Bounding Boxes & Rectangles ==========\n";
		Eigen::Matrix<double, 200, 3> testPointsMat;
		std::ifstream file("./testData/randomPoints.csv");
		std::string line;
		Eigen::Index row = 0;
		while (getline(file, line) && row < 200) {
			std::istringstream lineStream(line);
			std::string c1, c2, c3;
			if (std::getline(lineStream, c1, ',') &&
				std::getline(lineStream, c2, ',') &&
				std::getline(lineStream, c3)) {

				testPointsMat.row(row) << stod(c1), stod(c2), stod(c3);
				row++;
			}
		}
		{
			Eigen::Matrix<double, 8, 3> pcaBox = computePCAOBB3D(testPointsMat, true).first;
			Eigen::Matrix<double, 8, 3> mbBox = computeMinOBB3D(testPointsMat).first;

			std::cout << "\n--- PCA 3D Bounding Box ---\n";
			std::cout << pcaBox.format(Eigen::FullPrecision) << "\n";

			std::cout << "\n--- Minimum 3D Bounding Box ---\n";
			std::cout << mbBox.format(Eigen::FullPrecision) << "\n";

		}

		{
			Eigen::Matrix<double, 4, 2> pcaRect = computePCAOBB2D(testPointsMat.leftCols<2>(), true).first;
			Eigen::Matrix<double, 4, 2> mbRect = computeMinOBB2D(testPointsMat.leftCols<2>()).first;

			std::cout << "\n--- PCA 2D Bounding Rectangle (XY) ---\n";
			std::cout << pcaRect.format(Eigen::FullPrecision) << "\n";

			std::cout << "\n--- Minimum 2D Bounding Rectangle (XY) ---\n";
			std::cout << mbRect.format(Eigen::FullPrecision) << "\n";
		}
		std::cout << "========== End of Test: Bounding Boxes & Rectangles ==========\n";
	}

	{// test WorldToCameraImageCoords
		std::cout << "\n========== Test: WorldToCameraImageCoords ==========\n";
		Eigen::Matrix<double,8,3> testTarget;
		testTarget << 748793.66839939367, 2565000.6052587116, 48.909999999999997,
			748791.82177844818, 2565003.9931349289, 48.909999999999997,
			748795.20965466544, 2565005.8397558741, 48.909999999999997,
			748797.05627561093, 2565002.4518796569, 48.909999999999997,
			748793.66839939367, 2565000.6052587116, 50.909999999999997,
			748791.82177844818, 2565003.9931349289, 50.909999999999997,
			748795.20965466544, 2565005.8397558741, 50.909999999999997,
			748797.05627561093, 2565002.4518796569, 50.909999999999997;

		Eigen::RowVector3d cameraPosition(748780.10705481586, 2565024.6238759784, 57.51);
		Eigen::RowVector3d cameraDirPoint(748794.43902702956, 2565003.2225072929, 49.909999999999997);
		double FL = 32.0;
		std::cout << "Camera Position: " << cameraPosition.format(Eigen::FullPrecision) << "\n";
		std::cout << "Camera Look-at Point: " << cameraDirPoint.format(Eigen::FullPrecision) << "\n";
		std::cout << "Focal Length: " << FL << "\n";

		Eigen::Matrix<double, 8, 2> imagePoints = WorldToCameraImageCoords(testTarget, cameraDirPoint, cameraPosition, FL);
		std::cout << "\nProjected Image Coordinates:\n";
		std::cout << imagePoints.format(Eigen::FullPrecision) << "\n";
		std::cout << "========== End of Test: WorldToCameraImageCoords ==========\n";

	}


	return 0;
}