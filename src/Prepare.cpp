#include "Prepare.h"
#include <nlohmann/json.hpp>
using json = nlohmann::json;

//************************************************************************************************************************//
static json loadJsonFromFile(const fs::path& filePath) {
    std::ifstream file(filePath, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filePath.string());
    }
    json data;
    file >> data;
    return data;
}
//************************************************************************************************************************//

//************************************************************************************************************************//
static Eigen::MatrixX3d fromJson2Matrix3d(const json& j) {
    const size_t num_vertices = j.size();
    Eigen::MatrixX3d mat(num_vertices, 3);
    for (size_t i = 0; i < num_vertices; ++i) {
        const auto& vertex = j[i];
        mat(i, 0) = vertex[0].get<double>();
        mat(i, 1) = vertex[1].get<double>();
        mat(i, 2) = vertex[2].get<double>();
    }

    return mat;
}
//************************************************************************************************************************//

//************************************************************************************************************************//
static Eigen::RowVector3d fromJson2RowVector3d(const json& j) {
    return {
        j[0].get<double>(),
        j[1].get<double>(),
        j[2].get<double>()
    };
}
//************************************************************************************************************************//

//************************************************************************************************************************//
std::pair<std::unordered_map<size_t, SCData>, std::unordered_map<size_t, EFData>> getDataSet(const std::string& filepath) {
    json dataJson = loadJsonFromFile(filepath);

    std::unordered_map<size_t, SCData> SCSet;
    const auto& SCJson = dataJson["SCSet"];
    for (const auto& sc : SCJson) {
        SCData scData;
        size_t sc_id = sc["id"].get<size_t>();
        scData.polyhedronVertices = fromJson2Matrix3d(sc["polyhedronVertices"]);
        scData.obbBoxVertices = fromJson2Matrix3d(sc["obbBoxVertices"]);
        scData.obbBoxShrink1Vertices = fromJson2Matrix3d(sc["obbBoxShrink1Vertices"]);
        scData.obbBoxCenter = fromJson2RowVector3d(sc["obbBoxCenter"]);
        SCSet.emplace(sc_id, std::move(scData));
    }

    std::unordered_map<size_t, EFData> EFSet;
    const auto& EFJson = dataJson["EFSet"];
    for (const auto& ef : EFJson) {
        EFData efData;
        size_t ef_id = ef["id"].get<size_t>();
        efData.polyhedronVertices = fromJson2Matrix3d(ef["polyhedronVertices"]);
        efData.obbBoxVertices = fromJson2Matrix3d(ef["obbBoxVertices"]);
        efData.obbBoxExpand1Vertices = fromJson2Matrix3d(ef["obbBoxExpand1Vertices"]);
        efData.obbBoxCenter = fromJson2RowVector3d(ef["obbBoxCenter"]);
        EFSet.emplace(ef_id, std::move(efData));
    }

    return { SCSet, EFSet };
}
//************************************************************************************************************************//

