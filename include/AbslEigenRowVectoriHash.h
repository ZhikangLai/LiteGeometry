#ifndef ABSLEIGENROWVECTORIHASH_H
#define ABSLEIGENROWVECTORIHASH_H
#include "Utility.h"

#define XXH_INLINE_ALL
#include <xxhash.hpp>
#include <absl/container/flat_hash_set.h>

struct EigenRowVectoriEqual {
	bool operator()(const Eigen::RowVector2i& v1, const Eigen::RowVector2i& v2) const noexcept {
		return v1 == v2;
	}
	bool operator()(const Eigen::RowVector3i& v1, const Eigen::RowVector3i& v2) const noexcept {
		return v1 == v2;
	}
};

static constexpr uint64_t kSeed = 0x9E3779B97F4A7C15ULL;
struct EigenRowVectoriHash {
	size_t operator()(const Eigen::RowVector2i& v) const noexcept {
		return static_cast<size_t>(xxh::xxhash<64>(v.data(), sizeof(int) * 2, kSeed));
	}

	size_t operator()(const Eigen::RowVector3i& v) const {
		return static_cast<size_t>(xxh::xxhash<64>(v.data(), sizeof(int) * 3, kSeed));
	}
};


using RowVector2iSet = absl::flat_hash_set<Eigen::RowVector2i, EigenRowVectoriHash, EigenRowVectoriEqual>;
using RowVector3iSet = absl::flat_hash_set<Eigen::RowVector3i, EigenRowVectoriHash, EigenRowVectoriEqual>;

#endif