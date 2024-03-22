#include "MFGVRP_Solver.h"
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
int main() {
	std::string path = "instances/";

	for (const auto& entry : fs::directory_iterator(path)) {
		MFGVRP_Solver solver = MFGVRP_Solver();
		solver.LoadData(entry.path().string(), true);
	}
	return 0;
}
