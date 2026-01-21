#include<iostream>
#include <nlohmann/json.hpp>
#include "common/nfp_lib.hpp"
#include <fstream>

#if defined(_WIN32)
  #include <windows.h>
#elif defined(__APPLE__)
  #include <mach-o/dyld.h>
#else //Linux and Unix-like systems
  #include <unistd.h>
  #include <limits.h>
#endif


using json = nlohmann::json;

static bool read_all(std::istream& is, std::string& out) {
    std::ostringstream ss;
    ss << is.rdbuf();
    out = ss.str();
    return !out.empty() || !is.fail();
}

static json load_json_arg(const std::string& arg) {
    // 1) "-" => read JSON from stdin (large JSON support)
    if (arg == "-") {
        std::string text;
        if (!read_all(std::cin, text)) {
            throw std::runtime_error("Failed to read JSON from stdin");
        }
        return json::parse(text);
    }

    // 2) If it's an existing regular file => parse from file
    std::error_code ec;
    if (std::filesystem::is_regular_file(arg, ec)) {
        std::ifstream in(arg);
        if (!in) throw std::runtime_error("Failed to open input file: " + arg);
        json j;
        in >> j;
        return j;
    }

    // 3) Otherwise => treat as inline JSON text
    return json::parse(arg);
}
// ---------- Main ----------
int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr
            << "Usage:\n"
            << "  nfp_tool <input.json> <height> <width> [output.json]\n\n"
            << "Notes:\n"
            << "  - Assumes polygon entries are objects with VERTICES: [{x,y},...]\n"
            << "intended use with searchEdges=false.\n";
        return 1;
    }

    const std::string inputPath = argv[1];
    double height = 0.0;
    double width  = 0.0;

    try {
        height = std::stod(argv[2]);
        width  = std::stod(argv[3]);
    } catch (...) {
        std::cerr << "Invalid height/width numeric values.\n";
        return 1;
    }

    if (height <= 0.0 || width <= 0.0) {
        std::cerr << "Height and width must be > 0.\n";
        return 1;
    }

   json dataset;
    try {
        dataset = load_json_arg(argv[1]);
    } catch (const std::exception& e) {
        std::cerr << "Failed to load/parse JSON: " << e.what() << "\n";
        return 1;
    }    

    if (!dataset.is_object()) {
        std::cerr << "Expected top-level JSON object.\n";
        return 1;
    }

    try {
        dataset = NFPTool::processNFP(dataset, height, width);
    } catch (const std::exception& e) {
        std::cerr << "Error during NFP processing: " << e.what() << "\n";
        return 1;
    }

    // Output
    if (argc >= 5) {
        const std::string outputPath = argv[4];
        std::ofstream out(outputPath);
        if (!out) {
            std::cerr << "Failed to open output file: " << outputPath << "\n";
            return 1;
        }
        out << dataset.dump(2) << "\n";
    } 
    std::cout << dataset.dump(2) << "\n";

    return 0;
}