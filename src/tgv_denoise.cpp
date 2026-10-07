#include <fstream>
#include <string>
#include <stdexcept>
#include <chrono>
#include <vector>
#include <map>

#include "coord2d.hpp"
#include "json.hpp"
#include "Field.hpp"
#include "denoise.hpp"

///////////////////////////////////////////////////////////////////////////////////////////////////
// Create methods to convert string to model options
///////////////////////////////////////////////////////////////////////////////////////////////////
inline const std::map<std::string, VerboseOption> mapper_verbose_option {
    {"silent", VerboseOption::SILENT},
    {"verbose", VerboseOption::VERBOSE}
};

inline const std::map<std::string, ModelOption> mapper_model_option {
    {"TV", ModelOption::TV},
    {"TGV", ModelOption::TGV}
};

///////////////////////////////////////////////////////////////////////////////////////////////////
// Read json configuration file and store in structure
///////////////////////////////////////////////////////////////////////////////////////////////////
struct json_config {
    // Path to files
    std::string path_image;
    std::string path_output = "image_in/output.tiff";

    ModelOption option = ModelOption::TGV;

    // Denoising parameters
    double alpha0 = 1.0;
    double alpha1 = 2.0;
    double tau = 0.25;
    double sigma = 0.25;
    double lambda = 0.01;
    int niter = 1000;
    double convergence = 1e-5;

    // Verbose
    VerboseOption verbose = VerboseOption::SILENT;
};

json_config load_config(const std::string& filename) {
    // Open .json configuration file
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Could not open json configuration file: " + filename);
    nlohmann::json json;
    file >> json;

    // Create json configuration structure
    json_config config;

    // Set the image path
    config.path_image = json.at("path_image").get<std::string>();

    // Set output path
    if (json.contains("path_output"))
        config.path_output = json.at("path_output").get<std::string>();

    // Read model option
    if (json.contains("model_option")) {
        std::string str_option = json.at("model_option").get<std::string>();
        auto it = mapper_model_option.find(str_option);
        if (it != mapper_model_option.end())
            config.option = it->second;
    }

    // Set denoising parameters
    if (json.contains("parameters")) {
        const auto& parameters = json.at("parameters");

        if (parameters.contains("alpha0"))
            config.alpha0 = parameters.at("alpha0").get<double>();

        if (parameters.contains("alpha1"))
            config.alpha1 = parameters.at("alpha1").get<double>();

        if (parameters.contains("tau"))
            config.tau = parameters.at("tau").get<double>();

        if (parameters.contains("sigma"))
            config.sigma = parameters.at("sigma").get<double>();

        if (parameters.contains("lambda"))
            config.lambda = parameters.at("lambda").get<double>();

        if (parameters.contains("niter"))
            config.niter = parameters.at("niter").get<int>(); 

        if (parameters.contains("convergence"))
            config.convergence = parameters.at("convergence").get<double>();

        if (parameters.contains("verbose")) {
            std::string str_option = parameters.at("verbose").get<std::string>();
            auto it = mapper_verbose_option.find(str_option);
            if (it != mapper_verbose_option.end())
                config.verbose = it->second;
        }
    }

    // Done
    return config;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Main function
///////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {
    // Check input arguments
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " config.json\n";
        return 1;
    }

    // Load configuration
    json_config config = load_config(argv[1]);

    // Load image
    std::ifstream file(config.path_image.c_str(), std::ios::binary);
    if (!file)
        throw std::runtime_error("Could not open image file.");

    // Read the image dimensions
    double rows_d, cols_d;
    file.read(reinterpret_cast<char*>(&rows_d), sizeof(double));
    file.read(reinterpret_cast<char*>(&cols_d), sizeof(double));
    if (!file)
        throw std::runtime_error("Could not read image dimensions");

    const dim dimin(
        static_cast<std::size_t>(rows_d),
        static_cast<std::size_t>(cols_d)
    );

    // Read the image contents
    opticalflow::Image image_in(dimin);
    file.read(reinterpret_cast<char*>(image_in.get_field()), image_in.get_size() * sizeof(double));

    if (!file)
        throw std::runtime_error("Could not read image file.");

    // Denoise
    opticalflow::Image image_out(dimin);
    
    const auto start = std::chrono::steady_clock::now();
    if (config.option == ModelOption::TGV) {
        image_out = denoise::tgv_denoise(
            image_in,
            config.tau,
            config.sigma,
            config.lambda,
            config.alpha0,
            config.alpha1,
            config.niter,
            config.convergence,
            config.verbose
        );
    }
    else if (config.option == ModelOption::TV) {
        image_out = denoise::tv_denoise(
            image_in,
            config.tau,
            config.sigma,
            config.lambda,
            config.alpha0,
            config.niter,
            config.convergence,
            config.verbose
        );
    }
    const auto end = std::chrono::steady_clock::now();

    if (config.verbose == VerboseOption::VERBOSE) {
        const std::chrono::duration<double> elapsed = end-start;
        std::cout << "Runtime (s): " << elapsed.count() << "\n";
    }
        

    // Write to output file
    std::ofstream output(config.path_output.c_str(), std::ios::binary);
    if (!output)
        throw std::runtime_error("Could not open output file.");

    output.write(reinterpret_cast<const char*>(&rows_d), sizeof(double));
    output.write(reinterpret_cast<const char*>(&cols_d), sizeof(double));

    output.write(reinterpret_cast<const char*>(image_out.get_field()), image_out.get_size()*sizeof(double));
    if (!output)
        throw std::runtime_error("Could not write output file.");

    // Done
    return 0;
}