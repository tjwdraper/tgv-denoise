// #define cimg_display 0
// #include "CImg.h"

#include <fstream>
#include <string>
#include <stdexcept>
#include <chrono>
#include <vector>

#include "coord2d.hpp"
#include "json.hpp"
#include "Field.hpp"
#include "tv_denoise_2d.hpp"
#include "tgv_denoise_2d.hpp"

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
    int niter = 100;
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
    }

    // Done
    return config;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Read images using the cimg_library and convert to opticalflow::Image type from Field.hpp,
// which is the input for this ImageRegistration class implementation
///////////////////////////////////////////////////////////////////////////////////////////////////
// void convert_cimg_to_opticalflow(opticalflow::Image& image, cimg_library::CImg<double>& cimage) {
//     const dim dimin(cimage.width(), cimage.height());
//     const std::size_t size = dimin.x * dimin.y;

//     if (image.get_dimensions() != dimin)
//         throw std::runtime_error("In convert_cimg_to_opticalflow(Image&, cimg_library::CImg<double>&), dimensions of input and target have to equal");

//     // Convert cimage to grayscale, raw image. Average over 3 color channels
//     double* image_gs = new double[size];
//     double* cimage_rgb = cimage.data();

//     if (cimage.spectrum() == 1) {
//         for (std::size_t idx = 0; idx < size; ++idx) {
//             image_gs[idx] = cimage_rgb[idx];
//         }
//     }
//     else if (cimage.spectrum() == 3) {
//         for (std::size_t idx = 0; idx < size; ++idx) {
//             image_gs[idx] = (cimage_rgb[idx] + cimage_rgb[idx + size] + cimage_rgb[idx + 2*size]) / 3.0;
//         }
//     }
//     else {
//         throw std::runtime_error("In convert_cimg_to_opticalflow(Image&, cimg_library::CImg<double>&), number of channgels should be either 1 or 3.");
//     }

//     // Set data from opticalflow::Image target to raw data values
//     opticalflow::image::load_image(image_gs, image);

//     // Normalize intensities between zero and one
//     opticalflow::image::normalize(image);

//     // Free memory
//     delete[] image_gs;
// }

///////////////////////////////////////////////////////////////////////////////////////////////////
// Write images from opticalflow::Image to cimg_library for saving 
///////////////////////////////////////////////////////////////////////////////////////////////////
// void convert_opticalflow_to_cimg(cimg_library::CImg<double>& cimage, const opticalflow::Image& image) {
//     opticalflow::image::save_image(cimage.data(), image);
//     cimage *= 255;
// }

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

    opticalflow::Image image_in(dim(256, 256));
    file.read(reinterpret_cast<char*>(image_in.get_field()), image_in.get_size() * sizeof(double));

    if (!file)
        throw std::runtime_error("Could not read image file.");










    // std::cout << "Loading image...";
    // cimg_library::CImg<double> cimage_in(config.path_image.c_str());
    // std::cout << "Image loaded: " << cimage_in.width() << "x" << cimage_in.height() << std::endl;

    // Convert to opticalflow type:
    // std::cout << "Convert to opticalflow type...";
    // const dim dimin(256, 256);

    // opticalflow::Image image_in(dimin);
    // convert_cimg_to_opticalflow(image_in, cimage_in.data());

    // std::cout << "Complete!" << std::endl;

    // Denoise
    opticalflow::Image image_out(dim(256, 256));
    if (config.option == ModelOption::TGV) {
        image_out = tgv_denoise::denoise(
            image_in,
            config.tau,
            config.sigma,
            config.lambda,
            config.alpha0,
            config.alpha1,
            config.niter
        );
    }
    else if (config.option == ModelOption::TV) {
        image_out = tv_denoise::denoise(
            image_in,
            config.tau,
            config.sigma,
            config.lambda,
            config.alpha0,
            config.niter
        );
    }

    // Write to output file
    std::ofstream output(config.path_output.c_str(), std::ios::binary);
    if (!output)
        throw std::runtime_error("Could not open output file.");

    output.write(reinterpret_cast<const char*>(image_out.get_field()), image_out.get_size()*sizeof(double));
    if (!output)
        throw std::runtime_error("Could not write output file.");


    // std::cout << "Converting to CImg type...";
    // cimg_library::CImg<double> cimage_out(cimage_in.width(), cimage_in.height(), 1); // 1 channel, gray-scale output
    // convert_opticalflow_to_cimg(cimage_out, image_out);
    // std::cout << "Complete!" << std::endl;

    // std::cout << "Writing image...";
    // cimage_out.save(config.path_output.c_str());
    // std::cout << "Complete!" << std::endl;

    // Done
    return 0;
}