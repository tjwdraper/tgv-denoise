#ifndef _CONV_2D_HPP_
#define _CONV_2D_HPP_

#include "Field.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>

///////////////////////////////////////////////////////////////////////////////////////////////////
// Separable 2d convolution
///////////////////////////////////////////////////////////////////////////////////////////////////
class separable_conv2d {
    public:
        // Constructors and deconstructors
        separable_conv2d(const dim dimconv) : _dimconv(dimconv),
                                              _stepconv(1, dimconv.x),
                                              _sizeconv(dimconv.x*dimconv.y),
                                              _weights_i(new double[dimconv.x]),
                                              _weights_j(new double[dimconv.y]) {
            if (dimconv.x == 0 || dimconv.y == 0)
                throw std::runtime_error("In separable_conv2d(const dim), convolution dimensions must be non-zero.");

            if (dimconv.x % 2 == 0 || dimconv.y % 2 == 0)
                throw std::runtime_error("In separable_conv2d(const dim), convolution dimensions must be odd.");
        }
        virtual ~separable_conv2d() {
            delete[] _weights_i;
            delete[] _weights_j;
        }

        // Getters and setters
        dim get_dimensions() const {return _dimconv;}
        dim get_step() const {return _stepconv;}
        std::size_t get_size() const {return _sizeconv;}
        double get_weight_i(std::size_t i) const {
            if (i >= _dimconv.x)
                throw std::runtime_error("In double separable_conv2d::get_weigh_i(std::size_t) index out of bounds");
            return _weights_i[i];
        }
        double get_weight_j(std::size_t j) const {
            if (j >= _dimconv.x)
                throw std::runtime_error("In double separable_conv2d::get_weigh_j(std::size_t) index out of bounds");
            return _weights_j[j];
        }

        // Apply convolution
        template<typename T>
        void convolute(opticalflow::Field<T>& fin) const {
            const dim dimin = fin.get_dimensions();

            // Convolution centers
            const std::size_t center_i = _dimconv.x / 2;
            const std::size_t center_j = _dimconv.y / 2;

            // Apply convolution in y-direction
            opticalflow::Field<T> intermediate(dimin);

            for (std::size_t i = 0; i < dimin.x; ++i) {
                for (std::size_t j = 0; j < dimin.y; ++j) {
                    T val{};

                    for (std::size_t kj = 0; kj < _dimconv.y; ++kj) {
                        long jj = static_cast<long>(j) + static_cast<long>(kj) - static_cast<long>(center_j);

                        // Clip to boundary
                        if (jj < 0) {jj = 0;}
                        if (jj > dimin.y-1) {jj = dimin.y-1;}

                        val += fin.get_val(i, static_cast<std::size_t>(jj)) * _weights_j[kj];
                    }

                    intermediate.set_val(val, i, j);
                }
            }

            // Apply convolution in x-direction
            opticalflow::Field<T> result(dimin);
            
            for (std::size_t i = 0; i < dimin.x; ++i) {
                for (std::size_t j = 0; j < dimin.y; ++j) {
                    T val{};

                    for (std::size_t ki = 0; ki < _dimconv.x; ++ki) {
                        long ii = static_cast<long>(i) + static_cast<long>(ki) - static_cast<long>(center_i);

                        // Clip to boundary
                        if (ii < 0) {ii = 0;}
                        if (ii > dimin.x-1) {ii = dimin.x-1;}

                        val += intermediate.get_val(static_cast<std::size_t>(ii), j) * _weights_i[ki];
                    }

                    result.set_val(val, i, j);
                }
            }

            fin = std::move(result);

        }
        

    protected:
        void set_weight_i(double val, std::size_t i) {
            if (i >= _dimconv.x)
                throw std::runtime_error("In separable_conv2d::set_weight_i(double, std::size_t), index out of bounds");
            _weights_i[i] = val;
        }

        void set_weight_j(double val, std::size_t j) {
            if (j >= _dimconv.y)
                throw std::runtime_error("In separable_conv2d::set_weight_j(double, std::size_t), index out of bounds");
            _weights_j[j] = val;
        }

        void normalize_weights() {
            double sum_i(0.0);
            double sum_j(0.0);

            for (std::size_t i = 0; i < _dimconv.x; ++i)
                sum_i += _weights_i[i];
            for (std::size_t j = 0; j < _dimconv.y; ++j)
                sum_j += _weights_j[j];

            if (sum_i == 0.0 || sum_j == 0.0)
                throw std::runtime_error("In separable_conv2d::normalize_weights, sum of weights is zero.");

            for (std::size_t i = 0; i < _dimconv.x; ++i)
                _weights_i[i] /= sum_i;
            for (std::size_t j = 0; j < _dimconv.y; ++j)
                _weights_j[j] /= sum_j;
        }

        dim _dimconv;
        dim _stepconv;
        std::size_t _sizeconv;

        double* _weights_i;
        double* _weights_j;
};

class average_conv2d : public separable_conv2d {
    public:
        average_conv2d(const dim dimconv) : separable_conv2d(dimconv) {
            const double weight_i = 1.0 / static_cast<double>(dimconv.x);
            const double weight_j = 1.0 / static_cast<double>(dimconv.y);

            for (std::size_t i = 0; i < dimconv.x; ++i)
                set_weight_i(weight_i, i);
            for (std::size_t j = 0; j < dimconv.y; ++j)
                set_weight_j(weight_j, j);
        }
};

class gaussian_conv2d : public separable_conv2d {
    public:
        // Constructors and deconstructors
        gaussian_conv2d(const dim dimconv, const vector2d sigma) : separable_conv2d(dimconv) {
            if (sigma.x <= 0.0 || sigma.y <= 0.0)
                throw std::runtime_error("In gaussian_conv2d(const dim, const vector2d), sigma must be positive.");
            
            _sigma = sigma;

            const double center_i = static_cast<double>(dimconv.x-1)/2.0;
            const double center_j = static_cast<double>(dimconv.y-1)/2.0;

            for (std::size_t i = 0; i < dimconv.x; ++i) {
                const double x = static_cast<double>(i) - center_i;
                set_weight_i(std::exp(-x*x / (2*sigma.x*sigma.x)), i);
            }

            for (std::size_t j = 0; j < dimconv.y; ++j) {
                const double y = static_cast<double>(j) - center_j;
                set_weight_j(std::exp(-y*y / (2*sigma.y*sigma.y)), j);
            }

            normalize_weights();
        }
        gaussian_conv2d(const vector2d sigma) : gaussian_conv2d(dim(2.0*static_cast<std::size_t>(std::ceil(3.0 * sigma.x)) + 1,
                                                                    2.0*static_cast<std::size_t>(std::ceil(3.0 * sigma.y)) + 1),
                                                                    sigma) {}
        virtual ~gaussian_conv2d() {}

        //  Getters and setters
        vector2d get_sigma() const {return _sigma;}

    private:
        vector2d _sigma;
};

#endif
