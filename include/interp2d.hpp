#ifndef _INTERP_2D_HPP_
#define _INTERP_2D_HPP_

#include "Field.hpp"

#include <cmath>
#include <stdexcept>

namespace interp2d {
    // Template for interpolation
    template <typename T>
    T clamp_at_boundary(const opticalflow::Field<T>& field, long posx, long posy, const dim dim_field) {
        if (posx < 0) {posx = 0;}
        if (posy < 0) {posy = 0;}
        if (posx > dim_field.x-1) {posx = dim_field.x-1;}
        if (posy > dim_field.y-1) {posy = dim_field.y-1;}
        return field.get_val(posx, posy);
    }

    template <typename T>
    T interpolate_at_value(const opticalflow::Field<T>& field, const vector2d pos, const dim dim_field) {
        const double px = pos.x; const long dx = static_cast<long>(std::floor(px)); const double fx = px - static_cast<double>(dx);
        const double py = pos.y; const long dy = static_cast<long>(std::floor(py)); const double fy = py - static_cast<double>(dy);

        T val{};

        val += clamp_at_boundary(field, dx, dy, dim_field)*(1.0-fx)*(1.0-fy);
        val += clamp_at_boundary(field, dx+1, dy, dim_field)*fx*(1.0-fy);
        val += clamp_at_boundary(field, dx, dy+1, dim_field)*(1.0-fx)*fy;
        val += clamp_at_boundary(field, dx+1, dy+1, dim_field)*fx*fy;

        return val;
    }

    template <typename T>
    void interp2d(opticalflow::Field<T>& field_out, const opticalflow::Field<T>& field_in, const opticalflow::Motion& motion) {
        const dim dim_out = field_out.get_dimensions();
        const dim dim_in = field_in.get_dimensions();

        for (std::size_t i = 0; i < dim_out.x; ++i) {
            for (std::size_t j = 0; j < dim_out.y; ++j) {
                const vector2d p = vector2d(static_cast<double>(i), static_cast<double>(j)) + motion.get_val(i,j);

                T val = interpolate_at_value<T>(field_in, p, dim_in);

                field_out.set_val(val, i, j);
            }
        }
    }


    // Warp image with deformation vector field
    inline void warp2d(opticalflow::Image& image_out, const opticalflow::Image& image_in, const opticalflow::Motion& motion) {
        interp2d::interp2d<double>(image_out, image_in, motion);
    }

    // Accumulate deformation vector fields, i.e. calculate the result of 1 + u_acc <- (1 + u) \circ (1 + u_interp)
    inline void accumulate(opticalflow::Motion& motion_acc, const opticalflow::Motion& motion, const opticalflow::Motion& motion_interp) {
        interp2d::interp2d<vector2d>(motion_acc, motion, motion_interp);
        motion_acc += motion_interp;
    }
    inline void accumulate(opticalflow::Motion& motion, const opticalflow::Motion& motion_interp) {
        opticalflow::Motion motion_tmp(motion.get_dimensions());
        interp2d::accumulate(motion_tmp, motion, motion_interp);
        motion = std::move(motion_tmp);
    }

    // Resize
    inline void resize(opticalflow::Image& image_out, const opticalflow::Image& image_in) {
        const dim dim_out = image_out.get_dimensions();
        const dim dim_in = image_in.get_dimensions();

        const vector2d factor(
            static_cast<double>(dim_in.x) / static_cast<double>(dim_out.x),
            static_cast<double>(dim_in.y) / static_cast<double>(dim_out.y)
        );

        for (std::size_t i = 0; i < dim_out.x; ++i) {
            for (std::size_t j = 0; j < dim_out.y; ++j) {
                vector2d pos(
                    (static_cast<double>(i) + 0.5) * factor.x - 0.5,
                    (static_cast<double>(j) + 0.5) * factor.y - 0.5
                );

                image_out.set_val(
                    interpolate_at_value<double>(image_in, pos, dim_in),
                    i, j
                );
            }
        }
    }
    inline void resize(opticalflow::Motion& motion_out, const opticalflow::Motion& motion_in) {
        const dim dim_out = motion_out.get_dimensions();
        const dim dim_in = motion_in.get_dimensions();

        const vector2d factor(
            static_cast<double>(dim_in.x) / static_cast<double>(dim_out.x),
            static_cast<double>(dim_in.y) / static_cast<double>(dim_out.y)
        );

        for (std::size_t i = 0; i < dim_out.x; ++i) {
            for (std::size_t j = 0; j < dim_out.y; ++j) {
                vector2d pos(
                    (static_cast<double>(i) + 0.5) * factor.x - 0.5,
                    (static_cast<double>(j) + 0.5) * factor.y - 0.5
                );

                vector2d val = interpolate_at_value<vector2d>(motion_in, pos, dim_in);

                motion_out.set_val(vector2d(val.x / factor.x, val.y / factor.y), i, j);
            }
        }
    }
}

#endif