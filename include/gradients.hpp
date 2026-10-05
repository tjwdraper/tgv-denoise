#ifndef _GRADIENTS_H_
#define _GRADIENTS_H_

#include "coord2d.hpp"
#include "Field.hpp"

namespace gradients {
    // First order partial derivatives
    template <typename T>
    inline T partial_x(const opticalflow::Field<T>& field, const std::size_t i, const std::size_t j) {
        const dim dimin = field.get_dimensions();
        if (dimin.x < 2)
            throw std::runtime_error("In T gradients::partial_x(const Field<T>&, const std::size_t, const std::size_t), x-dimension must be at least 2.");

        if (i == 0)
            return field.get_val(i+1,j) - field.get_val(i,j); 
        else if (i == dimin.x-1)
            return field.get_val(i,j) - field.get_val(i-1, j);
        else
            return (field.get_val(i+1,j) - field.get_val(i-1,j)) / 2.0; 
    }

    template <typename T>
    inline T partial_y(const opticalflow::Field<T>& field, const std::size_t i, const std::size_t j) {
        const dim dimin = field.get_dimensions();
        if (dimin.y < 2)
            throw std::runtime_error("In T gradients::partial_y(const Field<T>&, const std::size_t, const std::size_t), y-dimension must be at least 2.");

        if (j == 0)
            return field.get_val(i,j+1) - field.get_val(i,j); 
        else if (j == dimin.y-1)
            return field.get_val(i,j) - field.get_val(i, j-1);
        else
            return (field.get_val(i,j+1) - field.get_val(i,j-1)) / 2.0; 
    }

    template <typename T>
    inline T partial_xx(const opticalflow::Field<T>& field, const std::size_t i, const std::size_t j) {
        const dim dimin = field.get_dimensions();
        if (dimin.x < 4)
            throw std::runtime_error("In T gradients::partial_xx(const Field<T>&, const std::size_t, const std::size_t), x-dimension must be at least 4");

        if (i == 0)
            return 2.0*field.get_val(i,j) - 5.0*field.get_val(i+1,j) + 4.0*field.get_val(i+2,j) - field.get_val(i+3,j);
        else if (i == dimin.x-1)
            return 2.0*field.get_val(i,j) - 5.0*field.get_val(i-1,j) + 4.0*field.get_val(i-2,j) - field.get_val(i-3,j);
        else
            return field.get_val(i+1,j) - 2.0*field.get_val(i,j) + field.get_val(i-1,j);
    }

    template <typename T>
    inline T partial_yy(const opticalflow::Field<T>& field, const std::size_t i, const std::size_t j) {
        const dim dimin = field.get_dimensions();
        if (dimin.y < 4)
            throw std::runtime_error("In T gradients::partial_yy(const Field<T>&, const std::size_t, const std::size_t), y-dimension must be at least 4");

        if (j == 0)
            return 2.0*field.get_val(i,j) - 5.0*field.get_val(i,j+1) + 4.0*field.get_val(i,j+2) - field.get_val(i,j+3);
        else if (j == dimin.y-1)
            return 2.0*field.get_val(i,j) - 5.0*field.get_val(i,j-1) + 4.0*field.get_val(i,j-2) - field.get_val(i,j-3);
        else
            return field.get_val(i,j+1) - 2.0*field.get_val(i,j) + field.get_val(i,j-1);
    }

    template <typename T>
    inline T partial_xy(const opticalflow::Field<T>& field, const std::size_t i, const std::size_t j) {
        const dim dimin = field.get_dimensions();
        if (dimin.x < 1 || dimin.y < 1)
            throw std::runtime_error("In T gradients::partial_xy(const Field<T>&, const std::size_t, const std::size_t), x and y dimensions must be at least 1");

        if (i == 0 || j == 0 || i == dimin.x-1 || j == dimin.y-1)
            return T{};

        return (field.get_val(i+1,j+1) - field.get_val(i+1,j-1) - field.get_val(i-1,j+1) + field.get_val(i-1,j-1)) / 4.0;
    }

    template <typename T>
    void horn_schunck_average(opticalflow::Field<T>& field_avg, const opticalflow::Field<T>& field) {
        // Check that input dimensions are OK
        if (field_avg.get_dimensions() != field.get_dimensions())
            throw std::runtime_error("Error in gradients::horn_schunck_average(opticalflow::Field<T>&, const opticalflow::Field<T>&): input dimensions have to be the same as target");

        // Get the dimensions of the image
        const dim dimin = field.get_dimensions();

        // Iterate over interior points
        for (std::size_t i = 1; i < dimin.x-1; ++i) {
            for (std::size_t j = 1; j < dimin.y-1; ++j) {
                field_avg.set_val(
                    0.25 * (field.get_val(i+1,j) 
                            + field.get_val(i-1,j) 
                            + field.get_val(i,j+1) 
                            + field.get_val(i,j-1)),
                    i,j
                );
            }
        }

        // Iterate over edges
        for (std::size_t j = 1; j < dimin.y-1; ++j) {
            field_avg.set_val(
                0.25 * (field.get_val(1,j) 
                        + field.get_val(0,j) 
                        + field.get_val(0,j+1) 
                        + field.get_val(0,j-1)),
                0, j
            );

            field_avg.set_val(
                0.25 * (field.get_val(dimin.x-2, j) 
                        + field.get_val(dimin.x-1, j) 
                        + field.get_val(dimin.x-1, j+1) 
                        + field.get_val(dimin.x-1, j-1)),
                dimin.x-1, j
            );
        }

        for (std::size_t i = 1; i < dimin.x-1; ++i) {
            field_avg.set_val(
                0.25 * (field.get_val(i,1) 
                        + field.get_val(i,0) 
                        + field.get_val(i+1,0) 
                        + field.get_val(i-1,0)),
                i,0
            );

            field_avg.set_val(
                0.25 * (field.get_val(i,dimin.y-2) 
                        + field.get_val(i,dimin.y-1) 
                        + field.get_val(i+1,dimin.y-1) 
                        + field.get_val(i-1,dimin.y-1)),
                i, dimin.y-1
            );
        }

        // Corner points
        field_avg.set_val(0.25 * (field.get_val(1,0) 
                                + field.get_val(0,0) 
                                + field.get_val(0,1) 
                                + field.get_val(0,0)),
                                0,0);
        field_avg.set_val(0.25 * (field.get_val(dimin.x-2,0) 
                                + field.get_val(dimin.x-1,0) 
                                + field.get_val(dimin.x-1,1) 
                                + field.get_val(dimin.x-1,0)),
                                dimin.x-1,0);
        field_avg.set_val(0.25 * (field.get_val(0,dimin.y-2) 
                                + field.get_val(0,dimin.y-1) 
                                + field.get_val(1,dimin.y-1) 
                                + field.get_val(0,dimin.y-1)),
                                0,dimin.y-1);
        field_avg.set_val(0.25 * (field.get_val(dimin.x-2,dimin.y-1) 
                                + field.get_val(dimin.x-1,dimin.y-1) 
                                + field.get_val(dimin.x-1,dimin.y-2) 
                                + field.get_val(dimin.x-1, dimin.y-1)),
                                dimin.x-1,dimin.y-1);
    }
    
    void gradient(opticalflow::Motion& dI, const opticalflow::Image& image) {
        // Check that input dimensions are OK
        if (image.get_dimensions() != dI.get_dimensions())
            throw std::runtime_error("Error in gradients::gradient(opticalflow::Motion&, const opticalflow::Image&, const opticalflow::Image&): input dimensions have to be the same as target");

        // Get the dimensions of the image
        const dim dimin = image.get_dimensions();

        for (std::size_t i = 0; i < dimin.x; i++) {
            for (std::size_t j = 0; j < dimin.y; j++) {
                dI.set_val(
                    vector2d(gradients::partial_x(image, i, j),
                             gradients::partial_y(image, i, j)),
                    i,j
                );
            }
        }
    }

    void jacobian(opticalflow::Image& image, const opticalflow::Motion& motion) {
        // Check that input dimensions are OK
        if (image.get_dimensions() != motion.get_dimensions())
            throw std::runtime_error("Error in gradients::jacobian(Image&, const opticalflow::Motion&): input dimensions have to be the same as target");

        // Get the dimensions of the image
        const dim dimin = image.get_dimensions();

        // Store the input in the object
        vector2d dudx, dudy;
        for (std::size_t i = 0; i < dimin.x; i++) {
            for (std::size_t j = 0; j < dimin.y; j++) {
                dudx = gradients::partial_x<vector2d>(motion, i, j);
                dudy = gradients::partial_y<vector2d>(motion, i, j);

                image.set_val(
                    (1.0 + dudx.x) * (1.0 + dudy.y) - dudx.y * dudy.x,
                    i, j
                );
            }
        }

        // Done
        return;

    }
}

#endif