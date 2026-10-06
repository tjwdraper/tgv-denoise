#ifndef _TV_DENOISE_2D_HPP_
#define _TV_DENOISE_2D_HPP_

#include "coord2d.hpp"
//#include "gradients.hpp"
#include "Field.hpp"
#include "ConfigurationOptions.hpp"

namespace tv_denoise {
    // Update operators
    void update_p(opticalflow::Image& px, opticalflow::Image& py,
                  const opticalflow::Image& ubar, 
                  double sigma) {
        const dim dimin = px.get_dimensions();

        for (std::size_t i = 0; i < dimin.x; ++i) {
            for (std::size_t j = 0; j < dimin.y; ++j) {
                px.set_val(px.get_val(i,j) + sigma * (opticalflow::gradients::partial_x_forward(ubar, i, j)),i,j);
                py.set_val(py.get_val(i,j) + sigma * (opticalflow::gradients::partial_y_forward(ubar, i, j)),i,j);
            }
        }
    }

    void update_u(opticalflow::Image& u,
                  const opticalflow::Image& px, const opticalflow::Image& py,
                  double tau) {
        const dim dimin = u.get_dimensions();        

        for (std::size_t i = 0; i < dimin.x; ++i) {
            for (std::size_t j = 0; j < dimin.y; ++j) {
                const double divp = opticalflow::gradients::partial_x_backwards(px,i,j) + opticalflow::gradients::partial_y_backwards(py,i,j);
                u.set_val(u.get_val(i,j) + tau*divp,i,j);
            }
        }

    }

    // Projection operators
    void proj_p(opticalflow::Image& px, opticalflow::Image& py, double alpha) {
        const dim dimin = px.get_dimensions();

        for (std::size_t i = 0; i < dimin.x; ++i) {
            for (std::size_t j = 0; j < dimin.y; ++j) {
                const double pxv = px.get_val(i,j);
                const double pyv = py.get_val(i,j);

                const double normp = std::sqrt(pxv*pxv + pyv*pyv) / alpha;

                if (normp > 1.0) {
                    px.set_val(px.get_val(i,j)/normp, i,j);
                    py.set_val(py.get_val(i,j)/normp, i,j);
                }
            }
        }
    }
    void prox(opticalflow::Image& u, const opticalflow::Image& f, double tau, double lambda) {
        const dim dimin = u.get_dimensions();

        for (std::size_t i = 0; i < dimin.x; ++i) {
            for (std::size_t j = 0; j < dimin.y; ++j) {
                u.set_val( (lambda*u.get_val(i,j) + tau * f.get_val(i,j)) / (lambda + tau), i, j);
            }
        }
    }

    // Primal-dual method TV-denoising algorithm
    opticalflow::Image denoise(const opticalflow::Image& f, double tau, double sigma, double lambda, double alpha0, int niter) {
        // Get the image dimensions
        const dim dimin = f.get_dimensions();

        // Initialize primal-dual variables
        opticalflow::Image u(f);
        opticalflow::Image ubar(f);

        opticalflow::Image px(dimin);
        opticalflow::Image py(dimin);

        // Tracking variable
        opticalflow::Image uold(dimin);

        // Primal-dual iterations
        for (int iter = 0; iter < niter; ++iter) {

            tv_denoise::update_p(px,py,ubar,sigma);
            tv_denoise::proj_p(px,py,alpha0);

            uold = u;

            tv_denoise::update_u(u,px,py,tau);
            tv_denoise::prox(u,f,tau,lambda);

            // Update ubar
            ubar = 2*u - uold;

            if (opticalflow::image::norm(u-uold)/opticalflow::image::norm(u) < 1e-6)
                break;

            // Update some norms:
            if (iter % 50 == 0) {
                std::cout << "iter " << iter
                    << " |u| = " << opticalflow::image::norm(u)
                    << " |ubar| = " << opticalflow::image::norm(ubar)
                    << std::endl;
            }
        }

        return u;
    }
}

#endif