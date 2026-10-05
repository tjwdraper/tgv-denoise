#ifndef _TGV_DENOISE_2D_HPP_
#define _TGV_DENOISE_2D_HPP_

#include "coord2d.hpp"
#include "gradients.hpp"
#include "Field.hpp"
#include "ConfigurationOptions.hpp"

namespace tgv_denoise {
    // Update operators
    void update_p(opticalflow::Image& px, opticalflow::Image& py,
                  const opticalflow::Image& ubar, 
                  const opticalflow::Image& vbarx, const opticalflow::Image& vbary,
                  double sigma) {
        const dim dimin = px.get_dimensions();

        for (std::size_t i = 0; i < dimin.x; ++i) {
            for (std::size_t j = 0; j < dimin.y; ++j) {
                px.set_val(px.get_val(i,j) + sigma * (gradients::partial_x_forward(ubar, i, j) - vbarx.get_val(i,j)),i,j);
                py.set_val(py.get_val(i,j) + sigma * (gradients::partial_y_forward(ubar, i, j) - vbary.get_val(i,j)),i,j);
            }
        }
    }

    void update_q(opticalflow::Image& qxx, opticalflow::Image& qyy, opticalflow::Image& qxy,
                  const opticalflow::Image& vbarx, const opticalflow::Image& vbary,
                  double sigma) {
        const dim dimin = qxx.get_dimensions();

        for (std::size_t i = 0; i < dimin.x; ++i) {
            for (std::size_t j = 0; j < dimin.y; ++j) {
                const double Evbarxx = gradients::partial_x_forward(vbarx, i, j);
                const double Evbaryy = gradients::partial_y_forward(vbary, i, j);
                const double Evbarxy = 0.5 * (gradients::partial_x_forward(vbary,i,j) + gradients::partial_y_forward(vbarx,i,j));

                qxx.set_val(qxx.get_val(i,j) + sigma * Evbarxx,i,j);
                qyy.set_val(qyy.get_val(i,j) + sigma * Evbaryy,i,j);
                qxy.set_val(qxy.get_val(i,j) + sigma * Evbarxy,i,j);
            }
        }
    }

    void update_u(opticalflow::Image& u,
                  const opticalflow::Image& px, const opticalflow::Image& py,
                  double tau) {
        const dim dimin = u.get_dimensions();        

        for (std::size_t i = 0; i < dimin.x; ++i) {
            for (std::size_t j = 0; j < dimin.y; ++j) {
                const double divp = gradients::partial_x_backwards(px,i,j) + gradients::partial_y_backwards(py,i,j);
                u.set_val(u.get_val(i,j) + tau*divp,i,j);
            }
        }

    }

    void update_v(opticalflow::Image& vx, opticalflow::Image& vy,
                  const opticalflow::Image& px, const opticalflow::Image& py,
                  const opticalflow::Image& qxx, const opticalflow::Image& qyy, const opticalflow::Image& qxy,
                  double tau) {
        const dim dimin = vx.get_dimensions();

        for (std::size_t i = 0; i < dimin.x; ++i) {
            for (std::size_t j = 0; j < dimin.y; ++j) {
                const double divqx = gradients::partial_x_backwards(qxx,i,j) + gradients::partial_y_backwards(qxy,i,j);
                const double divqy = gradients::partial_y_backwards(qyy,i,j) + gradients::partial_x_backwards(qxy,i,j);

                vx.set_val(vx.get_val(i,j) + tau * (px.get_val(i,j) + divqx), i,j);
                vy.set_val(vy.get_val(i,j) + tau * (py.get_val(i,j) + divqy), i,j);
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

    void proj_q(opticalflow::Image& qxx, opticalflow::Image& qyy, opticalflow::Image& qxy, double alpha) {
        const dim dimin = qxx.get_dimensions();

        for (std::size_t i = 0; i < dimin.x; ++i) {
            for (std::size_t j = 0; j < dimin.y; ++j) {
                const double qxxv = qxx.get_val(i,j);
                const double qyyv = qyy.get_val(i,j);
                const double qxyv = qxy.get_val(i,j);

                const double normq = std::sqrt(qxxv*qxxv + qyyv*qyyv + 2*qxyv*qxyv) / alpha;

                if (normq > 1.0) {
                    qxx.set_val(qxx.get_val(i,j)/normq, i,j);
                    qyy.set_val(qyy.get_val(i,j)/normq, i,j);
                    qxy.set_val(qxy.get_val(i,j)/normq, i,j);
                }
            }
        }
    }

    // Proximal operator
    void prox(opticalflow::Image& u, const opticalflow::Image& f, double tau, double lambda) {
        const dim dimin = u.get_dimensions();

        for (std::size_t i = 0; i < dimin.x; ++i) {
            for (std::size_t j = 0; j < dimin.y; ++j) {
                u.set_val( (lambda*u.get_val(i,j) + tau * f.get_val(i,j)) / (lambda + tau), i, j);
            }
        }
    }

    // The primal-dual TGV-denoising algorithm
    opticalflow::Image denoise(const opticalflow::Image& f, double tau, double sigma, double lambda, double alpha0, double alpha1, int niter) {
        if (tau <= 0.0)
            throw std::runtime_error("Tau has to be a positive scalar.");

        if (sigma <= 0) 
            throw std::runtime_error("Sigma has to be a positive scalar.");

        // Get image dimensions
        const dim dimin = f.get_dimensions();

        // Initialize Primal-Dual variables
        opticalflow::Image u(f);
        opticalflow::Image ubar(f);

        opticalflow::Image vx(dimin);
        opticalflow::Image vy(dimin);

        opticalflow::Image vbarx(dimin);
        opticalflow::Image vbary(dimin);

        opticalflow::Image px(dimin);
        opticalflow::Image py(dimin);

        opticalflow::Image qxx(dimin);
        opticalflow::Image qyy(dimin);
        opticalflow::Image qxy(dimin);

        // Tracking variables
        opticalflow::Image uold(dimin);
        opticalflow::Image vxold(dimin);
        opticalflow::Image vyold(dimin);

        // Primal-dual iterations
        for (int iter = 0; iter < niter; ++iter) {
            // Update dual variables
            tgv_denoise::update_p(px,py,ubar,vbarx,vbary,sigma);
            tgv_denoise::update_q(qxx, qyy, qxy, vbarx, vbary, sigma);

            tgv_denoise::proj_p(px,py,alpha1);
            tgv_denoise::proj_q(qxx,qyy,qxy,alpha0);

            // Track u
            uold = u;

            // Proximal operator
            tgv_denoise::update_u(u, px, py, tau);
            tgv_denoise::prox(u, f, tau, lambda);

            // Update ubar
            ubar = 2*u - uold;

            // Track v
            vxold = vx;
            vyold = vy;

            // Update v
            tgv_denoise::update_v(vx, vy, px, py, qxx, qyy, qxy, tau);

            // Update vbar
            vbarx = 2*vx - vxold;
            vbary = 2*vy - vyold;

            if (opticalflow::image::norm(u-uold)/opticalflow::image::norm(u) < 1e-6)
                break;

            // Update some norms:
            if (iter % 50 == 0) {
                std::cout << "iter " << iter
                    << " |u| = " << opticalflow::image::norm(u)
                    << " |ubar| = " << opticalflow::image::norm(ubar)
                    << " |v| = ("
                    << opticalflow::image::norm(vx) << ", "
                    << opticalflow::image::norm(vy) << ")"
                    << " |vbar| = ("
                    << opticalflow::image::norm(vbarx) << ", "
                    << opticalflow::image::norm(vbary) << ")"
                    << std::endl;
            }
        }

        return u;
    }
}

#endif