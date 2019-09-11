#define _USE_MATH_DEFINES
#include <cmath>
#undef NDEBUG
#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>
#include "geometry.h"
#include "model.h"

int main(int argc, char** argv) {
    if (argc<2) {
        std::cerr << "Usage: " << argv[0] << " obj/model.obj" << std::endl;
        return 1;
    }

    Model m(argv[1]);

    std::vector<bool> boundary_verts(m.nverts(), false);
    for (int i=0; i<m.nhalfedges(); i++) {
        if (m.opp(i)<0) {
            boundary_verts[m.from(i)] = true;
            boundary_verts[m.to  (i)] = true;
        }
    }

    std::vector<int> valency(m.nverts(), -1); // valency of interior vertices
    int nintv = 0;                            // number of interior vertices

    for (int v = 0; v < m.nverts(); v++) {
        if (boundary_verts[v]) continue;
        nintv++;

        { // compute valency
            int cnt = 0;
            int hfirst = m.first_halfedge(v);
            int hcur = hfirst;
            do {
                cnt++;
                hcur = m.opp(m.prev(hcur));
            } while (hcur != hfirst);
            valency[v] = cnt;
        }
    }


#if 0
    // Laplacian smoothing
    for (int iter = 0; iter < 1000; iter++) {
        std::cerr << "iter: " << iter << std::endl;

        for (int v = 0; v < m.nverts(); v++) {
            if (boundary_verts[v]) continue;
            int hfirst = m.first_halfedge(v);
            int hcur = hfirst;
            Vec3f p(0, 0, 0);
            do {
                p = p + m.point(m.to(hcur))/valency[v];
                hcur = m.opp(m.prev(hcur));
            } while (hcur != hfirst);
            m.point(v) = p;
        }

    }
#else
    // Winslow smoothing
    for (int iter = 0; iter < 10000; iter++) { // Gauss-Seidel iterations
        std::cerr << "iter: " << iter << std::endl;
        for (int v = 0; v < m.nverts(); v++) { // the variables of the system are the (non-boundary) vertex positions
            if (boundary_verts[v]) continue;
            std::vector<double> coeff(m.nverts(), 0); // coefficients for the current row of the linear system matrix

            int hfirst = m.first_halfedge(v);
            int cnt = 0;
            int hcur = hfirst;
            do { // iterate through all triangles incident to the vertex v (1-ring of v)
                std::vector<Vec3f> stencil(valency[v], Vec3f(0, 0, 0)); // virtual control volume for the 1-ring of the vertex v

                { // compute the virtual control volume, the vertices are placed uniformly on a circle
                    stencil[0] = Vec3f(1, 0, 0);
                    stencil[0] = Vec3f(-sqrtf(2), sqrtf(2), 0);
                    stencil[0] = (m.point(m.to(hfirst)) - m.point(m.from(hfirst))).normalize(); // in fact, the result should be independent from the stencil orientation; TODO: check it
                    for (int i = 1; i < valency[v]; i++) { // rotate the firt vertex to place the rest on the circle
                        double angle = i * 2. * M_PI / valency[v];
                        stencil[i] = Vec3f(stencil[0].x * cos(angle) - stencil[0].y * sin(angle), stencil[0].x * sin(angle) + stencil[0].y * cos(angle), 0);
                    }
                }

                Vec3f virtv1 = Vec3f(0, 0, 0);
                Vec3f virtv2 = stencil[cnt];
                Vec3f virtv3 = stencil[(cnt + 1) % valency[v]];
                Vec3f n1 = cross(virtv3 - virtv2, Vec3f(0, 0, 1)); // non-unit!
                Vec3f n2 = cross(virtv1 - virtv3, Vec3f(0, 0, 1));
                Vec3f n3 = cross(virtv2 - virtv1, Vec3f(0, 0, 1));
                Vec3f t = n1;

                int i1 = m.from(hcur); // i1=v
                int i2 = m.to(hcur);
                int i3 = m.to(m.next(hcur));

                Vec3f v1 = m.point(i1);
                Vec3f v2 = m.point(i2);
                Vec3f v3 = m.point(i3);

                double x_xi = n1.x * v1.x + n2.x * v2.x + n3.x * v3.x;
                double x_eta = n1.y * v1.x + n2.y * v2.x + n3.y * v3.x;
                double y_xi = n1.x * v1.y + n2.x * v2.y + n3.x * v3.y;
                double y_eta = n1.y * v1.y + n2.y * v2.y + n3.y * v3.y;
                double alpha = x_eta * x_eta + y_eta * y_eta;
                double beta = x_xi * x_eta + y_xi * y_eta;
                double gamma  = x_xi * x_xi + y_xi * y_xi;

                coeff[i1] += alpha * n1.x * t.x - beta * n1.y * t.x - beta * n1.x * t.y + gamma * n1.y * t.y;
                coeff[i2] += alpha * n2.x * t.x - beta * n2.y * t.x - beta * n2.x * t.y + gamma * n2.y * t.y;
                coeff[i3] += alpha * n3.x * t.x - beta * n3.y * t.x - beta * n3.x * t.y + gamma * n3.y * t.y;
                //coeff[i1] += alpha * n1.x * t.x - 2 * beta * n1.y * t.x + gamma * n1.y * t.y;
                //coeff[i2] += alpha * n2.x * t.x - 2 * beta * n2.y * t.x + gamma * n2.y * t.y;
                //coeff[i3] += alpha * n3.x * t.x - 2 * beta * n3.y * t.x + gamma * n3.y * t.y;
                cnt++;
                hcur = m.opp(m.prev(hcur));
            } while (hcur != hfirst);

            double lambda = .1; // under-relaxation
            m.point(v) = m.point(v) * (1 - lambda);
            for (int i = 0; i < m.nverts(); i++) {
                if (i == v) continue;
                m.point(v) = m.point(v) - m.point(i) * (lambda * coeff[i] / coeff[v]);
            }
        } // forall inner vertices
    } // system loop
#endif

    std::ofstream ofs("drop.obj");
    ofs << m << std::endl;
    ofs.close();
    return 0;
}

