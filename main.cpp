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
        if (m.opp(i)>=0) continue;
        boundary_verts[m.from(i)] = true;
        boundary_verts[m.to  (i)] = true;
    }

    std::vector<int> valency(m.nverts(), -1); // valency of interior vertices
    int max_valency = -1;

    for (int v = 0; v < m.nverts(); v++) {  // compute valencies over interior vertices
        if (boundary_verts[v]) continue;
        int cnt = 0;
        int hfirst = m.first_halfedge(v);
        int hcur = hfirst;
        do {
            cnt++;
            hcur = m.opp(m.prev(hcur));
        } while (hcur != hfirst);
        valency[v] = cnt;
        max_valency = std::max(max_valency, cnt);
    }

    std::vector<std::vector<vec3> > stencils(max_valency+1); // virtual control volumes for all possible valencies
    for (int i = 1; i <= max_valency; i++) {
        // compute the virtual control volume, the vertices are placed uniformly on the unit circle
        stencils[i].push_back(vec3(1, 0, 0));
        for (int j = 1; j < i; j++) {
            double angle = j * 2. * M_PI / i;
            stencils[i].push_back(vec3(cos(angle), sin(angle), 0));
        }
    }

    // Winslow smoothing
    for (int iter = 0; iter < 10000; iter++) { // Gauss-Seidel iterations
        std::cerr << "iter: " << iter << std::endl;
        for (int v = 0; v < m.nverts(); v++) { // the variables of the system are the (non-boundary) vertex positions
            if (boundary_verts[v]) continue;
            std::vector<double> coeff(1+valency[v], 0); // coefficients for the current row of the linear system matrix

            int hfirst = m.first_halfedge(v);
            int cnt = 0;
            int hcur = hfirst;
            do { // iterate through all triangles incident to the vertex v (1-ring of v)
                vec3 virtv1 = vec3(0, 0, 0);
                vec3 virtv2 = stencils[valency[v]][cnt];
                vec3 virtv3 = stencils[valency[v]][(cnt + 1) % valency[v]];
                vec3 n1 = cross(virtv3 - virtv2, vec3(0, 0, 1)); // non-unit!
                vec3 n2 = cross(virtv1 - virtv3, vec3(0, 0, 1));
                vec3 n3 = cross(virtv2 - virtv1, vec3(0, 0, 1));

                vec3 v1 = m.point(v);
                vec3 v2 = m.point(m.to(hcur));
                vec3 v3 = m.point(m.to(m.next(hcur)));

                double x_u = n1.x*v1.x + n2.x*v2.x + n3.x*v3.x;
                double x_v = n1.y*v1.x + n2.y*v2.x + n3.y*v3.x;
                double y_u = n1.x*v1.y + n2.x*v2.y + n3.x*v3.y;
                double y_v = n1.y*v1.y + n2.y*v2.y + n3.y*v3.y;

                double alpha = x_v*x_v + y_v*y_v;
                double beta  = x_u*x_v + y_u*y_v;
                double gamma = x_u*x_u + y_u*y_u;

                vec3 t = n1;
                int c1 = 1 + cnt, c2 = 1 + (cnt+1)%valency[v];
                coeff[0]  += alpha*n1.x*t.x - beta*n1.y*t.x - beta*n1.x*t.y + gamma*n1.y*t.y;
                coeff[c1] += alpha*n2.x*t.x - beta*n2.y*t.x - beta*n2.x*t.y + gamma*n2.y*t.y;
                coeff[c2] += alpha*n3.x*t.x - beta*n3.y*t.x - beta*n3.x*t.y + gamma*n3.y*t.y;

                cnt++;
                hcur = m.opp(m.prev(hcur));
            } while (hcur != hfirst);

            double lambda = .1; // under-relaxation
            m.point(v) = m.point(v) * (1 - lambda);
            cnt = 0;
            hcur = hfirst;
            do { // iterate through all triangles incident to the vertex v (1-ring of v)
                m.point(v) = m.point(v) - m.point(m.to(hcur)) * (lambda * coeff[cnt+1] / coeff[0]);
                cnt++;
                hcur = m.opp(m.prev(hcur));
            } while (hcur != hfirst);
        } // forall inner vertices
    } // system loop

    std::ofstream ofs("drop.obj");
    ofs << m << std::endl;
    ofs.close();
    return 0;
}

