#ifndef mesh_reader
#define mesh_reader

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;

namespace mesh {

    matrix<double> read_vertices() {
        /* index # is vertix id, contains x and y
        coordinate respectively of vertix, and
        a flag, 1 if boundary point 0 if not */
        std::ifstream infile("../triangle/pear.1.node");
        std::string line;
        double a, b, c, d;
        int line_nb = 0;
        std::getline(infile, line);
        std::istringstream iss(line);
        iss >> a >> b >> c;
        int nb_nodes = a;
        matrix<double> vertices (nb_nodes, 3);

        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            if (iss >> a >> b >> c >> d) {
                vertices(line_nb, 0) = b;
                vertices(line_nb, 1) = c;
                vertices(line_nb, 2) = d;
            }
            else {
                break;
            }
            line_nb += 1;
        }

        return vertices;
    }

    matrix<double> read_triangles(matrix<double> &vertices) {
        /* index # is triangle id, contains
        the 3 vertix id's of the corner nodes  and
        the area of the triangle, calculated using
        the vertices matrix */
        std::ifstream infile("../triangle/pear.1.ele");
        std::string line;
        int a, b, c, d;
        int line_nb = 0;
        std::getline(infile, line);
        std::istringstream iss(line);
        iss >> a >> b >> c;
        int nb_triangles = a;
        matrix<double> triangles (nb_triangles, 4);

        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            if (iss >> a >> b >> c >> d) {
                triangles(line_nb, 0) = b-1;
                triangles(line_nb, 1) = c-1;
                triangles(line_nb, 2) = d-1;
                triangles(line_nb, 3) = 0.5*fabs(vertices(b-1, 0)*(vertices(c-1, 1) - vertices(d-1, 1)) +
                    vertices(c-1, 0)*(vertices(d-1, 1) - vertices(b-1, 1)) +
                    vertices(d-1, 0)*(vertices(b-1, 1) - vertices(c-1, 1)));
            }
            else {
                break;
            }
            line_nb += 1;
        }

        return triangles;
    }

}

#endif
