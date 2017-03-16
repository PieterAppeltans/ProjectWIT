#ifndef mesh_reader
#define mesh_reader

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <Eigen/Dense>
using namespace Eigen;
typedef Matrix<int,Dynamic,Dynamic> MatrixXi;

namespace mesh {

    MatrixXd read_vertices(std::string location) {
        /* index # is vertix id, contains x and y
        coordinate respectively of vertix, and
        a flag, 1 if boundary point 0 if not */
        std::ifstream infile(location);
        std::string line;
        double a, b, c, d, line_nb = 0;
        std::getline(infile, line);
        std::istringstream iss(line);
        iss >> a >> b >> c;
        int nb_nodes = a;
        MatrixXd vertices(nb_nodes, 3);

        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            if (iss >> a >> b >> c >> d) {
                vertices(line_nb, 0) = b;
                vertices(line_nb, 1) = c;
                vertices(line_nb, 2) = d;
            }
            else {
                // std::cout << "ERROR: in .1.node file" << std::endl;
                break;
            }
            line_nb += 1;
        }

        return vertices;
    }

    MatrixXd read_triangles(MatrixXd &vertices,std::string location) {
        /* index # is triangle id, contains
        the 3 vertix id's of the corner nodes  and
        the area of the triangle, calculated using
        the vertices matrix */
        std::ifstream infile(location);
        std::string line;
        int a, b, c, d, line_nb = 0;
        std::getline(infile, line);
        std::istringstream iss(line);
        iss >> a >> b >> c;
        int nb_triangles = a;
        MatrixXd triangles(nb_triangles, 4);

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
                // std::cout << "ERROR: in .1.ele file" << std::endl;
                break;
            }
            line_nb += 1;
        }

        return triangles;
    }

    MatrixXi read_boundaries(MatrixXd &vertices,std::string location) {
        /* index # is triangle id, contains
        the 3 vertix id's of the corner nodes  and
        the area of the triangle, calculated using
        the vertices matrix */
        std::ifstream infile(location);
        std::string line;
        int a, b, c, d, count = 0;
        std::getline(infile, line);
        std::getline(infile, line);
        std::istringstream iss(line);
        iss >> a >> b;
        int nb_boundaries = a;
        MatrixXi boundaries(nb_boundaries, 2);

        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            if (iss >> a >> b >> c >> d) {
                if (vertices(b-1, 0) != 0 || vertices(c-1, 0) != 0) {
                  boundaries(count, 0) = b-1;
                  boundaries(count, 1) = c-1;
                  count += 1;
                }
            }
            else {
                // std::cout << "ERROR: in .1.poly file" << std::endl;
                break;
            }
        }

        MatrixXi boundaries_red(count, 2);
        boundaries_red = boundaries.block(0,0,count,2);
        return boundaries_red;
    }

}

#endif
