#include <fstream>
#include <vector>
#include <string>
#include <iomanip>      // dla std::setprecision
#include <iostream>
#include <sys/stat.h>   // mkdir
#include <sys/types.h>  // mode_t
#include "structs.h"

void write_to_vtu_file(int step,
                       const std::vector<Fem::Node>& nodes,
                       const std::vector<double>& temp,
                       const std::vector<Fem::Element>& elements)
{
    struct stat st;

    std::ostringstream fname;
    fname << "sol_" << step << ".vtu";
    const std::string path = fname.str();

    std::ofstream out(path);
    if (!out.is_open()) {
        std::perror(("Nie można otworzyć pliku " + path).c_str());
        return;
    }

    out << R"(<?xml version="1.0"?>)"
        << "\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
           "  <UnstructuredGrid>\n"
           "    <Piece NumberOfPoints=\"" << nodes.size()
        << "\" NumberOfCells=\"" << elements.size() << "\">\n";

    // --- 4) Points ---
    out << "      <Points>\n"
           "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    out << std::setprecision(6);
    for (auto& n : nodes) {
        out << "          " << n.x << " " << n.y << " 0\n";
    }
    out << "        </DataArray>\n"
           "      </Points>\n";

    // --- 5) Cells ---
    out << "      <Cells>\n"
           "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (auto& e : elements) {
        out << "          "
            << (e.node_ids[0] - 1) << " "
            << (e.node_ids[1] - 1) << " "
            << (e.node_ids[2] - 1) << " "
            << (e.node_ids[3] - 1) << "\n";
    }
    out << "        </DataArray>\n"
           "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int offset = 0;
    for (size_t i = 0; i < elements.size(); ++i) {
        offset += 4;  // Zmieniono z 3 na 4 (czworokąt ma 4 węzły)
        out << "          " << offset << "\n";
    }
    out << "        </DataArray>\n"
           "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t i = 0; i < elements.size(); ++i) {
        out << "          9\n";  // Zmieniono z 5 na 9 (typ VTK_QUAD)
    }
    out << "        </DataArray>\n"
           "      </Cells>\n";

    // --- 6) PointData: temperatura ---
    out << "      <PointData Scalars=\"Temperature\">\n"
           "        <DataArray type=\"Float32\" Name=\"Temperature\" format=\"ascii\">\n";
    for (double t : temp) {
        out << "          " << t << "\n";
    }
    out << "        </DataArray>\n"
           "      </PointData>\n";

    out << "    </Piece>\n"
           "  </UnstructuredGrid>\n"
           "</VTKFile>\n";

    out.close();
}
