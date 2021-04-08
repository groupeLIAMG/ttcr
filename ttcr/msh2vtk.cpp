//
//  msh2vtk.cpp
//  ttcr
//
//  Created by Bernard Giroux on 2014-10-18.
//  Copyright (c) 2014 Bernard Giroux. All rights reserved.
//

/*
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdio.h>

#include <cstdint>
#include <exception>
#include <iostream>
#include <set>
#include <sstream>


#include "MSHReader.h"
#include "Node2Dc.h"
#include "Node2Dn.h"
#include "Grid2Ducfm.h"
#include "Grid2Dunfm.h"
#include "Grid3Ducfm.h"
#include "Grid3Dunfm.h"
#include "Rcv2D.h"
#include "Rcv.h"

#include "msh2vtk_io.h"

using namespace std;
using namespace ttcr;

int main(int argc, char * argv[])
{

    input_parameters par;
    parse_input(argc, argv, par);


    string fname(par.mshFile);

    MSHReader reader(fname.c_str());

    if ( !reader.isValid() ) {
        cerr << "File " << fname << " invalid or not found\n";
        abort();
    }

    map<string, double> slownesses;

    int d = reader.get2Ddim();
    if ( d == 1 || d == 2 ) {

        if ( ttcr::verbose ) {
            cout << "Number of nodes: " << reader.getNumberOfNodes() << '\n';
            cout << "Number of elements: " << reader.getNumberOfElements() << '\n';
            cout << "Number of media: " << reader.getPhysicalNames(2).size() << '\n';
            cout << "Number of reflectors: " << reader.getPhysicalNames(1).size() << '\n';
        }

        vector<sxz<double> > nodes(reader.getNumberOfNodes());
        vector<triangleElem<uint32_t> > triangles(reader.getNumberOfTriangles());
        vector<double> slowness(reader.getNumberOfTriangles());

        reader.readNodes2D(nodes, d);
        reader.readTriangleElements(triangles);

        bool constCells = true;
        if ( !par.sloFile.empty() ) {
            std::ifstream fin(par.sloFile.c_str());
            if ( !fin ) {
                std::cout << "Error: cannot open file " << par.sloFile << std::endl;
                exit ( -1);
            }
            std::vector<double> tmp;
            double dtmp;
            fin >> dtmp;
            while ( fin ) {
                tmp.push_back( dtmp );
                fin >> dtmp;
            }
            fin.close();
            if ( tmp.size() != slowness.size() ) {
                if ( tmp.size() == nodes.size() ) {
                    slowness.resize( nodes.size() );
                    constCells = false;
                } else {
                    std::cerr << "Error: slowness file should contain " << slowness.size()
                    << " values.\nAborting." << std::endl;
                    abort();
                }
            }
            for ( size_t n=0; n<slowness.size(); ++n ) {
                slowness[n] = tmp[n];
            }

        } else if ( par.velFile == "" ) {
            for ( size_t n=0; n<reader.getPhysicalNames(2).size(); ++n ) {
                double s;

                cout << "  Input velocity in " << reader.getPhysicalNames(2)[n] << ": ";
                cin >> s;
                slownesses.insert( {reader.getPhysicalNames(2)[n], 1./s} );
            }
        } else {
            std::ifstream fin(par.velFile.c_str());
            if ( !fin ) {
                cerr << "Error: cannot open file " << par.velFile << endl;
                exit ( -1);
            }
            std::string line;
            while ( fin ) {
                getline( fin, line );
                if ( line.empty() ) continue;
                size_t i1 = line.find('"');
                size_t i2 = line.rfind('"');
                string name = line.substr(i1+1, i2-i1-1);
                istringstream sin( line.substr(i2+1, 100) );
                double val;
                sin >> val;
                slownesses.insert( {name, 1./val} );
            }
            fin.close();

            if ( ttcr::verbose ) {
                for ( size_t n=0; n<reader.getPhysicalNames(2).size(); ++n ) {
                    cout << "  Velocity for " << reader.getPhysicalNames(2)[n] << " is "
                    << 1./slownesses[ reader.getPhysicalNames(2)[n] ] << '\n';
                }
            }

            for ( size_t n=0; n<slowness.size(); ++n ) {
                slowness[n] = slownesses[ reader.getPhysicalNames(2)[triangles[n].physical_entity] ];
            }
        }


        if ( ttcr::verbose ) {
            cout << "Creating grid ... ";
            cout.flush();
        }
        Grid2D<double,uint32_t,sxz<double>> *g=nullptr;
        if ( constCells )
            g = new Grid2Ducfm<double,uint32_t, Node2Dc<double,uint32_t>,sxz<double>>(nodes, triangles);
        else
            g = new Grid2Dunfm<double,uint32_t, Node2Dn<double,uint32_t>,sxz<double>>(nodes, triangles);

        if ( ttcr::verbose ) {
            cout << "done.\n";
        }
        try {
            g->setSlowness(slowness);
        } catch (std::exception& e) {
            cerr << e.what() << endl;
            abort();
        }

        if ( par.rectilinear ) {
            double d[] = { par.d, par.d, par.d };
            if ( ttcr::verbose ) {
                cout << "Saving " << par.vtkFile << " ... ";
                cout.flush();
            }
            g->saveModelVTR(par.vtkFile, d, par.saveSlowness);
            if ( ttcr::verbose ) {
                cout << "done.\n";
            }
        } else {
            if ( ttcr::verbose ) {
                cout << "Saving " << par.vtkFile << " ... ";
                cout.flush();
            }
            g->saveModelVTU(par.vtkFile, par.saveSlowness, true);
            if ( ttcr::verbose ) {
                cout << "done.\n";
            }
        }

        if ( par.saveReflectors && reader.getPhysicalNames(1).size()>0 ) {
            vector<string> reflector_names = reader.getPhysicalNames(1);
            vector<int> indices = reader.getPhysicalIndices(1);

            if ( reflector_names.size() != indices.size() ) {
                cerr << "Error - definition of reflectors\n";
                exit(1);
            }

            std::vector<lineElem<uint32_t>> lines;
            reader.readLineElements(lines);

            for ( size_t ni=0; ni<indices.size(); ++ni ) {

                set<uint32_t> refl_ind;
                set<uint32_t>::iterator it;

                string fname = reflector_names[ni] + ".dat";
                Rcv2D<double> reflector(fname);

                for ( size_t nl=0; nl<lines.size(); ++nl ) {
                    if ( indices[ni] == lines[nl].physical_entity ) {

                        refl_ind.insert( lines[nl].i[0] );
                        refl_ind.insert( lines[nl].i[1] );

                    }
                }

                //				cout << "myset contains:";
                for (it=refl_ind.begin(); it!=refl_ind.end(); ++it) {
                    //					cout << ' ' << *it;
                    reflector.add_coord( nodes[*it] );
                }
                //				cout << '\n';
                if ( ttcr::verbose ) cout << "Saving reflector file " << fname << '\n';
                reflector.save_rcvfile();
            }
        }

        delete g;

    } else {

        if ( ttcr::verbose ) {
            cout << "Number of nodes: " << reader.getNumberOfNodes() << '\n';
            cout << "Number of elements: " << reader.getNumberOfElements() << '\n';
            cout << "Number of media: " << reader.getPhysicalNames(3).size() << '\n';
            cout << "Number of physical entities: " << reader.getPhysicalNames().size() << '\n';
        }

        vector<sxyz<double> > nodes(reader.getNumberOfNodes());
        vector<tetrahedronElem<uint32_t> > tetrahedra(reader.getNumberOfTetra());
        vector<double> slowness(reader.getNumberOfTetra());

        reader.readNodes3D(nodes);
        reader.readTetrahedronElements(tetrahedra);

        bool constCells = true;
        if ( !par.sloFile.empty() ) {
            std::ifstream fin(par.sloFile.c_str());
            if ( !fin ) {
                std::cout << "Error: cannot open file " << par.sloFile << std::endl;
                exit ( -1);
            }
            std::vector<double> tmp;
            double dtmp;
            fin >> dtmp;
            while ( fin ) {
                tmp.push_back( dtmp );
                fin >> dtmp;
            }
            fin.close();
            if ( tmp.size() != slowness.size() ) {
                if ( tmp.size() == nodes.size() ) {
                    slowness.resize( nodes.size() );
                    constCells = false;
                } else {
                    std::cerr << "Error: slowness file should contain " << slowness.size()
                    << " values.\nAborting." << std::endl;
                    abort();
                }
            }
            for ( size_t n=0; n<slowness.size(); ++n ) {
                slowness[n] = tmp[n];
            }

        } else if ( par.velFile == "" ) {
            for ( size_t n=0; n<reader.getPhysicalNames().size(); ++n ) {
                double s;

                cout << "  Input velocity in " << reader.getPhysicalNames()[n] << ": ";
                cin >> s;
                slownesses.insert( {reader.getPhysicalNames()[n], 1./s} );
            }
        } else {
            std::ifstream fin(par.velFile.c_str());
            if ( !fin ) {
                cerr << "Error: cannot open file " << par.velFile << endl;
                exit ( -1);
            }
            std::string line;
            while ( fin ) {
                getline( fin, line );
                if ( line.empty() ) continue;
                size_t i1 = line.find('"');
                size_t i2 = line.rfind('"');
                string name = line.substr(i1+1, i2-i1-1);
                istringstream sin( line.substr(i2+1, 100) );
                double val;
                sin >> val;
                slownesses.insert( {name, 1./val} );
            }
            fin.close();

            if ( ttcr::verbose ) {
                for ( size_t n=0; n<reader.getPhysicalNames().size(); ++n ) {
                    cout << "  Velocity for " << reader.getPhysicalNames()[n] << " is "
                    << 1./slownesses[ reader.getPhysicalNames()[n] ] << '\n';
                }
            }
            for ( size_t n=0; n<slowness.size(); ++n ) {
                slowness[n] = slownesses[ reader.getPhysicalNames()[tetrahedra[n].physical_entity] ];
            }
        }

        Grid3D<double, uint32_t> *g = nullptr;
        if ( ttcr::verbose ) {
            cout << "Creating grid ... ";
            cout.flush();
        }
        if ( constCells )
            g = new Grid3Ducfm<double, uint32_t>(nodes, tetrahedra, false, false, false, 1.e-5);
        else
            g = new Grid3Dunfm<double, uint32_t>(nodes, tetrahedra, false, false, false, 1.e-5);
        if ( ttcr::verbose ) {
            cout << "done.\n";
        }
        try {
            g->setSlowness(slowness);
        } catch (std::exception& e) {
            cerr << e.what() << endl;
            abort();
        }

        if ( par.rectilinear ) {
            double d[] = { par.d, par.d, par.d };
            if ( ttcr::verbose ) {
                cout << "Saving " << par.vtkFile << " ... ";
                cout.flush();
            }
            g->saveModelVTR(par.vtkFile, d, par.saveSlowness);
            if ( ttcr::verbose ) {
                cout << "done.\n";
            }
        } else {
            if ( ttcr::verbose ) {
                cout << "Saving " << par.vtkFile << " ... ";
                cout.flush();
            }
            g->saveModelVTU(par.vtkFile, par.saveSlowness, true);
            if ( ttcr::verbose ) {
                cout << "done.\n";
            }
        }

        if ( par.saveReflectors && reader.getPhysicalNames(2).size()>0 ) {
            vector<string> reflector_names = reader.getPhysicalNames(2);
            vector<int> indices = reader.getPhysicalIndices(2);

            if ( reflector_names.size() != indices.size() ) {
                cerr << "Error - definition of reflectors\n";
                exit(1);
            }

            std::vector<triangleElem<uint32_t>> triangles;
            reader.readTriangleElements(triangles);

            for ( size_t ni=0; ni<indices.size(); ++ni ) {

                set<uint32_t> refl_ind;
                set<uint32_t>::iterator it;

                string fname = reflector_names[ni] + ".dat";
                Rcv<double> reflector(fname);

                for ( size_t nl=0; nl<triangles.size(); ++nl ) {
                    if ( indices[ni] == triangles[nl].physical_entity ) {

                        refl_ind.insert( triangles[nl].i[0] );
                        refl_ind.insert( triangles[nl].i[1] );
                        refl_ind.insert( triangles[nl].i[2] );

                    }
                }

                for (it=refl_ind.begin(); it!=refl_ind.end(); ++it) {
                    reflector.add_coord( nodes[*it] );
                }

                if ( ttcr::verbose ) cout << "Saving reflector file " << fname << '\n';
                reflector.save_rcvfile();
            }

        }

        delete g;
    }
    return 0;
}
