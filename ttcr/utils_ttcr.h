//
//  utils_ttcr.h
//  ttcr_u
//
//  Created by Bernard Giroux on 2014-01-23.
//  Copyright (c) 2014 Bernard Giroux. All rights reserved.
//

#ifndef ttcr_u_utils_ttcr_h
#define ttcr_u_utils_ttcr_h

#include <chrono>
#include <set>
#include <string>
#include <vector>

#ifdef VTK
#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkRectilinearGrid.h>
#include <vtkSmartPointer.h>
#include <vtkXMLReader.h>
#include <vtkXMLRectilinearGridReader.h>
#endif

#include "Grid2Drc.h"
#include "Grid2Ducfm.h"
#include "Grid2Ducfs.h"
#include "Grid2Ducsp.h"
#include "Grid2Duifm.h"
#include "Grid2Duifs.h"
#include "Grid2Duisp.h"
#include "Grid3Drc.h"
#include "Grid3Dri.h"
#include "Grid3Ducfm.h"
#include "Grid3Ducfs.h"
#include "Grid3Ducsp.h"
#include "Grid3Duifm.h"
#include "Grid3Duifs.h"
#include "Grid3Duisp.h"
#include "Node2Dcsp.h"
#include "Node2Disp.h"
#include "MSHReader.h"
#ifdef VTK
#include "VTUReader.h"
#endif

#include "Rcv.h"
#include "Rcv2D.h"

#include "utils.h"

template<typename T>
std::string to_string( const T & Value )
{
    // utiliser un flux de sortie pour créer la chaîne
    std::ostringstream oss;
    // écrire la valeur dans le flux
    oss << Value;
    // renvoyer une string
    return oss.str();
}


#ifdef VTK
template<typename T>
Grid3Dr<T,uint32_t> *recti3D(const input_parameters &par, const size_t nt)
{
    Grid3Dr<T,uint32_t> *g = nullptr;
	vtkRectilinearGrid *dataSet;
    
    vtkSmartPointer<vtkXMLRectilinearGridReader> reader =
    vtkSmartPointer<vtkXMLRectilinearGridReader>::New();
    reader->SetFileName(par.modelfile.c_str());
    reader->Update();
    reader->GetOutput()->Register(reader);
    dataSet = reader->GetOutput();
	
    vtkIdType numberOfCells = dataSet->GetNumberOfCells();
    vtkIdType numberOfPoints = dataSet->GetNumberOfPoints();
    int nnodes[3], ncells[3];
    dataSet->GetDimensions(nnodes);
    ncells[0] = nnodes[0]-1;
    ncells[1] = nnodes[1]-1;
    ncells[2] = nnodes[2]-1;
    double xrange[2], yrange[2], zrange[2];
    dataSet->GetXCoordinates()->GetRange(xrange);
    dataSet->GetYCoordinates()->GetRange(yrange);
    dataSet->GetZCoordinates()->GetRange(zrange);
    double d[3];
    d[0] = (xrange[1]-xrange[0])/(nnodes[0]-1);
    d[1] = (yrange[1]-yrange[0])/(nnodes[1]-1);
    d[2] = (zrange[1]-zrange[0])/(nnodes[2]-1);
	
    if ( par.verbose ) {
        std::cout << "Reading model file " << par.modelfile
        << "\n  Rectilinear grid in file has"
        << "\n    " << nnodes[0]*nnodes[1]*nnodes[2] << " nodes"
        << "\n    " << ncells[0]*ncells[1]*ncells[2] << " cells"
        << "\n    (size: " << nnodes[0] << " x " << nnodes[1] << " x " << nnodes[2] << ')'
        << "\n  Dim\tmin\tmax\tinc.\t N. sec nodes"
        << "\n   X\t" << xrange[0] << '\t' << xrange[1] << '\t' << d[0] << "\t\t" << par.nn[0]
        << "\n   Y\t" << yrange[0] << '\t' << yrange[1] << '\t' << d[1] << "\t\t" << par.nn[1]
        << "\n   Z\t" << zrange[0] << '\t' << zrange[1] << '\t' << d[2] << "\t\t" << par.nn[2]
        << std::endl;
    }
    vtkPointData *pd = dataSet->GetPointData();
    vtkCellData *cd = dataSet->GetCellData();
	
	std::chrono::high_resolution_clock::time_point begin, end;
	
    bool foundSlowness = false;
	std::vector<T> slowness;
    if ( pd->HasArray("P-wave velocity") ||
		pd->HasArray("Velocity") ||
		pd->HasArray("Slowness") ) {
        for (int na = 0; na < pd->GetNumberOfArrays(); ++na) {
            if ( strcmp(pd->GetArrayName(na), "P-wave velocity")==0 ||
				strcmp(pd->GetArrayName(na), "Velocity")==0 ) {
                slowness.resize( numberOfPoints );
                for ( size_t k=0,n=0; k<nnodes[2]; ++k ) {
                    for ( size_t j=0; j<nnodes[1]; ++j ) {
                        for ( size_t i=0; i<nnodes[0]; ++i,++n ) {
                            slowness[n] = static_cast<T>(1./pd->GetArray(na)->GetTuple1(n));
                        }
                    }
                }
                foundSlowness = true;
				break;
			} else if ( strcmp(pd->GetArrayName(na), "Slowness")==0 ) {
				
				vtkSmartPointer<vtkDoubleArray> slo = vtkSmartPointer<vtkDoubleArray>::New();
				slo = vtkDoubleArray::SafeDownCast( cd->GetArray("Slowness") );
				
				if ( slo->GetSize() != dataSet->GetNumberOfPoints() ) {
					std::cerr << "Problem with Slowness data (wrong size)" << std::endl;
					return 0;
				}
				
				slowness.resize( slo->GetSize() );
				for ( size_t n=0; n<slo->GetSize(); ++n ) {
					slowness[n] = slo->GetComponent(n, 0);
				}
				foundSlowness = true;
				break;
			}
			
			
		}
		if ( foundSlowness ) {
			if ( par.verbose ) { std::cout << "Building grid (Grid3Dri) ... "; std::cout.flush(); }
			if ( par.time ) { begin = std::chrono::high_resolution_clock::now(); }
			g = new Grid3Dri<T, uint32_t>(ncells[0], ncells[1], ncells[2],
										  d[0], d[1], d[2],
										  xrange[0], yrange[0], zrange[0],
										  par.nn[0], par.nn[1], par.nn[2],
										  nt, par.inverseDistance);
			if ( par.time ) { end = std::chrono::high_resolution_clock::now(); }
			if ( par.verbose ) {
                std::cout << "done.\nTotal number of nodes: " << g->getNumberOfNodes()
                << "\nAssigning slowness at grid nodes ... ";
                std::cout.flush();
            }
			g->setSlowness( slowness );
			if ( par.verbose ) std::cout << "done.\n";
			if ( par.verbose && par.inverseDistance )
				std::cout << "  Inverse distance interpolation was used.\n";
			if ( par.time ) {
				std::cout.precision(12);
				std::cout << "Time to build grid: " << std::chrono::duration<double>(end-begin).count() << '\n';
			}
        } else {
			return 0;
		}
        
    } else if ( cd->HasArray("P-wave velocity") ||
			   cd->HasArray("Velocity") ||
			   cd->HasArray("Slowness") ) {
        for (int na = 0; na < cd->GetNumberOfArrays(); na++) {
            if ( strcmp(cd->GetArrayName(na), "P-wave velocity")==0 ||
				strcmp(cd->GetArrayName(na), "Velocity")==0 ) {
                slowness.resize( numberOfCells );
                for ( size_t k=0,n=0; k<ncells[2]; ++k ) {
                    for ( size_t j=0; j<ncells[1]; ++j ) {
                        for ( size_t i=0; i<ncells[0]; ++i,++n ) {
                            slowness[n] = static_cast<T>(1./cd->GetArray(na)->GetTuple1(n));
                        }
                    }
                }
                foundSlowness = true;
				break;
			} else if ( strcmp(cd->GetArrayName(na), "Slowness")==0 ) {
                
				vtkSmartPointer<vtkDoubleArray> slo = vtkSmartPointer<vtkDoubleArray>::New();
				slo = vtkDoubleArray::SafeDownCast( cd->GetArray("Slowness") );
				
				if ( slo->GetSize() != dataSet->GetNumberOfCells() ) {
					std::cerr << "Problem with Slowness data (wrong size)" << std::endl;
					return 0;
				}
				
				slowness.resize( slo->GetSize() );
				for ( size_t n=0; n<slo->GetSize(); ++n ) {
					slowness[n] = slo->GetComponent(n, 0);
				}
				foundSlowness = true;
				break;
			}
		}
		if ( foundSlowness ) {
			if ( par.verbose ) { std::cout << "Building grid (Grid3Drc) ... "; std::cout.flush(); }
			if ( par.time ) { begin = std::chrono::high_resolution_clock::now(); }
			g = new Grid3Drc<T, uint32_t>(ncells[0], ncells[1], ncells[2],
										  d[0], d[1], d[2],
										  xrange[0], yrange[0], zrange[0],
										  par.nn[0], par.nn[1], par.nn[2], nt);
			if ( par.time ) { end = std::chrono::high_resolution_clock::now(); }
			if ( par.verbose ) {
                std::cout << "done.\nTotal number of nodes: " << g->getNumberOfNodes()
                << "\nAssigning slowness at grid cells ... ";
            }
			g->setSlowness( slowness );
			if ( par.verbose ) std::cout << "done.\n";
			if ( par.time ) {
				std::cout.precision(12);
				std::cout << "Time to build grid: " << std::chrono::duration<double>(end-begin).count() << '\n';
			}
		} else {
			return 0;
		}
    }
    dataSet->Delete();
    return g;
};

template<typename T>
Grid3D<T, uint32_t> *unstruct3D_vtu(const input_parameters &par, const size_t nt)
{
	
    VTUReader reader( par.modelfile.c_str() );
    
    if ( !reader.isValid() ) {
        return 0;
    }
    
	if ( par.verbose ) {
		std::cout << "Reading model file " << par.modelfile << " ... ";
		std::cout.flush();
	}
	
    std::vector<sxyz<T>> nodes(reader.getNumberOfNodes());
	std::vector<tetrahedronElem<uint32_t>> tetrahedra(reader.getNumberOfElements());
    bool constCells = reader.isConstCell();
    
	std::vector<T> slowness;
    if ( constCells )
        slowness.resize(reader.getNumberOfElements());
    else
        slowness.resize(reader.getNumberOfNodes());
	
	reader.readNodes3D(nodes);
	reader.readTetrahedronElements(tetrahedra);
    reader.readSlowness(slowness, constCells);
    
    if ( par.verbose ) {
        std::cout << "  Unstructured mesh in file has"
        << "\n    " << nodes.size() << " nodes"
        << "\n    " << tetrahedra.size() << " cells";
		if ( constCells )
			std::cout << "\n  Mesh has cells of constant slowness";
		else
			std::cout << "\n  Mesh has slowness defined at nodes";
		std::cout << std::endl;
    }
	
	std::chrono::high_resolution_clock::time_point begin, end;
	Grid3D<T, uint32_t> *g;
    switch (par.method) {
        case SHORTEST_PATH:
        {
            if ( par.verbose ) {
                std::cout << "Creating grid using " << par.nn[0] << " secondary nodes ... ";
                std::cout.flush();
            }
            if ( par.time ) { begin = std::chrono::high_resolution_clock::now(); }
			if ( constCells )
				g = new Grid3Ducsp<T, uint32_t>(nodes, tetrahedra,par.nn[0], nt,
												par.verbose);
			else
				g = new Grid3Duisp<T, uint32_t>(nodes, tetrahedra,par.nn[0], nt,
												par.verbose);
            if ( par.time ) { end = std::chrono::high_resolution_clock::now(); }
            if ( par.verbose ) {
                std::cout << "done.\nTotal number of nodes: " << g->getNumberOfNodes()
                << "\n";
                std::cout.flush();
            }

            break;
        }
        case FAST_MARCHING:
        {
            if ( par.verbose ) {
                std::cout << "Creating grid ... ";
                std::cout.flush();
            }
            if ( par.time ) { begin = std::chrono::high_resolution_clock::now(); }
			if ( constCells )
				g = new Grid3Ducfm<T, uint32_t>(nodes, tetrahedra, nt);
			else
				g = new Grid3Duifm<T, uint32_t>(nodes, tetrahedra, nt);
            if ( par.time ) { end = std::chrono::high_resolution_clock::now(); }
            if ( par.verbose ) {
                std::cout << "done.\n";
                std::cout.flush();
            }

            break;
        }
        case FAST_SWEEPING:
        {
            if ( par.verbose ) {
                std::cout << "Creating grid ... ";
                std::cout.flush();
            }
            if ( par.time ) { begin = std::chrono::high_resolution_clock::now(); }
			if ( constCells )
				g = new Grid3Ducfs<T, uint32_t>(nodes, tetrahedra, par.epsilon,
												par.nitermax, nt);
			else
				g = new Grid3Duifs<T, uint32_t>(nodes, tetrahedra, par.epsilon,
												par.nitermax, nt);
            T xmin = g->getXmin();
            T xmax = g->getXmax();
            T ymin = g->getYmin();
            T ymax = g->getYmax();
            T zmin = g->getZmin();
            T zmax = g->getZmax();
            
            std::vector<sxyz<T>> ptsRef;
            ptsRef.push_back( {xmin, ymin, zmin} );
            ptsRef.push_back( {xmin, ymin, zmax} );
            ptsRef.push_back( {xmin, ymax, zmin} );
            ptsRef.push_back( {xmin, ymax, zmax} );
            ptsRef.push_back( {xmax, ymin, zmin} );
            ptsRef.push_back( {xmax, ymin, zmax} );
            ptsRef.push_back( {xmax, ymax, zmin} );
            ptsRef.push_back( {xmax, ymax, zmax} );
			if ( constCells )
				dynamic_cast<Grid3Ducfs<T, uint32_t>*>(g)->initOrdering( ptsRef, par.order );
			else
				dynamic_cast<Grid3Duifs<T, uint32_t>*>(g)->initOrdering( ptsRef, par.order );
            if ( par.time ) { end = std::chrono::high_resolution_clock::now(); }
            if ( par.verbose ) {
                std::cout << "done.\n";
                std::cout << "Initiated ordering with " << ptsRef.size() << " reference points and l-"
                << par.order << " metric\n";
                std::cout.flush();
            }
            
            break;
        }
        default:
            break;
    }
	if ( par.time ) {
		std::cout.precision(12);
		std::cout << "Time to build grid: " << std::chrono::duration<double>(end-begin).count() << '\n';
	}
    std::cout.flush();
	if ( par.verbose && par.method == SHORTEST_PATH ) {
		std::cout << "Interpolating slowness at secondary nodes ... ";
		std::cout.flush();
	}
	if ( par.time && par.method == SHORTEST_PATH ) {
		begin = std::chrono::high_resolution_clock::now();
	}
	g->setSlowness(slowness);
	if ( par.verbose && par.method == SHORTEST_PATH ) {
		std::cout << "done.\n";
		std::cout.flush();
	}
	if ( par.time && par.method == SHORTEST_PATH ) {
		end = std::chrono::high_resolution_clock::now();
		std::cout << "Time to interpolate slowness values: " << std::chrono::duration<double>(end-begin).count() << '\n';
	}
	
    return g;
}
#endif

template<typename T>
Grid3D<T, uint32_t> *unstruct3D(const input_parameters &par,
								std::vector<Rcv<T>> &reflectors,
								const size_t nt, const size_t ns)
{

	MSHReader reader( par.modelfile.c_str() );
    
    if ( !reader.isValid() ) {
        return nullptr;
    }
    
	if ( par.verbose ) {
		std::cout << "Reading model file " << par.modelfile << " ... ";
		std::cout.flush();
	}
	
    std::vector<sxyz<T>> nodes(reader.getNumberOfNodes());
	std::vector<tetrahedronElem<uint32_t>> tetrahedra(reader.getNumberOfElements());
	std::vector<T> slowness(reader.getNumberOfElements());
	
	reader.readNodes3D(nodes);
	reader.readTetrahedronElements(tetrahedra);
    if ( par.verbose ) std::cout << "done.\n";
	std::map<std::string, double> slownesses;
	
    bool constCells = true;
    if ( !par.slofile.empty() ) {
        
        std::ifstream fin(par.slofile.c_str());
        if ( !fin ) {
            std::cerr << "Error: cannot open file " << par.slofile << std::endl;
            exit ( -1);
        }
        std::vector<T> tmp;
        T dtmp;
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
        
    } else {
        std::ifstream fin(par.velfile.c_str());
        if ( !fin ) {
            std::cerr << "Error: cannot open file " << par.velfile << std::endl;
            exit ( -1);
        }
        std::string line;
        while ( fin ) {
            getline( fin, line );
            if ( line.empty() ) continue;
            size_t i1 = line.find('"');
            size_t i2 = line.rfind('"');
            std::string name = line.substr(i1+1, i2-i1-1);
            std::istringstream sin( line.substr(i2+1, 100) );
            double val;
            sin >> val;
            slownesses.insert( {name, 1./val} );
        }
        fin.close();
	
        if ( par.verbose ) {
            for ( size_t n=0; n<reader.getPhysicalNames(3).size(); ++n ) {
                std::cout << "  Velocity for " << reader.getPhysicalNames(3)[n] << " is "
                << 1./slownesses[ reader.getPhysicalNames(3)[n] ] << '\n';
            }
        }
	
        for ( size_t n=0; n<slowness.size(); ++n ) {
            slowness[n] = slownesses[ reader.getPhysicalNames(3)[tetrahedra[n].physical_entity] ];
        }
    }
	
    if ( par.verbose ) {
        std::cout << "  Unstructured mesh in file has"
        << "\n    " << nodes.size() << " nodes"
        << "\n    " << tetrahedra.size() << " cells";
		if ( constCells )
			std::cout << "\n  Mesh has cells of constant slowness";
		else
			std::cout << "\n  Mesh has slowness defined at nodes";
		std::cout << std::endl;
    }
	
	std::chrono::high_resolution_clock::time_point begin, end;
    Grid3D<T, uint32_t> *g = nullptr;
    switch (par.method) {
        case SHORTEST_PATH:
        {
            if ( par.verbose ) {
                std::cout << "Creating grid using " << par.nn[0] << " secondary nodes ... ";
                std::cout.flush();
            }
            if ( par.time ) { begin = std::chrono::high_resolution_clock::now(); }
			if ( constCells )
				g = new Grid3Ducsp<T, uint32_t>(nodes, tetrahedra,par.nn[0], nt,
												par.verbose);
			else
				g = new Grid3Duisp<T, uint32_t>(nodes, tetrahedra,par.nn[0], nt,
												par.verbose);
            if ( par.time ) { end = std::chrono::high_resolution_clock::now(); }
            if ( par.verbose ) {
                std::cout << "done.\nTotal number of nodes: " << g->getNumberOfNodes()
                << "\n";
                std::cout.flush();
            }
            
            break;
        }
        case FAST_MARCHING:
        {
            if ( par.verbose ) {
                std::cout << "Creating grid ... ";
                std::cout.flush();
            }
            if ( par.time ) { begin = std::chrono::high_resolution_clock::now(); }
			if ( constCells )
				g = new Grid3Ducfm<T, uint32_t>(nodes, tetrahedra, nt);
			else
				g = new Grid3Duifm<T, uint32_t>(nodes, tetrahedra, nt);
            if ( par.time ) { end = std::chrono::high_resolution_clock::now(); }
            if ( par.verbose ) {
                std::cout << "done.\n";
                std::cout.flush();
            }
            
            break;
        }
        case FAST_SWEEPING:
        {
            if ( par.verbose ) {
                std::cout << "Creating grid ... ";
                std::cout.flush();
            }
            if ( par.time ) { begin = std::chrono::high_resolution_clock::now(); }
			if ( constCells )
				g = new Grid3Ducfs<T, uint32_t>(nodes, tetrahedra, par.epsilon,
												par.nitermax, nt);
			else
				g = new Grid3Duifs<T, uint32_t>(nodes, tetrahedra, par.epsilon,
												par.nitermax, nt);
				
            T xmin = g->getXmin();
            T xmax = g->getXmax();
            T ymin = g->getYmin();
            T ymax = g->getYmax();
            T zmin = g->getZmin();
            T zmax = g->getZmax();
            
            std::vector<sxyz<T>> ptsRef;
            ptsRef.push_back( {xmin, ymin, zmin} );
            ptsRef.push_back( {xmin, ymin, zmax} );
            ptsRef.push_back( {xmin, ymax, zmin} );
            ptsRef.push_back( {xmin, ymax, zmax} );
            ptsRef.push_back( {xmax, ymin, zmin} );
            ptsRef.push_back( {xmax, ymin, zmax} );
            ptsRef.push_back( {xmax, ymax, zmin} );
            ptsRef.push_back( {xmax, ymax, zmax} );
			if ( constCells )
				dynamic_cast<Grid3Ducfs<T, uint32_t>*>(g)->initOrdering( ptsRef, par.order );
			else
				dynamic_cast<Grid3Duifs<T, uint32_t>*>(g)->initOrdering( ptsRef, par.order );
            if ( par.time ) { end = std::chrono::high_resolution_clock::now(); }
            if ( par.verbose ) {
                std::cout << "done.\n";
                std::cout << "Initiated ordering with " << ptsRef.size() << " reference points and l-"
                << par.order << " metric\n";
                std::cout.flush();
            }
            
            break;
        }
        default:
            break;
    }
	if ( par.time ) {
		std::cout.precision(12);
		std::cout << "Time to build grid: " << std::chrono::duration<double>(end-begin).count() << '\n';
	}
    std::cout.flush();
	if ( par.verbose && par.method == SHORTEST_PATH ) {
		std::cout << "Interpolating slowness at secondary nodes ... ";
		std::cout.flush();
	}
	if ( par.time && par.method == SHORTEST_PATH ) {
		begin = std::chrono::high_resolution_clock::now();
	}
	g->setSlowness(slowness);
	if ( par.verbose && par.method == SHORTEST_PATH ) {
		std::cout << "done.\n";
		std::cout.flush();
	}
	if ( par.time && par.method == SHORTEST_PATH ) {
		end = std::chrono::high_resolution_clock::now();
		std::cout << "Time to interpolate slowness values: " << std::chrono::duration<double>(end-begin).count() << '\n';
	}
	
	if ( par.processReflectors ) {
		buildReflectors(reader, nodes, ns, par.nn[0], reflectors);
	}
	
	if ( par.saveModelVTK ) {
#ifdef VTK
		std::string filename = par.modelfile;
		size_t i = filename.rfind(".msh");
		filename.replace(i, 4, ".vtu");
		
		if ( par.verbose ) std::cout << "Saving model in " << filename << " ... ";
		g->saveModelVTU(filename, false);
		if ( par.verbose ) std::cout << "done.\n";
#else
		std::cerr << "Error: program not compiled with VTK support" << std::endl;
		return nullptr;
#endif
	}
	
    return g;

}


#ifdef VTK
template<typename T>
Grid2Drc<T,uint32_t> *recti2Dc(const input_parameters &par, const size_t nt)
{
    Grid2Drc<T,uint32_t> *g = nullptr;
	vtkRectilinearGrid *dataSet;
    
    vtkSmartPointer<vtkXMLRectilinearGridReader> reader =
    vtkSmartPointer<vtkXMLRectilinearGridReader>::New();
    reader->SetFileName(par.modelfile.c_str());
    reader->Update();
    reader->GetOutput()->Register(reader);
    dataSet = reader->GetOutput();
    
    vtkIdType numberOfCells = dataSet->GetNumberOfCells();
    int nnodes[3], ncells[3];
    dataSet->GetDimensions(nnodes);
    ncells[0] = nnodes[0]-1;
    ncells[1] = nnodes[1]-1;
    ncells[2] = nnodes[2]-1;
    double xrange[2], yrange[2], zrange[2];
    dataSet->GetXCoordinates()->GetRange(xrange);
    dataSet->GetYCoordinates()->GetRange(yrange);
    dataSet->GetZCoordinates()->GetRange(zrange);
    double d[3];
    d[0] = (xrange[1]-xrange[0])/(nnodes[0]-1);
    d[1] = (yrange[1]-yrange[0])/(nnodes[1]-1);
    d[2] = (zrange[1]-zrange[0])/(nnodes[2]-1);
    
    if ( nnodes[1]>1 ) {
        std::cerr << "Error - model is not 2D\n";
        abort();
    }
	
    if ( par.verbose ) {
        std::cout << "Reading model file " << par.modelfile
        << "\n  Rectilinear grid in file has"
        << "\n    " << nnodes[0]*nnodes[1]*nnodes[2] << " nodes"
        << "\n    " << ncells[0]*ncells[1]*ncells[2] << " cells"
        << "\n    (size: " << nnodes[0] << " x " << nnodes[1] << " x " << nnodes[2] << ')'
        << "\n  Dim\tmin\tmax\tinc.\t N. sec nodes"
        << "\n   X\t" << xrange[0] << '\t' << xrange[1] << '\t' << d[0] << "\t\t" << par.nn[0]
        << "\n   Y\t" << yrange[0] << '\t' << yrange[1] << '\t' << d[1] << "\t\t" << par.nn[1]
        << "\n   Z\t" << zrange[0] << '\t' << zrange[1] << '\t' << d[2] << "\t\t" << par.nn[2]
        << std::endl;
    }
    vtkCellData *cd = dataSet->GetCellData();
    
	std::chrono::high_resolution_clock::time_point begin, end;
	
    bool foundSlowness = false;
	std::vector<T> slowness;
    if ( cd->HasArray("P-wave velocity") || cd->HasArray("Velocity") ||
		cd->HasArray("Slowness") ) {
        
        for (int na = 0; na < cd->GetNumberOfArrays(); na++) {
            if ( strcmp(cd->GetArrayName(na), "P-wave velocity")==0 ||
				strcmp(cd->GetArrayName(na), "Velocity")==0 ) {
                slowness.resize( numberOfCells );
                for ( size_t k=0,n=0; k<ncells[2]; ++k ) {
                    for ( size_t j=0; j<ncells[1]; ++j ) {
                        for ( size_t i=0; i<ncells[0]; ++i,++n ) {
                            slowness[n] = static_cast<T>(1./cd->GetArray(na)->GetTuple1(n));
                        }
                    }
                }
                foundSlowness = true;
				break;
			} else if ( strcmp(cd->GetArrayName(na), "Slowness")==0 ) {
                
				vtkSmartPointer<vtkDoubleArray> slo = vtkSmartPointer<vtkDoubleArray>::New();
				slo = vtkDoubleArray::SafeDownCast( cd->GetArray("Slowness") );
				
				if ( slo->GetSize() != dataSet->GetNumberOfCells() ) {
					std::cout << "Problem with Slowness data (wrong size)" << std::endl;
					return nullptr;
				}
				
				slowness.resize( slo->GetSize() );
				for ( size_t n=0; n<slo->GetSize(); ++n ) {
					slowness[n] = slo->GetComponent(n, 0);
				}
				foundSlowness = true;
				break;
			}
		}
		if ( foundSlowness ) {
			if ( par.verbose ) { cout << "Building grid (Grid2Drc) ... "; cout.flush(); }
			if ( par.time ) { begin = std::chrono::high_resolution_clock::now(); }
			g = new Grid2Drc<T, uint32_t>(ncells[0], ncells[2], d[0], d[2],
										  xrange[0], zrange[0],
                                          par.nn[0], par.nn[2], nt);
			if ( par.time ) { end = std::chrono::high_resolution_clock::now(); }
			if ( par.verbose ) {
                cout << "done.\nTotal number of nodes: " << g->getNumberOfNodes()
                << "\nAssigning slowness at grid cells ... ";
            }
			g->setSlowness( slowness );
			if ( par.verbose ) cout << "done.\n";
			if ( par.time ) {
				std::cout.precision(12);
				std::cout << "Time to build grid: " << std::chrono::duration<double>(end-begin).count() << '\n';
			}
		} else {
			return nullptr;
		}
    }
    dataSet->Delete();
    return g;
};

template<typename T>
Grid2D<T,uint32_t,sxz<T>> *unstruct2D_vtu(const input_parameters &par, const size_t nt)
{
    VTUReader reader( par.modelfile.c_str() );
    
    if ( !reader.isValid() ) {
        return nullptr;
    }
    
	if ( par.verbose ) {
		cout << "Reading model file " << par.modelfile << " ... ";
		cout.flush();
	}
	
    std::vector<sxz<T>> nodes(reader.getNumberOfNodes());
	std::vector<triangleElem<uint32_t>> triangles(reader.getNumberOfElements());
    bool constCells = reader.isConstCell();
    
	std::vector<T> slowness;
    if ( constCells )
        slowness.resize(reader.getNumberOfElements());
    else
        slowness.resize(reader.getNumberOfNodes());
	
	reader.readNodes2D(nodes);
	reader.readTriangleElements(triangles);
    reader.readSlowness(slowness, constCells);
    
    
    if ( par.verbose ) {
        cout << "  Unstructured mesh in file has"
        << "\n    " << nodes.size() << " nodes"
        << "\n    " << triangles.size() << " cells";
		if ( constCells )
			std::cout << "\n  Mesh has cells of constant slowness";
		else
			std::cout << "\n  Mesh has slowness defined at nodes";
		std::cout << std::endl;
    }
    
	std::chrono::high_resolution_clock::time_point begin, end;
    Grid2D<T, uint32_t,sxz<T>> *g=nullptr;
    switch (par.method) {
        case SHORTEST_PATH:
        {
			if ( par.verbose ) {
				cout << "Creating grid using " << par.nn[0] << " secondary nodes ... ";
				cout.flush();
			}
			if ( par.time ) { begin = std::chrono::high_resolution_clock::now(); }
			if ( constCells )
				g = new Grid2Ducsp<T, uint32_t, Node2Dcsp<T,uint32_t>,sxz<T>>(nodes, triangles, par.nn[0], nt);
			else
				g = new Grid2Duisp<T, uint32_t, Node2Disp<T,uint32_t>,sxz<T>>(nodes, triangles, par.nn[0], nt);
			if ( par.time ) { end = std::chrono::high_resolution_clock::now(); }
			if ( par.verbose ) {
				cout << "done.\nTotal number of nodes: " << g->getNumberOfNodes()
				<< "\n";
				cout.flush();
			}
			
            break;
        }
        case FAST_MARCHING:
        {
            if ( par.verbose ) {
                std::cout << "Creating grid ... ";
                std::cout.flush();
            }
            if ( par.time ) { begin = std::chrono::high_resolution_clock::now(); }
			if ( constCells )
				g = new Grid2Ducfm<T, uint32_t, Node2Dcsp<T,uint32_t>,sxz<T>>(nodes, triangles, nt);
			else
				g = new Grid2Duifm<T, uint32_t, Node2Disp<T,uint32_t>,sxz<T>>(nodes, triangles, nt);
            if ( par.time ) { end = std::chrono::high_resolution_clock::now(); }
            if ( par.verbose ) {
                std::cout << "done.\n";
                std::cout.flush();
            }
			
            break;
        }
        case FAST_SWEEPING:
        {
            if ( par.verbose ) {
                std::cout << "Creating grid ... ";
                std::cout.flush();
            }
            if ( par.time ) { begin = std::chrono::high_resolution_clock::now(); }
			if ( constCells )
				g = new Grid2Ducfs<T, uint32_t, Node2Dcsp<T,uint32_t>,sxz<T>>(nodes, triangles, par.epsilon,
																			  par.nitermax, nt);
			else
				g = new Grid2Duifs<T, uint32_t, Node2Disp<T,uint32_t>,sxz<T>>(nodes, triangles, par.epsilon,
																			  par.nitermax, nt);
            T xmin = g->getXmin();
            T xmax = g->getXmax();
            T zmin = g->getZmin();
            T zmax = g->getZmax();
            
            std::vector<sxz<T>> ptsRef;
            ptsRef.push_back( {xmin, zmin} );
            ptsRef.push_back( {xmin, zmax} );
            ptsRef.push_back( {xmax, zmin} );
            ptsRef.push_back( {xmax, zmax} );
			if ( constCells )
				dynamic_cast<Grid2Ducfs<T, uint32_t, Node2Dcsp<T,uint32_t>,sxz<T>>*>(g)->initOrdering( ptsRef, par.order );
			else
				dynamic_cast<Grid2Duifs<T, uint32_t, Node2Disp<T,uint32_t>,sxz<T>>*>(g)->initOrdering( ptsRef, par.order );
            if ( par.time ) { end = std::chrono::high_resolution_clock::now(); }
            if ( par.verbose ) {
                std::cout << "done.\n";
                std::cout << "Initiated ordering with " << ptsRef.size() << " reference points and l-"
                << par.order << " metric\n";
                std::cout.flush();
            }
            
            break;
        }
        default:
            break;
    }
	if ( par.time ) {
		cout.precision(12);
		cout << "Time to build grid: " << std::chrono::duration<double>(end-begin).count() << '\n';
	}
    cout.flush();
	if ( 1 == (g->setSlowness(slowness)) ) {
		delete g;
		return nullptr;
	}
    return g;
}
#endif

template<typename T>
Grid2D<T,uint32_t,sxz<T>> *unstruct2D(const input_parameters &par,
									  std::vector<Rcv2D<T>> &reflectors,
									  const size_t nt, const size_t ns)
{
    
    MSHReader reader( par.modelfile.c_str() );
    
    if ( !reader.isValid() ) {
        return nullptr;
    }
    
	if ( par.verbose ) {
		std::cout << "Reading model file " << par.modelfile << " ... ";
		std::cout.flush();
	}
	
    std::vector<sxz<T>> nodes(reader.getNumberOfNodes());
	std::vector<triangleElem<uint32_t>> triangles(reader.getNumberOfTriangles());
	std::vector<T> slowness(reader.getNumberOfTriangles());
	
	reader.readNodes2D(nodes);
	reader.readTriangleElements(triangles);
    if ( par.verbose ) std::cout << "done.\n";
	std::map<std::string, double> slownesses;

	bool constCells = true;
    if ( !par.slofile.empty() ) {
        
        std::ifstream fin(par.slofile.c_str());
        if ( !fin ) {
            std::cout << "Error: cannot open file " << par.slofile << std::endl;
            exit ( -1);
        }
        std::vector<T> tmp;
        T dtmp;
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
        
    } else {
        std::ifstream fin(par.velfile.c_str());
        if ( !fin ) {
            std::cerr << "Error: cannot open file " << par.velfile << std::endl;
            exit ( -1);
        }
        std::string line;
        while ( fin ) {
            getline( fin, line );
            if ( line.empty() ) continue;
            size_t i1 = line.find('"');
            size_t i2 = line.rfind('"');
            std::string name = line.substr(i1+1, i2-i1-1);
            std::istringstream sin( line.substr(i2+1, 100) );
            double val;
            sin >> val;
            slownesses.insert( {name, 1./val} );
        }
        fin.close();
	
        if ( par.verbose ) {
            for ( size_t n=0; n<reader.getPhysicalNames(2).size(); ++n ) {
                std::cout << "  Velocity for " << reader.getPhysicalNames(2)[n] << " is "
                << 1./slownesses[ reader.getPhysicalNames(2)[n] ] << '\n';
            }
        }

        for ( size_t n=0; n<slowness.size(); ++n ) {
            slowness[n] = slownesses[ reader.getPhysicalNames(2)[triangles[n].physical_entity] ];
        }
    }

    
    if ( par.verbose ) {
        std::cout << "  Unstructured mesh in file has"
        << "\n    " << nodes.size() << " nodes"
        << "\n    " << triangles.size() << " cells";
		if ( constCells )
			std::cout << "\n  Mesh has cells of constant slowness";
		else
			std::cout << "\n  Mesh has slowness defined at nodes";
		std::cout << std::endl;
    }
    
	std::chrono::high_resolution_clock::time_point begin, end;
    Grid2D<T,uint32_t,sxz<T>> *g=nullptr;
    switch (par.method) {
        case SHORTEST_PATH:
        {
			if ( par.verbose ) {
				cout << "Creating grid using " << par.nn[0] << " secondary nodes ... ";
				cout.flush();
			}
			if ( par.time ) { begin = std::chrono::high_resolution_clock::now(); }
			if ( constCells )
				g = new Grid2Ducsp<T, uint32_t, Node2Dcsp<T,uint32_t>,sxz<T>>(nodes, triangles, par.nn[0], nt);
			else
				g = new Grid2Duisp<T, uint32_t, Node2Disp<T,uint32_t>,sxz<T>>(nodes, triangles, par.nn[0], nt);
			if ( par.time ) { end = std::chrono::high_resolution_clock::now(); }
			if ( par.verbose ) {
				cout << "done.\nTotal number of nodes: " << g->getNumberOfNodes()
				<< "\n";
				cout.flush();
			}
			
            break;
        }
        case FAST_MARCHING:
        {
            if ( par.verbose ) {
                std::cout << "Creating grid ... ";
                std::cout.flush();
            }
            if ( par.time ) { begin = std::chrono::high_resolution_clock::now(); }
			if ( constCells )
				g = new Grid2Ducfm<T, uint32_t, Node2Dcsp<T,uint32_t>,sxz<T>>(nodes, triangles, nt);
			else
				g = new Grid2Duifm<T, uint32_t, Node2Disp<T,uint32_t>,sxz<T>>(nodes, triangles, nt);
            if ( par.time ) { end = std::chrono::high_resolution_clock::now(); }
            if ( par.verbose ) {
                std::cout << "done.\n";
                std::cout.flush();
            }
			
            break;
        }
        case FAST_SWEEPING:
        {
            if ( par.verbose ) {
                std::cout << "Creating grid ... ";
                std::cout.flush();
            }
            if ( par.time ) { begin = std::chrono::high_resolution_clock::now(); }
			if ( constCells )
				g = new Grid2Ducfs<T, uint32_t, Node2Dcsp<T,uint32_t>,sxz<T>>(nodes, triangles, par.epsilon,
																			  par.nitermax, nt);
			else
				g = new Grid2Duifs<T, uint32_t, Node2Disp<T,uint32_t>,sxz<T>>(nodes, triangles, par.epsilon,
																			  par.nitermax, nt);

            T xmin = g->getXmin();
            T xmax = g->getXmax();
            T zmin = g->getZmin();
            T zmax = g->getZmax();
            
            std::vector<sxz<T>> ptsRef;
            ptsRef.push_back( {xmin, zmin} );
            ptsRef.push_back( {xmin, zmax} );
            ptsRef.push_back( {xmax, zmin} );
            ptsRef.push_back( {xmax, zmax} );
			if ( constCells )
				dynamic_cast<Grid2Ducfs<T, uint32_t, Node2Dcsp<T,uint32_t>,sxz<T>>*>(g)->initOrdering( ptsRef, par.order );
			else
				dynamic_cast<Grid2Duifs<T, uint32_t, Node2Disp<T,uint32_t>,sxz<T>>*>(g)->initOrdering( ptsRef, par.order );
			if ( par.time ) { end = std::chrono::high_resolution_clock::now(); }
            if ( par.verbose ) {
                std::cout << "done.\n";
                std::cout << "Initiated ordering with " << ptsRef.size() << " reference points and l-"
                << par.order << " metric\n";
                std::cout.flush();
            }
            
            break;
        }
        default:
            break;
    }
	if ( par.time ) {
		std::cout.precision(12);
		std::cout << "Time to build grid: " << std::chrono::duration<double>(end-begin).count() << '\n';
	}
    std::cout.flush();
	g->setSlowness(slowness);
    
	if ( par.processReflectors ) {
		std::vector<std::string> reflector_names = reader.getPhysicalNames(1);
		std::vector<int> indices = reader.getPhysicalIndices(1);
	
		if ( reflector_names.size() != indices.size() ) {
			std::cerr << "Error - definition of reflectors\n";
			exit(1);
		}

		std::vector<lineElem<uint32_t>> lines;
		reader.readLineElements(lines);

		sxz<T> pt1, pt2, pt3, d;
		int nsecondary = par.nn[0];
	
		for ( size_t ni=0; ni<indices.size(); ++ni ) {
		
			reflectors.push_back( reflector_names[ni] );
		
			std::set<sxz<T>> refl_pts;  // use set to avoid duplicate points
			typename std::set<sxz<T>>::iterator it;
		
			for ( size_t nl=0; nl<lines.size(); ++nl ) {
				if ( indices[ni] == lines[nl].physical_entity ) {
					pt1 = nodes[ lines[nl].i[0] ];
					pt2 = nodes[ lines[nl].i[1] ];
				
					d.x = (pt2.x-pt1.x)/(nsecondary+1);
					d.z = (pt2.z-pt1.z)/(nsecondary+1);

					refl_pts.insert( pt1 );
					for ( size_t n2=0; n2<nsecondary; ++n2 ) {
						pt3.x = pt1.x+(1+n2)*d.x;
						pt3.z = pt1.z+(1+n2)*d.z;
						refl_pts.insert( pt3 );
					}
					refl_pts.insert( pt2 );
				}
			}
		
			for (it=refl_pts.begin(); it!=refl_pts.end(); ++it) {
				reflectors.back().add_coord( *it );
			}
			reflectors.back().init_tt( ns );
		}
	}
	if ( par.saveModelVTK ) {
#ifdef VTK
		std::string filename = par.modelfile;
		size_t i = filename.rfind(".msh");
		filename.replace(i, 4, ".vtu");
		
		if ( par.verbose ) std::cout << "Saving model in " << filename << " ... ";
		g->saveModelVTU(filename, false);
		if ( par.verbose ) std::cout << "done.\n";
#else
		std::cerr << "Error: Program not compiled with VTK support" << std::endl;
		return nullptr;
#endif
	}
	
    return g;
}

template<typename T>
Grid2D<T, uint32_t, sxyz<T>> *unstruct2Ds_vtu(const input_parameters &par, const size_t nt)
{
    VTUReader reader( par.modelfile.c_str() );
    
    if ( !reader.isValid() ) {
        return nullptr;
    }
    
	if ( par.verbose ) {
		cout << "Reading model file " << par.modelfile << " ... ";
		cout.flush();
	}
	
    std::vector<sxyz<T>> nodes(reader.getNumberOfNodes());
	std::vector<triangleElem<uint32_t>> triangles(reader.getNumberOfElements());
    
    bool constCells = reader.isConstCell();
    
	std::vector<T> slowness;
    if ( constCells )
        slowness.resize(reader.getNumberOfElements());
    else
        slowness.resize(reader.getNumberOfNodes());
        
	reader.readNodes3D(nodes);
	reader.readTriangleElements(triangles);
    reader.readSlowness(slowness, constCells);
    	
    if ( par.verbose ) {
        cout << "  Unstructured mesh in file has"
        << "\n    " << nodes.size() << " nodes"
        << "\n    " << triangles.size() << " cells";
		if ( constCells )
			std::cout << "\n  Mesh has cells of constant slowness";
		else
			std::cout << "\n  Mesh has slowness defined at nodes";
		std::cout << std::endl;
    }
    
	std::chrono::high_resolution_clock::time_point begin, end;
    Grid2D<T, uint32_t, sxyz<T>> *g=nullptr;
	
	if ( par.verbose ) {
		cout << "Creating grid using " << par.nn[0] << " secondary nodes ... ";
		cout.flush();
	}
	if ( par.time ) { begin = std::chrono::high_resolution_clock::now(); }
    if ( constCells )
        g = new Grid2Ducsp<T, uint32_t, Node3Dcsp<T,uint32_t>, sxyz<T>>(nodes, triangles, par.nn[0], nt);
    else
        g = new Grid2Duisp<T, uint32_t, Node3Disp<T,uint32_t>, sxyz<T>>(nodes, triangles, par.nn[0], nt);
	if ( par.time ) { end = std::chrono::high_resolution_clock::now(); }
	if ( par.verbose ) {
		cout << "done.\nTotal number of nodes: " << g->getNumberOfNodes()
		<< "\n";
		cout.flush();
	}
			
	if ( par.time ) {
		cout.precision(12);
		cout << "Time to build grid: " << std::chrono::duration<double>(end-begin).count() << '\n';
	}
    cout.flush();
	if ( 1 == (g->setSlowness(slowness)) ) {
		delete g;
		return nullptr;
	}
    return g;

}



template<typename T>
Grid2D<T, uint32_t, sxyz<T>> *unstruct2Ds(const input_parameters &par,
                                          const size_t nt, const size_t ns)
{
    
    MSHReader reader( par.modelfile.c_str() );
    
    if ( !reader.isValid() ) {
        return nullptr;
    }
    
	if ( par.verbose ) {
		std::cout << "Reading model file " << par.modelfile << " ... ";
		std::cout.flush();
	}
	
    std::vector<sxyz<T>> nodes(reader.getNumberOfNodes());
	std::vector<triangleElem<uint32_t>> triangles(reader.getNumberOfTriangles());
	std::vector<T> slowness(reader.getNumberOfTriangles());
	
	reader.readNodes3D(nodes);
	reader.readTriangleElements(triangles);
    if ( par.verbose ) std::cout << "done.\n";
	std::map<std::string, double> slownesses;
	
    bool constCells = true;
    if ( !par.slofile.empty() ) {
        
        std::ifstream fin(par.slofile.c_str());
        if ( !fin ) {
            std::cout << "Error: cannot open file " << par.slofile << std::endl;
            exit ( -1);
        }
        std::vector<T> tmp;
        T dtmp;
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
        
    } else {
        std::ifstream fin(par.velfile.c_str());
        if ( !fin ) {
            std::cerr << "Error: cannot open file " << par.velfile << std::endl;
            exit ( -1);
        }
        std::string line;
        while ( fin ) {
            getline( fin, line );
            if ( line.empty() ) continue;
            size_t i1 = line.find('"');
            size_t i2 = line.rfind('"');
            std::string name = line.substr(i1+1, i2-i1-1);
            std::istringstream sin( line.substr(i2+1, 100) );
            double val;
            sin >> val;
            slownesses.insert( {name, 1./val} );
        }
        fin.close();
		
        if ( par.verbose ) {
            for ( size_t n=0; n<reader.getPhysicalNames(2).size(); ++n ) {
                std::cout << "  Velocity for " << reader.getPhysicalNames(2)[n] << " is "
                << 1./slownesses[ reader.getPhysicalNames(2)[n] ] << '\n';
            }
        }
		
        for ( size_t n=0; n<slowness.size(); ++n ) {
            slowness[n] = slownesses[ reader.getPhysicalNames(2)[triangles[n].physical_entity] ];
        }
    }
	    
    if ( par.verbose ) {
        std::cout << "  Unstructured mesh in file has"
        << "\n    " << nodes.size() << " nodes"
        << "\n    " << triangles.size() << " cells";
		if ( constCells )
			std::cout << "\n  Mesh has cells of constant slowness";
		else
			std::cout << "\n  Mesh has slowness defined at nodes";
		std::cout << std::endl;
    }
    
	std::chrono::high_resolution_clock::time_point begin, end;
    Grid2D<T, uint32_t, sxyz<T>> *g=nullptr;
	
	if ( par.verbose ) {
		cout << "Creating grid using " << par.nn[0] << " secondary nodes ... ";
		cout.flush();
	}
	if ( par.time ) { begin = std::chrono::high_resolution_clock::now(); }
    if ( constCells )
        g = new Grid2Ducsp<T, uint32_t, Node3Dcsp<T,uint32_t>, sxyz<T>>(nodes, triangles, par.nn[0], nt);
    else
        g = new Grid2Duisp<T, uint32_t, Node3Disp<T,uint32_t>, sxyz<T>>(nodes, triangles, par.nn[0], nt);
	if ( par.time ) { end = std::chrono::high_resolution_clock::now(); }
	if ( par.verbose ) {
		cout << "done.\nTotal number of nodes: " << g->getNumberOfNodes()
		<< "\n";
		cout.flush();
	}
			
	if ( par.time ) {
		std::cout.precision(12);
		std::cout << "Time to build grid: " << std::chrono::duration<double>(end-begin).count() << '\n';
	}
    std::cout.flush();
	g->setSlowness(slowness);
    
	if ( par.saveModelVTK ) {
#ifdef VTK
		std::string filename = par.modelfile;
		size_t i = filename.rfind(".msh");
		filename.replace(i, 4, ".vtu");
		
		if ( par.verbose ) std::cout << "Saving model in " << filename << " ... ";
		g->saveModelVTU(filename, false);
		if ( par.verbose ) std::cout << "done.\n";
#else
		std::cerr << "Error: Program not compiled with VTK support" << std::endl;
		return nullptr;
#endif
	}
	
    return g;
}

#endif
