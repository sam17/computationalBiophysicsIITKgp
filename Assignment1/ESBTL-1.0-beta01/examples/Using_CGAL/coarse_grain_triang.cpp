// Copyright (c) 2009-2010  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
//This file is part of ESBTL.
//
//ESBTL is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//ESBTL is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with ESBTL.  If not, see <http://www.gnu.org/licenses/>.
//
//
//Additional permission under GNU GPL version 3 section 7
//
//If you modify this Library, or any covered work, by linking or
//combining it with CGAL (or a modified version of that library), the
//licensors of this Library grant you additional permission to convey
//the resulting work. Corresponding Source for a non-source form of
//such a combination shall include the source code for the parts of CGAL
//used as well as that of the covered work. 
//
//
//
// Author(s)     :  Sébastien Loriot



#include <iostream>
#include <cstring>
#include <list>
#include <vector>


#include <ESBTL/constants.h>
#include <ESBTL/coarse_grain.h>
#include <ESBTL/coarse_creators.h>
#include <ESBTL/PDB.h>
#include <ESBTL/line_selectors.h>
#include <ESBTL/builder.h>
#include <ESBTL/line_reader.h>
#include <ESBTL/coarse_classifier.h>
#include <ESBTL/selected_atom_iterator.h>
#include <ESBTL/occupancy_handlers.h>
#include <ESBTL/global_functions.h>
#include <ESBTL/CGAL/EPIC_kernel_with_atom.h>

//alternative for selection
#include <boost/iterator/filter_iterator.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>

typedef ESBTL::CGAL::EPIC_kernel_with_coarse_atom                       Kernel;
typedef ESBTL::CGAL::System_with_coarse_grain                           My_system;
typedef ESBTL::Accept_none_occupancy_policy<ESBTL::PDB::Line_format<> > Accept_none_occupancy_policy;

typedef CGAL::Delaunay_triangulation_3<Kernel>                       Delaunay;

int main(int argc, char** argv){
  
  if (argc != 2 ){
    std::cerr << "Please provide a filename"  << std::endl;  
    return EXIT_FAILURE;
  }

  ESBTL::PDB_line_selector_two_systems sel;
  
  std::vector<My_system> systems;
  //indicate the number of system to be created (+1 for solvant manually added)
  unsigned nb_systems=sel.max_nb_systems()+1;  
  
  ESBTL::All_atom_system_builder<My_system> builder(systems,nb_systems);
  
  ESBTL::Coarse_creator_closest_to_barycenter<My_system::Residue,Kernel::FT> creator;
  
  if (ESBTL::read_a_pdb_file(argv[1],sel,builder,Accept_none_occupancy_policy())){
    
    assert(systems.size()==3);
    //generate fake solvant
    std::vector<Kernel::Point_3> solvant;
    solvant.reserve(1000);
    for (double i=0;i<10;++i)
      for (double j=0;j<10;++j)
        for (double k=0;k<10;++k)
          solvant.push_back(Kernel::Point_3(i,j,k));
    insert_coarse_atoms(solvant.begin(),solvant.end(),systems[2],1);
    
    
    for (My_system::Models_iterator it_model=systems[0].models_begin();it_model!=systems[0].models_end();++it_model){
      unsigned nb_atom=0;
      unsigned nb_residue=0;
      unsigned nb_coarse_created=0;
      {
        My_system::Model& model=*it_model;
        
        //create coarse atoms        
        for (My_system::Model::Residues_iterator it_res=model.residues_begin();it_res!=model.residues_end();++it_res){
          nb_coarse_created+=it_res->create_coarse_atoms(creator);
          ++nb_residue;
        }
      }
      
      
      
      const My_system::Model& model=*it_model;

      typedef ESBTL::Generic_classifier<ESBTL::Radius_of_coarse_atom<double,My_system::Residue::Coarse_atom> > T_Coarse_atom_classifier;
      T_Coarse_atom_classifier classifier;
      
      for (My_system::Model::Atoms_const_iterator it_atm=model.atoms_begin();it_atm!=model.atoms_end();++it_atm){
        ++nb_atom;
      }
      std::cout << "nb atoms " << nb_atom << std::endl;
      std::cout << "Nb res "<< nb_residue << std::endl;
      
      typedef ESBTL::Coarse_atoms_iterators<My_system::Model>::const_iterator Coarse_atoms_iterator;
      
      unsigned nb_coarse=0;
      for (Coarse_atoms_iterator itc=coarse_atoms_begin(model);itc!=coarse_atoms_end(model);++itc){
        ++nb_coarse;
        Kernel::FT radius=classifier.get_properties(*itc).value();
        assert(radius > 0);
      }
      const My_system::Model& solvant=*systems[2].models_begin();
      unsigned nb_fake_sol=0;
      for (Coarse_atoms_iterator itc=coarse_atoms_begin(solvant);itc!=coarse_atoms_end(solvant);++itc){
        ++nb_fake_sol;
      }
      assert( nb_coarse_created==nb_coarse );
      std::cout << "nb coarse " << nb_coarse << std::endl;
      std::cout << "nb_fake_sol " << nb_fake_sol << std::endl;
      
      Delaunay dt(coarse_atoms_begin(model),coarse_atoms_end(model));
      dt.insert(coarse_atoms_begin(solvant),coarse_atoms_end(solvant));
      
      for (Delaunay::Finite_vertices_iterator it=dt.finite_vertices_begin();it!=dt.finite_vertices_end();++it)
        assert(it->point().index()==0);
      
      std::cout << "Delaunay " << dt.number_of_vertices()<< std::endl;
      
      write_to_cgo("/tmp/cgo.py",coarse_atoms_begin(model),coarse_atoms_end(model),classifier,ESBTL::Color_of_atom<My_system::Residue::Coarse_atom>());
      std::cout << "CGO file created for first model of system 0 in /tmp/cgo.py\n";
    }
  }
  else
    return EXIT_FAILURE;
  
  return EXIT_SUCCESS;
}

