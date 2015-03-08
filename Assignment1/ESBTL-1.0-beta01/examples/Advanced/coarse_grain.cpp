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
#include <stdlib.h>
#include <cstring>
#include <list>
#include <vector>

#include <ESBTL/default.h>
#include <ESBTL/coarse_grain.h>
#include <ESBTL/coarse_creators.h>
#include <ESBTL/global_functions.h>
#include <ESBTL/selected_atom_iterator.h>
#include <ESBTL/atom_selectors.h>

//alternative for selection
#include <boost/iterator/filter_iterator.hpp>

struct Residue_filter_ALA{
  std::string name;
  Residue_filter_ALA():name("ALA"){}
  template <class Atom>
  bool operator()(const Atom& atom){
    return atom.residue().residue_name()==name;
  }
};

typedef ESBTL::Default_system_with_coarse_grain System;

typedef ESBTL::Accept_none_occupancy_policy<ESBTL::PDB::Line_format<> > Accept_none_occupancy_policy;


int main(int argc, char** argv){
  
  if (argc != 2 ){
    std::cerr << "Please provide a filename"  << std::endl;  
    return EXIT_FAILURE;
  }

  ESBTL::PDB_line_selector_two_systems sel;
  
  std::vector<System> systems;
  //indicate the number of system to be created (+1 for the solvant manually added)
  unsigned nb_systems=sel.max_nb_systems() + 1;

  
  ESBTL::All_atom_system_builder<System> builder(systems,nb_systems);
  
  //class responsible for building a coarse grain model
  ESBTL::Coarse_creator_two_barycenters<System::Residue> creator;
  
  
  if (ESBTL::read_a_pdb_file(argv[1],sel,builder,Accept_none_occupancy_policy())){
    assert(systems.size()==3);
    //generate fake solvant
    std::vector<ESBTL::Point_3> solvant;
    solvant.reserve(1000);
    for (double i=0;i<10;++i)
      for (double j=0;j<10;++j)
        for (double k=0;k<10;++k)
          solvant.push_back(ESBTL::Point_3(i,j,k));
    //insert this fake solvant into systems 2.
    insert_coarse_atoms(solvant.begin(),solvant.end(),systems[2],1);
    
    
    for (System::Models_iterator it_model=systems[0].models_begin();it_model!=systems[0].models_end();++it_model){
      unsigned nb_atom=0;
      unsigned nb_residue=0;
      unsigned nb_coarse_created=0;
      {
        System::Model& model=*it_model;

        //create coarse atoms        
        for (System::Model::Residues_iterator it_res=model.residues_begin();it_res!=model.residues_end();++it_res){
          nb_coarse_created+=it_res->create_coarse_atoms(creator);
          ++nb_residue;
        }
      }
      
      
      
      const System::Model& model=*it_model;

      for (System::Model::Atoms_const_iterator it_atm=model.atoms_begin();it_atm!=model.atoms_end();++it_atm){
        ++nb_atom;
      }
      std::cout << "nb atoms " << nb_atom << std::endl;
      std::cout << "Nb res "<< nb_residue << std::endl;
      
      typedef ESBTL::Coarse_atoms_iterators<System::Model>::const_iterator Coarse_atoms_iterator;
      
      unsigned nb_coarse=0;
      //~ std::cout << "from pymol.cgo import *\n"; //commented for development
      //~ std::cout << "from pymol import cmd\n";   //commented for development
      //~ std::cout << "obj=[COLOR, 1.0, 0, 0\n";   //commented for development
      for (Coarse_atoms_iterator itc=coarse_atoms_begin(model);itc!=coarse_atoms_end(model);++itc){
        ++nb_coarse;
        //~ std::cout << ",SPHERE, " <<itc->x() << ", " << itc->y() << ", " << itc->z() << ", " << 3 << std::endl; //commented for development
      }
      const System::Model& solvant=*systems[2].models_begin();
      unsigned nb_fake_sol=0;
      for (Coarse_atoms_iterator itc=coarse_atoms_begin(solvant);itc!=coarse_atoms_end(solvant);++itc){
        ++nb_fake_sol;
        //~ std::cout << ",SPHERE, " <<itc->x() << ", " << itc->y() << ", " << itc->z() << ", " << 1 << std::endl; //commented for development
      }
      //~ std::cout << "]\ncmd.load_cgo(obj,\'cgo01\')\n"; //commented for development
      assert( nb_coarse_created==nb_coarse );
      std::cout << "nb coarse " << nb_coarse << std::endl;
      std::cout << "nb_fake_sol " << nb_fake_sol << std::endl;
      
      //first possibility to iterate over alanine atoms
      typedef ESBTL::Selected_atom_iterator<System::Model,ESBTL::Select_by_resname,true>  Restrict_iterator;
      unsigned nb_ala1=0;
      for (Restrict_iterator itr=ESBTL::make_selected_atom_iterator(model.atoms_begin(),ESBTL::Select_by_resname("ALA"));
           itr!=ESBTL::make_selected_atom_iterator<ESBTL::Select_by_resname>(model.atoms_end());++itr
          )
            ++nb_ala1;
      std::cout << "nb ALA atoms " << nb_ala1 << std::endl;
      
      //second possibility to iterate over alanine atoms
      typedef boost::filter_iterator<Residue_filter_ALA,System::Model::Atoms_const_iterator> Filter_iterator;
      unsigned nb_ala2=0;
      for (Filter_iterator itr = boost::make_filter_iterator<Residue_filter_ALA>(model.atoms_begin(),model.atoms_end());
                           itr!= boost::make_filter_iterator<Residue_filter_ALA>(model.atoms_end(),model.atoms_end());
                        ++itr
          )
            ++nb_ala2;
      assert(nb_ala1==nb_ala2);
            
    }
  }
  else
    return EXIT_FAILURE;
  
  return EXIT_SUCCESS;
}

