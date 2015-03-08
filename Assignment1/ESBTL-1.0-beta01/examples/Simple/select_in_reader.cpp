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



#include <ESBTL/default.h>
#include <ESBTL/selected_atom_iterator.h>
#include <ESBTL/atom_selectors.h>

//discard atoms with no alternate location identification and with occupancy != 1
typedef ESBTL::Accept_none_occupancy_policy<ESBTL::PDB::Line_format<> > Accept_none_occupancy_policy;

int main(int argc, char** argv){
  
  if (argc != 3 ){
    std::cerr << "Please provide a filename and a residue name"  << std::endl;  
    return EXIT_FAILURE;
  }
  
  //creates one system with all atoms matching the criteria
  std::string resname(argv[2]);
  ESBTL::Select_by_resname selector(resname);
  ESBTL::Generic_line_selector<ESBTL::Select_by_resname> sel(selector);
  
  std::vector<ESBTL::Default_system> systems;
  //declare the building that will contruct the system from the pdb file.
  
  std::cout << sel.max_nb_systems() << std::endl;
  
  ESBTL::All_atom_system_builder<ESBTL::Default_system> builder(systems,sel.max_nb_systems());
  
  //read the pdb file
  if (ESBTL::read_a_pdb_file(argv[1],sel,builder,Accept_none_occupancy_policy())){
    if ( systems.empty() || systems[0].has_no_model() ){
      std::cerr << "No atoms found" << std::endl;
      return EXIT_FAILURE;
    }
    //consider the first model of the system
    const ESBTL::Default_system::Model& model=*systems[0].models_begin();
    
    std::cout << model.number_of_atoms() << " atoms in residue(s) named " << resname << std::endl;
  }
  else
    return EXIT_FAILURE;
  
  return EXIT_SUCCESS;
}

