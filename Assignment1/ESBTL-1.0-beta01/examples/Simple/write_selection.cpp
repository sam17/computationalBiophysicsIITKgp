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
#include <ESBTL/atom_classifier.h>
#include <ESBTL/selected_atom_iterator.h>

#include <iostream>

typedef ESBTL::Accept_none_occupancy_policy<ESBTL::PDB::Line_format<> > Accept_none_occupancy_policy;

namespace ESBTL{
  std::ostream&
  operator<<(std::ostream& os, const ESBTL::Default_system::Atom& atm)
  {
    os << ESBTL::PDB::get_atom_pdb_format(atm);
    return os;
  }
}

int main(int argc, char** argv){
  
  if (argc != 3 ){
    std::cerr << "Please provide a filename and the chains to select"  << std::endl;  
    return EXIT_FAILURE;
  }
  
  std::list<std::string> chains_to_select;
  std::string chains(argv[2]);
  chains_to_select.push_back(chains);
  
  //discard hydrogens and water
  ESBTL::PDB_line_selector_chain sel(chains_to_select.begin(),chains_to_select.end(),false,false);
  
  std::vector<ESBTL::Default_system> systems;
  
  ESBTL::All_atom_system_builder<ESBTL::Default_system> builder(systems,sel.max_nb_systems());

  if (ESBTL::read_a_pdb_file(argv[1],sel,builder,Accept_none_occupancy_policy())){
    const ESBTL::Default_system::Model& model=*systems[0].models_begin();
    
    std::copy(model.atoms_begin(),model.atoms_end(),std::ostream_iterator<ESBTL::Default_system::Atom>(std::cout,"\n"));
  }
  
  return 0;
}

