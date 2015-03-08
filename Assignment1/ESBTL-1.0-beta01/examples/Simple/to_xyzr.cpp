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
#include <ESBTL/weighted_atom_iterator.h>


typedef ESBTL::Generic_classifier<ESBTL::Radius_of_atom<double,ESBTL::Default_system::Atom> >                                  T_Atom_classifier;
typedef ESBTL::Accept_none_occupancy_policy<ESBTL::PDB::Line_format<> > Accept_none_occupancy_policy;

int main(int argc, char** argv){
  
  if (argc != 2 ){
    std::cerr << "Please provide a filename"  << std::endl;  
    return EXIT_FAILURE;
  }
  
  ESBTL::PDB_line_selector_two_systems sel;
  
  std::vector<ESBTL::Default_system> systems;
  ESBTL::All_atom_system_builder<ESBTL::Default_system> builder(systems,sel.max_nb_systems());
  
  
  
  T_Atom_classifier atom_classifier;
  
  if (ESBTL::read_a_pdb_file(argv[1],sel,builder,Accept_none_occupancy_policy())){
    if ( systems.empty() || systems[0].has_no_model() ){
      std::cerr << "No atoms found" << std::endl;
      return EXIT_FAILURE;
    }
    //Consider only the first model of the first system
    const ESBTL::Default_system::Model& model=* systems[0].models_begin();
    for (ESBTL::Default_system::Model::Atoms_const_iterator it_atm=model.atoms_begin();it_atm!=model.atoms_end();++it_atm){
      std::cout << it_atm->x() << " " << it_atm->y() << " " << it_atm->z() << " " ;
      std::cout << atom_classifier.get_properties(*it_atm).value() << "\n";
    }
  }
  else
    return EXIT_FAILURE;
  
  return EXIT_SUCCESS;
}

