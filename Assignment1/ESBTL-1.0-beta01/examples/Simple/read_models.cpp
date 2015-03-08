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

//discard atoms with no alternate location identification and with occupancy != 1
typedef ESBTL::Accept_none_occupancy_policy<ESBTL::PDB::Line_format<> > Accept_none_occupancy_policy;

int main(int argc, char** argv){
  
  if (argc != 2 ){
    std::cerr << "Please provide a filename"  << std::endl;  
    return EXIT_FAILURE;
  }
  
  //creates two systems of heavy atoms, one containing non-water molecules, the second containing water molecules
  ESBTL::PDB_line_selector_two_systems sel;
  
  std::vector<ESBTL::Default_system> systems;
  //declare the building that will contruct the systems from the pdb file.
  ESBTL::All_atom_system_builder<ESBTL::Default_system> builder(systems,sel.max_nb_systems());
  
  //read the pdb file
  if (ESBTL::read_a_pdb_file(argv[1],sel,builder,Accept_none_occupancy_policy())){
    if ( systems.empty() || systems[0].has_no_model() ){
      std::cerr << "No atoms found" << std::endl;
      return EXIT_FAILURE;
    }
    //Iterate over all models in the first system
    for (ESBTL::Default_system::Models_iterator 
          it_model=systems[0].models_begin();
          it_model!=systems[0].models_end();
          ++it_model)
    {
      const ESBTL::Default_system::Model& model=*it_model;
      unsigned nb_atom=0,nb_hetatm=0;
      //iterate over all atoms of a model
      for (ESBTL::Default_system::Model::Atoms_const_iterator it_atm=model.atoms_begin();it_atm!=model.atoms_end();++it_atm){
        if (it_atm->is_hetatm()) ++nb_hetatm;
        else ++nb_atom;
      }
      std::cout << "Model " << model.model_number() << " contains "
                << nb_atom << " atoms and " << nb_hetatm << " hetero-atoms within "
                << model.number_of_residues() << " residues.\n";
    }
  }
  else
    return EXIT_FAILURE;
  
  return EXIT_SUCCESS;
}

