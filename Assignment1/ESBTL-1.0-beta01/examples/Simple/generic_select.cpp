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
#include <ESBTL/combine_boolean_operator.h>

//discard atoms with no alternate location identification and with occupancy != 1
typedef ESBTL::Accept_none_occupancy_policy<ESBTL::PDB::Line_format<> > Accept_none_occupancy_policy;

int main(int argc, char** argv){
  
  if (argc != 3 ){
    std::cerr << "Please provide a filename and a residue name"  << std::endl;  
    return EXIT_FAILURE;
  }

  //defines a criterion: consider all non-hydrogen atoms of all but HOH residues.
  typedef ESBTL::Not_functor< ESBTL::Or_functors<ESBTL::Select_by_element,ESBTL::Select_by_resname> > Criterion;
  Criterion criterion( 
    ESBTL::Or_functors<ESBTL::Select_by_element,ESBTL::Select_by_resname>( ESBTL::Select_by_element("H"),
                                                                           ESBTL::Select_by_resname("HOH") 
    ) 
  );
  
  //creates one system with all atoms matching the criterion
  ESBTL::Generic_line_selector<Criterion> sel(criterion);
  
  std::vector<ESBTL::Default_system> systems;
  //declare the building that will contruct the system from the pdb file.
  ESBTL::All_atom_system_builder<ESBTL::Default_system> builder(systems,sel.max_nb_systems());
  
  //read the pdb file
  if (ESBTL::read_a_pdb_file(argv[1],sel,builder,Accept_none_occupancy_policy())){
    if ( systems.empty() || systems[0].has_no_model() ){
      std::cerr << "No atoms found" << std::endl;
      return EXIT_FAILURE;
    }
    //consider the first model of the system
    const ESBTL::Default_system::Model& model=*systems[0].models_begin();

    
    //Criterion to select all carbon alpha atoms of a residues with a given name
    typedef ESBTL::And_functors<ESBTL::Select_by_resname, ESBTL::Select_by_atmname> Criterion_subselect;
    //iterator restricted to atoms of a given residue type
    typedef ESBTL::Selected_atom_iterator<ESBTL::Default_system::Model,Criterion_subselect,true>  Restrict_iterator;
    
    
    unsigned nb_atoms=0;
    std::string resname(argv[2]);
    //declare a function object which operator() returns true if an atom belongs
    //to a residue named 'resname' and if its atom name is "CA"
    Criterion_subselect sel_ca_resn(ESBTL::Select_by_resname(resname),ESBTL::Select_by_atmname("CA"));
    for (Restrict_iterator itr=ESBTL::make_selected_atom_iterator(model.atoms_begin(),sel_ca_resn);
                           itr!=ESBTL::make_selected_atom_iterator(model.atoms_end(),sel_ca_resn);
                           ++itr
        )
    {
      ++nb_atoms;
    }
    std::cout << nb_atoms << " residue(s) named " << resname << std::endl;
  }
  else
    return EXIT_FAILURE;
  
  return EXIT_SUCCESS;
}

