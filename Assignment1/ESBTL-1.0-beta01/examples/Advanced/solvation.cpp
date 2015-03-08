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
#include <ESBTL/fcc_lattice.h>
#include <ESBTL/grid_of_cubes.h>

#include <iostream>



typedef ESBTL::Accept_none_occupancy_policy<ESBTL::PDB::Line_format<> > Accept_none_occupancy_policy;


namespace internal{
  template <class Grid,class Clash_functor,class Iterator>
  inline
  bool do_have_a_clash(const Clash_functor& clash,typename Grid::Cube_unit& cube,Iterator obj){
    for (typename Grid::In_cube_iterator it=cube.begin();it!=cube.end();++it)
      if (clash(* (*it),*obj)) return true;
      return false;
    return false;
  }
}

struct Simple_clash{
  template <class Point1,class Point2>
  bool operator()(const Point1& p1, const Point2& p2) const {
    return ESBTL::squared_distance(p1,p2) <= ESBTL::square(4);
  }
};

template <class Point, class Iterator1,class Iterator2,class Output_iterator,class Clash_functor>
void report_steric_clashes( Iterator1 test_begin,Iterator1 test_end,
                            Iterator2 system_begin,Iterator2 system_end,
                            Clash_functor clash,
                            Output_iterator out)
{
  typedef ESBTL::Grid_of_cubes<ESBTL::Traits_for_grid<Point,ESBTL::Default_system::Model::Atoms_const_iterator> > Grid;
  Grid grid(system_begin,system_end);
  
  std::cout << grid.nb_element() << std::endl;

  for (Iterator1 it=test_begin;it!=test_end;++it){
    if ( grid.is_outside_grid(it) ) continue;
    bool in_conflict=false;
    for (typename Grid::neighbor_iterator it_n=grid.first_neighbor(grid.locate_cube(it));it_n!=grid.nend();++it_n){
      if ( internal::do_have_a_clash<Grid>(clash,*it_n,it) ){
        in_conflict=true;
        break;
      }
    }
    if (!in_conflict){
      typename Grid::Cube_unit* cube_ptr=grid.get_cube(it);
      if (cube_ptr!=NULL)
        in_conflict=internal::do_have_a_clash<Grid>(clash,*cube_ptr,it);
    }
    if (in_conflict)
      *out++=it;
  }
}


namespace ESBTL{
  std::ostream&
  operator<<(std::ostream& os, const Point_3& p)
  {
    os << ", SPHERE, " <<  p.x() << " , " << p.y() << " , " << p.z() << " , " << 1 <<"\n";
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
  
  
  
  typedef ESBTL::Generic_classifier<ESBTL::Name_and_radius_of_atom<double,ESBTL::Default_system::Atom> > T_Atom_classifier;
  T_Atom_classifier atom_classif;
  
  std::list<ESBTL::Point_3> solvant;
  
  if (ESBTL::read_a_pdb_file(argv[1],sel,builder,Accept_none_occupancy_policy())){
    const ESBTL::Default_system::Model& model=*systems[0].models_begin();
    
    std::pair<ESBTL::Point_3,double> cube=ESBTL::bounding_cube<ESBTL::Point_3>(model.atoms_begin(),model.atoms_end());
    ESBTL::fcc_lattice<ESBTL::Point_3>(cube,1.4,10,std::back_inserter(solvant));
    
    std::list<std::list<ESBTL::Point_3>::iterator > in_conflict;
    //find water molecule in conflict with the model
    report_steric_clashes<ESBTL::Point_3>( solvant.begin(),solvant.end(),model.atoms_begin(),model.atoms_end(),Simple_clash(),std::back_inserter(in_conflict) );
    
    std::cout << "Before removing " << solvant.size() << std::endl;
    
    for( std::list<std::list<ESBTL::Point_3>::iterator >::iterator it=in_conflict.begin();it!=in_conflict.end();++it ){
      solvant.erase(*it);
    }
    
    std::cout << "After removing " << solvant.size() << std::endl;
    
    unsigned nb_atoms=0;
    for (ESBTL::Default_system::Model::Atoms_const_iterator it_atm=model.atoms_begin();it_atm!=model.atoms_end();++it_atm){
      T_Atom_classifier::Properties type=atom_classif.get_properties(*it_atm);
      ++nb_atoms;
    }
    
    
    
    std::cout << nb_atoms << std::endl;
    
    std::cout << "from pymol.cgo import *\nfrom pymol import cmd\nobj=[COLOR, 1.0, 0, 0," << 1 << "\n";
    std::copy(solvant.begin(),solvant.end(),std::ostream_iterator<ESBTL::Point_3>(std::cout,""));
    std::cout <<  "]\ncmd.load_cgo(obj,'cgo01')" <<std::endl;
  }
  
}

