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
#include <ESBTL/atom_classifier.h>
#include <ESBTL/selected_atom_iterator.h>

struct Residue_filter{
  std::string name;
  Residue_filter(){};
  Residue_filter(const std::string& str):name(str){}
  template <class Atom>
  bool operator()(const Atom& atom){
    return atom.residue().residue_name()==name;
  }
};

typedef ESBTL::Accept_none_occupancy_policy<ESBTL::PDB::Line_format<> > Accept_none_occupancy_policy;

int main(int argc, char** argv){
  
  if (argc != 2 ){
    std::cerr << "Please provide a filename"  << std::endl;  
    return EXIT_FAILURE;
  }
  
  std::list<std::string> chains_to_select;
  chains_to_select.push_back("A");
  chains_to_select.push_back("C");
  
  bool keep_water=true;
  bool keep_remaining_chains=true;
  
  //discard hydrogens and water
  ESBTL::PDB_line_selector_chain sel(chains_to_select.begin(),chains_to_select.end(),keep_water,keep_remaining_chains);
  
  std::vector<ESBTL::Default_system> systems;
  ESBTL::All_atom_system_builder<ESBTL::Default_system> builder(systems,sel.max_nb_systems());
  
  
  
  typedef ESBTL::Generic_classifier<ESBTL::Name_and_radius_of_atom<double,ESBTL::Default_system::Atom> > T_Atom_classifier;
  T_Atom_classifier atom_classif;
  
  
  if (ESBTL::read_a_pdb_file(argv[1],sel,builder,Accept_none_occupancy_policy())){
    for (std::vector<ESBTL::Default_system>::const_iterator it_sys=systems.begin();it_sys!=systems.end();++it_sys){
      for (ESBTL::Default_system::Models_const_iterator it_model=it_sys->models_begin();it_model!=it_sys->models_end();++it_model){
        const ESBTL::Default_system::Model& model=*it_model;
        unsigned nb_atm_1=0;
        unsigned nb_atm_2=0;
        unsigned nb_atm_3=0;
        unsigned nb_hetatm=0;
        unsigned nb_chain=0;
        unsigned nb_residue_1=0;
        unsigned nb_residue_2=0;
        for (ESBTL::Default_system::Model::Chains_const_iterator it_ch=model.chains_begin();it_ch!=model.chains_end();++it_ch){
          ++nb_chain;
          for (ESBTL::Default_system::Chain::Residues_const_iterator it_res=it_ch->residues_begin();it_res!=it_ch->residues_end();++it_res){
            ++nb_residue_1;
            for (ESBTL::Default_system::Residue::Atoms_const_iterator it_atm=it_res->atoms_begin();it_atm!=it_res->atoms_end();++it_atm){
              assert(&(it_atm->residue())==&(*it_res));
              ++nb_atm_1;
            }
          }
        }
        
        for (ESBTL::Default_system::Model::Residues_const_iterator it_res=model.residues_begin();it_res!=model.residues_end();++it_res){
          ++nb_residue_2;
            for (ESBTL::Default_system::Residue::Atoms_const_iterator it_atm=it_res->atoms_begin();it_atm!=it_res->atoms_end();++it_atm){
              assert(&(it_atm->residue())==&(*it_res));
              ++nb_atm_2;
            }
        }
        
        for (ESBTL::Default_system::Model::Atoms_const_iterator it_atm=model.atoms_begin();it_atm!=model.atoms_end();++it_atm){
          T_Atom_classifier::Properties type=atom_classif.get_properties(*it_atm);
          assert( &(model.system()) == &(*it_sys) );
          
          ESBTL::Default_system::Atom atm=*it_atm;
          //test operator==
          assert(atm.atom_name()+atm.residue_name()==it_atm->atom_name()+it_atm->residue_name());
          //test operator=
          atm=*model.atoms_begin();
          assert(atm.atom_name()==model.atoms_begin()->atom_name());
          
          if (it_atm->is_hetatm())
            ++nb_hetatm;
          else
            ++nb_atm_3;
        }
        
        assert(nb_hetatm+nb_atm_3==nb_atm_1);
        assert(nb_atm_2==nb_atm_1);
        assert(nb_residue_1==nb_residue_2);

        

        std::cout << "(nb_atm,nb_hetatm) = (" <<nb_atm_3<<","<<nb_hetatm<<")" << std::endl;
        std::cout << "nb res " <<  nb_residue_1 << std::endl;
        std::cout << "nb chains " <<nb_chain << std::endl;
        
        typedef ESBTL::Selected_atom_iterator<ESBTL::Default_system::Model,Residue_filter,true>  Restrict_iterator;
        unsigned nb_ala=0;
        for (Restrict_iterator itr=ESBTL::make_selected_atom_iterator(model.atoms_begin(),Residue_filter("ALA"));
             itr!=ESBTL::make_selected_atom_iterator<Residue_filter>(model.atoms_end());++itr)
        {
          ++nb_ala;
        }
        std::cout << "nb ALA atoms " << nb_ala << std::endl;
        
      }
    }
  }
  else
    return EXIT_FAILURE;
  
  return EXIT_SUCCESS;
}

