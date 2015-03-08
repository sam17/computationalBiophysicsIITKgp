#include <iostream>
#include <ESBTL/default.h>
#include <ESBTL/atom_classifier.h>
#include <ESBTL/weighted_atom_iterator.h>
#include <algorithm>
#include <list>
#include <cmath>
#include <limits>

#define PI 3.14159265

typedef ESBTL::Accept_none_occupancy_policy<ESBTL::PDB::Line_format<> > Accept_none_occupancy_policy;

//Class to define Direction Vectors
template< class T >
class D3Vector {
 
template< class U >
friend std::ostream & operator<<( std::ostream & , const D3Vector<U> & ) ;   
 
public :
   D3Vector( T a , T b , T c ) {
      x = a ;
      y = b ;
      z = c ;
   }

    D3Vector& operator=(D3Vector const& f)
    {
	this->x = f.x;
	this->y = f.y;
	this->z = f.z;
	return *this;
    }

    D3Vector& operator/(T f)
    {
	this->x = this->x/f;
	this->y = this->y/f;
	this->z = this->z/f;
	return *this;
    }

   T dotproduct ( const D3Vector & rhs )
   {
      T scalar = x * rhs.x + y * rhs.y + z * rhs.z ;
      return scalar ;
   }

   T absolute ()
   {
       T absolute = sqrt((this->x*this->x) + (this->y*this->y) + (this->z*this->z));
       return absolute;
   }
	    
	    
   D3Vector crossproduct ( const D3Vector & rhs )
   {
      T a = y * rhs.z - z * rhs.y ;
      T b = z * rhs.x - x * rhs.z ;
      T c = x * rhs.y - y * rhs.x ;
      D3Vector product( a , b , c ) ;
      return product ;
   }
 
   D3Vector triplevec( D3Vector & a , D3Vector & b ) {
      return crossproduct ( a.crossproduct( b ) ) ;
   }
 
   T triplescal( D3Vector & a, D3Vector & b )
   {
      return dotproduct( a.crossproduct( b ) ) ;
   }

   T x , y , z ;  
} ;
 
template< class T > std::ostream & operator<< ( std::ostream & os , const D3Vector<T> & vec )
{
   os << "( "  << vec.x << " ,  " << vec.y << " ,  " << vec.z << " )" ;
   return os ;
}

//Function to calculate dihedral angles from three vectors
template < class T> T CalculateDihedralAngle(D3Vector<T> a,D3Vector<T> b,D3Vector<T> c)
{
    T angle;
    D3Vector<T> d = b.crossproduct(a);
    d = d/d.absolute();
    D3Vector<T> e = c.crossproduct(b);
    e = e/e.absolute();
    b = b/b.absolute();
    D3Vector<T> m = b.crossproduct(d);

    T y,x;
    x = d.dotproduct(e);
    y = e.dotproduct(m);
    angle = atan2(y,x);
    return angle;
}

int main(int argc, char** argv)
{
    // Check for 2 arguments
    if (argc != 2 )
    {
	std::cerr << "Please provide a filename"  << std::endl;  
	return EXIT_FAILURE;
    }
    
    char *chain = "A";

    //System Building functions from PDB file
    ESBTL::PDB_line_selector sel;
    std::vector<ESBTL::Default_system> systems;
    ESBTL::All_atom_system_builder<ESBTL::Default_system> builder(systems,sel.max_nb_systems());
  
    
    if (ESBTL::read_a_pdb_file(argv[1],sel,builder,Accept_none_occupancy_policy()))
    {
	if ( systems.empty() || systems[0].has_no_model() )
	{
	    std::cerr << "No atoms found" << std::endl;
	    return EXIT_FAILURE;
	}

	const ESBTL::Default_system::Model& model=* systems[0].models_begin();

	D3Vector<double> atom(0.0,0.0,0.0);
	std::vector< D3Vector<double> > vectorList;
	std::vector< D3Vector<double> > atomList;

	//Loop to propagate list of atoms of backbone only and of specified chain
	for (ESBTL::Default_system::Model::Atoms_const_iterator it_atm=model.atoms_begin();it_atm!=model.atoms_end();++it_atm)
	{
	    if((it_atm->atom_name() =="CA" || it_atm->atom_name() =="N" || it_atm->atom_name() =="C")&& (it_atm->chain_identifier()==chain[0]))
	    {
		atom.x = it_atm->x();
		atom.y = it_atm->y();
		atom.z = it_atm->z();
		atomList.push_back(atom);
	    }
	}
	
	//Loop to generate vectors from atom x,y,z
	for(int i = 1;i<atomList.size();i++)
	{
	    atom.x = atomList[i].x-atomList[i-1].x;
	    atom.y = atomList[i].y-atomList[i-1].y;
	    atom.z = atomList[i].z-atomList[i-1].z;
	    vectorList.push_back(atom);
	}

	std::cout<<"PHI \t PSI \t OMEGA"<<std::endl;

	//Loop to calculate Dihedral Angles
	for(int i=0;i<vectorList.size();i=i+3)
	{
	    double phi=std::numeric_limits<double>::quiet_NaN(),psi = std::numeric_limits<double>::quiet_NaN(),omega=std::numeric_limits<double>::quiet_NaN();

	    if(i>0 && i<vectorList.size()-1)
		phi = CalculateDihedralAngle(vectorList[i-1],vectorList[i],vectorList[i+1])* 180.0 / PI;
	    if(i<vectorList.size()-2)
		psi = CalculateDihedralAngle(vectorList[i],vectorList[i+1],vectorList[i+2])* 180.0 / PI;
	    if(i<vectorList.size()-3)
		omega = CalculateDihedralAngle(vectorList[i+1],vectorList[i+2],vectorList[i+3])* 180.0 / PI;
	    
	    std::cout<<phi<<" \t "<<psi<<" \t "<<omega<<std::endl;
	}    
    }
    else
	return EXIT_FAILURE;
  
  return EXIT_SUCCESS;
}
