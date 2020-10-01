#include <iostream>
#include <cmath>
#include <numeric>
#include <fstream>
#include <algorithm>

#include "field2d.h"
#include "weno.h"
#include "rk.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/math/constants/constants.hpp>
namespace math = boost::math::constants;
using namespace boost::numeric;

/*
TODO :

* mesurer l'ordre (être certain de l'ordre en temps et en espace)
* chercher une autre source de bug qu'une erreur d'implémentation des schémas...

*/

auto
pacman_factory ( double r , double alpha , double value=1.0 ) {
  return [=]( double x , double y ) {
    if ( ( x*x + y*y < r*r ) && ( (x<0.0) || (std::abs(y)>alpha*x) ) ) {
      return value;
    }
    return 0.0;
  };
}
#define SQ(X) ((X)*(X))
auto
gauss_factory ( double X0 , double Y0 , double tx , double ty ) {
  return [=]( double x , double y ) {
    return std::exp(-SQ(x-X0)/tx - SQ(y-Y0)/ty);
  };
}


int
main(int,char**)
{
  std::size_t N = 100;
  double xmax = 5.0;
  std::size_t Nx=N,Ny=N;
  field2d<double> f(boost::extents[Nx][Ny]);
  field2d<double> df(boost::extents[Nx][Ny]);
  f.range.x_min = -xmax; f.range.x_max = xmax;
  f.range.y_min = -xmax; f.range.y_max = xmax;
  f.compute_steps();
  df.range = f.range;
  df.steps = f.steps;

  auto pacman = pacman_factory(1.0,0.5);
  auto g = gauss_factory(0.0,1.0,0.75,0.25);
  for ( auto i=0u ; i<f.size(0) ; ++i ) {
    for ( auto j=0u ; j<f.size(1) ; ++j ) {
      // f[i][j] = g(f.x(i),f.y(j));
       f[i][j] = pacman(f.x(i),f.y(j));
      // f[i][j] = std::cos(2.*math::pi<double>()*f.x(i)/xmax)*std::sin(2.*math::pi<double>()*f.y(j)/xmax);
      //f[i][j] = std::cos(2.*math::pi<double>()*f.x(i)/xmax);
    }
  }

  ublas::vector<double> x(Nx),y(Ny);
  std::generate( x.begin() , x.end() , [&,count=0] () mutable { return f.x(count++); } );
  std::generate( y.begin() , y.end() , [&,count=0] () mutable { return f.y(count++); } );

  auto save_f = [&](std::string filename) {
    std::ofstream of(filename);
    of << f << std::endl;
    of.close();
  };

  save_f("finit.dat");

  double dt= 1.3*f.steps.dx/(2.0*xmax);
  double Tf= .5*math::pi<double>();
  double current_time = 0.0;

  auto Lij = [&](double tn , const field2d<double> & u , std::size_t i, std::size_t j){
    //return - (weno2d::weno_x( -y[j] ,u,i,j) + weno2d::weno_y(  x[i] ,u,i,j));
    return - (weno2d::weno_x( -u.y(j) ,u,i,j) + weno2d::weno_y(  u.x(i) ,u,i,j));
    //return - ( weno2d::weno_x( -2.0 ,u,i,j) + weno2d::weno_y( 0.0 ,u,i,j) );
  };

  auto Lij_x = [](double tn , const field2d<double> & u , std::size_t i, std::size_t j){
    return -weno2d::weno_x( -u.y(j) ,u,i,j);
  };
  auto Lij_y = [](double tn , const field2d<double> & u , std::size_t i, std::size_t j){
    return  -weno2d::weno_y( u.x(j) ,u,i,j);
  };

  std::cout << "dt: "<< dt << "\tNiter: " << std::floor(Tf/dt) << "\n";
  std::size_t i_iter = 0;
  while ( i_iter*dt < Tf ) {
    std::cout << current_time << " \r" << std::flush;

    //f = rk33( Lij , current_time , f , dt );

    f = Lie::phi1( current_time , f , 0.5*dt );
    f = Lie::phi2( current_time , f , dt );
    f = Lie::phi1( current_time , f , 0.5*dt );

    ++i_iter;
    current_time += dt;
  }
  std::cout << current_time << std::endl;

  save_f("fend.dat");

  return 0;
}