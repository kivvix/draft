#include <iostream>
#include <cmath>
#include <numeric>
#include <fstream>
#include <algorithm>
#include <iomanip>

#include "field2d.h"
#include "weno.h"
#include "rk.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/math/constants/constants.hpp>
namespace math = boost::math::constants;
using namespace boost::numeric;

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
#undef SQ

int
main(int,char**)
{
  std::size_t N = 100;
  double xmax = 5.0;
  std::size_t Nx=N,Ny=N;
  field2d<double> f(boost::extents[Nx][Ny]);
  f.range.x_min = -xmax; f.range.x_max = xmax;
  f.range.y_min = -xmax; f.range.y_max = xmax;
  f.compute_steps();

  field2d<double> g(boost::extents[Nx][Ny]);
  g.range = f.range;
  g.steps = f.steps;

  auto pacman = pacman_factory(1.0,0.5);
  auto gauss = gauss_factory(0.0,1.0,0.75,0.25);
  for ( auto i=0u ; i<f.size(0) ; ++i ) {
    for ( auto j=0u ; j<f.size(1) ; ++j ) {
       //f[i][j] = gauss(f.x(i),f.y(j));
       f[i][j] = pacman(f.x(i),f.y(j));
      // f[i][j] = std::cos(2.*math::pi<double>()*f.x(i)/xmax)*std::sin(2.*math::pi<double>()*f.y(j)/xmax);
      //f[i][j] = std::cos(2.*math::pi<double>()*f.x(i)/xmax);

       g[i][j] = f[i][j];
    }
  }
  ublas::vector<double> E(2,0); E[0] = 0.1; E[1] = 0.1;

  //double dt= 1.3*f.steps.dx/(std::abs(E[0])+std::abs(E[1])+2.0*xmax);
  double dt= 1.3*f.steps.dx/(std::abs(E[0])+std::abs(E[1]));
  double Tf= 1.0; //.5*math::pi<double>();


  auto save_f = [](field2d<double> const& f,std::string filename) {
    std::ofstream of(filename);
    of << f << std::endl;
    of.close();
  };

  save_f(f,"finit.dat");
  save_f(g,"ginit.dat");



  auto Lij_f = [&](double tn , const field2d<double> & u , std::size_t i, std::size_t j){
    return - (weno2d::weno_x( E[0]+u.y(j) ,u,i,j) + weno2d::weno_y(  E[1]-u.x(i) ,u,i,j));
  };
  auto Lij_g = [&](double tn , const field2d<double> & u , std::size_t i , std::size_t j){
    return - ( weno2d::weno_x( E[0]*std::cos(tn) - E[1]*std::sin(tn) ,u,i,j)
             + weno2d::weno_y( E[0]*std::sin(tn) + E[1]*std::cos(tn) ,u,i,j) );
  };

  std::cout << "dt: "<< dt << "\tNiter: " << std::floor(Tf/dt) << "\n";
  double current_time = 0.0;
  std::size_t i_iter = 0;
  while ( current_time < Tf ) {
    std::cout << "[" << std::setw(10) << current_time << "] \r" << std::flush;

    f = rk33( Lij_f , current_time , f , dt );

    g = rk33( Lij_g , current_time , g , dt );

    ++i_iter;
    current_time += dt;
    if ( current_time+dt > Tf ) { dt = Tf - current_time; }
  }
  std::cout << "[" << std::setw(10) << current_time << "]" << std::endl;

  save_f(f,"fend.dat");
  save_f(g,"gend.dat");

  f = rotation(g,current_time);

  save_f(f,"ffend.dat");

  return 0;
}
