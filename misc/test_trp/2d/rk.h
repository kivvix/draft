#ifndef _RK_H
#define _RK_H

#include "field2d.h"
#include "lagrange5.h"

#include <boost/multi_array.hpp>

template <typename _T, typename L_t>
auto
rk33 ( L_t Lij , _T tn , const field2d<_T> & un , _T dt )
{
  field2d<_T> u1(boost::extents[un.size(0)][un.size(1)]);
  field2d<_T> u2(boost::extents[un.size(0)][un.size(1)]);
  field2d<_T> un1(boost::extents[un.size(0)][un.size(1)]);
  u1.range  = un.range; u1.steps  = un.steps;
  u2.range  = un.range; u2.steps  = un.steps;
  un1.range = un.range; un1.steps = un.steps;

  for ( auto i=0u ; i<un.size(0) ; ++i ) {
    for ( auto j=0u ; j<un.size(1) ; ++j ) {
      u1[i][j] = un[i][j] + dt*Lij(tn,un,i,j);
    }
  }
  for ( auto i=0u ; i<un.size(0) ; ++i ) {
    for ( auto j=0u ; j<un.size(1) ; ++j ) {
      u2[i][j] = 0.75*un[i][j] + 0.25*u1[i][j] + 0.25*dt*Lij(tn+dt,u1,i,j);
    }
  }
  for ( auto i=0u ; i<un.size(0) ; ++i ) {
    for ( auto j=0u ; j<un.size(1) ; ++j ) {
      un1[i][j] = (1./3.)*un[i][j] + (2./3.)*u2[i][j] + (2./3.)*dt*Lij(tn+0.5*dt,u2,i,j);
    }
  }
  return un1;
}


namespace Lie {
  template <typename _T>
  auto
  phi1 ( _T tn , const field2d<_T> & un , _T dt )
  {
    field2d<_T> un1(boost::extents[un.size(0)][un.size(1)]);
    un1.range = un.range; un1.steps = un.steps;
    std::size_t Nx = un.size(0);

    for ( auto i=0u ; i<un.size(0) ; ++i ) {
      for ( auto j=0u ; j<un.size(1) ; ++j ) {
        _T xstar = un.x(i) + dt*un.y(j);
        int istar = std::ceil((xstar - un.range.x_min)/un.steps.dx);

        auto N = lagrange5::generator(
            un[(istar-3)%Nx][j],
            un[(istar-2)%Nx][j],
            un[(istar-1)%Nx][j],
            un[(istar  )%Nx][j],
            un[(istar+1)%Nx][j],
            un[(istar+2)%Nx][j],
            un.steps.dx , istar*un.steps.dx + un.range.x_min
          );
        un1[i][j] = N(xstar);
      }
    }

    return un1;
  }

  template <typename _T>
  auto
  phi2 ( _T tn , const field2d<_T> & un , _T dt )
  {
    field2d<_T> un1(boost::extents[un.size(0)][un.size(1)]);
    un1.range = un.range; un1.steps = un.steps;
    std::size_t Ny = un.size(1);

    for ( auto i=0u ; i<un.size(0) ; ++i ) {
      for ( auto j=0u ; j<un.size(1) ; ++j ) {
        _T ystar = un.y(j) - dt*un.x(i);
        int jstar = std::ceil((ystar - un.range.y_min)/un.steps.dy);

        auto N = lagrange5::generator(
            un[i][(jstar-3)%Ny],
            un[i][(jstar-2)%Ny],
            un[i][(jstar-1)%Ny],
            un[i][(jstar  )%Ny],
            un[i][(jstar+1)%Ny],
            un[i][(jstar+2)%Ny],
            un.steps.dy , jstar*un.steps.dy + un.range.y_min
          );
        un1[i][j] = N(ystar);
      }
    }

    return un1;
  }
}

template <typename _T>
auto
rotation ( const field2d<_T> & h , _T t ) {
  field2d<_T> hs(boost::extents[h.size(0)][h.size(1)]);
  field2d<_T> hss(boost::extents[h.size(0)][h.size(1)]);
  field2d<_T> hr(boost::extents[h.size(0)][h.size(1)]);
  hr.range = h.range; hr.steps = h.steps;

  const std::size_t Nx = h.size(0);
  const std::size_t Ny = h.size(1);

  _T dt1 = std::tan(0.5*t);
  for ( auto i=0u ; i<Nx ; ++i ) {
    for ( auto j=0u ; j<Ny ; ++j ) {
      _T xstar = h.x(i) - dt1*h.y(j);
      int istar = std::ceil((xstar - h.range.x_min)/h.steps.dx);

      auto N = lagrange5::generator(
          h[(istar-3)%Nx][j],
          h[(istar-2)%Nx][j],
          h[(istar-1)%Nx][j],
          h[(istar  )%Nx][j],
          h[(istar+1)%Nx][j],
          h[(istar+2)%Nx][j],
          h.steps.dx , istar*h.steps.dx + h.range.x_min
        );
      hs[i][j] = N(xstar);
    }
  }

  _T dt2 = std::sin(t);
  for ( auto i=0u ; i<Nx ; ++i ) {
    for ( auto j=0u ; j<Ny ; ++j ) {
      _T ystar = h.y(j) + dt2*h.x(i);
      int jstar = std::ceil((ystar - h.range.y_min)/h.steps.dy);

      auto N = lagrange5::generator(
          hs[i][(jstar-3)%Ny],
          hs[i][(jstar-2)%Ny],
          hs[i][(jstar-1)%Ny],
          hs[i][(jstar  )%Ny],
          hs[i][(jstar+1)%Ny],
          hs[i][(jstar+2)%Ny],
          h.steps.dy , jstar*h.steps.dy + h.range.y_min
        );
      hss[i][j] = N(ystar);
    }
  }

  _T dt3 = std::tan(0.5*t);
  for ( auto i=0u ; i<Nx ; ++i ) {
    for ( auto j=0u ; j<Ny ; ++j ) {
      _T xstar = h.x(i) - dt3*h.y(j);
      int istar = std::ceil((xstar - h.range.x_min)/h.steps.dx);

      auto N = lagrange5::generator(
          hss[(istar-3)%Nx][j],
          hss[(istar-2)%Nx][j],
          hss[(istar-1)%Nx][j],
          hss[(istar  )%Nx][j],
          hss[(istar+1)%Nx][j],
          hss[(istar+2)%Nx][j],
          h.steps.dx , istar*h.steps.dx + h.range.x_min
        );
      hr[i][j] = N(xstar);
    }
  }

  return hr;
}

#endif
