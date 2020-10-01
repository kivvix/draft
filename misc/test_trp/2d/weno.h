#ifndef _WENO_H
#define _WENO_H

#include <algorithm>
#include <iostream>
#include <iterator>
#include <utility>
#include <type_traits>
#include <tuple>

#include "field2d.h"


namespace weno {
  #define SQ(X) ((X)*(X))

  template < typename _T >
  auto
  local_flux ( _T fim2 , _T fim1 , _T fi , _T fip1 , _T fip2 , _T fip3 )
  {
    const static _T epsi = 1e-6;

    // NB pour optimisation minuscule : 2 termes sont calculés 2 fois :
    // * `(13./12.)*SQ( f_it->template get<1>() - 2.*f_it->template get<2>() + f_it->template get<3>() )` dans `w0p` et `w2m`
    // * `(13./12.)*SQ( f_it->template get<2>() - 2.*f_it->template get<3>() + f_it->template get<4>() )` dans `w2p` et `w1m`

    // $f_{i,k+1/2}^+$
    _T w0p = (13./12.)*SQ( fim2 - 2.*fim1 + fi   ) + 0.25*SQ(    fim2 - 4.*fim1 + 3.*fi   );
    _T w1p = (13./12.)*SQ( fim1 - 2.*fi   + fip1 ) + 0.25*SQ(    fim1           -    fip1 );
    _T w2p = (13./12.)*SQ( fi   - 2.*fip1 + fip2 ) + 0.25*SQ( 3.*fi   - 4.*fip1 +    fip2 );

    w0p = 0.1/(SQ( epsi + w0p )); w1p = 0.6/(SQ( epsi + w1p )); w2p = 0.3/(SQ( epsi + w2p ));

    _T sum_wp = w0p+w1p+w2p;
    w0p /= sum_wp; w1p /= sum_wp; w2p /= sum_wp;

    _T fikp12p = w0p*( (2./6.)*fim2 - (7./6.)*fim1 + (11./6.)*fi   )
               + w1p*(-(1./6.)*fim1 + (5./6.)*fi   +  (2./6.)*fip1 )
               + w2p*( (2./6.)*fi   + (5./6.)*fip1 -  (1./6.)*fip2 );

    // $f_{i,k+1/2}^-$
    _T w0m = (13./12.)*SQ( fip1 - 2.*fip2 + fip3 ) + 0.25*SQ( 3.*fip1 - 4.*fip2 +    fip3 );
    _T w1m = (13./12.)*SQ( fi   - 2.*fip1 + fip2 ) + 0.25*SQ(    fi             -    fip2 );
    _T w2m = (13./12.)*SQ( fim1 - 2.*fi   + fip1 ) + 0.25*SQ(    fim1 - 4.*fi   + 3.*fip1 );

    w0m = 0.1/SQ( epsi + w0m ); w1m = 0.6/SQ( epsi + w1m ); w2m = 0.3/SQ( epsi + w2m );

    _T sum_wm = w0m+w1m+w2m;
    w0m /= sum_wm; w1m /= sum_wm; w2m /= sum_wm;

    _T fikp12m = w2m*(-(1./6.)*fim1 + (5./6.)*fi   + (2./6.)*fip1 )
               + w1m*( (2./6.)*fi   + (5./6.)*fip1 - (1./6.)*fip2 )
               + w0m*((11./6.)*fip1 - (7./6.)*fip2 + (2./6.)*fip3 );
    
    // return a pair of fluxes, one loop on the data structure for two fluxes
    return std::make_pair(fikp12p,fikp12m);
  }

  #undef SQ
}

namespace weno2d {

  template < typename _T >
  _T
  du_weno ( _T velocity , _T uim3 , _T uim2 , _T uim1 , _T ui , _T uip1 , _T uip2 , _T uip3 , _T dx ) {
    _T vp = std::max(velocity,0.);
    _T vm = std::min(velocity,0.);

    auto fip12 = weno::local_flux(uim2,uim1,ui,uip1,uip2,uip3);
    auto fim12 = weno::local_flux(uim3,uim2,uim1,ui,uip1,uip2);

    return ( vp*(fip12.first - fim12.first) + vm*(fip12.second - fim12.second) )/dx;
  }

  std::tuple<std::size_t,std::size_t,std::size_t,std::size_t,std::size_t,std::size_t,std::size_t>
  periodic_index ( const std::size_t N , const std::size_t k ) {
    std::size_t km3 = static_cast<std::size_t>(k-3), km2 = static_cast<std::size_t>(k-2), km1 = static_cast<std::size_t>(k-1),
                kp1 = static_cast<std::size_t>(k+1), kp2 = static_cast<std::size_t>(k+2), kp3 = static_cast<std::size_t>(k+3);
    if ( k < 3u ) {
      km3 = static_cast<std::size_t>(( k-3 +N)%N); km2 = static_cast<std::size_t>(( k-2 +N)%N); km1 = static_cast<std::size_t>(( k-1 +N)%N);
    }
    if ( k >= N-3 ) {
      kp1 = static_cast<std::size_t>(( k+1 )%N); kp2 = static_cast<std::size_t>(( k+2 )%N); kp3 = static_cast<std::size_t>(( k+3 )%N);
    }
    return std::make_tuple(km3,km2,km1,k,kp1,kp2,kp3);
  }

  template < typename _T >
  _T
  weno_x ( _T velocity , const field2d<_T> & u , std::size_t i , std::size_t j ) {
    const std::size_t Nx = u.size(0);
    std::size_t im3, im2, im1, ip1, ip2, ip3;
    std::tie(im3,im2,im1,i,ip1,ip2,ip3) = periodic_index(Nx,i);

    return du_weno( velocity , u[im3][j] , u[im2][j] , u[im1][j] , u[i][j] , u[ip1][j] , u[ip2][j] , u[ip3][j] , u.steps.dx );
  }

  template < typename _T >
  _T
  weno_y ( _T velocity , const field2d<_T> & u , std::size_t i , std::size_t j ) {
    const std::size_t Ny = u.size(1);
    std::size_t jm3, jm2, jm1, jp1, jp2, jp3;
    std::tie(jm3,jm2,jm1,j,jp1,jp2,jp3) = periodic_index(Ny,j);

    return du_weno( velocity , u[i][jm3] , u[i][jm2] , u[i][jm1] , u[i][j] , u[i][jp1] , u[i][jp2] , u[i][jp3] , u.steps.dy );
  }

}

#endif
