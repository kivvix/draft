#ifndef _FIELD2D_H
#define _FIELD2D_H

#include <algorithm>
#include <array>
#include <fstream>
#include <tuple>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/multi_array.hpp>

#include <boost/iterator/zip_iterator.hpp>

using namespace boost::numeric;

template < typename _T >
struct field2d: public boost::multi_array<_T,2>
{
  typedef typename boost::multi_array< _T , 2 >::value_type             value_type;
  typedef typename boost::multi_array< _T , 2 >::reference              reference;
  typedef typename boost::multi_array< _T , 2 >::const_reference        const_reference;
  typedef typename boost::multi_array< _T , 2 >::iterator               iterator;
  typedef typename boost::multi_array< _T , 2 >::const_iterator         const_iterator;
  typedef typename boost::multi_array< _T , 2 >::reverse_iterator       reverse_iterator;
  typedef typename boost::multi_array< _T , 2 >::element                element;
  typedef typename boost::multi_array< _T , 2 >::size_type              size_type; 
  typedef typename boost::multi_array< _T , 2 >::difference_type        difference_type; 
  typedef typename boost::multi_array< _T , 2 >::index                  index; 
  typedef typename boost::multi_array< _T , 2 >::extent_range           extent_range;
  typedef typename boost::multi_array< _T , 2 >::const_reverse_iterator const_reverse_iterator;

  typedef typename boost::multi_array< _T , 2 >::index_gen   index_gen;
  typedef typename boost::multi_array< _T , 2 >::index_range index_range;

  typedef std::array< index,2 > index_array;
  
  struct const_array_view {
    typedef boost::detail::multi_array::const_multi_array_view<_T,2> type;
  };
  struct array_view {
    typedef boost::detail::multi_array::multi_array_view<_T,2> type;
  };

  using boost::multi_array<_T,2>::multi_array;

  struct step_t
  {
    _T dx;
    _T dy;
  } steps;
  struct range_t
  {
    _T x_min, x_max;
    _T y_min, y_max;
  } range;

  ~field2d () {;}

  using boost::multi_array<_T,2>::shape;

  inline size_type
  size ( size_type n ) const
  { return shape()[n]; }

  _T
  x ( size_type i ) const
  { return i*steps.dx + range.x_min; }
  _T
  y ( size_type j ) const
  { return j*steps.dy + range.y_min; }

  void
  compute_steps ()
  {
    steps.dx = (range.x_max - range.x_min)/size(0);
    steps.dy = (range.y_max - range.y_min)/size(1);
  }
};

template < typename _T >
std::ostream &
operator << ( std::ostream & os , const field2d<_T> & f ) {
  for ( typename field2d<_T>::size_type i=0u ; i<f.size(0) ; ++i ) {
    for ( typename field2d<_T>::size_type j=0u ; j<f.size(0) ; ++j ) {
      os << f.x(i) << " " << f.y(j) << " " << f[i][j] << "\n";
    }
    os << "\n";
  }
  return os;
}


#endif
