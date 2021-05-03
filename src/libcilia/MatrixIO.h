#include "Matrix.h"
#include "Vector.h"
#include <iostream>

namespace Numeric_lib {

//-----------------------------------------------------------------------------

template <class T>
std::ostream &operator<<(std::ostream &os, const Matrix<T> &v) {
  Index i = 0;
  for (; i < v.dim1() - 1; ++i) {
    os << v(i) << " ";
  }
  os << v(i) << endl;
  return os;
}

//-----------------------------------------------------------------------------

template <class T>
std::ostream &operator<<(std::ostream &os, const Matrix<T, 2> &m) {
  for (Index i = 0; i < m.dim1(); ++i)
    os << m[i];
  os << endl;
  return os;
}

template <class T>
std::ostream &operator<<(std::ostream &os, const Matrix<T, 3> &m) {
  for (Index i = 0; i < m.dim1(); ++i)
    os << m[i];
  return os;
}

template <class T>
std::ostream &operator<<(std::ostream &os, const Matrix<T, 4> &m) {
  for (Index i = 0; i < m.dim1(); ++i)
    os << m[i];
  return os;
}
//-----------------------------------------------------------------------------

template <class T> std::istream &operator>>(std::istream &is, Matrix<T> &v) {

  for (Index i = 0; i < v.dim1(); ++i)
    is >> v(i);

  return is;
}

//-----------------------------------------------------------------------------

template <class T> std::istream &operator>>(std::istream &is, Matrix<T, 2> &m) {
  for (Index i = 0; i < m.dim1(); ++i) {
    Matrix<T, 1> tmp(m.dim2());
    is >> tmp;
    m[i] = tmp;
  }

  return is;
}

//-----------------------------------------------------------------------------

template <class T> std::istream &operator>>(std::istream &is, Matrix<T, 3> &m) {

  for (Index i = 0; i < m.dim1(); ++i) {
    Matrix<T, 2> tmp(m.dim2(), m.dim3());
    is >> tmp;
    m[i] = tmp;
  }

  return is;
}

//=============== vector stuff ==================
template <class T> std::istream &operator>>(std::istream &is, Vector3<T> &v) {

  for (Index i = 0; i < 3; ++i)
    is >> v[i];

  return is;
}
template <class T>
std::ostream &operator<<(std::ostream &os, const Vector3<T> &v) {
  Index i = 0;
  for (; i < v.size() - 1; ++i) {
    os << v[i] << " ";
  }
  os << v[i] << endl;
  return os;
}

} // namespace Numeric_lib
