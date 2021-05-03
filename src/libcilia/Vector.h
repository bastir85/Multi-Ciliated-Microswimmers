#ifndef VECTOR_LIB_H
#define VECTOR_LIB_H

#include "Matrix.h"
#include <algorithm>
#include <iostream>
#include <string>

using namespace std;
// using namespace Numeric_lib;
namespace Numeric_lib {

// Faster class to do 3d vector operations!
template <class T> class Vector3 {

protected:
  T elem[3]; // vector? no: we couldn't easily provide a vector for a slice
public:
  Vector3()
  // matrix of n elements (default initialized)
  {
    // std::cerr << "new[" << n << "]->" << elem << "\n";
  }

  Vector3(double x, double y, double z) {
    elem[0] = x;
    elem[1] = y;
    elem[2] = z;
  }

  // copy constructor: let the base do the copy:
  Vector3(const Matrix<T, 1> &a) {
    a.range_check(2);
    elem[0] = a[0];
    elem[1] = a[1];
    elem[2] = a[2];
  }

  ~Vector3() {}

  T *data() { return elem; }
  const T *data() const { return elem; }
  static int const size() { return 3; };

  template <class F> void base_apply(F f) {
    for (Index i = 0; i < size(); ++i)
      f(elem[i]);
  }
  template <class F> void base_apply(F f, const Matrix_base<T> &mat) {
    for (Index i = 0; i < size(); ++i)
      f(elem[i], mat.data()[i]);
  }
  template <class F> void base_apply(F f, const Vector3<T> &vec) {
    for (Index i = 0; i < size(); ++i)
      f(elem[i], vec.data()[i]);
  }
  template <class F> void base_apply(F f, const T &c) {
    for (Index i = 0; i < size(); ++i)
      f(elem[i], c);
  }

  Vector3 &operator=(const Matrix<T, 1> &a)
  // copy assignment: let the base do the copy
  {
    a.range_check(2);
    elem[0] = a[0];
    elem[1] = a[1];
    elem[2] = a[2];
    return *this;
  }

  // slicing (the same as subscripting for 1D matrixs):
  T &operator[](int n) { return elem[n]; }
  const T &operator[](int n) const { return elem[n]; }

  // element-wise operations:
  template <class F> Vector3 &apply(F f) {
    this->base_apply(f);
    return *this;
  }
  template <class F> Vector3 &apply(F f, const T &c) {
    this->base_apply(f, c);
    return *this;
  }

  Vector3 &operator=(const T &c) {
    this->base_apply(Assign<T>(), c);
    return *this;
  }

  Vector3 &operator*=(const T &c) {
    this->base_apply(Mul_assign<T>(), c);
    return *this;
  }
  Vector3 &operator/=(const T &c) {
    this->base_apply(Div_assign<T>(), c);
    return *this;
  }
  Vector3 &operator%=(const T &c) {
    this->base_apply(Mod_assign<T>(), c);
    return *this;
  }
  Vector3 &operator+=(const T &c) {
    this->base_apply(Add_assign<T>(), c);
    return *this;
  }
  Vector3 &operator-=(const T &c) {
    this->base_apply(Minus_assign<T>(), c);
    return *this;
  }

  Vector3 &operator&=(const T &c) {
    this->base_apply(And_assign<T>(), c);
    return *this;
  }
  Vector3 &operator|=(const T &c) {
    this->base_apply(Or_assign<T>(), c);
    return *this;
  }
  Vector3 &operator^=(const T &c) {
    this->base_apply(Xor_assign<T>(), c);
    return *this;
  }

  // Vector3d& operator+=(Matrix& a) { this->base_apply(Add_assign<T>(), a);
  // return *this; }
  Vector3 &operator+=(Matrix<T> a) {
    this->base_apply(Add_assign<T>(), a);
    return *this;
  }
  Vector3 &operator-=(Matrix<T> a) {
    this->base_apply(Minus_assign<T>(), a);
    return *this;
  }

  Vector3 &operator+=(Vector3 a) {
    this->base_apply(Add_assign<T>(), a);
    return *this;
  }
  Vector3 &operator-=(Vector3 a) {
    this->base_apply(Minus_assign<T>(), a);
    return *this;
  }
};

typedef Vector3<double> Vector3d;

//-----------------------------------------------------------------------------

/*Vector3d operator*(const Matrix<double,1>& m, const double& c) {Vector3d r(m);
return r*=c;} Vector3d operator/(const Matrix<double,1>& m, const double& c)
{Vector3d r(m); return r/=c;} Vector3d operator+(const Matrix<double,1>& m,
const double& c) {Vector3d r(m); return r+=c;} Vector3d operator-(const
Matrix<double,1>& m, const double& c) {Vector3d r(m); return r-=c;}


Vector3d operator+(const Vector3d& m, const double& c) {Vector3d r(m); return
r+=c;} Vector3d operator-(const Vector3d& m, const double& c) {Vector3d r(m);
return r-=c;}
*/

template <class T> Vector3<T> operator*(const Vector3<T> &m, const double &c) {
  Vector3<T> r(m);
  return r *= c;
}
template <class T> Vector3<T> operator*(const double &c, const Vector3<T> &m) {
  Vector3<T> r(m);
  return r *= c;
}
template <class T> Vector3<T> operator/(const Vector3<T> &m, const double &c) {
  Vector3<T> r(m);
  return r /= c;
}
template <class T> Vector3<T> operator+(const Vector3<T> &m, const double &c) {
  Vector3<T> r(m);
  return r += c;
}
template <class T> Vector3<T> operator-(const Vector3<T> &m, const double &c) {
  Vector3<T> r(m);
  return r -= c;
}

template <class T>
Vector3<T> operator+(const Vector3<T> &m, const Vector3<T> &m2) {
  Vector3<T> r(m);
  r.base_apply(Add_assign<double>(), m2);
  return r;
}
template <class T>
Vector3<T> operator-(const Vector3<T> &m, const Vector3<T> &m2) {
  Vector3<T> r(m);
  r.base_apply(Minus_assign<double>(), m2);
  return r;
}

template <class T>
Vector3<T> operator+(const Matrix<T, 1> &m, const Vector3<T> &m2) {
  Vector3<T> r(m);
  r.base_apply(Add_assign<double>(), m2);
  return r;
}
template <class T>
Vector3<T> operator-(const Matrix<T, 1> &m, const Vector3<T> &m2) {
  Vector3<T> r(m);
  r.base_apply(Minus_assign<double>(), m2);
  return r;
}
// commutative
template <class T>
Vector3<T> operator+(const Vector3<T> &m2, const Matrix<T, 1> &m) {
  Vector3<T> r(m);
  r.base_apply(Add_assign<double>(), m2);
  return r;
}
template <class T>
Vector3<T> operator-(const Vector3<T> &m2, const Matrix<T, 1> &m) {
  Vector3<T> r(m);
  r.base_apply(Minus_assign<double>(), m2);
  return r;
}

template <class T>
Vector3<T> operator/(const Vector3<T> &m2, const Matrix<T, 1> &m) {
  Vector3<T> r(m2);
  r.base_apply(Div_assign<double>(), m);
  return r;
}

template <class T>
Vector3<T> operator+(const Matrix<T, 1> &m, const Matrix<T, 1> &m2) {
  Vector3<T> r(m);
  r.base_apply(Add_assign<double>(), m2);
  return r;
}
template <class T>
Vector3<T> operator-(const Matrix<T, 1> &m, const Matrix<T, 1> &m2) {
  Vector3<T> r(m);
  r.base_apply(Minus_assign<double>(), m2);
  return r;
}

template <class T> T dot_product(const Vector3<T> &a, const Vector3<T> &b) {
  T sum = 0;
  for (Index i = 0; i < a.size(); ++i)
    sum += a[i] * b[i];
  return sum;
}

template <class T> T dot_product(const Matrix<T, 1> &a, const Vector3<T> &b) {
  T sum = 0;
  for (Index i = 0; i < a.size(); ++i)
    sum += a[i] * b[i];
  return sum;
}

template <class T>
void add_to_matrix(Matrix<T, 1> &&mat, const Vector3<T> &vec, const double c) {
  mat.data()[0] += vec[0] * c;
  mat.data()[1] += vec[1] * c;
  mat.data()[2] += vec[2] * c;
}

template <class T> void to_matrix(Matrix<T, 1> &&mat, const Vector3<T> &vec) {
  mat.data()[0] = vec[0];
  mat.data()[1] = vec[1];
  mat.data()[2] = vec[2];
}

template <class T>
void to_matrix(Matrix<T, 1> &&mat, const Vector3<T> &vec, const double c) {
  mat.data()[0] = vec[0] * c;
  mat.data()[1] = vec[1] * c;
  mat.data()[2] = vec[2] * c;
}

} // namespace Numeric_lib
#endif
