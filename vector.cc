#include <cmath>

#include "geometry.hh"

namespace Geometry {

Vector2D::Vector2D() {
}

Vector2D::Vector2D(double x, double y) {
   v_[0] = x;
   v_[1] = y; 
}

Vector2D &
Vector2D::operator+=(const Vector2D &v) {
  v_[0] += v.v_[0];
  v_[1] += v.v_[1];
  return *this;
}

Vector2D &
Vector2D::operator-=(const Vector2D &v) {
  v_[0] -= v.v_[0];
  v_[1] -= v.v_[1];
  return *this;
}

Vector2D &
Vector2D::operator*=(double x) {
  v_[0] *= x;
  v_[1] *= x;
  return *this;
}

Vector2D &
Vector2D::operator/=(double x) {
  v_[0] /= x;
  v_[1] /= x;
  return *this;
}

double &
Vector2D::operator[](size_t i) {
  return v_[i];
}

const double &
Vector2D::operator[](size_t i) const {
  return v_[i];
}

Vector2D
Vector2D::operator-() const {
  return Vector2D(-v_[0], -v_[1]);
}

Vector2D
Vector2D::operator+(const Vector2D &v) const {
  return Vector2D(v_[0] + v.v_[0], v_[1] + v.v_[1]);
}

Vector2D
Vector2D::operator-(const Vector2D &v) const {
  return Vector2D(v_[0] - v.v_[0], v_[1] - v.v_[1]);
}

double
Vector2D::operator*(const Vector2D &v) const {
  return v_[0] * v.v_[0] + v_[1] * v.v_[1];
}

Vector2D
Vector2D::operator*(double x) const {
  return Vector2D(v_[0] * x, v_[1] * x);
}

Vector2D
Vector2D::operator/(double x) const {
  return Vector2D(v_[0] / x, v_[1] / x);
}

double Vector2D::norm() const {
  return std::sqrt(normSqr());
}

double
Vector2D::normSqr() const {
  return operator*(*this);
}

Vector2D &
Vector2D::normalize() {
  return operator*=(1.0 / norm());
}

Vector3D::Vector3D() {
}

Vector3D::Vector3D(double x, double y, double z) {
	v_[0] = x;
    v_[1] = y;
	v_[2] = z;
}

Vector3D &
Vector3D::operator+=(const Vector3D &v) {
  v_[0] += v.v_[0];
  v_[1] += v.v_[1];
  v_[2] += v.v_[2];
  return *this;
}

Vector3D &
Vector3D::operator-=(const Vector3D &v) {
  v_[0] -= v.v_[0];
  v_[1] -= v.v_[1];
  v_[2] -= v.v_[2];
  return *this;
}

Vector3D &
Vector3D::operator*=(double x) {
  v_[0] *= x;
  v_[1] *= x;
  v_[2] *= x;
  return *this;
}

Vector3D &
Vector3D::operator/=(double x) {
  v_[0] /= x;
  v_[1] /= x;
  v_[2] /= x;
  return *this;
}

double &
Vector3D::operator[](size_t i) {
  return v_[i];
}

const double &
Vector3D::operator[](size_t i) const {
  return v_[i];
}

Vector3D
Vector3D::operator-() const {
  return Vector3D(-v_[0], -v_[1], -v_[2]);
}

Vector3D
Vector3D::operator+(const Vector3D &v) const {
  return Vector3D(v_[0] + v.v_[0], v_[1] + v.v_[1], v_[2] + v.v_[2]);
}

Vector3D
Vector3D::operator-(const Vector3D &v) const {
  return Vector3D(v_[0] - v.v_[0], v_[1] - v.v_[1], v_[2] - v.v_[2]);
}

Vector3D
Vector3D::operator^(const Vector3D &v) const {
  return Vector3D(v_[1] * v.v_[2] - v_[2] * v.v_[1],
                  v_[2] * v.v_[0] - v_[0] * v.v_[2],
                  v_[0] * v.v_[1] - v_[1] * v.v_[0]);
}

double
Vector3D::operator*(const Vector3D &v) const {
  return v_[0] * v.v_[0] + v_[1] * v.v_[1] + v_[2] * v.v_[2];
}

Vector3D
Vector3D::operator*(double x) const {
  return Vector3D(v_[0] * x, v_[1] * x, v_[2] * x);
}

Vector3D
Vector3D::operator/(double x) const {
  return Vector3D(v_[0] / x, v_[1] / x, v_[2] / x);
}

double Vector3D::norm() const {
  return std::sqrt(normSqr());
}

double
Vector3D::normSqr() const {
  return operator*(*this);
}

Vector3D &
Vector3D::normalize() {
  return operator*=(1.0 / norm());
}

} // namespace Geometry
