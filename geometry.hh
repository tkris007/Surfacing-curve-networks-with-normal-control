#pragma once

#define NO_SURFACE_FIT

#include <array>
#include <list>
#include <memory>
#include <vector>


namespace Geometry
{

const double epsilon = 1.0e-8;

class Vector2D;
class Vector3D;
class BSCurve;
/*
using Point2D = Vector2D;
using Point3D = Vector3D;
using DoubleVector = std::vector<double>;
using Vector2DVector = std::vector<Vector2D>;
using VectorVector = std::vector<Vector3D>;
using Point2DVector = std::vector<Point2D>;
using PointVector = std::vector<Point3D>;
using CurveVector = std::vector<std::shared_ptr<BSCurve>>;
*/

class Vector2D {
public:
	// Constructors
	Vector2D();
	Vector2D ( double x, double y );

	// Assignments
	Vector2D &operator+= ( const Vector2D &v );
	Vector2D &operator-= ( const Vector2D &v );
	Vector2D &operator*= ( double x );
	Vector2D &operator/= ( double x );

	// Coordinates
	double &operator[] ( size_t i );
	const double &operator[] ( size_t i ) const;

	// Arithmetic
	Vector2D operator-() const;
	Vector2D operator+ ( const Vector2D &v ) const;
	Vector2D operator- ( const Vector2D &v ) const;
	double operator* ( const Vector2D &v ) const;
	Vector2D operator* ( double x ) const;
	Vector2D operator/ ( double x ) const;

	// Other
	double norm() const;
	double normSqr() const;
	Vector2D &normalize();

private:
	std::array<double, 2> v_;
};

class Vector3D {
public:
	// Constructors
	Vector3D();
	Vector3D ( double x, double y, double z );

	// Assignments
	Vector3D &operator+= ( const Vector3D &v );
	Vector3D &operator-= ( const Vector3D &v );
	Vector3D &operator*= ( double x );
	Vector3D &operator/= ( double x );

	// Coordinates
	double &operator[] ( size_t i );
	const double &operator[] ( size_t i ) const;

	// Arithmetic
	Vector3D operator-() const;
	Vector3D operator+ ( const Vector3D &v ) const;
	Vector3D operator- ( const Vector3D &v ) const;
	Vector3D operator^ ( const Vector3D &v ) const;
	double operator* ( const Vector3D &v ) const;
	Vector3D operator* ( double x ) const;
	Vector3D operator/ ( double x ) const;

	// Other
	double norm() const;
	double normSqr() const;
	Vector3D &normalize();

private:
	std::array<double, 3> v_;
};

class Matrix3x3 {
public:
	// Special matrices
	static Matrix3x3 identity();
	static Matrix3x3 rotation ( const Vector3D &axis, double angle );

	// Arithmetic
	Matrix3x3 operator+ ( const Matrix3x3 &m ) const;
	Matrix3x3 &operator+= ( const Matrix3x3 &m );
	Matrix3x3 operator* ( double x ) const;
	Matrix3x3 &operator*= ( double x );
	Vector3D operator* ( const Vector3D &v ) const;
	Matrix3x3 operator* ( const Matrix3x3 &m ) const;
	Matrix3x3 &operator*= ( const Matrix3x3 &m );

private:
	std::array<double, 9> m_;
};

class BSCurve {
public:
	// Constructors
	BSCurve();
	BSCurve ( const std::vector<Vector3D> &cpts );
	BSCurve ( size_t degree, const std::vector<double> &knots, const std::vector<Vector3D> &cpts );

	// Evaluation
	Vector3D eval ( double u ) const;
	Vector3D eval ( double u, size_t nr_der, std::vector<Vector3D> &der ) const;

	// Parameterization
	void reverse();
	void normalize();

	// Other
	double arcLength ( double from, double to ) const;
	std::vector<Vector3D> getControlPoints();
	std::vector<Vector3D> getControlPoints() const;
	size_t getDegree();

private:
	size_t findSpan ( double u ) const;
	void basisFunctions ( size_t i, double u, std::vector<double> &coeff ) const;
	void basisFunctionsAll ( size_t i, double u, std::vector<std::vector<double> > &coeff ) const;
	void derivativeControlPoints ( size_t d, size_t r1, size_t r2, std::vector<std::vector<Vector3D> > &dcp ) const;

	size_t p_;
	size_t n_;
	std::vector<double> knots_;
	std::vector<Vector3D> cp_;
};



} // namespace Geometry
