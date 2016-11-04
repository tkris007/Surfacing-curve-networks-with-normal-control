#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <QDebug>

#include <QtGui/QKeyEvent>

#include "MyViewer.h"

#ifdef _WIN32
#define GL_CLAMP_TO_EDGE 0x812F
#define GL_BGRA 0x80E1
#endif

using qglviewer::Vec;

MyViewer::MyViewer ( QWidget *parent ) :
	QGLViewer ( parent ),
	mean_min ( 0.0 ), mean_max ( 0.0 ), cutoff_ratio ( 0.05 ),
	show_control_points ( true ), show_solid ( true ), show_wireframe ( false ), coloring ( COLOR_PLAIN )
{
	setSelectRegionWidth ( 10 );
	setSelectRegionHeight ( 10 );
	axes.shown = false;
}

MyViewer::~MyViewer()
{
	glDeleteTextures ( 1, &isophote_texture );
}

void MyViewer::updateMeanMinMax()
{
	size_t n = mesh.n_vertices();
	if ( n == 0 )
	{ return; }

	std::vector<double> mean;
	mean.reserve ( n );
	for ( MyMesh::VertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i )
	{ mean.push_back ( mesh.data ( *i ).mean ); }

	std::sort ( mean.begin(), mean.end() );
	size_t k = ( double ) n * cutoff_ratio;
	mean_min = std::min ( mean[k - 1], 0.0 );
	mean_max = std::max ( mean[n - k], 0.0 );
}

void MyViewer::updateMeanCurvature ( bool update_min_max )
{
	for ( MyMesh::ConstFaceIter i = mesh.faces_begin(), ie = mesh.faces_end(); i != ie; ++i )
	{
		MyMesh::HalfedgeHandle h1 = mesh.halfedge_handle ( *i );
		MyMesh::HalfedgeHandle h2 = mesh.next_halfedge_handle ( h1 );
		mesh.data ( *i ).area = ( halfedgeVector ( h1 ) % halfedgeVector ( h2 ) ).norm() / 2.0;
	}

	// Compute triangle strip areas
	for ( MyMesh::VertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i )
	{
		mesh.data ( *i ).area = 0;
		mesh.data ( *i ).mean = 0;
		for ( MyMesh::ConstVertexFaceIter j ( mesh, *i ); j.is_valid(); ++j )
		{ mesh.data ( *i ).area += mesh.data ( *j ).area; }
	}

	/*for ( MyMesh::VertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i )
	{
	double E, F, G, L, M, N, H;
	size_t const n = degree[0], m = degree[1];

	std::vector<double> coeff_Su, coeff_v, coeff_Sv, coeff_u, coeff_Suu, coeff_Svv;

	double u = mesh.data ( *i ).u;
	double v = mesh.data ( *i ).v;

	bernsteinAll ( n - 1, u, coeff_Su );
	bernsteinAll ( m, v, coeff_v );
	bernsteinAll ( n , u, coeff_u );
	bernsteinAll ( m - 1, v, coeff_Sv );
	bernsteinAll ( n - 2, u, coeff_Suu );
	bernsteinAll ( m - 2, v, coeff_Svv );


	Vec Su ( 0.0, 0.0, 0.0 );
	Vec Sv ( 0.0, 0.0, 0.0 );

	Vec Suu ( 0.0, 0.0, 0.0 );
	Vec Svv ( 0.0, 0.0, 0.0 );
	Vec Suv ( 0.0, 0.0, 0.0 );

	Vec normal;

	for ( size_t k = 0; k <= n - 1; ++k )
	for ( size_t l = 0; l <= m; ++l )
	{

	size_t const index1 = k * ( m + 1 ) + l;
	size_t const index2 = ( k + 1 ) * ( m + 1 ) + l;
	Su += ( control_points[index2] - control_points[index1] ) * coeff_Su[k] * coeff_v[l];
	}
	Su *= n;

	for ( size_t k = 0; k <= n; ++k )
	for ( size_t l = 0; l <= m - 1; ++l )
	{

	size_t const index1 = k * ( m + 1 ) + l;
	size_t const index2 = ( k ) * ( m + 1 ) + l + 1;
	Sv += ( control_points[index2] - control_points[index1] ) * coeff_u[k] * coeff_Sv[l];
	}
	Sv *= m;

	for ( size_t k = 0; k <= n - 2; ++k )
	for ( size_t l = 0; l <= m; ++l )
	{

	size_t const index1 = k * ( m + 1 ) + l;
	size_t const index2 = ( k + 1 ) * ( m + 1 ) + l;
	size_t const index3 = ( k + 2 ) * ( m + 1 ) + l;
	Suu += ( control_points[index3] - 2.0 * control_points[index2] +  control_points[index1] ) * coeff_Suu[k] * coeff_v[l];
	}

	Suu *= n * ( n - 1 );

	for ( size_t k = 0; k <= n ; ++k )
	for ( size_t l = 0; l <= m - 2; ++l )
	{

	size_t const index1 = k * ( m + 1 ) + l;
	size_t const index2 = k * ( m + 1 ) + l + 1;
	size_t const index3 = k * ( m + 1 ) + l + 2;
	Svv += ( control_points[index3] - 2.0 * control_points[index2] + control_points[index1] ) * coeff_u[k] * coeff_Svv[l];
	}

	Svv *= m * ( m - 1 );

	for ( size_t k = 0; k <= n - 1; ++k )
	for ( size_t l = 0; l <= m - 1; ++l )
	{

	size_t const index1 = k * ( m + 1 ) + l;
	size_t const index2 = ( k + 1 ) * ( m + 1 ) + l + 1;

	Suv += ( control_points[index2] + control_points[index1] ) * coeff_Su[k] * coeff_Sv[l];
	}

	Suv *= m * n;
	normal = cross ( Su, Sv );
	normal.normalize();
	E = Su * Su;
	F = Su * Sv;
	G = Sv * Sv;

	L = normal * Suu;
	M = normal * Suv;
	N = normal * Svv;

	H = ( N * E - 2 * M * F + L * G ) / ( 2 * ( E * G - F * F ) );
	mesh.data ( *i ).mean = H;
	}*/

	// Compute mean values using normal difference angles
	for ( MyMesh::VertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i )
	{
		for ( MyMesh::ConstVertexEdgeIter j ( mesh, *i ); j.is_valid(); ++j )
		{
			double angle;
			MyMesh::HalfedgeHandle h1 = mesh.halfedge_handle ( *j, 0 );
			MyMesh::HalfedgeHandle h2 = mesh.halfedge_handle ( *j, 1 );
			Vector v = halfedgeVector ( h1 );
			if ( mesh.is_boundary ( h1 ) || mesh.is_boundary ( h2 ) )
			{ angle = 0.0; }
			else
			{
				Vector n1 = mesh.normal ( mesh.face_handle ( h1 ) );
				Vector n2 = mesh.normal ( mesh.face_handle ( h2 ) );
				angle = acos ( std::min ( std::max ( n1 | n2, -1.0f ), 1.0f ) );
				angle *= ( ( n1 % n2 ) | v ) >= 0.0 ? 1.0 : -1.0;
			}
			mesh.data ( *i ).mean += angle * v.norm();
		}
		mesh.data ( *i ).mean *= 3.0 / 4.0 / mesh.data ( *i ).area;
	}

	if ( update_min_max )
	{ updateMeanMinMax(); }
}

void MyViewer::meanMapColor ( double d, double *color ) const
{
	if ( d <= mean_min )
	{
		color[0] = 0.0;
		color[1] = 0.0;
		color[2] = 1.0;
	}
	else if ( d >= mean_max )
	{
		color[0] = 1.0;
		color[1] = 0.0;
		color[2] = 0.0;
	}
	else if ( d < 0 )
	{
		double alpha = d / mean_min;
		color[0] = 0.0;
		color[1] = 1.0 - alpha;
		color[2] = alpha;
	}
	else
	{
		double alpha = d / mean_max;
		color[0] = alpha;
		color[1] = 1.0 - alpha;
		color[2] = 0;
	}
}

//bool MyViewer::openBezier ( std::string const &filename )
//{
//	size_t n, m;
//	try
//	{
//		std::ifstream f ( filename.c_str() );
//		f >> n >> m;
//		degree[0] = n++; degree[1] = m++;
//		control_points.resize ( n * m );
//		for ( size_t i = 0; i < n; ++i )
//			for ( size_t j = 0; j < m; ++j )
//			{
//				size_t const index = i * m + j;
//				f >> control_points[index][0] >> control_points[index][1] >> control_points[index][2];
//			}
//		f.close();
//	}
//	catch ( std::ifstream::failure )
//	{
//		return false;
//	}
//
//	generateMesh();
//	mesh.request_face_normals(); mesh.request_vertex_normals();
//	mesh.update_face_normals();  mesh.update_vertex_normals();
//
//	updateMeanCurvature();
//
//	// Set camera on the model
//	MyMesh::Point box_min, box_max;
//	box_min = box_max = mesh.point ( *mesh.vertices_begin() );
//	for ( MyMesh::ConstVertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i )
//	{
//		box_min.minimize ( mesh.point ( *i ) );
//		box_max.maximize ( mesh.point ( *i ) );
//	}
//	camera()->setSceneBoundingBox ( Vec ( box_min[0], box_min[1], box_min[2] ),
//	                                Vec ( box_max[0], box_max[1], box_max[2] ) );
//	camera()->showEntireScene();
//
//	setSelectedName ( -1 );
//	axes.shown = false;
//
//	updateGL();
//	return true;
//}

//bool MyViewer::saveBezier ( std::string const & filename )
//{
//	size_t n = degree[0], m = degree[1];
//	try
//	{
//		std::ofstream f ( filename.c_str() );
//		f << n++ << " " << m++ << std::endl;
//		for ( size_t i = 0; i < n; ++i )
//			for ( size_t j = 0; j < m; ++j )
//			{
//				size_t const index = i * m + j;
//				f << control_points[index][0] << " " << control_points[index][1] << " " << control_points[index][2] << std::endl;
//			}
//		f.close();
//		return true;
//	}
//	catch ( std::ofstream::failure )
//	{
//		return false;
//	}
//}

std::istringstream MyViewer::nextLine(std::ifstream &file) {

	std::string line;
	std::getline(file, line); 
	// Ignore empty lines
	while(!file.eof() && line.length() == 0) {

		std::getline(file, line); 
	}
	std::istringstream ss(line);

	QString qstr = QString::fromStdString(line);
	qDebug() << qstr;

	return ss;
}

void MyViewer::readBSCurve(std::ifstream &file) {

	// Read degree
	int degree;
	nextLine(file) >> degree;

	// Read knots
	std::vector<double> knots;
	int numberOfKnots;
	nextLine(file) >> numberOfKnots;
	std::istringstream ss = nextLine(file);
	for(int i = 0; i < numberOfKnots; i++) {

		double knot;
		ss >> knot;
		knots.push_back(knot);
	}

	// Read control points
	std::vector<Geometry::Vector3D> cpts;
	int numberOfCpts;
	nextLine(file) >> numberOfCpts;
	for(int i = 0; i < numberOfCpts; i++) {

		std::istringstream ss = nextLine(file);
		double x; 
		ss >> x;
		double y; 
		ss >> y;
		double z; 
		ss >> z;
		cpts.push_back(Geometry::Vector3D(x,y,z));
	}

	// Create new BSCurve and add it to bsCurves
	bsCurves.push_back(Geometry::BSCurve ( degree, knots, cpts ));
}

bool MyViewer::openBSpline ( std::string const &filename ) {

	try {

		std::ifstream file ( filename.c_str() );
		if (file.good())
		{
			int numberOfCurves;
			nextLine(file) >>  numberOfCurves;
			// ... you now get a number ...
			qDebug() << numberOfCurves;

			for(int i = 0; i < numberOfCurves; i++)
			{
				readBSCurve(file);
			}
		}

		return true;

	} catch (std::ifstream::failure) {

		return false;
	}
}

bool MyViewer::saveBSpline ( std::string const& filename ) {

	return true;
}


void MyViewer::init()
{
	glLightModeli ( GL_LIGHT_MODEL_TWO_SIDE, 1 );
	QImage img ( ":/isophotes.png" );
	glGenTextures ( 1, &isophote_texture );
	glBindTexture ( GL_TEXTURE_2D, isophote_texture );
	glTexParameteri ( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
	glTexParameteri ( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
	glTexParameteri ( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE );
	glTexParameteri ( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE );
	glTexImage2D ( GL_TEXTURE_2D, 0, GL_RGBA8, img.width(), img.height(), 0, GL_BGRA,
		GL_UNSIGNED_BYTE, img.convertToFormat ( QImage::Format_ARGB32 ).bits() );
}

void MyViewer::draw()
{
	if ( mesh.n_vertices() == 0 )
	{ return; }

	if ( show_control_points )
	{
		glDisable ( GL_LIGHTING );
		glLineWidth ( 3.0 );
		glColor3d ( 0.3, 0.3, 1.0 );
		size_t const m = degree[0] + 1;
		for ( size_t k = 0; k < 2; ++k )
			for ( size_t i = 0; i <= degree[k]; ++i )
			{
				glBegin ( GL_LINE_STRIP );
				for ( size_t j = 0; j <= degree[1 - k]; ++j )
				{
					size_t const index = k ? j * m + i : i * m + j;
					Vec const &p = control_points[index];
					glVertex3d ( p[0], p[1], p[2] );
				}
				glEnd();
			}
			glLineWidth ( 1.0 );
			glPointSize ( 8.0 );
			glColor3d ( 1.0, 0.0, 1.0 );
			glBegin ( GL_POINTS );
			for ( size_t i = 0, ie = control_points.size(); i < ie; ++i )
			{
				Vec const &p = control_points[i];
				glVertex3d ( p[0], p[1], p[2] );
			}
			glEnd();
			glPointSize ( 1.0 );
			glEnable ( GL_LIGHTING );
	}

	if ( !show_solid && show_wireframe )
	{ glPolygonMode ( GL_FRONT_AND_BACK, GL_LINE ); }
	else
	{ glPolygonMode ( GL_FRONT_AND_BACK, GL_FILL ); }

	glEnable ( GL_POLYGON_OFFSET_FILL );
	glPolygonOffset ( 1, 1 );

	std::vector<double> color ( 3, 1.0 );
	if ( show_solid || show_wireframe )
	{
		if ( coloring == COLOR_PLAIN )
		{ glColor3dv ( &color[0] ); }
		else if ( coloring == COLOR_ISOPHOTES )
		{
			glBindTexture ( GL_TEXTURE_2D, isophote_texture );
			glTexEnvf ( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL );
			glEnable ( GL_TEXTURE_2D );
			glTexGeni ( GL_S, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP );
			glTexGeni ( GL_T, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP );
			glEnable ( GL_TEXTURE_GEN_S );
			glEnable ( GL_TEXTURE_GEN_T );
		}
		for ( MyMesh::ConstFaceIter i = mesh.faces_begin(), ie = mesh.faces_end(); i != ie; ++i )
		{
			glBegin ( GL_POLYGON );
			for ( MyMesh::ConstFaceVertexIter j ( mesh, *i ); j.is_valid(); ++j )
			{
				if ( coloring == COLOR_MEAN )
				{
					meanMapColor ( mesh.data ( *j ).mean, &color[0] );
					glColor3dv ( &color[0] );
				}
				glNormal3fv ( mesh.normal ( *j ).data() );
				glVertex3fv ( mesh.point ( *j ).data() );
			}
			glEnd();
		}
		if ( coloring == COLOR_ISOPHOTES )
		{
			glDisable ( GL_TEXTURE_GEN_S );
			glDisable ( GL_TEXTURE_GEN_T );
			glDisable ( GL_TEXTURE_2D );
			glTexEnvf ( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
		}
	}

	if ( show_solid && show_wireframe )
	{
		glPolygonMode ( GL_FRONT, GL_LINE );
		glColor3d ( 0.0, 0.0, 0.0 );
		glDisable ( GL_LIGHTING );
		for ( MyMesh::ConstFaceIter i = mesh.faces_begin(), ie = mesh.faces_end(); i != ie; ++i )
		{
			glBegin ( GL_POLYGON );
			for ( MyMesh::ConstFaceVertexIter j ( mesh, *i ); j.is_valid(); ++j )
			{ glVertex3fv ( mesh.point ( *j ).data() ); }
			glEnd();
		}
		glEnable ( GL_LIGHTING );
	}

	if ( axes.shown )
	{ drawAxes(); }
}

void MyViewer::drawAxes() const
{
	Vec const p ( axes.position[0], axes.position[1], axes.position[2] );
	glColor3f ( 1.0, 0.0, 0.0 );
	drawArrow ( p, p + Vec ( axes.size, 0.0, 0.0 ), axes.size / 50.0 );
	glColor3f ( 0.0, 1.0, 0.0 );
	drawArrow ( p, p + Vec ( 0.0, axes.size, 0.0 ), axes.size / 50.0 );
	glColor3f ( 0.0, 0.0, 1.0 );
	drawArrow ( p, p + Vec ( 0.0, 0.0, axes.size ), axes.size / 50.0 );
	glEnd();
}

void MyViewer::drawWithNames()
{
	if ( axes.shown )
	{ drawAxesWithNames(); }
	else
	{
		if ( !show_control_points )
		{ return; }

		for ( size_t i = 0, ie = control_points.size(); i < ie; ++i )
		{
			Vec const &p = control_points[i];
			glPushName ( static_cast<GLuint> ( i ) );
			glRasterPos3f ( p[0], p[1], p[2] );
			glPopName();
		}
	}
}

void MyViewer::drawAxesWithNames() const
{
	Vec const p ( axes.position[0], axes.position[1], axes.position[2] );
	glPushName ( 0 );
	drawArrow ( p, p + Vec ( axes.size, 0.0, 0.0 ), axes.size / 50.0 );
	glPopName();
	glPushName ( 1 );
	drawArrow ( p, p + Vec ( 0.0, axes.size, 0.0 ), axes.size / 50.0 );
	glPopName();
	glPushName ( 2 );
	drawArrow ( p, p + Vec ( 0.0, 0.0, axes.size ), axes.size / 50.0 );
	glPopName();
}

void MyViewer::postSelection ( const QPoint &p )
{
	int sel = selectedName();
	if ( sel == -1 )
	{
		axes.shown = false;
		return;
	}

	if ( axes.shown )
	{
		axes.selected_axis = sel;
		bool found;
		axes.grabbed_pos = camera()->pointUnderPixel ( p, found );
		axes.original_pos[0] = axes.position[0];
		axes.original_pos[1] = axes.position[1];
		axes.original_pos[2] = axes.position[2];
		if ( !found )
		{ axes.shown = false; }
	}
	else
	{
		selected = sel;
		axes.position[0] = control_points[sel][0];
		axes.position[1] = control_points[sel][1];
		axes.position[2] = control_points[sel][2];
		Vec const pos ( axes.position[0], axes.position[1], axes.position[2] );
		double const depth = camera()->projectedCoordinatesOf ( pos ) [2];
		Vec const q1 = camera()->unprojectedCoordinatesOf ( Vec ( 0.0, 0.0, depth ) );
		Vec const q2 = camera()->unprojectedCoordinatesOf ( Vec ( width(), height(), depth ) );
		axes.size = ( q1 - q2 ).norm() / 10.0;
		axes.shown = true;
		axes.selected_axis = -1;
	}
}

void MyViewer::keyPressEvent ( QKeyEvent *e )
{
	if ( e->modifiers() == Qt::NoModifier )
		switch ( e->key() )
	{
		case Qt::Key_P:
			coloring = COLOR_PLAIN;
			updateGL();
			break;
		case Qt::Key_M:
			coloring = COLOR_MEAN;
			updateGL();
			break;
		case Qt::Key_I:
			coloring = COLOR_ISOPHOTES;
			updateGL();
			break;
		case Qt::Key_C:
			show_control_points = !show_control_points;
			updateGL();
			break;
		case Qt::Key_S:
			show_solid = !show_solid;
			updateGL();
			break;
		case Qt::Key_W:
			show_wireframe = !show_wireframe;
			updateGL();
			break;
		case Qt::Key_E:
			increaseDegree();
			updateGL();
			break;
		default:
			QGLViewer::keyPressEvent ( e );
	}
	else
	{ QGLViewer::keyPressEvent ( e ); }
}

void MyViewer::increaseDegree()
{
	size_t n = degree[0] + 2, m = degree[1] + 1;
	++degree[0];
	++degree[1];
	std::vector<qglviewer::Vec> tmpControlPoints = control_points;
	control_points.resize ( ( n ) * ( m  ) );
	tmpControlPoints.resize ( ( n  ) * ( m  ) );

	for ( size_t j = 0; j < m; ++j )
	{
		control_points[j] = tmpControlPoints[j];
		for ( size_t i = 1; i < n - 1; ++i )
		{
			size_t const index = i * m + j;
			control_points[index] = tmpControlPoints[index - m] + ( tmpControlPoints[index] - tmpControlPoints[index - m] ) / n * ( n - i );
		}
		control_points[ ( n - 1 ) *m + j] = tmpControlPoints[ ( n - 2 ) * m + j];

	}
	tmpControlPoints = control_points;
	++m;
	control_points.resize ( n * m );
	for ( size_t i = 0; i < n; ++i )
	{
		control_points[i * m] = tmpControlPoints[i * ( m - 1 )];
		for ( size_t j = 1; j < m - 1; ++j )
		{
			size_t const index = i * m + j;
			control_points[index] = tmpControlPoints[i * ( m - 1 ) + j - 1] + ( tmpControlPoints[i * ( m - 1 ) + j] - tmpControlPoints[i * ( m - 1 ) + j - 1] ) / m * ( m - j );
		}
		control_points[ ( i + 1 ) *m - 1] = tmpControlPoints[ ( i + 1 ) * ( m - 1 ) - 1];
	}

	generateMesh();
	mesh.request_face_normals(); mesh.request_vertex_normals();
	mesh.update_face_normals();  mesh.update_vertex_normals();

	updateMeanCurvature();
}

Vec MyViewer::intersectLines ( Vec const &ap, Vec const &ad, Vec const &bp, Vec const &bd ) const
{
	// always returns a point on the (ap, ad) line
	double a = ad * ad, b = ad * bd, c = bd * bd;
	double d = ad * ( ap - bp ), e = bd * ( ap - bp );
	if ( a * c - b * b < 1.0e-7 )
	{ return ap; }
	double s = ( b * e - c * d ) / ( a * c - b * b );
	return ap + s * ad;
}

void MyViewer::bernsteinAll ( size_t n, double u, std::vector<double> &coeff )
{
	coeff.clear(); coeff.reserve ( n + 1 );
	coeff.push_back ( 1.0 );
	double const u1 = 1.0 - u;
	for ( size_t j = 1; j <= n; ++j )
	{
		double saved = 0.0;
		for ( size_t k = 0; k < j; ++k )
		{
			double const tmp = coeff[k];
			coeff[k] = saved + tmp * u1;
			saved = tmp * u;
		}
		coeff.push_back ( saved );
	}
}

void MyViewer::generateMesh()
{
	size_t const resolution = 30;

	mesh.clear();
	std::vector<MyMesh::VertexHandle> handles, tri;
	size_t const n = degree[0], m = degree[1];

	std::vector<double> coeff_u, coeff_v;
	for ( size_t i = 0; i < resolution; ++i )
	{
		double u = ( double ) i / ( double ) ( resolution - 1 );
		bernsteinAll ( n, u, coeff_u );
		for ( size_t j = 0; j < resolution; ++j )
		{
			double v = ( double ) j / ( double ) ( resolution - 1 );
			bernsteinAll ( m, v, coeff_v );
			Vec p ( 0.0, 0.0, 0.0 );
			for ( size_t k = 0; k <= n; ++k )
				for ( size_t l = 0; l <= m; ++l )
				{
					size_t const index = k * ( m + 1 ) + l;
					p += control_points[index] * coeff_u[k] * coeff_v[l];
				}
				handles.push_back ( mesh.add_vertex ( MyMesh::Point ( p[0], p[1], p[2] ) ) );
				mesh.data ( handles[handles.size() - 1] ).u = u;
				mesh.data ( handles[handles.size() - 1] ).v = v;
		}
	}
	for ( size_t i = 0; i < resolution - 1; ++i )
		for ( size_t j = 0; j < resolution - 1; ++j )
		{
			tri.clear();
			tri.push_back ( handles[i * resolution + j] );
			tri.push_back ( handles[i * resolution + j + 1] );
			tri.push_back ( handles[ ( i + 1 ) * resolution + j] );
			mesh.add_face ( tri );
			tri.clear();
			tri.push_back ( handles[ ( i + 1 ) * resolution + j] );
			tri.push_back ( handles[i * resolution + j + 1] );
			tri.push_back ( handles[ ( i + 1 ) * resolution + j + 1] );
			mesh.add_face ( tri );
		}
}

void MyViewer::mouseMoveEvent ( QMouseEvent *e )
{
	if ( axes.shown && axes.selected_axis >= 0 &&
		e->modifiers() & Qt::ShiftModifier && e->buttons() & Qt::LeftButton )
	{
		Vec axis = Vec ( axes.selected_axis == 0, axes.selected_axis == 1, axes.selected_axis == 2 );
		Vec from, dir;
		camera()->convertClickToLine ( e->pos(), from, dir );
		Vec p = intersectLines ( axes.grabbed_pos, axis, from, dir );
		float d = ( p - axes.grabbed_pos ) * axis;
		axes.position[axes.selected_axis] = axes.original_pos[axes.selected_axis] + d;
		control_points[selected] = Vec ( axes.position[0], axes.position[1], axes.position[2] );
		generateMesh();
		mesh.request_face_normals(); mesh.request_vertex_normals();
		mesh.update_face_normals();  mesh.update_vertex_normals();
		updateMeanCurvature();
		updateGL();
	}
	else
	{ QGLViewer::mouseMoveEvent ( e ); }
}

QString MyViewer::helpString() const
{
	QString text ( "<h2>Sample Framework</h2>"
		"<p>This is a minimal framework for 3D mesh manipulation, which can be "
		"extended and used as a base for various projects, for example "
		"prototypes for fairing algorithms, or even displaying/modifying "
		"parametric surfaces, etc.</p>"
		"<p>The following hotkeys are available:</p>"
		"<ul>"
		"<li>&nbsp;P: Set plain map (no coloring)</li>"
		"<li>&nbsp;M: Set mean curvature map</li>"
		"<li>&nbsp;I: Set isophote line map</li>"
		"<li>&nbsp;C: Toggle control polygon visualization</li>"
		"<li>&nbsp;S: Toggle solid (filled polygon) visualization</li>"
		"<li>&nbsp;W: Toggle wireframe visualization</li>"
		"</ul>"
		"<p>There is also a simple selection and movement interface, enabled "
		"only when control points are displayed: a control point can be selected "
		"by shift-clicking, and it can be moved by shift-dragging one of the "
		"displayed axes.</p>"
		"<p>Note that libQGLViewer is furnished with a lot of useful features, "
		"such as storing/loading view positions, or saving screenshots. "
		"OpenMesh also has a nice collection of tools for mesh manipulation: "
		"decimation, subdivision, smoothing, etc. These can provide "
		"good comparisons to the methods you implement.</p>"
		"<p>This software can be used as a sample GUI base for handling "
		"parametric or procedural surfaces, as well. The power of "
		"Qt and libQGLViewer makes it easy to set up a prototype application. "
		"Feel free to modify and explore!</p>"
		"<p align=\"right\">Peter Salvi</p>" );
	return text;
}
