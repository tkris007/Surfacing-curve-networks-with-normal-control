#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cctype>
#include <QDebug>

#include <QtGui/QKeyEvent>

#include "MyViewer.h"

#ifdef _WIN32
#define GL_CLAMP_TO_EDGE 0x812F
#define GL_BGRA 0x80E1
#endif

using qglviewer::Vec;
using Geometry::Vector3D;

double cot ( double radian )
{
	return  cos ( radian ) / sin ( radian );
}

MyViewer::MyViewer ( QWidget *parent ) :
	QGLViewer ( parent ),
	mean_min ( 0.0 ), mean_max ( 0.0 ), cutoff_ratio ( 0.05 ),
	show_control_points ( true ), show_solid ( true ), show_wireframe ( false ),
	show_curves ( true ), show_normals ( false ),
	coloring ( COLOR_PLAIN ), normalSize ( 10 )
{
	setSelectRegionWidth ( 10 );
	setSelectRegionHeight ( 10 );
	axes.shown = false;
	step = 0.01f;
	sampling = 10;
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

std::string MyViewer::nextLine ( std::ifstream &file )
{

	std::string line;
	std::getline ( file, line );
	// Ignore empty lines
	while ( !file.eof() && line.length() == 0 )
	{

		std::getline ( file, line );
	}

	QString qstr = QString::fromStdString ( line );
	// qDebug() << qstr;

	return line;
}

void MyViewer::readBSCurve ( std::ifstream &file )
{

	// Read degree
	int degree;
	std::stringstream ss ( nextLine ( file ) );
	ss >> degree;

	// Read knots
	std::vector<double> knots;
	int numberOfKnots;
	ss.str ( nextLine ( file ) );
	ss.seekg ( 0 );
	ss >> numberOfKnots;
	for ( int i = 0; i < numberOfKnots; i++ )
	{

		double knot;
		ss >> knot;
		knots.push_back ( knot );
	}

	// Read control points
	std::vector<Geometry::Vector3D> cpts;
	int numberOfCpts;
	ss.str ( nextLine ( file ) );
	ss.seekg ( 0 );
	ss >> numberOfCpts;
	for ( int i = 0; i < numberOfCpts; i++ )
	{
		double x;
		ss >> x;
		double y;
		ss >> y;
		double z;
		ss >> z;
		cpts.push_back ( Geometry::Vector3D ( x, y, z ) );
	}

	// Create new BSCurve and add it to bsCurves
	std::shared_ptr<Geometry::BSCurve> p_BSCurve ( new Geometry::BSCurve ( degree, knots, cpts ) );
	p_BSCurve->normalize();
	//p_BSCurve->reverse();
	bsCurves.push_back ( p_BSCurve );

}

bool MyViewer::openBSpline ( std::string const &filename )
{
	try
	{
		bsCurves.clear();
		pointsOnPlain.clear();
		normals.clear();
		std::ifstream file ( filename.c_str() );
		if ( file.good() )
		{
			int numberOfCurves;
			std::stringstream ( nextLine ( file ) ) >>  numberOfCurves;

			for ( int i = 0; i < numberOfCurves; i++ )
			{
				readBSCurve ( file );
			}
		}
		file.close();


		generateMesh();
		mesh.request_face_normals(); mesh.request_vertex_normals();
		mesh.update_face_normals();  mesh.update_vertex_normals();
		updateMeanCurvature();

		focusOnMesh();

		setSelectedName ( -1 );
		axes.shown = false;

		updateGL();
		return true;
	}
	catch ( std::ifstream::failure )
	{

		return false;
	}


}

bool MyViewer::saveBSpline ( std::string const& filename )
{

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
	if ( bsCurves.size() == 0 )
	{ return; }
	if ( show_curves )
	{ drawCurves(); }

	if ( show_control_points )
	{drawControlPoints();	}

	if ( show_normals )
	{
		drawNormals();
	}
	if ( !show_solid && show_wireframe )
	{ glPolygonMode ( GL_FRONT_AND_BACK, GL_LINE ); }
	else
	{ glPolygonMode ( GL_FRONT_AND_BACK, GL_FILL ); }

	glEnable ( GL_POLYGON_OFFSET_FILL );
	glPolygonOffset ( 1, 1 );

	std::vector<double> color ( 3, 1.0 );
	if ( !mesh.n_vertices() == 0 )
	{
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
	}
	if ( axes.shown )
	{ drawAxes(); }
}

void MyViewer::drawCurves() const
{
	glDisable ( GL_LIGHTING );
	glLineWidth ( 3.0 );
	glColor3d ( 0.3, 0.3, 1.0 );
	glBegin ( GL_LINE_STRIP );
	for ( size_t i = 0; i < bsCurves.size(); ++i )
		for ( float t = 0; t < 1; t += step )
		{
			Vector3D const &p = bsCurves[i]->eval ( t );
			glVertex3d ( p[0], p[1], p[2] );
		}
	glEnd();
	glLineWidth ( 1.0 );
	glEnable ( GL_LIGHTING );
}

void MyViewer::drawControlPoints() const
{
	glDisable ( GL_LIGHTING );
	glPointSize ( 8.0 );
	glColor3d ( 1.0, 0.0, 1.0 );
	glBegin ( GL_POINTS );
	for ( size_t i = 0; i < bsCurves.size(); ++i )
	{
		std::vector<Vector3D> controlPoints = bsCurves[i]->getControlPoints();
		size_t numberOfControlPoints = controlPoints.size();
		for ( size_t k = 0; k < numberOfControlPoints; ++k )
		{
			Vector3D const &p = controlPoints[k];
			glVertex3d ( p[0], p[1], p[2] );
		}

	}
	glColor3d ( 1.0, 1.0, 0.0 );
	for ( auto p : pointsOnPlain )
	{
		glVertex3d ( p[0], p[1], p[2] );
	}
	glVertex3d ( plainPoint[0], plainPoint[1], plainPoint[2] );
	glEnd();
	glPointSize ( 1.0 );
	glColor3d ( 0.0, 1.0, 0.0 );
	glLineWidth ( 3.0 );
	glBegin ( GL_LINE_STRIP );
	for ( size_t i = 0; i < bsCurves.size(); ++i )
	{
		std::vector<Vector3D> controlPoints = bsCurves[i]->getControlPoints();
		size_t numberOfControlPoints = controlPoints.size();
		for ( size_t k = 0; k < numberOfControlPoints; ++k )
		{
			Vector3D const &p = controlPoints[k];
			glVertex3d ( p[0], p[1], p[2] );
		}

	}
	glEnd();
	glLineWidth ( 1.0 );
	glPointSize ( 1.0 );
	glEnable ( GL_LIGHTING );

}

void MyViewer::drawNormals() const
{
	for ( size_t i = 0; i < bsCurves.size(); ++i )
	{
		int j = 0;
		for ( size_t t = 0; t < sampling; t += 1 )
		{

			/*Vector3D p1 = rmf.eval ( t ) ;*/
			Vector3D p2 = bsCurves[i]->eval ( ( double ) t / sampling );
			Vec const arrowEndPoint = normals.at ( i ) [t]; //.at ( i ) [j * 10];
			//arrowEndPoint.normalize();
			Vec const & arrowStartPoint = Vec ( p2[0], p2[1], p2[2] );

			glColor3f ( 1.0, 0.0, 0.0 );
			drawArrow ( arrowStartPoint, ( arrowStartPoint + arrowEndPoint * normalSize ) , normalSize / 50.0 );


			++j;

		}
	}
	Vec const arrowEndPoint = plainNormal;
	//arrowEndPoint.normalize();
	Vec const &arrowStartPoint = plainPoint;

	glColor3f ( 0.0, 1.0, 1.0 );
	drawArrow ( arrowStartPoint, ( arrowStartPoint + arrowEndPoint * normalSize ), normalSize / 50.0 );
	glEnd();
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
		for ( size_t i = 0; i < bsCurves.size(); ++i )
		{
			std::vector<Vector3D> controlPoints = bsCurves[i]->getControlPoints();
			size_t numberOfControlPoints = controlPoints.size();
			for ( size_t k = 0; k < numberOfControlPoints; ++k )
			{
				Vector3D const &p = controlPoints[k];
				glPushName ( static_cast<GLuint> ( i * 1000 + k ) );
				glRasterPos3f ( p[0], p[1], p[2] );
				glPopName();
			}
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

void MyViewer::postSelection ( const QPoint & p )
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
		axes.position[0] = bsCurves[sel / 1000]->getControlPoints() [sel % 1000][0];
		axes.position[1] = bsCurves[sel / 1000]->getControlPoints() [sel % 1000][1];
		axes.position[2] = bsCurves[sel / 1000]->getControlPoints() [sel % 1000][2];
		Vec const pos ( axes.position[0], axes.position[1], axes.position[2] );

		double const depth = camera()->projectedCoordinatesOf ( pos ) [2];
		Vec const q1 = camera()->unprojectedCoordinatesOf ( Vec ( 0.0, 0.0, depth ) );
		Vec const q2 = camera()->unprojectedCoordinatesOf ( Vec ( width(), height(), depth ) );
		axes.size = ( q1 - q2 ).norm() / 10.0;
		axes.shown = true;
		axes.selected_axis = -1;
	}
}

void MyViewer::keyPressEvent ( QKeyEvent * e )
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
			//increaseDegree();
			updateGL();
			break;
		case Qt::Key_N:
			show_normals = !show_normals;
			updateGL();
			break;
		case Qt::Key_B:
			show_curves = !show_curves;
			updateGL();
			break;
		case Qt::Key_G:
			generateMesh();
			mesh.request_face_normals(); mesh.request_vertex_normals();
			mesh.update_face_normals();  mesh.update_vertex_normals();
			updateMeanCurvature();
			focusOnMesh();

			updateGL();
			break;
		case Qt::Key_8:
			normalSize += 0.5;
			updateGL();
			break;
		case  Qt::Key_9:
			normalSize -= 0.5;
			updateGL();
			break;

		default:
			QGLViewer::keyPressEvent ( e );
		}
	else
	{ QGLViewer::keyPressEvent ( e ); }
}

void MyViewer::calculateNormals ( size_t sampling )
{
	normals.clear();
	normalSamples.clear();
	Transfinite::RMF rmf;

	for ( size_t i = 0; i < bsCurves.size(); ++i )
	{

		rmf.setCurve ( bsCurves[i] );
		std::vector<Geometry::Vector3D> derstart1;
		std::vector<Geometry::Vector3D> derstart2;
		std::vector<Geometry::Vector3D> derend1;
		std::vector<Geometry::Vector3D> derend2;
		if ( i == 0 )
		{
			bsCurves[0]->eval ( 0, 1, derstart1 );
			bsCurves[bsCurves.size() - 1]->eval ( 1, 1, derstart2 );

			bsCurves[0]->eval ( 1, 1, derend1 );
			bsCurves[1]->eval ( 0, 1, derend2 );

		}

		else if ( i == bsCurves.size() - 1 )
		{
			bsCurves[bsCurves.size() - 1]->eval ( 0, 1, derstart1 );
			bsCurves[bsCurves.size() - 2]->eval ( 1, 1, derstart2 );

			bsCurves[bsCurves.size() - 1]->eval ( 1, 1, derend1 );
			bsCurves[0]->eval ( 0, 1, derend2 );

		}
		else
		{
			bsCurves[i]->eval ( 0, 1, derstart1 );
			bsCurves[i - 1]->eval ( 1, 1, derstart2 );

			bsCurves[i]->eval ( 1, 1, derend1 );
			bsCurves[i + 1]->eval ( 0, 1, derend2 );

		}
		rmf.setStart ( ( derstart1[1] ^ derstart2[1] ).normalize() );
		rmf.setEnd ( - ( derend1[1] ^ derend2[1] ).normalize() );
		rmf.update();
		std::vector<Vec>  normalsOfThisCurve;
		for ( size_t t = 0; t < sampling; ++t )
		{
			Vector3D p = rmf.eval ( ( double ) t / sampling );
			normalSamples.push_back ( Vector ( p[0], p[1], p[2] ) );
			normalsOfThisCurve.push_back ( Vec ( p[0], p[1], p[2] ) );
		}
		normals[i] = normalsOfThisCurve ;
	}
}

void MyViewer::calculatePlain()
{
	pointsOnPlain.clear();
	Vector3D sum;
	float j = 0;
	for ( auto i : bsCurves )
	{

		for ( size_t t = 0; t < sampling; ++t )
		{
			sum = sum + i->eval ( ( double ) t / sampling );

			++j;
		}
	}
	sum = sum / j;
	plainPoint = Vec ( sum[0], sum[1], sum[2] );

	calculateNormals ( sampling );
	sum = Vector3D ( 0, 0, 0 );
	j = 0;
	for ( auto curve : normals )
		for ( auto normal : curve.second )
		{
			sum[0] = sum[0] + normal[0];
			sum[1] = sum[1] + normal[1];
			sum[2] = sum[2] + normal[2];
			++j;
		}
	sum = sum / j;
	plainNormal = Vec ( sum[0], sum[1], sum[2] );
	plainNormal.normalize();
	calculate2DPoints();
}

void MyViewer::calculate2DPoints (  )
{
	pointsOnPlain.clear();
	samples.clear();
	for ( auto i : bsCurves )
	{
		for ( size_t t = 0; t < sampling; ++t )
		{
			Vector3D tmp = i->eval ( ( double ) t / sampling );
			Vec q = Vec ( tmp[0], tmp[1], tmp[2] );
			samples.push_back ( Vector ( tmp[0], tmp[1], tmp[2] ) );
			float l = ( q - plainPoint ) * plainNormal ;
			Vec result = q - l * plainNormal;
			pointsOnPlain.push_back ( result );
		}
	}
}

Vec MyViewer::intersectLines ( Vec const & ap, Vec const & ad, Vec const & bp, Vec const & bd ) const
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
	calculatePlain();
	Vec w1 ( 1, 0, 0 ), w2 ( 0, 1, 0 ), u, v;
	u = plainNormal ^ w1;
	v = plainNormal ^ w2;
	if ( v.norm() > u.norm() )
	{ u = v; }
	u.normalize();
	v = u ^ plainNormal;
	v.normalize();

	std::vector<double> points;
	for ( auto i : pointsOnPlain )
	{
		points.push_back ( ( ( i - plainPoint ) *u ) );
		points.push_back ( ( ( i - plainPoint ) *v ) );

	}

	int n = static_cast<int> ( points.size() / 2 );
//
//// Input segments : just a closed polygon
	std::vector<int> segments ( n * 2 );
	for ( int i = 0; i < n; ++i )
	{
		segments[2 * i] = i;
		segments[2 * i + 1] = i + 1;
	}
	segments[2 * n - 1] = 0;

// Setup input data structure
	struct triangulateio in, out;
	in.pointlist = &points[0];
	in.numberofpoints = n;
	in.numberofpointattributes = 0;
	in.pointmarkerlist = NULL;
	in.segmentlist = &segments[0];
	in.numberofsegments = n;
	in.segmentmarkerlist = NULL;
	in.numberofholes = 0;
	in.numberofregions = 0;

// Setup output data structure
	out.pointlist = NULL;
	out.pointattributelist = NULL;
	out.pointmarkerlist = NULL;
	out.trianglelist = NULL;
	out.triangleattributelist = NULL;
	out.segmentlist = NULL;
	out.segmentmarkerlist = NULL;

// Call the library function [with maximum triangle area = 10]
	triangulate ( ( char * ) "pa10qzQY", &in, &out, ( struct triangulateio * ) NULL );

	mesh.clear();
	std::vector<MyMesh::VertexHandle> handles, tri;

	for ( int i = 0; i < out.numberofpoints; ++i )
	{
		Vec p = plainPoint + u * out.pointlist[2 * i] + v * out.pointlist[2 * i + 1];
		handles.push_back ( mesh.add_vertex ( MyMesh::Point ( p[0], p[1], p[2] ) ) );
	}

	for ( int i = 0; i < out.numberoftriangles; ++i )
	{
		tri.clear();
		for ( int j = 0; j < 3; ++j )
		{ tri.push_back ( handles[out.trianglelist[3 * i + j]] ); }
		mesh.add_face ( tri );
	}
	mesh.request_vertex_normals();

	int boundaryVertexIndex = 0;
	for ( MyMesh::VertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i )
	{
		if ( mesh.is_boundary ( *i ) )
		{
			mesh.set_point ( *i, samples[boundaryVertexIndex] );
			mesh.set_normal ( *i, normalSamples[boundaryVertexIndex] );
			++boundaryVertexIndex;
		}
		else
		{
			mesh.set_normal ( *i, Vector ( plainNormal[0], plainNormal[1], plainNormal[2] ) );
		}
	}
	calculateWeights();
	modifyMesh();
}

void MyViewer::calculateWeights()
{
	Eigen::MatrixXd M ( mesh.n_vertices(), mesh.n_vertices() );
	Eigen::MatrixXd Ls ( mesh.n_vertices(), mesh.n_vertices() );
	Eigen::MatrixXd L ( mesh.n_vertices(), mesh.n_vertices() );
	Eigen::MatrixXd Vcs ( mesh.n_vertices(), 3 );
	Eigen::MatrixXd Ncs ( mesh.n_vertices(), 3 );
	// halfEdge oppositeAngles
	for ( MyMesh::HalfedgeIter i = mesh.halfedges_begin(), ie = mesh.halfedges_end(); i != ie; ++i )
	{
		if ( mesh.is_boundary ( *i ) )
		{
			MyMesh::HalfedgeHandle h1 = mesh.opposite_halfedge_handle ( *i );
			MyMesh::HalfedgeHandle h2 = mesh.next_halfedge_handle ( h1 );
			h1 = mesh.next_halfedge_handle ( h2 );
			mesh.data ( *i ).oppositeAngle =
			    acos ( -halfedgeVector ( h1 ).normalize() | halfedgeVector ( h2 ).normalize() );
		}
		else
		{
			MyMesh::HalfedgeHandle h1 = mesh.next_halfedge_handle ( *i );
			MyMesh::HalfedgeHandle h2 = mesh.next_halfedge_handle ( h1 );
			h1 = mesh.opposite_halfedge_handle ( h1 );
			mesh.data ( *i ).oppositeAngle =
			    acos ( halfedgeVector ( h1 ).normalize() | halfedgeVector ( h2 ).normalize() );
			qDebug() << "angles";
			qDebug() << mesh.data ( *i ).oppositeAngle;
		}
	}
	//Edge weights
	for ( MyMesh::EdgeIter i = mesh.edges_begin(), ie = mesh.edges_end(); i != ie; ++i )
	{
		MyMesh::HalfedgeHandle h1 = mesh.halfedge_handle ( *i, 0 );
		MyMesh::HalfedgeHandle h2 = mesh.halfedge_handle ( *i, 1 );
		mesh.data ( *i ).weight = ( cot ( mesh.data ( h1 ).oppositeAngle ) + cot ( mesh.data ( h2 ).oppositeAngle ) ) / 2.0;
		qDebug() << "edge weigths";
		qDebug() << mesh.data ( *i ).weight;
	}

	//Voronoi Area and M matrix calculation
	int index = 0;
	for ( MyMesh::VertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i )
	{
		for ( MyMesh::ConstVertexEdgeIter j ( mesh, *i ); j.is_valid(); ++j )
		{
			mesh.data ( *i ).VoronoiArea += mesh.calc_edge_sqr_length ( *j ) * mesh.data ( *j ).weight / 4.0;
		}
		qDebug() << mesh.data ( *i ).VoronoiArea;
		M ( index, index ) = 1.0 / mesh.data ( *i ).VoronoiArea;
		++index;
	}

	// calculate M and Ls matrices
	for ( int i = 0; i < mesh.n_vertices(); i++ )
	{
		Vcs ( i, 0 ) = mesh.point ( mesh.vertex_handle ( i ) ) [0];
		Vcs ( i, 1 ) = mesh.point ( mesh.vertex_handle ( i ) ) [1];
		Vcs ( i, 2 ) = mesh.point ( mesh.vertex_handle ( i ) ) [2];

		Ncs ( i, 0 ) = mesh.normal ( mesh.vertex_handle ( i ) ) [0];
		Ncs ( i, 1 ) = mesh.normal ( mesh.vertex_handle ( i ) ) [1];
		Ncs ( i, 2 ) = mesh.normal ( mesh.vertex_handle ( i ) ) [2];
		for ( int j = 0; j < mesh.n_vertices(); j++ )
		{
			Ls ( i, j ) = 0;
			if ( i == j )
			{
				MyMesh::VertexHandle vi = mesh.vertex_handle ( i );
				double sum = 0;
				for ( MyMesh::ConstVertexEdgeIter k ( mesh, vi ); k; ++k )
				{
					sum += mesh.data ( k ).weight;
				}
				Ls ( i, j ) = -sum;
			}
			else
			{
				MyMesh::VertexHandle vi = mesh.vertex_handle ( i );
				MyMesh::VertexHandle vj = mesh.vertex_handle ( j );
				MyMesh::HalfedgeHandle he = mesh.find_halfedge ( vi, vj );
				if ( he.is_valid() )
				{
					Ls ( i, j ) = mesh.data ( mesh.edge_handle ( he ) ).weight;
				}

			}
		}
	}

	//calculate L matrice
	L = M.inverse() * Ls;
	L = L * L;

	Eigen::MatrixXd bNorm = Eigen::MatrixXd::Zero ( mesh.n_vertices() - samples.size(), 3 );

	Eigen::MatrixXd bVert ( mesh.n_vertices() - samples.size(), 3 );
	bVert = Eigen::MatrixXd::Zero ( bVert.rows(), bVert.cols() );

	Eigen::MatrixXd A = Eigen::MatrixXd::Zero ( mesh.n_vertices() - samples.size(), mesh.n_vertices() - samples.size() );

	int BrowIndex = 0;
	int BcolIndex;
	int ArowIndex = 0;
	int AcolIndex;
	for ( int i = 0; i < mesh.n_vertices(); i++ )
	{

		AcolIndex = 0;
		if ( !mesh.is_boundary ( mesh.vertex_handle ( i ) ) )
		{
			for ( int k = 0; k < mesh.n_vertices(); k++ )
			{
				BcolIndex = 0;
				for ( int j = 0; j < 3; j++ )
				{
					if ( mesh.is_boundary ( mesh.vertex_handle ( k ) ) )
					{
						bVert ( BrowIndex, BcolIndex ) -= L ( i, k ) * Vcs ( k, j );

						bNorm ( BrowIndex, BcolIndex ) -= L ( i, k ) * Ncs ( k, j );
					}
					++BcolIndex;

				}

				if ( !mesh.is_boundary ( mesh.vertex_handle ( k ) ) )
				{
					A ( ArowIndex, AcolIndex ) = L ( i, k );

					++AcolIndex;
				}

			}
			++BrowIndex;
			++ArowIndex;
		}
	}

	newCoords = A.colPivHouseholderQr().solve ( bVert );
	newNormals = A.colPivHouseholderQr().solve ( bNorm );
}

void MyViewer::modifyMesh()
{
	int boundaryVertexIndex = 0;
	int innerVertexIndex = 0;

	for ( MyMesh::VertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i )
	{
		if ( mesh.is_boundary ( *i ) )
		{
			mesh.set_point ( *i, samples[boundaryVertexIndex] );
			mesh.set_normal ( *i, normalSamples[boundaryVertexIndex] );
			++boundaryVertexIndex;
		}
		else
		{
			mesh.set_point ( *i, Vector ( newCoords ( innerVertexIndex, 0 ),
			                              newCoords ( innerVertexIndex, 1 ),
			                              newCoords ( innerVertexIndex, 2 ) ) );
			mesh.set_normal ( *i, Vector ( newNormals ( innerVertexIndex, 0 ),
			                               newNormals ( innerVertexIndex, 1 ),
			                               newNormals ( innerVertexIndex, 2 ) ) );
			++innerVertexIndex;
		}
	}
}

void MyViewer::mouseMoveEvent ( QMouseEvent * e )
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

		bsCurves[selected / 1000]->getControlPoints() [selected % 1000][0] = Vec ( axes.position[0], axes.position[1], axes.position[2] ) [0];
		bsCurves[selected / 1000]->getControlPoints() [selected % 1000][1] = Vec ( axes.position[0], axes.position[1], axes.position[2] ) [1];
		bsCurves[selected / 1000]->getControlPoints() [selected % 1000][2] = Vec ( axes.position[0], axes.position[1], axes.position[2] ) [2];


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

void MyViewer::focusOnMesh()
{
	MyMesh::Point box_min, box_max;
	box_min = box_max = mesh.point ( *mesh.vertices_begin() );
	for ( MyMesh::ConstVertexIter i = mesh.vertices_begin(), ie = mesh.vertices_end(); i != ie; ++i )
	{
		box_min.minimize ( mesh.point ( *i ) );
		box_max.maximize ( mesh.point ( *i ) );
	}
	camera()->setSceneBoundingBox ( Vec ( box_min[0], box_min[1], box_min[2] ),
	                                Vec ( box_max[0], box_max[1], box_max[2] ) );
	camera()->showEntireScene();
}
