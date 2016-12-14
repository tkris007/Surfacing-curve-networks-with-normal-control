#pragma once



#include <string>

#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <vector>

#define ANSI_DECLARATORS
#define REAL double
#define VOID void

extern "C" {
#include "triangle.h"
}

#undef REAL
#undef VOID
#undef ANSI_DECLARATORS

#include <QGLViewer/qglviewer.h>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include "rmf.hh"
#include <Eigen\Eigen>

using qglviewer::Vec;

class MyViewer : public QGLViewer {
	Q_OBJECT

public:
	MyViewer ( QWidget *parent );
	virtual ~MyViewer();

	inline double getCutoffRatio() const;
	inline void setCutoffRatio ( double ratio );
	inline double getMeanMin() const;
	inline void setMeanMin ( double min );
	inline double getMeanMax() const;
	inline void setMeanMax ( double max );

	bool openBSpline ( std::string const &filename );
	bool saveBSpline ( std::string const& filename );

signals:
	void startComputation ( QString message );
	void midComputation ( int percent );
	void endComputation();

protected:
	virtual void init();
	virtual void draw();
	virtual void drawWithNames();
	virtual void postSelection ( const QPoint &p );
	virtual void keyPressEvent ( QKeyEvent *e );
	virtual void mouseMoveEvent ( QMouseEvent *e );
	virtual QString helpString() const;

private:
	struct MyTraits : public OpenMesh::DefaultTraits
	{
		VertexTraits
		{
			double area;              // total area of the surrounding triangles
			double mean;              // approximated mean curvature
			double VoronoiArea;
		};
		FaceTraits
		{
			double area;
		};
		EdgeTraits
		{
			double weight;
		};
		HalfedgeTraits
		{
			double oppositeAngle;
		};
	};
	typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;
	typedef MyMesh::Point Vector;

	inline Vector halfedgeVector ( MyMesh::HalfedgeHandle const &h ) const;
	void updateMeanMinMax();
	void updateMeanCurvature ( bool update_min_max = true );
	void meanMapColor ( double d, double *color ) const;

	void drawAxes() const;
	void drawAxesWithNames() const;
	void drawCurves() const;
	void drawControlPoints() const;
	void drawNormals() const;
	Vec intersectLines ( Vec const &ap, Vec const &ad, Vec const &bp, Vec const &bd ) const;
	static void bernsteinAll ( size_t n, double u, std::vector<double> &coeff );
	void generateMesh();

	void calculateNormals ( size_t _step );
	void calculatePlain();
	void calculate2DPoints (  );
	void calculateWeights();

	std::vector<Vec> control_points;

	float step;
	size_t sampling;
	Vec plainPoint;
	Vec plainNormal;
	std::vector<std::shared_ptr<Geometry::BSCurve>> bsCurves;
	std::map<size_t, std::vector<Vec> > normals;
	std::vector<Vec> pointsOnPlain;
	MyMesh mesh;
	double mean_min, mean_max;
	double cutoff_ratio;
	bool show_control_points, show_solid, show_wireframe, show_curves, show_normals;

	enum { COLOR_PLAIN, COLOR_MEAN, COLOR_ISOPHOTES } coloring;
	GLuint isophote_texture;
	size_t selected;
	float normalSize;
	struct ModificationAxes
	{
		bool shown;
		float size;
		int selected_axis;
		GLfloat position[3];
		Vec grabbed_pos, original_pos;
	} axes;

	/**
	 * Returns the next line from file as a string
	 */
	std::string nextLine ( std::ifstream &file );

	/**
	 * Reads 1 curve from file (file input stream)
	 */
	void readBSCurve ( std::ifstream &file );
};

#include "MyViewer.hpp"
