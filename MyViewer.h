#pragma once

#include <string>

#include <QGLViewer/qglviewer.h>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include "rmf.hh"

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

	//void increaseDegree();


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
			double u;
			double v;
		};
		FaceTraits
		{
			double area;
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

	size_t degree[2];
	std::vector<Vec> control_points;

	std::vector<std::shared_ptr<Geometry::BSCurve>> bsCurves;
	MyMesh mesh;
	double mean_min, mean_max;
	double cutoff_ratio;
	bool show_control_points, show_solid, show_wireframe, show_curves, show_normals;

	enum { COLOR_PLAIN, COLOR_MEAN, COLOR_ISOPHOTES } coloring;
	GLuint isophote_texture;
	size_t selected;
	struct ModificationAxes
	{
		bool shown;
		float size;
		int selected_axis;
		GLfloat position[3];
		Vec grabbed_pos, original_pos;
	} axes;

	/**
	 * Returns the next line from file as a stringstream
	 */
	std::istringstream nextLine(std::ifstream &file);

	/** 
	 * Reads 1 curve from file (file input stream)
	 *
	 * Format must be as follows:
	 *
	 * 3                   # degree
	 * 9                   # length of knots (vector)
	 * 0 0 0 0 0.5 1 1 1 1 # data in knots
	 * 5                   # number of control points (cpts)
	 * 0 0 1               # coordinates of 1st control point (cpts[0])
	 * 0 0 2               # coordinates of 2nd control point (cpts[1])
	 * 0 2 1               # ...
	 * 2 0 1               # 
	 * 2 2 1               # coordinates of last control point
	 */
	void readBSCurve(std::ifstream &file);
};

#include "MyViewer.hpp"
