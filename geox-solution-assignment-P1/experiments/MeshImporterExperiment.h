//---------------------------------------------------------------------------
#ifndef MeshImporterExperimentH
#define MeshImporterExperimentH
//---------------------------------------------------------------------------
#include "Experiment.h"
#include "LinearAlgebra.h"
#include "BasicGLViewer3D.h"
#include "BasicGLViewer3D.h"
#include "ObjectListProperty.h"
#include "TriangleMesh.h"
//---------------------------------------------------------------------------


class TriangleMesh;
class SimpleGLMeshMaterial;
class Raytracer;
class Raytracer;

///
/// This is an example experiment that demonstrates importing and rendering triangle meshes.
/// It supports raytracing as well now.
///
class MeshImporterExperiment : public Experiment, public RenderableObject {
	GEOX_CLASS(MeshImporterExperiment)
 private:
	BasicGLViewer3D* viewer;
	vector<TriangleMesh*> meshes;
	SimpleGLMeshMaterial *renderer;
	Raytracer *raytracer;
	card32 rtWidth;
	card32 rtHeight;
	float32 groundPlane;
	Camera oldCam;
	ViewFrustum oldViewFrustum;
	bool freeze;
	bool renderMesh;
	bool renderKDTree;
	bool renderRay;
	card32 numLeafTriangles;
	//IMPLEMENT_OBJECT_LIST_ACCESS_STL(meshes, TriangleMesh)

	void clearScene();
	void rebuildKDTree();

 public:                                       

	MeshImporterExperiment();

	virtual void renderGL();
	void importMesh();
	void traceImage();
	virtual QWidget *createViewer();

	~MeshImporterExperiment();
};


#endif                                         
