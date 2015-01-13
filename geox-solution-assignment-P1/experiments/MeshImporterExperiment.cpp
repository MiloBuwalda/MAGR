//---------------------------------------------------------------------------
#include "stdafx.h"
//---------------------------------------------------------------------------
#include "MeshImporterExperiment.h"
//---------------------------------------------------------------------------
#include "Properties.h"
#include "GLGeometryViewer3D.h"
#include "GeoXOutput.h"   
#include "MeshImporterSMFOBJ.h"
#include "SimpleGLMeshMaterial.h"
#include "TriangleMesh.h"
#include "ObjectViewsTable.h"
#include "ObjectListProperty.h"
#include "raytracing/Raytracer.h"
#include "raytracing/RaytracingMaterial.h"
//---------------------------------------------------------------------------

IMPLEMENT_GEOX_CLASS( MeshImporterExperiment ,0)
{
   BEGIN_CLASS_INIT( MeshImporterExperiment );
//	ADD_OBJECT_LIST_PROP(meshes, 0, TriangleMesh::getClass());
	ADD_OBJECT_PROP(renderer, 0, SimpleGLMeshMaterial::getClass(), true)
	ADD_OBJECT_PROP(raytracer, 0, Raytracer::getClass(), true)
	ADD_CARD32_PROP(rtWidth, 0)
	ADD_CARD32_PROP(rtHeight, 0)
	ADD_CARD32_PROP(numLeafTriangles, 0)
	ADD_FLOAT32_PROP(groundPlane, 0)
	ADD_BOOLEAN_PROP(freeze, 0)
	ADD_BOOLEAN_PROP(renderMesh, 0)
	ADD_BOOLEAN_PROP(renderKDTree, 0)
	ADD_BOOLEAN_PROP(renderRay, 0)
	ADD_NOARGS_METHOD(MeshImporterExperiment::importMesh)
   ADD_NOARGS_METHOD(MeshImporterExperiment::traceImage)
	ADD_NOARGS_METHOD(MeshImporterExperiment::rebuildKDTree)	
	ADD_NOARGS_METHOD(MeshImporterExperiment::clearScene)	
}

QWidget *MeshImporterExperiment::createViewer() {
   viewer = new BasicGLViewer3D();
	viewer->addRenderCallback(this);
   return viewer;
}

MeshImporterExperiment::MeshImporterExperiment() 
{
	viewer = NULL;
	renderer = new SimpleGLMeshMaterial();
	raytracer = new Raytracer();
	rtWidth = 320;
	rtHeight = 240;
	groundPlane = -1;
	freeze = false;
	renderMesh = true;
	renderKDTree = true;
	renderRay = true;
	numLeafTriangles = 1;
}

MeshImporterExperiment::~MeshImporterExperiment()
{
	for (size_t i=0; i<meshes.size(); i++) delete meshes[i];
	delete renderer;
}



void MeshImporterExperiment::importMesh() 
{
	
	QString filename = QFileDialog::getOpenFileName(viewer, "Open OBJ File...", QString(), "3D Objects (*.obj *.smf)", 0, 0 );
	if (filename != "") {
		MeshImporterSMFOBJ importer;
		importer.performImport(qString2STLString(filename));
		TriangleMesh *mesh = importer.createTriangleMesh();
		if (mesh) {
			meshes.push_back(mesh);
			raytracer->setNumLeafTriangles(numLeafTriangles);
			raytracer->addObject((TriangleMesh*)mesh->copy(), new RaytracingMaterial());
		}
	}

	// display changes
	ObjectViewsTable::update(this);
	viewer->refresh();
}

void MeshImporterExperiment::renderGL()
{
	if (renderer) {
		if (renderMesh) for (size_t i=0; i<meshes.size(); i++) renderer->draw(meshes[i]);
		if (renderKDTree) raytracer->renderKDTreeGL();
		if (!freeze) {
			if (renderRay) raytracer->renderCenterRayGL(rtWidth, rtHeight, viewer->getCamera(), viewer->getViewFrustum());
			oldCam = viewer->getCamera();
			oldViewFrustum = viewer->getViewFrustum();
		} else {
			if (renderRay) raytracer->renderCenterRayGL(rtWidth, rtHeight, oldCam, oldViewFrustum);
		}
	}
}

void MeshImporterExperiment::traceImage()
{
	QImage img;
	raytracer->render(img, rtWidth, rtHeight, viewer->getCamera(), viewer->getViewFrustum(), groundPlane, viewer);
	viewer->viewImage(img);
}

void MeshImporterExperiment::rebuildKDTree()
{
	raytracer->setNumLeafTriangles(numLeafTriangles);
	raytracer->rebuildKDTree();
}

void MeshImporterExperiment::clearScene()
{
	raytracer->clearObjects();
}

