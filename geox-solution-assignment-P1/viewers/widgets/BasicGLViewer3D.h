//---------------------------------------------------------------------------
#ifndef GLGeometryViewer3DImprovedH
#define GLGeometryViewer3DImprovedH
//---------------------------------------------------------------------------
#include <QtWidgets\QWidget.h>
#include <QtWidgets\QLabel.h>
#include "moc\ui_BasicGLViewer3D.h"
#include "VisTypes.h"
#include "ViewerImgUpdateCallBack.h"
//---------------------------------------------------------------------------


/// This interface describes objects that output OpenGL commands
class RenderableObject { // interface
 public:
	 virtual void renderGL() = 0;
};


class BasicGLViewer3D : public QWidget, public ViewerImgUpdateCallBack {
   Q_OBJECT
 private:
   Ui_GLGeometryViewer3DImproved ui;
	vector<RenderableObject*> renderables;

	void drawAxes();


 private slots:
	void on_btnResetCamera_pressed();
	void on_btnOrtho_clicked(bool);
	void on_btnDrawAxes_clicked(bool);
	void on_btnYAxis_clicked(bool);
	void on_btnScaleImg_clicked(bool);
	
	void sltWidgetContentChanged();
	void sltCustomRenderGL();

signals:
	void viewerContentChanged();

 public:
   BasicGLViewer3D(QWidget * parent=NULL, Qt::WindowFlags f=0);

	// connect to a viewer
	void addRenderCallback(RenderableObject *obj) {renderables.push_back(obj);}

	// user interface appearance
	void      setBackgroundColor(Vector3f color){ui.glFrame->setBackgroundColor(color);};
	void      setBackgroundColor(float red, float green, float blue){ui.glFrame->setBackgroundColor(red, green,blue);};
	void      refresh() {ui.glFrame->updateGL();}

	/// displays an QImage immediately (processes event loop)
	void viewImage(QImage &img);

	const Camera &getCamera() {return ui.glFrame->getCamera();}
	const ViewFrustum &getViewFrustum() {return ui.glFrame->getViewFrustum();}

	~BasicGLViewer3D();
};





#endif
