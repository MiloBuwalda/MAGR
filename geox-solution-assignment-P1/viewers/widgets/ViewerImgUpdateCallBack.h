//---------------------------------------------------------------------------
#ifndef ViewerImgUpdateCallBackH
#define ViewerImgUpdateCallBackH
//---------------------------------------------------------------------------
#include <QtGui\qimage.h>
//---------------------------------------------------------------------------

class ViewerImgUpdateCallBack {
 public:
	virtual void viewImage(QImage &currentImg) = 0;
};


#endif