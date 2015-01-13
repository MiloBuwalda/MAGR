/********************************************************************************
** Form generated from reading UI file 'BasicGLViewer3D.ui'
**
** Created: Tue 4. Mar 17:21:40 2014
**      by: Qt User Interface Compiler version 4.8.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_BASICGLVIEWER3D_H
#define UI_BASICGLVIEWER3D_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QFrame>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QSpacerItem>
#include <QtGui/QTabWidget>
#include <QtGui/QToolButton>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "basicgeoxglwidget3d.h"

QT_BEGIN_NAMESPACE

class Ui_GLGeometryViewer3DImproved
{
public:
    QHBoxLayout *hboxLayout;
    QFrame *frame_2;
    QVBoxLayout *vboxLayout;
    QToolButton *btnResetCamera;
    QToolButton *btnDrawAxes;
    QToolButton *btnOrtho;
    QToolButton *btnYAxis;
    QToolButton *btnScaleImg;
    QSpacerItem *spacerItem;
    QTabWidget *tabWidget;
    QWidget *tab_3;
    QHBoxLayout *horizontalLayout;
    BasicGeoXGLWidget3D *glFrame;
    QWidget *tab_4;
    QHBoxLayout *horizontalLayout_2;
    QLabel *lbImage;

    void setupUi(QWidget *GLGeometryViewer3DImproved)
    {
        if (GLGeometryViewer3DImproved->objectName().isEmpty())
            GLGeometryViewer3DImproved->setObjectName(QString::fromUtf8("GLGeometryViewer3DImproved"));
        GLGeometryViewer3DImproved->resize(785, 531);
        GLGeometryViewer3DImproved->setAutoFillBackground(true);
        hboxLayout = new QHBoxLayout(GLGeometryViewer3DImproved);
        hboxLayout->setSpacing(7);
        hboxLayout->setContentsMargins(0, 0, 0, 0);
        hboxLayout->setObjectName(QString::fromUtf8("hboxLayout"));
        frame_2 = new QFrame(GLGeometryViewer3DImproved);
        frame_2->setObjectName(QString::fromUtf8("frame_2"));
        frame_2->setMinimumSize(QSize(54, 16));
        frame_2->setMaximumSize(QSize(54, 16777215));
        frame_2->setFrameShape(QFrame::StyledPanel);
        frame_2->setFrameShadow(QFrame::Raised);
        vboxLayout = new QVBoxLayout(frame_2);
#ifndef Q_OS_MAC
        vboxLayout->setSpacing(6);
#endif
        vboxLayout->setContentsMargins(3, 3, 3, 3);
        vboxLayout->setObjectName(QString::fromUtf8("vboxLayout"));
        btnResetCamera = new QToolButton(frame_2);
        btnResetCamera->setObjectName(QString::fromUtf8("btnResetCamera"));
        btnResetCamera->setMinimumSize(QSize(46, 58));
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/images/Camera.png"), QSize(), QIcon::Normal, QIcon::Off);
        btnResetCamera->setIcon(icon);
        btnResetCamera->setIconSize(QSize(32, 32));
        btnResetCamera->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);

        vboxLayout->addWidget(btnResetCamera);

        btnDrawAxes = new QToolButton(frame_2);
        btnDrawAxes->setObjectName(QString::fromUtf8("btnDrawAxes"));
        btnDrawAxes->setMinimumSize(QSize(46, 58));
        btnDrawAxes->setCheckable(true);
        btnDrawAxes->setChecked(true);
        btnDrawAxes->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);

        vboxLayout->addWidget(btnDrawAxes);

        btnOrtho = new QToolButton(frame_2);
        btnOrtho->setObjectName(QString::fromUtf8("btnOrtho"));
        btnOrtho->setMinimumSize(QSize(46, 58));
        QIcon icon1;
        icon1.addFile(QString::fromUtf8(":/images/wirecube.png"), QSize(), QIcon::Normal, QIcon::Off);
        btnOrtho->setIcon(icon1);
        btnOrtho->setCheckable(true);
        btnOrtho->setChecked(false);
        btnOrtho->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);

        vboxLayout->addWidget(btnOrtho);

        btnYAxis = new QToolButton(frame_2);
        btnYAxis->setObjectName(QString::fromUtf8("btnYAxis"));
        btnYAxis->setMinimumSize(QSize(46, 58));
        btnYAxis->setCheckable(true);
        btnYAxis->setChecked(false);
        btnYAxis->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);

        vboxLayout->addWidget(btnYAxis);

        btnScaleImg = new QToolButton(frame_2);
        btnScaleImg->setObjectName(QString::fromUtf8("btnScaleImg"));
        btnScaleImg->setMinimumSize(QSize(46, 58));
        btnScaleImg->setCheckable(true);
        btnScaleImg->setChecked(true);
        btnScaleImg->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);

        vboxLayout->addWidget(btnScaleImg);

        spacerItem = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        vboxLayout->addItem(spacerItem);


        hboxLayout->addWidget(frame_2);

        tabWidget = new QTabWidget(GLGeometryViewer3DImproved);
        tabWidget->setObjectName(QString::fromUtf8("tabWidget"));
        tab_3 = new QWidget();
        tab_3->setObjectName(QString::fromUtf8("tab_3"));
        horizontalLayout = new QHBoxLayout(tab_3);
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        glFrame = new BasicGeoXGLWidget3D(tab_3);
        glFrame->setObjectName(QString::fromUtf8("glFrame"));
        QSizePolicy sizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(glFrame->sizePolicy().hasHeightForWidth());
        glFrame->setSizePolicy(sizePolicy);
        glFrame->setMinimumSize(QSize(0, 0));
        glFrame->setMouseTracking(false);

        horizontalLayout->addWidget(glFrame);

        tabWidget->addTab(tab_3, QString());
        tab_4 = new QWidget();
        tab_4->setObjectName(QString::fromUtf8("tab_4"));
        horizontalLayout_2 = new QHBoxLayout(tab_4);
        horizontalLayout_2->setSpacing(0);
        horizontalLayout_2->setContentsMargins(0, 0, 0, 0);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        lbImage = new QLabel(tab_4);
        lbImage->setObjectName(QString::fromUtf8("lbImage"));
        lbImage->setScaledContents(true);
        lbImage->setAlignment(Qt::AlignCenter);

        horizontalLayout_2->addWidget(lbImage);

        tabWidget->addTab(tab_4, QString());

        hboxLayout->addWidget(tabWidget);


        retranslateUi(GLGeometryViewer3DImproved);

        tabWidget->setCurrentIndex(1);


        QMetaObject::connectSlotsByName(GLGeometryViewer3DImproved);
    } // setupUi

    void retranslateUi(QWidget *GLGeometryViewer3DImproved)
    {
        GLGeometryViewer3DImproved->setWindowTitle(QApplication::translate("GLGeometryViewer3DImproved", "Form", 0, QApplication::UnicodeUTF8));
        btnResetCamera->setText(QApplication::translate("GLGeometryViewer3DImproved", "Reset", 0, QApplication::UnicodeUTF8));
        btnDrawAxes->setText(QApplication::translate("GLGeometryViewer3DImproved", "Show\n"
"Axes", 0, QApplication::UnicodeUTF8));
        btnOrtho->setText(QApplication::translate("GLGeometryViewer3DImproved", "Orhto", 0, QApplication::UnicodeUTF8));
        btnYAxis->setText(QApplication::translate("GLGeometryViewer3DImproved", "Rot.\n"
"Y-Axis", 0, QApplication::UnicodeUTF8));
        btnScaleImg->setText(QApplication::translate("GLGeometryViewer3DImproved", "Scale\n"
"Image", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tab_3), QApplication::translate("GLGeometryViewer3DImproved", "OpenGL View", 0, QApplication::UnicodeUTF8));
        lbImage->setText(QApplication::translate("GLGeometryViewer3DImproved", "no image", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tab_4), QApplication::translate("GLGeometryViewer3DImproved", "Image View", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class GLGeometryViewer3DImproved: public Ui_GLGeometryViewer3DImproved {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_BASICGLVIEWER3D_H
