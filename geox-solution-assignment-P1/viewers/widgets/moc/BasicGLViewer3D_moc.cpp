/****************************************************************************
** Meta object code from reading C++ file 'BasicGLViewer3D.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.4.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../BasicGLViewer3D.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'BasicGLViewer3D.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.4.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_BasicGLViewer3D_t {
    QByteArrayData data[10];
    char stringdata[192];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_BasicGLViewer3D_t, stringdata) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_BasicGLViewer3D_t qt_meta_stringdata_BasicGLViewer3D = {
    {
QT_MOC_LITERAL(0, 0, 15), // "BasicGLViewer3D"
QT_MOC_LITERAL(1, 16, 20), // "viewerContentChanged"
QT_MOC_LITERAL(2, 37, 0), // ""
QT_MOC_LITERAL(3, 38, 25), // "on_btnResetCamera_pressed"
QT_MOC_LITERAL(4, 64, 19), // "on_btnOrtho_clicked"
QT_MOC_LITERAL(5, 84, 22), // "on_btnDrawAxes_clicked"
QT_MOC_LITERAL(6, 107, 19), // "on_btnYAxis_clicked"
QT_MOC_LITERAL(7, 127, 22), // "on_btnScaleImg_clicked"
QT_MOC_LITERAL(8, 150, 23), // "sltWidgetContentChanged"
QT_MOC_LITERAL(9, 174, 17) // "sltCustomRenderGL"

    },
    "BasicGLViewer3D\0viewerContentChanged\0"
    "\0on_btnResetCamera_pressed\0"
    "on_btnOrtho_clicked\0on_btnDrawAxes_clicked\0"
    "on_btnYAxis_clicked\0on_btnScaleImg_clicked\0"
    "sltWidgetContentChanged\0sltCustomRenderGL"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_BasicGLViewer3D[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       8,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    0,   54,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       3,    0,   55,    2, 0x08 /* Private */,
       4,    1,   56,    2, 0x08 /* Private */,
       5,    1,   59,    2, 0x08 /* Private */,
       6,    1,   62,    2, 0x08 /* Private */,
       7,    1,   65,    2, 0x08 /* Private */,
       8,    0,   68,    2, 0x08 /* Private */,
       9,    0,   69,    2, 0x08 /* Private */,

 // signals: parameters
    QMetaType::Void,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void BasicGLViewer3D::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        BasicGLViewer3D *_t = static_cast<BasicGLViewer3D *>(_o);
        switch (_id) {
        case 0: _t->viewerContentChanged(); break;
        case 1: _t->on_btnResetCamera_pressed(); break;
        case 2: _t->on_btnOrtho_clicked((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 3: _t->on_btnDrawAxes_clicked((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 4: _t->on_btnYAxis_clicked((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 5: _t->on_btnScaleImg_clicked((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 6: _t->sltWidgetContentChanged(); break;
        case 7: _t->sltCustomRenderGL(); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (BasicGLViewer3D::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&BasicGLViewer3D::viewerContentChanged)) {
                *result = 0;
            }
        }
    }
}

const QMetaObject BasicGLViewer3D::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_BasicGLViewer3D.data,
      qt_meta_data_BasicGLViewer3D,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *BasicGLViewer3D::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *BasicGLViewer3D::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_BasicGLViewer3D.stringdata))
        return static_cast<void*>(const_cast< BasicGLViewer3D*>(this));
    if (!strcmp(_clname, "ViewerImgUpdateCallBack"))
        return static_cast< ViewerImgUpdateCallBack*>(const_cast< BasicGLViewer3D*>(this));
    return QWidget::qt_metacast(_clname);
}

int BasicGLViewer3D::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 8)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 8;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 8)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 8;
    }
    return _id;
}

// SIGNAL 0
void BasicGLViewer3D::viewerContentChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 0, Q_NULLPTR);
}
QT_END_MOC_NAMESPACE
