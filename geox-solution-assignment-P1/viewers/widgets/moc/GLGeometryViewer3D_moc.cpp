/****************************************************************************
** Meta object code from reading C++ file 'GLGeometryViewer3D.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.4.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../GLGeometryViewer3D.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'GLGeometryViewer3D.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.4.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_GLGeometryViewer3D_t {
    QByteArrayData data[13];
    char stringdata[268];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_GLGeometryViewer3D_t, stringdata) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_GLGeometryViewer3D_t qt_meta_stringdata_GLGeometryViewer3D = {
    {
QT_MOC_LITERAL(0, 0, 18), // "GLGeometryViewer3D"
QT_MOC_LITERAL(1, 19, 20), // "viewerContentChanged"
QT_MOC_LITERAL(2, 40, 0), // ""
QT_MOC_LITERAL(3, 41, 20), // "on_btnCamera_pressed"
QT_MOC_LITERAL(4, 62, 18), // "on_btnPick_pressed"
QT_MOC_LITERAL(5, 81, 19), // "on_btnClear_pressed"
QT_MOC_LITERAL(6, 101, 19), // "on_btnLight_pressed"
QT_MOC_LITERAL(7, 121, 22), // "on_btnFillMode_pressed"
QT_MOC_LITERAL(8, 144, 22), // "on_btnDrawAxis_pressed"
QT_MOC_LITERAL(9, 167, 21), // "on_btnHandles_pressed"
QT_MOC_LITERAL(10, 189, 25), // "on_btnResetCamera_pressed"
QT_MOC_LITERAL(11, 215, 28), // "on_btnVertexLighting_pressed"
QT_MOC_LITERAL(12, 244, 23) // "sltWidgetContentChanged"

    },
    "GLGeometryViewer3D\0viewerContentChanged\0"
    "\0on_btnCamera_pressed\0on_btnPick_pressed\0"
    "on_btnClear_pressed\0on_btnLight_pressed\0"
    "on_btnFillMode_pressed\0on_btnDrawAxis_pressed\0"
    "on_btnHandles_pressed\0on_btnResetCamera_pressed\0"
    "on_btnVertexLighting_pressed\0"
    "sltWidgetContentChanged"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_GLGeometryViewer3D[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      11,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    0,   69,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       3,    0,   70,    2, 0x08 /* Private */,
       4,    0,   71,    2, 0x08 /* Private */,
       5,    0,   72,    2, 0x08 /* Private */,
       6,    0,   73,    2, 0x08 /* Private */,
       7,    0,   74,    2, 0x08 /* Private */,
       8,    0,   75,    2, 0x08 /* Private */,
       9,    0,   76,    2, 0x08 /* Private */,
      10,    0,   77,    2, 0x08 /* Private */,
      11,    0,   78,    2, 0x08 /* Private */,
      12,    0,   79,    2, 0x08 /* Private */,

 // signals: parameters
    QMetaType::Void,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void GLGeometryViewer3D::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        GLGeometryViewer3D *_t = static_cast<GLGeometryViewer3D *>(_o);
        switch (_id) {
        case 0: _t->viewerContentChanged(); break;
        case 1: _t->on_btnCamera_pressed(); break;
        case 2: _t->on_btnPick_pressed(); break;
        case 3: _t->on_btnClear_pressed(); break;
        case 4: _t->on_btnLight_pressed(); break;
        case 5: _t->on_btnFillMode_pressed(); break;
        case 6: _t->on_btnDrawAxis_pressed(); break;
        case 7: _t->on_btnHandles_pressed(); break;
        case 8: _t->on_btnResetCamera_pressed(); break;
        case 9: _t->on_btnVertexLighting_pressed(); break;
        case 10: _t->sltWidgetContentChanged(); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (GLGeometryViewer3D::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&GLGeometryViewer3D::viewerContentChanged)) {
                *result = 0;
            }
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject GLGeometryViewer3D::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_GLGeometryViewer3D.data,
      qt_meta_data_GLGeometryViewer3D,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *GLGeometryViewer3D::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *GLGeometryViewer3D::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_GLGeometryViewer3D.stringdata))
        return static_cast<void*>(const_cast< GLGeometryViewer3D*>(this));
    return QWidget::qt_metacast(_clname);
}

int GLGeometryViewer3D::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 11)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 11;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 11)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 11;
    }
    return _id;
}

// SIGNAL 0
void GLGeometryViewer3D::viewerContentChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 0, Q_NULLPTR);
}
QT_END_MOC_NAMESPACE
