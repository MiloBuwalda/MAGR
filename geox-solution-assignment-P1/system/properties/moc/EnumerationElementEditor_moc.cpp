/****************************************************************************
** Meta object code from reading C++ file 'EnumerationElementEditor.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.4.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../EnumerationElementEditor.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'EnumerationElementEditor.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.4.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_EnumerationElementEditor_t {
    QByteArrayData data[4];
    char stringdata[44];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_EnumerationElementEditor_t, stringdata) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_EnumerationElementEditor_t qt_meta_stringdata_EnumerationElementEditor = {
    {
QT_MOC_LITERAL(0, 0, 24), // "EnumerationElementEditor"
QT_MOC_LITERAL(1, 25, 15), // "sltIndexChanged"
QT_MOC_LITERAL(2, 41, 0), // ""
QT_MOC_LITERAL(3, 42, 1) // "i"

    },
    "EnumerationElementEditor\0sltIndexChanged\0"
    "\0i"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_EnumerationElementEditor[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    1,   19,    2, 0x0a /* Public */,

 // slots: parameters
    QMetaType::Void, QMetaType::Int,    3,

       0        // eod
};

void EnumerationElementEditor::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        EnumerationElementEditor *_t = static_cast<EnumerationElementEditor *>(_o);
        switch (_id) {
        case 0: _t->sltIndexChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObject EnumerationElementEditor::staticMetaObject = {
    { &ElementEditor::staticMetaObject, qt_meta_stringdata_EnumerationElementEditor.data,
      qt_meta_data_EnumerationElementEditor,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *EnumerationElementEditor::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *EnumerationElementEditor::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_EnumerationElementEditor.stringdata))
        return static_cast<void*>(const_cast< EnumerationElementEditor*>(this));
    return ElementEditor::qt_metacast(_clname);
}

int EnumerationElementEditor::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = ElementEditor::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 1)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 1;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 1)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 1;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
