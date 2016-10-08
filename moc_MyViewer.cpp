/****************************************************************************
** Meta object code from reading C++ file 'MyViewer.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.3.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "MyViewer.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MyViewer.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.3.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_MyViewer_t {
    QByteArrayData data[7];
    char stringdata[73];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_MyViewer_t, stringdata) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_MyViewer_t qt_meta_stringdata_MyViewer = {
    {
QT_MOC_LITERAL(0, 0, 8),
QT_MOC_LITERAL(1, 9, 16),
QT_MOC_LITERAL(2, 26, 0),
QT_MOC_LITERAL(3, 27, 7),
QT_MOC_LITERAL(4, 35, 14),
QT_MOC_LITERAL(5, 50, 7),
QT_MOC_LITERAL(6, 58, 14)
    },
    "MyViewer\0startComputation\0\0message\0"
    "midComputation\0percent\0endComputation"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_MyViewer[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       3,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   29,    2, 0x06 /* Public */,
       4,    1,   32,    2, 0x06 /* Public */,
       6,    0,   35,    2, 0x06 /* Public */,

 // signals: parameters
    QMetaType::Void, QMetaType::QString,    3,
    QMetaType::Void, QMetaType::Int,    5,
    QMetaType::Void,

       0        // eod
};

void MyViewer::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        MyViewer *_t = static_cast<MyViewer *>(_o);
        switch (_id) {
        case 0: _t->startComputation((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 1: _t->midComputation((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: _t->endComputation(); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (MyViewer::*_t)(QString );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&MyViewer::startComputation)) {
                *result = 0;
            }
        }
        {
            typedef void (MyViewer::*_t)(int );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&MyViewer::midComputation)) {
                *result = 1;
            }
        }
        {
            typedef void (MyViewer::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&MyViewer::endComputation)) {
                *result = 2;
            }
        }
    }
}

const QMetaObject MyViewer::staticMetaObject = {
    { &QGLViewer::staticMetaObject, qt_meta_stringdata_MyViewer.data,
      qt_meta_data_MyViewer,  qt_static_metacall, 0, 0}
};


const QMetaObject *MyViewer::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *MyViewer::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MyViewer.stringdata))
        return static_cast<void*>(const_cast< MyViewer*>(this));
    return QGLViewer::qt_metacast(_clname);
}

int MyViewer::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QGLViewer::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 3)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 3;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 3)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 3;
    }
    return _id;
}

// SIGNAL 0
void MyViewer::startComputation(QString _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void MyViewer::midComputation(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void MyViewer::endComputation()
{
    QMetaObject::activate(this, &staticMetaObject, 2, 0);
}
QT_END_MOC_NAMESPACE
