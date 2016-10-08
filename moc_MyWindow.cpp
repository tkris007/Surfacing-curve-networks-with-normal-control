/****************************************************************************
** Meta object code from reading C++ file 'MyWindow.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.3.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "MyWindow.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MyWindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.3.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_MyWindow_t {
    QByteArrayData data[10];
    char stringdata[97];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_MyWindow_t, stringdata) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_MyWindow_t qt_meta_stringdata_MyWindow = {
    {
QT_MOC_LITERAL(0, 0, 8),
QT_MOC_LITERAL(1, 9, 4),
QT_MOC_LITERAL(2, 14, 0),
QT_MOC_LITERAL(3, 15, 9),
QT_MOC_LITERAL(4, 25, 8),
QT_MOC_LITERAL(5, 34, 16),
QT_MOC_LITERAL(6, 51, 7),
QT_MOC_LITERAL(7, 59, 14),
QT_MOC_LITERAL(8, 74, 7),
QT_MOC_LITERAL(9, 82, 14)
    },
    "MyWindow\0open\0\0setCutoff\0setRange\0"
    "startComputation\0message\0midComputation\0"
    "percent\0endComputation"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_MyWindow[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       6,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   44,    2, 0x08 /* Private */,
       3,    0,   45,    2, 0x08 /* Private */,
       4,    0,   46,    2, 0x08 /* Private */,
       5,    1,   47,    2, 0x08 /* Private */,
       7,    1,   50,    2, 0x08 /* Private */,
       9,    0,   53,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QString,    6,
    QMetaType::Void, QMetaType::Int,    8,
    QMetaType::Void,

       0        // eod
};

void MyWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        MyWindow *_t = static_cast<MyWindow *>(_o);
        switch (_id) {
        case 0: _t->open(); break;
        case 1: _t->setCutoff(); break;
        case 2: _t->setRange(); break;
        case 3: _t->startComputation((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 4: _t->midComputation((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 5: _t->endComputation(); break;
        default: ;
        }
    }
}

const QMetaObject MyWindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_MyWindow.data,
      qt_meta_data_MyWindow,  qt_static_metacall, 0, 0}
};


const QMetaObject *MyWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *MyWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MyWindow.stringdata))
        return static_cast<void*>(const_cast< MyWindow*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int MyWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 6)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 6;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 6)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 6;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
