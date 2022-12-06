/****************************************************************************
** Meta object code from reading C++ file 'mainwindow.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.12.8)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "include/mainwindow.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'mainwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.12.8. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_AnimatedButton_t {
    QByteArrayData data[5];
    char stringdata0[61];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_AnimatedButton_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_AnimatedButton_t qt_meta_stringdata_AnimatedButton = {
    {
QT_MOC_LITERAL(0, 0, 14), // "AnimatedButton"
QT_MOC_LITERAL(1, 15, 15), // "updateAnimation"
QT_MOC_LITERAL(2, 31, 0), // ""
QT_MOC_LITERAL(3, 32, 14), // "startAnimation"
QT_MOC_LITERAL(4, 47, 13) // "stopAnimation"

    },
    "AnimatedButton\0updateAnimation\0\0"
    "startAnimation\0stopAnimation"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_AnimatedButton[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    1,   29,    2, 0x08 /* Private */,
       3,    0,   32,    2, 0x0a /* Public */,
       4,    0,   33,    2, 0x0a /* Public */,

 // slots: parameters
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void AnimatedButton::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<AnimatedButton *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->updateAnimation((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: _t->startAnimation(); break;
        case 2: _t->stopAnimation(); break;
        default: ;
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject AnimatedButton::staticMetaObject = { {
    &QPushButton::staticMetaObject,
    qt_meta_stringdata_AnimatedButton.data,
    qt_meta_data_AnimatedButton,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *AnimatedButton::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *AnimatedButton::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_AnimatedButton.stringdata0))
        return static_cast<void*>(this);
    return QPushButton::qt_metacast(_clname);
}

int AnimatedButton::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QPushButton::qt_metacall(_c, _id, _a);
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
struct qt_meta_stringdata_C_TableViewDelegate_t {
    QByteArrayData data[1];
    char stringdata0[20];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_C_TableViewDelegate_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_C_TableViewDelegate_t qt_meta_stringdata_C_TableViewDelegate = {
    {
QT_MOC_LITERAL(0, 0, 19) // "C_TableViewDelegate"

    },
    "C_TableViewDelegate"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_C_TableViewDelegate[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
       0,    0, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

       0        // eod
};

void C_TableViewDelegate::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    Q_UNUSED(_o);
    Q_UNUSED(_id);
    Q_UNUSED(_c);
    Q_UNUSED(_a);
}

QT_INIT_METAOBJECT const QMetaObject C_TableViewDelegate::staticMetaObject = { {
    &QItemDelegate::staticMetaObject,
    qt_meta_stringdata_C_TableViewDelegate.data,
    qt_meta_data_C_TableViewDelegate,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *C_TableViewDelegate::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *C_TableViewDelegate::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_C_TableViewDelegate.stringdata0))
        return static_cast<void*>(this);
    return QItemDelegate::qt_metacast(_clname);
}

int C_TableViewDelegate::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QItemDelegate::qt_metacall(_c, _id, _a);
    return _id;
}
struct qt_meta_stringdata_MainWindow_t {
    QByteArrayData data[97];
    char stringdata0[1372];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_MainWindow_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_MainWindow_t qt_meta_stringdata_MainWindow = {
    {
QT_MOC_LITERAL(0, 0, 10), // "MainWindow"
QT_MOC_LITERAL(1, 11, 15), // "progress_append"
QT_MOC_LITERAL(2, 27, 0), // ""
QT_MOC_LITERAL(3, 28, 16), // "progress_replace"
QT_MOC_LITERAL(4, 45, 25), // "material3DLocationChanged"
QT_MOC_LITERAL(5, 71, 4), // "open"
QT_MOC_LITERAL(6, 76, 4), // "save"
QT_MOC_LITERAL(7, 81, 6), // "saveAs"
QT_MOC_LITERAL(8, 88, 10), // "saveScreen"
QT_MOC_LITERAL(9, 99, 11), // "importGoCad"
QT_MOC_LITERAL(10, 111, 11), // "exportVTU3D"
QT_MOC_LITERAL(11, 123, 11), // "exportVTU2D"
QT_MOC_LITERAL(12, 135, 12), // "exportFeFlow"
QT_MOC_LITERAL(13, 148, 9), // "exportOGS"
QT_MOC_LITERAL(14, 158, 9), // "exportTIN"
QT_MOC_LITERAL(15, 168, 12), // "exportCOMSOL"
QT_MOC_LITERAL(16, 181, 12), // "exportABAQUS"
QT_MOC_LITERAL(17, 194, 12), // "exportEXODUS"
QT_MOC_LITERAL(18, 207, 7), // "addUnit"
QT_MOC_LITERAL(19, 215, 8), // "addFault"
QT_MOC_LITERAL(20, 224, 9), // "addBorder"
QT_MOC_LITERAL(21, 234, 7), // "addWell"
QT_MOC_LITERAL(22, 242, 13), // "deleteSurface"
QT_MOC_LITERAL(23, 256, 5), // "about"
QT_MOC_LITERAL(24, 262, 5), // "reset"
QT_MOC_LITERAL(25, 268, 8), // "viewAxis"
QT_MOC_LITERAL(26, 277, 14), // "FillNameCombos"
QT_MOC_LITERAL(27, 292, 12), // "setUShowGBox"
QT_MOC_LITERAL(28, 305, 10), // "setUMeshes"
QT_MOC_LITERAL(29, 316, 12), // "setFShowGBox"
QT_MOC_LITERAL(30, 329, 10), // "setFMeshes"
QT_MOC_LITERAL(31, 340, 12), // "setBShowGBox"
QT_MOC_LITERAL(32, 353, 10), // "setBMeshes"
QT_MOC_LITERAL(33, 364, 12), // "setWShowGBox"
QT_MOC_LITERAL(34, 377, 10), // "setWMeshes"
QT_MOC_LITERAL(35, 388, 12), // "setMShowGBox"
QT_MOC_LITERAL(36, 401, 10), // "setMMeshes"
QT_MOC_LITERAL(37, 412, 12), // "setSShowGBox"
QT_MOC_LITERAL(38, 425, 22), // "interpolationSetMethod"
QT_MOC_LITERAL(39, 448, 17), // "interpolationFill"
QT_MOC_LITERAL(40, 466, 19), // "material3dFillValue"
QT_MOC_LITERAL(41, 486, 33), // "material3dSetLocationFromDSpi..."
QT_MOC_LITERAL(42, 520, 31), // "material3dSetLocationFromSlider"
QT_MOC_LITERAL(43, 552, 25), // "ExportRotationAngelUpdate"
QT_MOC_LITERAL(44, 578, 14), // "refinementFill"
QT_MOC_LITERAL(45, 593, 16), // "refinementUpdate"
QT_MOC_LITERAL(46, 610, 14), // "material1dFill"
QT_MOC_LITERAL(47, 625, 16), // "material1dUpdate"
QT_MOC_LITERAL(48, 642, 16), // "QListWidgetItem*"
QT_MOC_LITERAL(49, 659, 14), // "material2dFill"
QT_MOC_LITERAL(50, 674, 16), // "material2dUpdate"
QT_MOC_LITERAL(51, 691, 14), // "material3dFill"
QT_MOC_LITERAL(52, 706, 22), // "material3dFillLocation"
QT_MOC_LITERAL(53, 729, 8), // "Material"
QT_MOC_LITERAL(54, 738, 13), // "material3dAdd"
QT_MOC_LITERAL(55, 752, 16), // "material3dRemove"
QT_MOC_LITERAL(56, 769, 18), // "material3dIndexing"
QT_MOC_LITERAL(57, 788, 21), // "material3dLocationAdd"
QT_MOC_LITERAL(58, 810, 24), // "material3dLocationRemove"
QT_MOC_LITERAL(59, 835, 17), // "selectionSetFlags"
QT_MOC_LITERAL(60, 853, 9), // "selection"
QT_MOC_LITERAL(61, 863, 13), // "findSelection"
QT_MOC_LITERAL(62, 877, 21), // "preMeshGradientUpdate"
QT_MOC_LITERAL(63, 899, 18), // "meshGradientUpdate"
QT_MOC_LITERAL(64, 918, 17), // "settingCallByMenu"
QT_MOC_LITERAL(65, 936, 17), // "settingCallByTree"
QT_MOC_LITERAL(66, 954, 16), // "QTreeWidgetItem*"
QT_MOC_LITERAL(67, 971, 11), // "settingShow"
QT_MOC_LITERAL(68, 983, 5), // "phase"
QT_MOC_LITERAL(69, 989, 7), // "setting"
QT_MOC_LITERAL(70, 997, 7), // "preMesh"
QT_MOC_LITERAL(71, 1005, 4), // "Mesh"
QT_MOC_LITERAL(72, 1010, 7), // "readPLC"
QT_MOC_LITERAL(73, 1018, 7), // "meshPLC"
QT_MOC_LITERAL(74, 1026, 22), // "applyMaterialSelection"
QT_MOC_LITERAL(75, 1049, 9), // "fillTable"
QT_MOC_LITERAL(76, 1059, 20), // "onHighlightingChange"
QT_MOC_LITERAL(77, 1080, 12), // "tableClicked"
QT_MOC_LITERAL(78, 1093, 17), // "QTableWidgetItem*"
QT_MOC_LITERAL(79, 1111, 4), // "item"
QT_MOC_LITERAL(80, 1116, 14), // "clearErrorInfo"
QT_MOC_LITERAL(81, 1131, 15), // "updateErrorInfo"
QT_MOC_LITERAL(82, 1147, 7), // "message"
QT_MOC_LITERAL(83, 1155, 14), // "fillErrorTable"
QT_MOC_LITERAL(84, 1170, 17), // "errorTableClicked"
QT_MOC_LITERAL(85, 1188, 17), // "clearErrorMarkers"
QT_MOC_LITERAL(86, 1206, 12), // "FinishedRead"
QT_MOC_LITERAL(87, 1219, 21), // "threadFinishedPreMesh"
QT_MOC_LITERAL(88, 1241, 18), // "threadFinishedMesh"
QT_MOC_LITERAL(89, 1260, 12), // "callMakeMats"
QT_MOC_LITERAL(90, 1273, 12), // "callMakeTets"
QT_MOC_LITERAL(91, 1286, 10), // "updateInfo"
QT_MOC_LITERAL(92, 1297, 7), // "replace"
QT_MOC_LITERAL(93, 1305, 17), // "printErrorMessage"
QT_MOC_LITERAL(94, 1323, 22), // "onRotationCenterChange"
QT_MOC_LITERAL(95, 1346, 5), // "index"
QT_MOC_LITERAL(96, 1352, 19) // "onSensitivityChange"

    },
    "MainWindow\0progress_append\0\0"
    "progress_replace\0material3DLocationChanged\0"
    "open\0save\0saveAs\0saveScreen\0importGoCad\0"
    "exportVTU3D\0exportVTU2D\0exportFeFlow\0"
    "exportOGS\0exportTIN\0exportCOMSOL\0"
    "exportABAQUS\0exportEXODUS\0addUnit\0"
    "addFault\0addBorder\0addWell\0deleteSurface\0"
    "about\0reset\0viewAxis\0FillNameCombos\0"
    "setUShowGBox\0setUMeshes\0setFShowGBox\0"
    "setFMeshes\0setBShowGBox\0setBMeshes\0"
    "setWShowGBox\0setWMeshes\0setMShowGBox\0"
    "setMMeshes\0setSShowGBox\0interpolationSetMethod\0"
    "interpolationFill\0material3dFillValue\0"
    "material3dSetLocationFromDSpinBox\0"
    "material3dSetLocationFromSlider\0"
    "ExportRotationAngelUpdate\0refinementFill\0"
    "refinementUpdate\0material1dFill\0"
    "material1dUpdate\0QListWidgetItem*\0"
    "material2dFill\0material2dUpdate\0"
    "material3dFill\0material3dFillLocation\0"
    "Material\0material3dAdd\0material3dRemove\0"
    "material3dIndexing\0material3dLocationAdd\0"
    "material3dLocationRemove\0selectionSetFlags\0"
    "selection\0findSelection\0preMeshGradientUpdate\0"
    "meshGradientUpdate\0settingCallByMenu\0"
    "settingCallByTree\0QTreeWidgetItem*\0"
    "settingShow\0phase\0setting\0preMesh\0"
    "Mesh\0readPLC\0meshPLC\0applyMaterialSelection\0"
    "fillTable\0onHighlightingChange\0"
    "tableClicked\0QTableWidgetItem*\0item\0"
    "clearErrorInfo\0updateErrorInfo\0message\0"
    "fillErrorTable\0errorTableClicked\0"
    "clearErrorMarkers\0FinishedRead\0"
    "threadFinishedPreMesh\0threadFinishedMesh\0"
    "callMakeMats\0callMakeTets\0updateInfo\0"
    "replace\0printErrorMessage\0"
    "onRotationCenterChange\0index\0"
    "onSensitivityChange"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_MainWindow[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
      86,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       3,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,  444,    2, 0x06 /* Public */,
       3,    1,  447,    2, 0x06 /* Public */,
       4,    2,  450,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       5,    0,  455,    2, 0x08 /* Private */,
       6,    0,  456,    2, 0x08 /* Private */,
       7,    0,  457,    2, 0x08 /* Private */,
       8,    0,  458,    2, 0x08 /* Private */,
       9,    0,  459,    2, 0x08 /* Private */,
      10,    0,  460,    2, 0x08 /* Private */,
      11,    0,  461,    2, 0x08 /* Private */,
      12,    0,  462,    2, 0x08 /* Private */,
      13,    0,  463,    2, 0x08 /* Private */,
      14,    0,  464,    2, 0x08 /* Private */,
      15,    0,  465,    2, 0x08 /* Private */,
      16,    0,  466,    2, 0x08 /* Private */,
      17,    0,  467,    2, 0x08 /* Private */,
      18,    0,  468,    2, 0x08 /* Private */,
      19,    0,  469,    2, 0x08 /* Private */,
      20,    0,  470,    2, 0x08 /* Private */,
      21,    0,  471,    2, 0x08 /* Private */,
      22,    0,  472,    2, 0x08 /* Private */,
      23,    0,  473,    2, 0x08 /* Private */,
      24,    0,  474,    2, 0x08 /* Private */,
      25,    0,  475,    2, 0x08 /* Private */,
      26,    0,  476,    2, 0x08 /* Private */,
      27,    1,  477,    2, 0x08 /* Private */,
      28,    0,  480,    2, 0x08 /* Private */,
      29,    1,  481,    2, 0x08 /* Private */,
      30,    0,  484,    2, 0x08 /* Private */,
      31,    1,  485,    2, 0x08 /* Private */,
      32,    0,  488,    2, 0x08 /* Private */,
      33,    1,  489,    2, 0x08 /* Private */,
      34,    0,  492,    2, 0x08 /* Private */,
      35,    1,  493,    2, 0x08 /* Private */,
      36,    0,  496,    2, 0x08 /* Private */,
      37,    1,  497,    2, 0x08 /* Private */,
      38,    1,  500,    2, 0x08 /* Private */,
      39,    0,  503,    2, 0x08 /* Private */,
      40,    1,  504,    2, 0x08 /* Private */,
      41,    1,  507,    2, 0x08 /* Private */,
      42,    1,  510,    2, 0x08 /* Private */,
      43,    1,  513,    2, 0x08 /* Private */,
      44,    0,  516,    2, 0x08 /* Private */,
      45,    1,  517,    2, 0x08 /* Private */,
      46,    0,  520,    2, 0x08 /* Private */,
      47,    1,  521,    2, 0x08 /* Private */,
      49,    0,  524,    2, 0x08 /* Private */,
      50,    1,  525,    2, 0x08 /* Private */,
      51,    0,  528,    2, 0x08 /* Private */,
      52,    1,  529,    2, 0x08 /* Private */,
      54,    0,  532,    2, 0x08 /* Private */,
      55,    0,  533,    2, 0x08 /* Private */,
      56,    0,  534,    2, 0x08 /* Private */,
      57,    0,  535,    2, 0x08 /* Private */,
      58,    0,  536,    2, 0x08 /* Private */,
      59,    1,  537,    2, 0x08 /* Private */,
      60,    1,  540,    2, 0x08 /* Private */,
      61,    3,  543,    2, 0x08 /* Private */,
      62,    1,  550,    2, 0x08 /* Private */,
      63,    1,  553,    2, 0x08 /* Private */,
      64,    0,  556,    2, 0x08 /* Private */,
      65,    2,  557,    2, 0x08 /* Private */,
      67,    2,  562,    2, 0x08 /* Private */,
      70,    0,  567,    2, 0x08 /* Private */,
      71,    0,  568,    2, 0x08 /* Private */,
      72,    0,  569,    2, 0x08 /* Private */,
      73,    0,  570,    2, 0x08 /* Private */,
      74,    0,  571,    2, 0x08 /* Private */,
      75,    0,  572,    2, 0x08 /* Private */,
      76,    0,  573,    2, 0x08 /* Private */,
      77,    1,  574,    2, 0x08 /* Private */,
      80,    0,  577,    2, 0x08 /* Private */,
      81,    1,  578,    2, 0x08 /* Private */,
      83,    0,  581,    2, 0x08 /* Private */,
      84,    1,  582,    2, 0x08 /* Private */,
      85,    0,  585,    2, 0x08 /* Private */,
      86,    0,  586,    2, 0x08 /* Private */,
      87,    0,  587,    2, 0x08 /* Private */,
      88,    0,  588,    2, 0x08 /* Private */,
      89,    0,  589,    2, 0x08 /* Private */,
      90,    0,  590,    2, 0x08 /* Private */,
      91,    0,  591,    2, 0x08 /* Private */,
      92,    1,  592,    2, 0x08 /* Private */,
      93,    1,  595,    2, 0x08 /* Private */,
      94,    1,  598,    2, 0x08 /* Private */,
      96,    0,  601,    2, 0x08 /* Private */,

 // signals: parameters
    QMetaType::Void, QMetaType::QString,    2,
    QMetaType::Void, QMetaType::QString,    2,
    QMetaType::Void, QMetaType::QString, QMetaType::Double,    2,    2,

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
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QString,    2,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QString,    2,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QString,    2,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QString,    2,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QString,    2,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QString,    2,
    QMetaType::Void, QMetaType::QString,    2,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Double,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Double,    2,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Double,    2,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 48,    2,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 48,    2,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,   53,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::UChar, QMetaType::UChar, QMetaType::UChar,    2,    2,    2,
    QMetaType::Void, QMetaType::Double,    2,
    QMetaType::Void, QMetaType::Double,    2,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 66, QMetaType::Int,    2,    2,
    QMetaType::Void, QMetaType::QString, QMetaType::QString,   68,   69,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 78,   79,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QString,   82,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 78,   79,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QString,    2,
    QMetaType::Void, QMetaType::QString,   82,
    QMetaType::Void, QMetaType::Int,   95,
    QMetaType::Void,

       0        // eod
};

void MainWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<MainWindow *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->progress_append((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 1: _t->progress_replace((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 2: _t->material3DLocationChanged((*reinterpret_cast< QString(*)>(_a[1])),(*reinterpret_cast< double(*)>(_a[2]))); break;
        case 3: _t->open(); break;
        case 4: _t->save(); break;
        case 5: _t->saveAs(); break;
        case 6: _t->saveScreen(); break;
        case 7: _t->importGoCad(); break;
        case 8: _t->exportVTU3D(); break;
        case 9: _t->exportVTU2D(); break;
        case 10: _t->exportFeFlow(); break;
        case 11: _t->exportOGS(); break;
        case 12: _t->exportTIN(); break;
        case 13: _t->exportCOMSOL(); break;
        case 14: _t->exportABAQUS(); break;
        case 15: _t->exportEXODUS(); break;
        case 16: _t->addUnit(); break;
        case 17: _t->addFault(); break;
        case 18: _t->addBorder(); break;
        case 19: _t->addWell(); break;
        case 20: _t->deleteSurface(); break;
        case 21: _t->about(); break;
        case 22: _t->reset(); break;
        case 23: _t->viewAxis(); break;
        case 24: _t->FillNameCombos(); break;
        case 25: _t->setUShowGBox((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 26: _t->setUMeshes(); break;
        case 27: _t->setFShowGBox((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 28: _t->setFMeshes(); break;
        case 29: _t->setBShowGBox((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 30: _t->setBMeshes(); break;
        case 31: _t->setWShowGBox((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 32: _t->setWMeshes(); break;
        case 33: _t->setMShowGBox((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 34: _t->setMMeshes(); break;
        case 35: _t->setSShowGBox((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 36: _t->interpolationSetMethod((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 37: _t->interpolationFill(); break;
        case 38: _t->material3dFillValue((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 39: _t->material3dSetLocationFromDSpinBox((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 40: _t->material3dSetLocationFromSlider((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 41: _t->ExportRotationAngelUpdate((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 42: _t->refinementFill(); break;
        case 43: _t->refinementUpdate((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 44: _t->material1dFill(); break;
        case 45: _t->material1dUpdate((*reinterpret_cast< QListWidgetItem*(*)>(_a[1]))); break;
        case 46: _t->material2dFill(); break;
        case 47: _t->material2dUpdate((*reinterpret_cast< QListWidgetItem*(*)>(_a[1]))); break;
        case 48: _t->material3dFill(); break;
        case 49: _t->material3dFillLocation((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 50: _t->material3dAdd(); break;
        case 51: _t->material3dRemove(); break;
        case 52: _t->material3dIndexing(); break;
        case 53: _t->material3dLocationAdd(); break;
        case 54: _t->material3dLocationRemove(); break;
        case 55: _t->selectionSetFlags((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 56: _t->selection((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 57: _t->findSelection((*reinterpret_cast< unsigned char(*)>(_a[1])),(*reinterpret_cast< unsigned char(*)>(_a[2])),(*reinterpret_cast< unsigned char(*)>(_a[3]))); break;
        case 58: _t->preMeshGradientUpdate((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 59: _t->meshGradientUpdate((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 60: _t->settingCallByMenu(); break;
        case 61: _t->settingCallByTree((*reinterpret_cast< QTreeWidgetItem*(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 62: _t->settingShow((*reinterpret_cast< QString(*)>(_a[1])),(*reinterpret_cast< QString(*)>(_a[2]))); break;
        case 63: _t->preMesh(); break;
        case 64: _t->Mesh(); break;
        case 65: _t->readPLC(); break;
        case 66: _t->meshPLC(); break;
        case 67: _t->applyMaterialSelection(); break;
        case 68: _t->fillTable(); break;
        case 69: _t->onHighlightingChange(); break;
        case 70: _t->tableClicked((*reinterpret_cast< QTableWidgetItem*(*)>(_a[1]))); break;
        case 71: _t->clearErrorInfo(); break;
        case 72: _t->updateErrorInfo((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 73: _t->fillErrorTable(); break;
        case 74: _t->errorTableClicked((*reinterpret_cast< QTableWidgetItem*(*)>(_a[1]))); break;
        case 75: _t->clearErrorMarkers(); break;
        case 76: _t->FinishedRead(); break;
        case 77: _t->threadFinishedPreMesh(); break;
        case 78: _t->threadFinishedMesh(); break;
        case 79: _t->callMakeMats(); break;
        case 80: _t->callMakeTets(); break;
        case 81: _t->updateInfo(); break;
        case 82: _t->replace((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 83: _t->printErrorMessage((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 84: _t->onRotationCenterChange((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 85: _t->onSensitivityChange(); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        {
            using _t = void (MainWindow::*)(QString );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&MainWindow::progress_append)) {
                *result = 0;
                return;
            }
        }
        {
            using _t = void (MainWindow::*)(QString );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&MainWindow::progress_replace)) {
                *result = 1;
                return;
            }
        }
        {
            using _t = void (MainWindow::*)(QString , double );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&MainWindow::material3DLocationChanged)) {
                *result = 2;
                return;
            }
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject MainWindow::staticMetaObject = { {
    &QMainWindow::staticMetaObject,
    qt_meta_stringdata_MainWindow.data,
    qt_meta_data_MainWindow,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *MainWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *MainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_MainWindow.stringdata0))
        return static_cast<void*>(this);
    return QMainWindow::qt_metacast(_clname);
}

int MainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 86)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 86;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 86)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 86;
    }
    return _id;
}

// SIGNAL 0
void MainWindow::progress_append(QString _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void MainWindow::progress_replace(QString _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void MainWindow::material3DLocationChanged(QString _t1, double _t2)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}
struct qt_meta_stringdata_C_Thread_t {
    QByteArrayData data[1];
    char stringdata0[9];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_C_Thread_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_C_Thread_t qt_meta_stringdata_C_Thread = {
    {
QT_MOC_LITERAL(0, 0, 8) // "C_Thread"

    },
    "C_Thread"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_C_Thread[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
       0,    0, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

       0        // eod
};

void C_Thread::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    Q_UNUSED(_o);
    Q_UNUSED(_id);
    Q_UNUSED(_c);
    Q_UNUSED(_a);
}

QT_INIT_METAOBJECT const QMetaObject C_Thread::staticMetaObject = { {
    &QThread::staticMetaObject,
    qt_meta_stringdata_C_Thread.data,
    qt_meta_data_C_Thread,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *C_Thread::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *C_Thread::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_C_Thread.stringdata0))
        return static_cast<void*>(this);
    return QThread::qt_metacast(_clname);
}

int C_Thread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QThread::qt_metacall(_c, _id, _a);
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
