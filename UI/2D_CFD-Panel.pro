QT       += core gui opengl openglwidgets

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    Boundary.cpp \
    Domain.cpp \
    Mesh.cpp \
    MeshList.cpp \
    boundaryList.cpp \
    comboboxdelegate.cpp \
    domainlist.cpp \
    main.cpp \
    mainwindow.cpp \
    newopenglwidget.cpp

HEADERS += \
    Boundary.h \
    Domain.h \
    Mesh.h \
    MeshList.h \
    boundaryList.h \
    comboboxdelegate.h \
    domainlist.h \
    mainwindow.h \
    newopenglwidget.h

FORMS += \
    mainwindow.ui

LIBS += -lglut -lGLU

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
