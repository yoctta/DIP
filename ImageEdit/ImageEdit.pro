#-------------------------------------------------
#
# Project created by QtCreator 2016-01-03T16:39:12
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = ImageEdit
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
        expcode.cpp \

HEADERS  += mainwindow.h \
    expcode.h

FORMS    += \
    mainwindow.ui

RC_FILE += ImageEdit.rc
