#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QImage>
#include <QTextCodec>
#include "math.h"
#include "stdlib.h"


namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow *ui;
    QImage *winimg;

    int winimgheight;
    int winimgwidth;
    int** winimgpixelmat;
    int** framepixelmat;

	bool geometricflag;

    bool winimgtopixelmat();
    void pixelmattowinimg();
    void initframepixelmat();
    void framepixelmattoshow();
    bool winimgshow();
    bool winimgpixelmatshow();

private slots:
    bool on_slotOpenImage_triggered();
    void on_slotSaveImage_triggered();
    void on_slotFreshImage_triggered();
    void on_slotExit_triggered();
    void on_slotHomepage_triggered();
    void on_slotDemo_triggered();
    void on_slotInfo_triggered();

	bool on_slotMidptvalueBtn_triggered();
    bool on_slotMaxvalueBtn_triggered();
    bool on_slotMinvalueBtn_triggered();
    bool on_slotAvgvalueBtn_triggered();
    bool on_slotVarvalueBtn_triggered();
    bool on_slotHistogramBtn_triggered();
	bool on_slotBinaryBtn_triggered();
    bool on_slotHistequalBtn_triggered();
    bool on_slotGraystretchBtn_triggered();
    bool on_slotCentreBtn_triggered();
    bool on_slotRotationBtn_triggered();
    bool on_slotMoveBtn_triggered();
    bool on_slotScalingBtn_triggered();
    bool on_slotDFTBtn_triggered();
    bool on_slotDCTBtn_triggered();
    bool on_slotHaarBtn_triggered();
    bool on_slotWalshBtn_triggered();
    bool on_slotMedianfitBtn_triggered();
    bool on_slotAveragefitBtn_triggered();
    bool on_slotLowpassfitBtn_triggered();
    bool on_slotSobelBtn_triggered();
    bool on_slotLaplaceBtn_triggered();
    bool on_slotHighpassfitBtn_triggered();
    bool on_slotRandomnoiseBtn_triggered();
    bool on_slotImpulsenoiseBtn_triggered();
    bool on_slotInversefitBtn_triggered();
    bool on_slotWienerfitBtn_triggered();

    bool on_slotJpegBtn_triggered();
};





#endif // MAINWINDOW_H
