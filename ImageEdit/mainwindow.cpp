#include "mainwindow.h"
#include "expcode.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QDesktopServices>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    winimg = new QImage();

	geometricflag = true;

    QMenu *parametermenu = new QMenu(this);
    
	QAction *MidptvalueBtn = parametermenu->addAction(QStringLiteral("中心点灰度(示例)"));
    QAction *MaxvalueBtn = parametermenu->addAction(QStringLiteral("最大值"));
    QAction *MinvalueBtn = parametermenu->addAction(QStringLiteral("最小值"));
    QAction *AvgvalueBtn = parametermenu->addAction(QStringLiteral("平均值"));
    QAction *VarvalueBtn = parametermenu->addAction(QStringLiteral("方差"));
    QAction *HistogramBtn = parametermenu->addAction(QStringLiteral("直方图统计"));

    QMenu *enhancemenu = new QMenu(this);
    QMenu *contrastmenu = enhancemenu->addMenu(QStringLiteral("对比度增强"));
	QAction *BinaryBtn = contrastmenu->addAction(QStringLiteral("二值化(示例)"));
    QAction *HistequalBtn = contrastmenu->addAction(QStringLiteral("直方图均衡"));
    QAction *GraystretchBtn = contrastmenu->addAction(QStringLiteral("灰度拉伸"));
    QMenu *smoothmenu = enhancemenu->addMenu(QStringLiteral("图像平滑"));
    QAction *MedianfitBtn = smoothmenu->addAction(QStringLiteral("中值滤波"));
    QAction *AveragefitBtn = smoothmenu->addAction(QStringLiteral("均值滤波"));
    QAction *LowpassfitBtn = smoothmenu->addAction(QStringLiteral("理想低通滤波"));
    QMenu *sharpenmenu = enhancemenu->addMenu(QStringLiteral("图像锐化"));
    QAction *SobelBtn = sharpenmenu->addAction(QStringLiteral("Sobel算子"));
    QAction *LaplaceBtn = sharpenmenu->addAction(QStringLiteral("Laplace算子"));
    QAction *HighpassfitBtn = sharpenmenu->addAction(QStringLiteral("理想高通滤波"));

    QMenu *transformmenu = new QMenu(this);
    QMenu *geometrymenu = transformmenu->addMenu(QStringLiteral("几何变换"));
    QAction *CentreBtn = geometrymenu->addAction(QStringLiteral("平移至中心(示例)"));
    QAction *RotationBtn = geometrymenu->addAction(QStringLiteral("旋转"));
    QAction *MoveBtn = geometrymenu->addAction(QStringLiteral("平移"));
    QAction *ScalingBtn = geometrymenu->addAction(QStringLiteral("缩放"));
    QMenu *linearmenu = transformmenu->addMenu(QStringLiteral("线性变换"));
    QAction *DFTBtn = linearmenu->addAction(QStringLiteral("DFT"));
    QAction *DCTBtn = linearmenu->addAction(QStringLiteral("DCT"));
    QAction *HaarBtn = linearmenu->addAction(QStringLiteral("HAAR"));
    QAction *WalshBtn = linearmenu->addAction(QStringLiteral("WALSH"));

    QMenu *recovermenu = new QMenu(this);
    QAction *RandomnoiseBtn = recovermenu->addAction(QStringLiteral("生成随机噪声"));
    QAction *ImpulsenoiseBtn = recovermenu->addAction(QStringLiteral("生成椒盐噪声"));
    QAction *InversefitBtn = recovermenu->addAction(QStringLiteral("逆滤波复原"));
    QAction *WienerfitBtn = recovermenu->addAction(QStringLiteral("维纳滤波复原"));

    QMenu *codingmenu = new QMenu(this);
    QAction *JpegBtn = codingmenu->addAction(QStringLiteral("JPEG压缩(示例)"));

    ui->ParameterBtn->setMenu(parametermenu);
    ui->EnhanceBtn->setMenu(enhancemenu);
    ui->TransformBtn->setMenu(transformmenu);
    ui->RecoverBtn->setMenu(recovermenu);
    ui->CodingBtn->setMenu(codingmenu);

    connect(ui->OpenImage,SIGNAL(triggered()),this,SLOT(on_slotOpenImage_triggered()));
    connect(ui->SaveImage,SIGNAL(triggered()),this,SLOT(on_slotSaveImage_triggered()));
    connect(ui->FreshImage,SIGNAL(triggered()),this,SLOT(on_slotFreshImage_triggered()));
    connect(ui->ActionExit,SIGNAL(triggered()),this,SLOT(on_slotExit_triggered()));
    connect(ui->Actionhomepage,SIGNAL(triggered()),this,SLOT(on_slotHomepage_triggered()));
    connect(ui->ActionDemo,SIGNAL(triggered()),this,SLOT(on_slotDemo_triggered()));
    connect(ui->ActionInfo,SIGNAL(triggered()),this,SLOT(on_slotInfo_triggered()));

    connect(MidptvalueBtn,SIGNAL(triggered()), this, SLOT(on_slotMidptvalueBtn_triggered()));
    connect(MaxvalueBtn,SIGNAL(triggered()),this,SLOT(on_slotMaxvalueBtn_triggered()));
    connect(MinvalueBtn,SIGNAL(triggered()),this,SLOT(on_slotMinvalueBtn_triggered()));
    connect(AvgvalueBtn,SIGNAL(triggered()),this,SLOT(on_slotAvgvalueBtn_triggered()));
    connect(VarvalueBtn,SIGNAL(triggered()),this,SLOT(on_slotVarvalueBtn_triggered()));
    connect(HistogramBtn,SIGNAL(triggered()),this,SLOT(on_slotHistogramBtn_triggered()));
	connect(BinaryBtn,SIGNAL(triggered()),this,SLOT(on_slotBinaryBtn_triggered()));
    connect(HistequalBtn,SIGNAL(triggered()),this,SLOT(on_slotHistequalBtn_triggered()));
    connect(GraystretchBtn,SIGNAL(triggered()),this,SLOT(on_slotGraystretchBtn_triggered()));
    connect(CentreBtn,SIGNAL(triggered()),this,SLOT(on_slotCentreBtn_triggered()));
    connect(RotationBtn,SIGNAL(triggered()),this,SLOT(on_slotRotationBtn_triggered()));
    connect(MoveBtn,SIGNAL(triggered()),this,SLOT(on_slotMoveBtn_triggered()));
    connect(ScalingBtn,SIGNAL(triggered()),this,SLOT(on_slotScalingBtn_triggered()));
    connect(DFTBtn,SIGNAL(triggered()),this,SLOT(on_slotDFTBtn_triggered()));
    connect(DCTBtn,SIGNAL(triggered()),this,SLOT(on_slotDCTBtn_triggered()));
    connect(WalshBtn,SIGNAL(triggered()),this,SLOT(on_slotWalshBtn_triggered()));
    connect(HaarBtn,SIGNAL(triggered()),this,SLOT(on_slotHaarBtn_triggered()));
    connect(MedianfitBtn,SIGNAL(triggered()),this,SLOT(on_slotMedianfitBtn_triggered()));
    connect(AveragefitBtn,SIGNAL(triggered()),this,SLOT(on_slotAveragefitBtn_triggered()));
    connect(LowpassfitBtn,SIGNAL(triggered()),this,SLOT(on_slotLowpassfitBtn_triggered()));
    connect(SobelBtn,SIGNAL(triggered()),this,SLOT(on_slotSobelBtn_triggered()));
    connect(LaplaceBtn,SIGNAL(triggered()),this,SLOT(on_slotLaplaceBtn_triggered()));
    connect(HighpassfitBtn,SIGNAL(triggered()),this,SLOT(on_slotHighpassfitBtn_triggered()));
    connect(RandomnoiseBtn,SIGNAL(triggered()),this,SLOT(on_slotRandomnoiseBtn_triggered()));
    connect(ImpulsenoiseBtn,SIGNAL(triggered()),this,SLOT(on_slotImpulsenoiseBtn_triggered()));
    connect(InversefitBtn,SIGNAL(triggered()),this,SLOT(on_slotInversefitBtn_triggered()));
    connect(WienerfitBtn,SIGNAL(triggered()),this,SLOT(on_slotWienerfitBtn_triggered()));

    connect(JpegBtn,SIGNAL(triggered()),this,SLOT(on_slotJpegBtn_triggered()));
}

MainWindow::~MainWindow()
{
    delete ui;
}

bool MainWindow::on_slotOpenImage_triggered()
{
    QString filename = QFileDialog::getOpenFileName(this,tr("Open File Name"),".",tr("Image File(*.bmp *.jpg)"));
    if(filename!="")
    {
        if(winimg->load(filename))
        {
            if(winimg->isNull())
            {
                int ret = QMessageBox::warning(this,QStringLiteral("错误"),QStringLiteral("无法打开图片"),QMessageBox::Abort);
                if(ret == QMessageBox::Abort)
                    return false;
            }
            QGraphicsScene *scene = new QGraphicsScene;
            scene->addPixmap(QPixmap::fromImage(*winimg));
            ui->graphicsView->setAlignment(Qt::AlignLeft | Qt::AlignTop);
            ui->graphicsView->setScene(scene);
            ui->graphicsView->show();

			geometricflag = true;
            return true;
        }
        else    return false;
    }
    else    return false;
}

void MainWindow::on_slotSaveImage_triggered()
{
    QString filename = QFileDialog::getSaveFileName(this,tr("Save File"),tr("filename"),tr("*.bmp"));
    if(filename!="")
    {
        winimg->save(filename + ".bmp");
    }
}

void MainWindow::on_slotFreshImage_triggered()
{
    winimg = new QImage(256, 256, QImage::Format_Indexed8);
    winimg->fill(255);
    winimgshow();
}

void MainWindow::on_slotExit_triggered()
{
    qApp->quit();
}

void MainWindow::on_slotHomepage_triggered()
{
    QString weburl = "http://image.ustc.edu.cn";
    QDesktopServices::openUrl(QUrl(weburl));
}

void MainWindow::on_slotDemo_triggered()
{
    QString weburl = "http://image.ustc.edu.cn/course/dip/DIP%E5%AE%9E%E9%AA%8C%E8%AF%B4%E6%98%8E.pdf";
    QDesktopServices::openUrl(QUrl(weburl));
}

void MainWindow::on_slotInfo_triggered()
{
    QMessageBox msg(QMessageBox::NoIcon, QStringLiteral("平台信息"), QStringLiteral("数字图像处理实验平台ver1.0\nupdated by QingKe\n20160321"));
    msg.setIconPixmap(QPixmap("ImageEdit.ico"));
    msg.exec();
}

bool MainWindow::on_slotMidptvalueBtn_triggered()
{
	if(winimgtopixelmat() == false)
	{
		return false;
	}
	int ans = midptvalue(winimgpixelmat, winimgheight, winimgwidth);
	QMessageBox::about(this, QStringLiteral(" "), QStringLiteral("中心点灰度值：%1").arg(ans));
	return true;
}

bool MainWindow::on_slotMaxvalueBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
    int ans = maxvalue(winimgpixelmat, winimgheight, winimgwidth);
    QMessageBox::about(this, QStringLiteral(" "), QStringLiteral("最大值：%1").arg(ans));
    return true;
}


bool MainWindow::on_slotMinvalueBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
    int ans = minvalue(winimgpixelmat, winimgheight, winimgwidth);
    QMessageBox::about(this, QStringLiteral(" "), QStringLiteral("最小值：%1").arg(ans));
    return true;
}

bool MainWindow::on_slotAvgvalueBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
    float ans = avgvalue(winimgpixelmat, winimgheight, winimgwidth);
    QMessageBox::about(this, QStringLiteral(" "), QStringLiteral("平均值：%1").arg(ans));
    return true;
}

bool MainWindow::on_slotVarvalueBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
    float ans = varvalue(winimgpixelmat, winimgheight, winimgwidth);
    QMessageBox::about(this, QStringLiteral(" "), QStringLiteral("方差：%1").arg(ans));
    return true;
}

bool MainWindow::on_slotHistogramBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
    int* histarr = new int[256];
    histarr = histogram(winimgpixelmat, winimgheight, winimgwidth);
    if(histarr == NULL)
    {
        return false;
    }
    int histmax = 0;
    for(int i = 0; i < 256; i++)
    {
        if(histmax < histarr[i])
        {
            histmax = histarr[i];
        }
    }
    float magratio = 512.0 / float(histmax);

    winimgheight = 512;
    winimgwidth = 512;
    winimgpixelmat = new int*[winimgheight];
    for(int i = 0; i < winimgheight; i++)
    {
        winimgpixelmat[i] = new int[winimgwidth];
        for(int j = 0; j < winimgwidth; j++)
        {
            winimgpixelmat[i][j] = 255;
        }
    }

    for(int j = 0; j < 256; j++)
    {
        for(int i = 512 - 1; i > 512 - histarr[j] * magratio - 1 && i >= 0; i--)
        {
            winimgpixelmat[i][j * 2] = 0;
            winimgpixelmat[i][j * 2 + 1] = 0;
        }
    }
    if(winimgpixelmatshow() == false)
    {
        return false;
    }
    return true;
}

bool MainWindow::on_slotBinaryBtn_triggered()
{
	if(winimgtopixelmat() == false)
	{
		return false;
	}
	winimgpixelmat = binaryimg(winimgpixelmat, winimgheight, winimgwidth);
	if(winimgpixelmatshow() == false)
	{
		return false;
	}
	return true;
}

bool MainWindow::on_slotHistequalBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
    winimgpixelmat = histogramequ(winimgpixelmat, winimgheight, winimgwidth);
    if(winimgpixelmatshow() == false)
    {
        return false;
    }
    return true;
}


bool MainWindow::on_slotGraystretchBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
    winimgpixelmat = graystretch(winimgpixelmat, winimgheight, winimgwidth);
    if(winimgpixelmatshow() == false)
    {
        return false;
    }
    return true;
}

bool MainWindow::on_slotCentreBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
	if(!geometricflag)
	{
		int ret = QMessageBox::warning(this,QStringLiteral("错误"),QStringLiteral("请重新打开新的图片重置显示框"),QMessageBox::Abort);
		if(ret == QMessageBox::Abort)
			return false;
	}
    initframepixelmat();
    framepixelmat = centralize(framepixelmat, winimgpixelmat, winimgheight, winimgwidth);
    if(framepixelmat == NULL)
    {
        return false;
    }
	geometricflag = false;
    framepixelmattoshow();
    return true;
}

bool MainWindow::on_slotRotationBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
	if(!geometricflag)
	{
		int ret = QMessageBox::warning(this,QStringLiteral("错误"),QStringLiteral("请重新打开新的图片重置显示框"),QMessageBox::Abort);
		if(ret == QMessageBox::Abort)
			return false;
	}
    initframepixelmat();
    framepixelmat = rotation(framepixelmat, winimgpixelmat, winimgheight, winimgwidth);
    if(framepixelmat == NULL)
    {
        return false;
    }
	geometricflag = false;
    framepixelmattoshow();
    return true;
}


bool MainWindow::on_slotMoveBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
	if(!geometricflag)
	{
		int ret = QMessageBox::warning(this,QStringLiteral("错误"),QStringLiteral("请重新打开新的图片重置显示框"),QMessageBox::Abort);
		if(ret == QMessageBox::Abort)
			return false;
	}
    initframepixelmat();
    framepixelmat = moveimage(framepixelmat, winimgpixelmat, winimgheight, winimgwidth);
    if(framepixelmat == NULL)
    {
        return false;
    }
	geometricflag = false;
    framepixelmattoshow();
    return true;
}


bool MainWindow::on_slotScalingBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
	if(!geometricflag)
	{
		int ret = QMessageBox::warning(this,QStringLiteral("错误"),QStringLiteral("请重新打开新的图片重置显示框"),QMessageBox::Abort);
		if(ret == QMessageBox::Abort)
			return false;
	}
    initframepixelmat();
    framepixelmat = scaling(framepixelmat, winimgpixelmat, winimgheight, winimgwidth);
    if(framepixelmat == NULL)
    {
        return false;
    }
	geometricflag = false;
    framepixelmattoshow();
    return true;
}


bool MainWindow::on_slotDFTBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
    winimgpixelmat = DFT(winimgpixelmat, winimgheight, winimgwidth);
    if(winimgpixelmatshow() == false)
    {
        return false;
    }
    return true;
}


bool MainWindow::on_slotDCTBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
    winimgpixelmat = DCT(winimgpixelmat, winimgheight, winimgwidth);
    if(winimgpixelmatshow() == false)
    {
        return false;
    }
    return true;
}


bool MainWindow::on_slotWalshBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
    winimgpixelmat = walsh(winimgpixelmat, winimgheight, winimgwidth);
    if(winimgpixelmatshow() == false)
    {
        return false;
    }
    return true;
}


bool MainWindow::on_slotHaarBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
    winimgpixelmat = haar(winimgpixelmat, winimgheight, winimgwidth);
    if(winimgpixelmatshow() == false)
    {
        return false;
    }
    return true;
}


bool MainWindow::on_slotMedianfitBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
    winimgpixelmat = medianfit(winimgpixelmat, winimgheight, winimgwidth);
    if(winimgpixelmatshow() == false)
    {
        return false;
    }
    return true;
}


bool MainWindow::on_slotAveragefitBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
    winimgpixelmat = averagefit(winimgpixelmat, winimgheight, winimgwidth);
    if(winimgpixelmatshow() == false)
    {
        return false;
    }
    return true;
}

bool MainWindow::on_slotLowpassfitBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
    winimgpixelmat = lowpassfit(winimgpixelmat, winimgheight, winimgwidth);
    if(winimgpixelmatshow() == false)
    {
        return false;
    }
    return true;
}


bool MainWindow::on_slotSobelBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
    winimgpixelmat = sobel(winimgpixelmat, winimgheight, winimgwidth);
    if(winimgpixelmatshow() == false)
    {
        return false;
    }
    return true;
}


bool MainWindow::on_slotLaplaceBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
    winimgpixelmat = laplace(winimgpixelmat, winimgheight, winimgwidth);
    if(winimgpixelmatshow() == false)
    {
        return false;
    }
    return true;
}

bool MainWindow::on_slotHighpassfitBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
    winimgpixelmat = highpassfit(winimgpixelmat, winimgheight, winimgwidth);
    if(winimgpixelmatshow() == false)
    {
        return false;
    }
    return true;
}

bool MainWindow::on_slotRandomnoiseBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
    winimgpixelmat = randomnoise(winimgpixelmat, winimgheight, winimgwidth);
    if(winimgpixelmatshow() == false)
    {
        return false;
    }
    return true;
}


bool MainWindow::on_slotImpulsenoiseBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
    winimgpixelmat = impulsenoise(winimgpixelmat, winimgheight, winimgwidth);
    if(winimgpixelmatshow() == false)
    {
        return false;
    }
    return true;
}


bool MainWindow::on_slotInversefitBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
    winimgpixelmat = inversefit(winimgpixelmat, winimgheight, winimgwidth);
    if(winimgpixelmatshow() == false)
    {
        return false;
    }
    return true;
}


bool MainWindow::on_slotWienerfitBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
    winimgpixelmat = wienerfit(winimgpixelmat, winimgheight, winimgwidth);
    if(winimgpixelmatshow() == false)
    {
        return false;
    }
    return true;
}


bool MainWindow::on_slotJpegBtn_triggered()
{
    if(winimgtopixelmat() == false)
    {
        return false;
    }
    winimgpixelmat = jpeg(winimgpixelmat, winimgheight, winimgwidth);
    if(winimgpixelmatshow() == false)
    {
        return false;
    }
    QMessageBox::about(this, QStringLiteral(" "), QStringLiteral("压缩成功并解压"));
    return true;
}


bool MainWindow::winimgtopixelmat()
{
    if(winimg->isNull())
    {
        return false;
    }
    uchar* winimgpt;
    winimgheight = winimg->height();
    winimgwidth = winimg->width();
    winimgpixelmat = new int*[winimgheight];
    for(int i = 0; i < winimgheight; i++)
    {
        winimgpixelmat[i] = new int[winimgwidth];
        winimgpt = winimg->scanLine(i);
        for(int j = 0; j < winimgwidth; j++)
        {
            winimgpixelmat[i][j] = int(winimgpt[j]);
        }
    }
    return true;
}

void MainWindow::pixelmattowinimg()
{
    winimg = new QImage(winimgheight, winimgwidth, QImage::Format_Indexed8);
    int* dataptr;
    uchar* imgptr;
    for(int i = 0; i < winimgheight; i++)
    {
        dataptr = winimgpixelmat[i];
        imgptr = winimg->scanLine(i);
        for(int j = 0; j < winimgwidth; j++)
        {
            imgptr[j] = uchar(dataptr[j]);
        }
    }
}


void MainWindow::initframepixelmat()
{
    framepixelmat = new int*[FRAME_WIDTH];
    for(int i = 0; i < FRAME_WIDTH; i++)
    {
        framepixelmat[i] = new int[FRAME_WIDTH];
        for(int j = 0; j < FRAME_WIDTH; j++)
        {
            framepixelmat[i][j] = 255;
        }
    }
}

void MainWindow::framepixelmattoshow()
{
    int* dataptr;
    uchar* winimgpt;
    winimg = new QImage(FRAME_WIDTH, FRAME_WIDTH, QImage::Format_Indexed8);
    for(int i = 0; i < FRAME_WIDTH; i++)
    {
        dataptr = framepixelmat[i];
        winimgpt = winimg->scanLine(i);
        for(int j = 0; j < FRAME_WIDTH; j++)
        {
            winimgpt[j] = uchar(dataptr[j]);
        }
    }
    winimgshow();
}

bool MainWindow::winimgshow()
{
    if(winimg->isNull())
    {
        return false;
    }
    QGraphicsScene *scene = new QGraphicsScene;
    scene->addPixmap(QPixmap::fromImage(*winimg));
    ui->graphicsView->setAlignment(Qt::AlignLeft | Qt::AlignTop);
    ui->graphicsView->setScene(scene);
    ui->graphicsView->show();

    return true;
}

bool MainWindow::winimgpixelmatshow()
{
    if(winimgpixelmat == NULL)
    {
        return false;
    }
    pixelmattowinimg();
    winimgshow();
    return true;
}


