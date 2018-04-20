#include "expcode.h"
#include <iostream>
#include<algorithm>
#include<numeric>
#include<vector>
#include<random>
#include<fft.h>
using namespace std;

//示例: 求图像中心点像素的灰度值
int midptvalue(int** pixelmat, int mheight, int mwidth)
{
	//pixelmat为指向图像像素值二维数组指针, 指向图像左上角的像素, 像素值类型为int;
	//mheight为图像的高度(行), mwidth为图像的宽度(列);
	int middlerow = mheight / 2 - 1;
	int middlecol = mwidth / 2 - 1;
	return pixelmat[middlerow][middlecol];
}

//求最大值
int maxvalue(int** pixelmat, int mheight, int mwidth)
{
	int x = pixelmat[0][0];
	for (int a = 0; a < mheight; ++a)
		for (int b = 0; b < mwidth; ++b)
			if (pixelmat[a][b] > x)
				x = pixelmat[a][b];
    return x;
}

//求最小值
int minvalue(int** pixelmat, int mheight, int mwidth)
{
	int x = pixelmat[0][0];
	for (int a = 0; a < mheight; ++a)
		for (int b = 0; b < mwidth; ++b)
			if (pixelmat[a][b] < x)
				x = pixelmat[a][b];
	return x;
}

//求平均值
float avgvalue(int** pixelmat, int mheight, int mwidth)
{
	float x = 0;
	for (int a = 0; a < mheight; ++a)
		for (int b = 0; b < mwidth; ++b)
			x += pixelmat[a][b];
	return x/mheight/mwidth;
}

//求方差
float varvalue(int** pixelmat, int mheight, int mwidth)
{
	float aver = avgvalue(pixelmat, mheight, mwidth);
	float x = 0;
	for (int a = 0; a < mheight; ++a)
		for (int b = 0; b < mwidth; ++b)
			x += (pixelmat[a][b]-aver)*(pixelmat[a][b] - aver);
	return x / mheight / mwidth;
}

//统计直方图, 返回长度为256的1维数组
int* histogram(int** pixelmat, int mheight, int mwidth)
{
	int* arr = new int[256]();
	for (int a = 0; a < mheight; ++a)
		for (int b = 0; b < mwidth; ++b)
			++arr[pixelmat[a][b]];
	//注意:函数内分配数组必须使用动态分配;
	return arr;
}

//示例,将灰度图转化为二值图像,返回处理后的图像
int** binaryimg(int** pixelmat, int mheight, int mwidth)
{
	for(int i = 0; i < mheight; i++)
	{
		for(int j = 0; j < mwidth; j++)
		{
			//从左上角开始遍历整幅图像, 实现二值化;

			pixelmat[i][j] = pixelmat[i][j] > 128 ? 255 : 0;
		}
	}
	//根据实验要求返回对应的值;
	return pixelmat;
}

//直方图均衡, 返回处理后的图像
int** histogramequ(int** pixelmat, int mheight, int mwidth)
{
	int* arr = histogram(pixelmat, mheight, mwidth);
	int* brr = new int[256]();
	for (int i = 0; i < 256; ++i) {
		arr[i] += i ? arr[i - 1] : 0;
		brr[i] = 255 * arr[i] / mheight / mwidth;
	}
	for (int a = 0; a < mheight; ++a)
		for (int b = 0; b < mwidth; ++b)
			pixelmat[a][b] = brr[pixelmat[a][b]];
    return pixelmat;
}

//灰度拉伸, 返回处理后的图像
int** graystretch(int** pixelmat, int mheight, int mwidth)
{
	int min = minvalue(pixelmat, mheight, mwidth), max = maxvalue(pixelmat, mheight, mwidth);
	for (int i = 0; i < mheight; ++i)
		for (int j = 0; j < mwidth; ++j)
			pixelmat[i][j] = (pixelmat[i][j] - min) * 255 / (max - min);
    return pixelmat;
}

template<class T>
int** somthing(int** src, int mh, int mw,T func) {
	int** a = new int*[mh];
	for (auto i = 0; i < mh;++i) {
		a[i] = new int[mw];
	}
	for(auto i=1;i<mh-1;++i)
		for (auto j = 1; j < mw - 1; ++j) {
			int range[3][3];
			for (auto ii = 0; ii < 3; ++ii)
				for (auto jj = 0; jj < 3; ++jj)
					range[ii][jj] = src[i + ii - 1][j + jj - 1];
			a[i - 1][j - 1] = func(range);
		}
	return a;
}
//中值滤波, 返回处理后的图像
int** medianfit(int** pixelmat, int mheight, int mwidth)
{
	auto func = [](int(&x)[3][3]) {std::sort(*x, *x + 9); return x[1][1]; };
    return somthing(pixelmat,mheight,mwidth,func);
}

//均值滤波, 返回处理后的图像
int** averagefit(int** pixelmat, int mheight, int mwidth)
{
	auto func = [](int(&x)[3][3]) {return std::accumulate(*x, *x + 9, 0) / 9; };
	return somthing(pixelmat, mheight, mwidth, func);
}


//理想低通滤波, 返回处理后的图像
int** lowpassfit(int** pixelmat, int mheight, int mwidth)
{
	auto pma = aloc<complex<float> >(mheight, mwidth);
	for (int i = 0; i < mheight; ++i)
		for (int j = 0; j < mwidth; ++j)
			pma[i][j] = complex<float>(pixelmat[i][j]);
	auto ptr = fft2d(pma, mheight, mwidth);
	for(int i=0;i<mheight;++i)
		for (int j = 0; j < mwidth; ++j)
		{
			int x = j + 1 > mwidth - j ? mwidth - j : j + 1;
			int y = i + 1 > mheight - i ? mheight - i : i + 1;
			if (x*x + y*y > 25 * 25)
				ptr[i][j] = 0;
		}
	ptr = ifft2d(ptr, mheight, mwidth);
	for (int i = 0; i < mheight; ++i)
		for (int j = 0; j < mwidth; ++j)
			pixelmat[i][j] = abs(ptr[i][j]);
	del(ptr, mheight);
    return pixelmat;
}

//sobel算子, 返回处理后的图像
int** sobel(int** pixelmat, int mheight, int mwidth)
{
	auto sobel1 = [](int(&x)[3][3]) {return(-x[0][0] - 2 * x[0][1] - x[0][2] + x[2][0] + 2*x[2][1] + x[2][2]); };
	auto sobel2 = [](int(&x)[3][3]) {return(-x[0][0] - 2 * x[1][0] - x[2][0] + x[0][2] + 2*x[1][2] + x[2][2]); };
	auto sobel = [sobel1,sobel2](int(&x)[3][3]) {return (sobel1(x) + sobel2(x))/4; };
    return somthing(pixelmat, mheight, mwidth, sobel);
}

//laplace算子, 返回处理后的图像
int** laplace(int** pixelmat, int mheight, int mwidth)
{
	auto func = [](int(&x)[3][3]) {return(4 * x[1][1] - x[0][1] - x[1][0] - x[1][2] - x[2][1]); };
	return somthing(pixelmat, mheight, mwidth, func);
}

//理想高通滤波, 返回处理后的图像
int** highpassfit(int** pixelmat, int mheight, int mwidth)
{
	auto pma = aloc<complex<float> >(mheight, mwidth);
	for (int i = 0; i < mheight; ++i)
		for (int j = 0; j < mwidth; ++j)
			pma[i][j] = complex<float>(pixelmat[i][j]);
	auto ptr = fft2d(pma, mheight, mwidth);
	for (int i = 0; i<mheight; ++i)
		for (int j = 0; j < mwidth; ++j)
		{
			int x = j + 1 > mwidth - j ? mwidth - j : j + 1;
			int y = i + 1 > mheight - i ? mheight - i : i + 1;
			if (x*x + y*y < 10 * 10)
				ptr[i][j] = 0;
		}
	ptr = ifft2d(ptr, mheight, mwidth);
	for (int i = 0; i < mheight; ++i)
		for (int j = 0; j < mwidth; ++j)
			pixelmat[i][j] = abs(ptr[i][j]);
	del(ptr, mheight);
	return pixelmat;
}

//示例, 将图像平移到显示区域的中心
int** centralize(int** framemat, int** pixelmat, int mheight, int mwidth)
{
	//framemat为指向显示区域(画板)的二维数组指针, 大小为FRAME_HEIGHT x FRAMEWIDTH = 800 x 800
	int xpt = (FRAME_HEIGHT - mheight) / 2;
	int ypt = (FRAME_WIDTH - mwidth) / 2;
	for(int i = 0; i < mheight; i++)
	{
		for(int j = 0; j < mwidth; j++)
		{
			framemat[i + xpt][j + ypt] = pixelmat[i][j];
		}
	}
	return framemat;
}
template<class F>
void insert(int** frame, int** pix, int mh, int mw, F func) {
	for(int i=0;i<FRAME_HEIGHT;++i)
		for (int j = 0; j < FRAME_WIDTH; ++j)
		{
			pair<float, float> p = func(i, j);
			float x = p.first, y = p.second;
			if (x<0 || x>=mh-1 || y<0 || y>=mw-1)
				frame[i][j] = 0;
			else {
				int a = x, b = y;
				int c = a + 1, d = b + 1;
				float d1 = x - a, d2 = c - x, d3 = y - b, d4 = d - y;
				frame[i][j] = int(d1*d3*pix[c][d] + d1*d4*pix[c][b] + d2*d3*pix[a][d] + d2*d4*pix[a][b]);
			}
		}
}
//旋转图像, 返回显示区域(画板)指针
int** rotation(int** framemat, int** pixelmat, int mheight, int mwidth)
{
	float angel = PI / 6;
	insert(framemat, pixelmat, mheight, mwidth, [=](int a, int b) {return make_pair(cos(angel)*(a - (FRAME_HEIGHT-mheight) / 2) - sin(angel)*(b - (FRAME_WIDTH-mwidth) / 2),sin(angel)*(a - (FRAME_HEIGHT-mheight) / 2) + cos(angel)*(b - (FRAME_WIDTH-mwidth) / 2)); });
	return framemat;
}

//平移图像, 返回显示区域(画板)指针
int** moveimage(int** framemat, int** pixelmat, int mheight, int mwidth)
{
	int dx = 100, dy = 200;
	insert(framemat, pixelmat, mheight, mwidth, [=](int a, int b) {return make_pair(float(a-dx),float(b-dy)); });
	return framemat;
}

//缩放图像, 返回显示区域(画板)指针
int** scaling(int** framemat, int** pixelmat, int mheight, int mwidth)
{
	float rate = 3;
	insert(framemat, pixelmat, mheight, mwidth, [=](int a, int b) {return make_pair(float(a)/rate,float(b)/rate); });
	return framemat;
}

//DFT变换, 返回处理后的图像, 注意缩放到0~255的整型
int** DFT(int** pixelmat, int mheight, int mwidth)
{
	auto pma = aloc<complex<float> >(mheight, mwidth);
	for (int i = 0; i < mheight; ++i)
		for (int j = 0; j < mwidth; ++j)
			pma[i][j] = complex<float>(pixelmat[i][j]);
	pma = fft2d(pma, mheight, mwidth);
	for (int i = 0; i < mheight; ++i)
		for (int j = 0; j < mwidth; ++j) {
			pixelmat[i][j] = abs(pma[i][j]);
			if (pixelmat[i][j] > 255)
				pixelmat[i][j] = 255;
		}

	del(pma, mheight);
	return pixelmat;
}


//DCT变换, 返回处理后的图像
int** DCT(int** pixelmat, int mheight, int mwidth)
{
	auto ptr = aloc<float>(mheight, mwidth);
	for (int i = 0; i < mheight; ++i)
		for (int j = 0; j < mwidth; ++j)
			ptr[i][j] = float(pixelmat[i][j]);
	auto ptr2 = dct2d(ptr, mheight, mwidth);
	for (int i = 0; i < mheight; ++i)
		for (int j = 0; j < mwidth; ++j)
		{
			pixelmat[i][j] = abs(ptr2[i][j]);
			if (pixelmat[i][j] > 255)
				pixelmat[i][j] = 255;
		}
    return pixelmat;
}


inline int g(int x, int u) {
	return __popcnt(x&BitReverseTable256[u])&1?-1:1;
}

int* walsh1d(int* vec) {
	auto ptr = new int[256];
	for (int i = 0; i < 256; ++i) {
		ptr[i] = 0;
		for (int j = 0; j < 256; ++j)
			ptr[i] += vec[j] * g(j, i);
	}
	return ptr;
}


trans2d<int> walsh2d(REUSE_MAT|SKIP_FI, [](int* a1, int a2, int a3) {return walsh1d(a1); }, [](int* a1, int a2, int a3) {return walsh1d(a1); },
	[](int** a1, int a2, int a3, int a4, int a5) {return a1[a2][a3] / 256; });

//walsh变换, 返回处理后的图像
int** walsh(int** pixelmat, int mheight, int mwidth)
{
	int** xp = walsh2d(pixelmat, mheight, mwidth);
	for(int i=0;i<mheight;++i)
		for (int j = 0; j < mwidth; ++j)
		{
			if (xp[i][j] > 255)
				xp[i][j] = 255;
			if (xp[i][j] < 0)
				xp[i][j] = 0;
		}
    return xp;
}

float* haar_core(float* vec, int N) {
	if (N == 1)
		return vec;
	int mid = N / 2;
	auto ptr = new float[N];
	for (int i = 0; i < N; ++i)
		ptr[i] = vec[i];
	for (int i = 0; i < mid; ++i) {
		vec[i] = ptr[2 * i] + ptr[2 * i + 1];
		vec[mid + i] = ptr[2 * i] - ptr[2 * i + 1];
	}
	delete[] ptr;
	haar_core(vec, mid);
	return vec;
}

float* haar1d(float* vec, int N) {
	auto ptr = new float[N];
	haar_core(vec, N);
	int M = __popcnt(N - 1);
	int reg = 1;
	for (int i = 0; i < N; ++i) {
		if (i >= (1 << reg))
			++reg;
		ptr[i] = vec[i]*pow(sqrt(2), reg - 1)/sqrt(N);
	}
	return ptr;
}


trans2d<float> haar2d(REUSE_MAT | SKIP_FI | SKIP_FO, [](float* a1, int a2, int a3) {return haar1d(a1, a2); },\
	[](float* a1, int a2, int a3) {return haar1d(a1, a2); });

//haar变换, 返回处理后的图像
int** haar(int** pixelmat, int mheight, int mwidth)
{
	auto ptr = aloc<float>(mheight, mwidth);
	for (int i = 0; i < mheight; ++i)
		for (int j = 0; j < mwidth; ++j)
			ptr[i][j] = pixelmat[i][j];
	ptr = haar2d(ptr, mheight, mwidth);
	for (int i = 0; i < mheight; ++i)
		for (int j = 0; j < mwidth; ++j)
			pixelmat[i][j] = ptr[i][j];
    return pixelmat;
}

//生成随机噪声, 返回处理后的图像;
int** randomnoise(int** pixelmat, int mheight, int mwidth)
{
	random_device r;
	default_random_engine e1(r());
	std::normal_distribution<> normal_dist(0 , 10);
	for(int i=0;i<mheight;++i)
		for (int j = 0; j < mwidth; ++j)
		{
			pixelmat[i][j] += round(normal_dist(e1));
			if (pixelmat[i][j] > 255)
				pixelmat[i][j] = 255;
			if (pixelmat[i][j] < 0)
				pixelmat[i][j] = 0;
		}
    return pixelmat;
}

//生成椒盐噪声, 返回处理后的图像
int** impulsenoise(int** pixelmat, int mheight, int mwidth)
{
	random_device r;
	default_random_engine e1(r());
	std::uniform_int_distribution<int> u1(0,mheight-1);
	std::uniform_int_distribution<int> u2(0, mwidth-1);
	for (int i = 0; i < 50; ++i) {
		pixelmat[u1(e1)][u2(e1)] = 255;
	}
	for (int i = 0; i < 50; ++i) {
		pixelmat[u1(e1)][u2(e1)] = 0;
	}
    return pixelmat;
}

//逆滤波复原
int** inversefit(int** pixelmat, int mheight, int mwidth)
{
    return NULL;
}

//维纳滤波
int** wienerfit(int** pixelmat, int mheight, int mwidth)
{
    return NULL;
}


//示例: JPEG压缩及解压缩
int** jpeg(int** pixelmat, int mheight, int mwidth)
{
    return NULL;
}



