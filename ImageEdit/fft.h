
/*PB15000119 赵汉卿  数字信号处理实验二 快速傅里叶变换*/
#pragma once
#include<complex>
#include<cmath>
#include<functional>
#include<cstdlib>
#include<atomic>
#include<thread>
#include<vector>
#include<Windows.h>
using namespace std;
#define PI 3.14159265
const int workers = 8;//threads must be 2 4 8 16 or 32
const int loga = 3;  // loga=log2(workers)
static const unsigned char BitReverseTable256[] =
{
	0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0, 0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0,
	0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8, 0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8,
	0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4, 0x64, 0xE4, 0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4,
	0x0C, 0x8C, 0x4C, 0xCC, 0x2C, 0xAC, 0x6C, 0xEC, 0x1C, 0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC,
	0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2, 0x12, 0x92, 0x52, 0xD2, 0x32, 0xB2, 0x72, 0xF2,
	0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA, 0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A, 0xFA,
	0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6, 0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6,
	0x0E, 0x8E, 0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE, 0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE,
	0x01, 0x81, 0x41, 0xC1, 0x21, 0xA1, 0x61, 0xE1, 0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1,
	0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9, 0x19, 0x99, 0x59, 0xD9, 0x39, 0xB9, 0x79, 0xF9,
	0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5, 0x15, 0x95, 0x55, 0xD5, 0x35, 0xB5, 0x75, 0xF5,
	0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD, 0x6D, 0xED, 0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD, 0x7D, 0xFD,
	0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3, 0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3,
	0x0B, 0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB, 0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB,
	0x07, 0x87, 0x47, 0xC7, 0x27, 0xA7, 0x67, 0xE7, 0x17, 0x97, 0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7,
	0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F, 0xEF, 0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF
};

////////////////////////////////////////表初始化
static unsigned int brt0[256], brt1[256], brt2[256], brt3[256];
static bool inited = false;
void fftinit() {
	if (inited == true)
		return;
	inited = true;
	for (int iii = 0; iii < 256; ++iii) {
		brt0[iii] = BitReverseTable256[iii];
		brt1[iii] = brt0[iii] << 8;
		brt2[iii] = brt1[iii] << 8;
		brt3[iii] = brt2[iii] << 8;
	}
}
/////////////////////////////////////////线程同步

class racer {
public:
	racer(int num) :c1(0), c2(0), lim(num) {}
	pair<int, int> race() {
		int ap1, ap2;
		while (c2 != 0)
			this_thread::yield();
		ap1 = ++c1;
		if (ap1 == lim) {
			ap2 = --c1;
			c2 = 1;
		}
		else {
			while (c2 != 1)
				this_thread::yield();
			ap2 = --c1;
			if (ap2 == 0)
				c2 = 0;
		}
		return make_pair(ap1 - 1, ap2);
	}
	atomic<int> c1, c2;
	const int lim;
};

//单线程位逆序函数

complex<float>* rev(complex<float>* vec, int N) {
	fftinit();
	complex<float>* pt = new complex<float>[1 << N];
	for (int32_t t = 0; t < 1 << N; ++t) {
		uint32_t c = 0;
		c += brt3[t & 0xFF];
		c += brt2[(t >> 8) & 0xFF];
		c += brt1[(t >> 16) & 0xFF];
		c += brt0[(t >> 24) & 0xFF];
		pt[t] = vec[c >> (32 - N)];
	}
	return pt;
}
///////////////////////////针对内存访问优化
void twork(int num, int N, complex<float>* vec, complex<float>* pt) {
	for (uint32_t t = num*(1 << (N - loga)); t <(num + 1)*(1 << (N - loga)); ++t) {
		uint32_t c = 0;
		c += brt3[t & 0xFF];
		c += brt2[(t >> 8) & 0xFF];
		c += brt1[(t >> 16) & 0xFF];
		c += brt0[(t >> 24) & 0xFF];
		pt[t] = vec[c >> (32 - N)];
	}
}
//////////////////////////////////////////


/////////////////////////多线程位逆序
complex<float>* rev2(complex<float>* vec, int N) {
	fftinit();
	complex<float>* pt = new complex<float>[1 << N];
	vector<thread> vct;
	for (int y = 0; y < workers; ++y) {
		vct.push_back(thread(bind(&twork, y, N, vec, pt)));
	}
	for (int y = 0; y < workers; ++y)
		vct[y].join();
	return pt;
}


template<class F>
complex<float>* cztfft(complex<float>* vec, int M, F func);
////////////////////////////单线程FFT
complex<float>* simplefft(complex<float>* vec_, int M) {
	int N = __popcnt(M - 1);
	if (1 << N != M)
		return cztfft(vec_, M, simplefft);
	if (N == 1) {
		complex<float> temp = vec_[0] - vec_[1];
		vec_[0] += vec_[1];
		vec_[1] = temp;
		return vec_;
	}
	complex<float> *k = new complex<float>[N];
	complex<float> ba = complex<float>(cos(2 * PI / M), -sin(2 * PI / M));
	k[0] = ba;
	for (int it = 1; it < N; ++it)
		k[it] = k[it - 1] * k[it - 1];
	complex<float>* vec = rev(vec_, N);
	for (int u = 1; u <= N; ++u) {
		int x = 1 << u, y = 1 << (N - u);
		complex<float>* ini = vec;
		complex<float> base0 = k[N - u];
		for (int j = 0; j < y; ++j, ini += x) {
			int d = 0;
			complex<float> base = 1;
			for (int b = x >> 1; b < x; ++b, ++d) {
				ini[b] *= base;
				base *= base0;
				complex<float> temp = ini[d] - ini[b];
				ini[d] += ini[b];
				ini[b] = temp;
			}
		}
	}
	delete[] k;
	return vec;
}


/////////////////////多线程FFT工作线程
void worker(int num, int N, complex<float>* vec, complex<float>* k, racer& rc) {
	for (int layer = 1; layer <= N; ++layer) {

		int batch_size = 1 << layer, batch_num = 1 << (N - layer);
		complex<float> base0 = k[N - layer];
		int batch_ever = batch_num >> loga;
		if (batch_ever == 0)
			batch_ever = 1;
		int stride = batch_size / 2;
		int batch_in = (1 << (N - 1 - loga)) / batch_ever;
		int b_offset = 0, b_num = 0;
		b_num = (num*batch_num) / workers;
		b_offset = ((num*batch_num) % workers)*batch_size / 2 / workers;
		complex<float>* ini = vec + batch_size*b_num + b_offset;
		complex<float> basex = pow(base0, b_offset);
		if (workers > batch_num) {
			rc.race();
		}
		for (int j = 0; j<batch_ever; ++j, ini += batch_size) {
			int d = 0, b = stride; complex<float> base = basex;
			for (; d < batch_in; ++b, ++d) {
				ini[b] *= base;
				base *= base0;
				complex<float> temp = ini[d] - ini[b];
				ini[d] += ini[b];
				ini[b] = temp;
			}
		}
	}
}


/////////////////////多线程FFT函数
complex<float>* fft2(complex<float>* vec_, int M) {
	if (M <= 64)
		return simplefft(vec_, M);
	int N = __popcnt(M - 1);
	if (1 << N != M)
		return cztfft(vec_, M, fft2);
	racer rc(workers);
	complex<float> *k = new complex<float>[N];
	complex<float> ba = complex<float>(cos(2 * PI / M), -sin(2 * PI / M));
	k[0] = ba;
	for (int it = 1; it < N; ++it)
		k[it] = k[it - 1] * k[it - 1];
	complex<float>* vec = rev2(vec_, N);
	vector<thread> vot;
	for (int th = 0; th < workers; ++th)
		vot.push_back(thread(bind(&worker, th, N, vec, k, ref(rc))));
	for (int th = 0; th < workers; ++th)
		vot[th].join();
	delete[] k;
	return vec;
}

template<class F>
complex<float>* ifft(complex<float>* vec, int M, F func) {
	complex<float>* fk = func(vec, M);
	float m = M;
	fk[0] /= m;
	complex<float> temp;
	for (int i = 1; i < m / 2 - 0.1; ++i) {
		temp = fk[i] / m;
		fk[i] = fk[M - i] / m;
		fk[M - i] = temp;
	}
	return fk;
}

template<class F>
complex<float>* cztfft(complex<float>* vec, int M, F func) {
	int L = 2;
	while (L < M)
		L *= 2;
	L *= 2;
	auto g = new complex<float>[L];
	auto h = new complex<float>[L];
	complex<float> G(cos(PI / M), sin(PI / M));
	complex<float> aux = 1;
	h[0] = 1;
	g[0] = vec[0];
	for (int i = 1; i < M; ++i) {
		h[i] = h[i - 1] * aux*aux*G;
		aux *= G;
		h[L - i] = h[i];
		g[i] = vec[i] / h[i];
	}
	auto hf = func(h, L);
	auto gf = func(g, L);
	for (int i = 0; i < L; ++i)
		hf[i] *= gf[i];
	auto rtf = ifft(hf, L, func);
	delete[] hf;
	delete[] gf;
	delete[] g;
	auto rt = new complex<float>[M];
	for (int i = 0; i < M; ++i)
		rt[i] = rtf[i] / (h[i]);
	delete[] h;
	delete[] rtf;
	return rt;
}

template<class T>
T** aloc(int mh, int mw) {
	auto ptr = new T*[mh];
	for (int i = 0; i < mh; ++i)
		ptr[i] = new T[mw];
	return ptr;
}

template<class T>
void del(T** ptr, int mh) {
	for (int i = 0; i < mh; ++i)
		delete[] ptr[i];
	delete[] ptr;
}

pair < complex<float>*, complex<float>* >  fft2in1(float* A, float* B, int M) {
	auto fpa = new complex<float>[M];
	auto fpb = new complex<float>[M];
	for (int i = 0; i < M; ++i)
		fpa[i] = complex<float>(A[i], B[i]);
	auto fpc = simplefft(fpa, M);
	for (int i = 1; i < M; ++i) {
		fpa[i] = 0.5f*(fpc[i] + conj(fpc[M - i]));
		fpb[i] = 0.5f*(fpc[i] - conj(fpc[M - i]))*complex<float>(0, -1);
	}
	fpa[0] = fpc[0].real();
	fpb[0] = fpc[0].imag();
	delete[] fpc;
	return make_pair(fpa, fpb);
}


float* simpledct(float* vec, int N) {
	float* ptr = new float[N];
	for (int i = 0; i < N; ++i) {
		ptr[i] = 0;
		for (int j = 0; j < N; ++j)
			ptr[i] += vec[j] * cos((2 * j + 1)*i*PI / 2 / N);
	}
	ptr[0] *= sqrt(1.0 / N);
	for (int i = 1; i < N; ++i)
		ptr[i] *= sqrt(2.0 / N);
	return ptr;
}



const int REUSE_MAT = 1, REUSE_ROW = 2, REUSE_COL = 4, SKIP_FI = 8, SKIP_FO = 16, SKIP_COL = 32, SKIP_ROW = 64, NEW_MAT = 128;
template<class T, class T2>
T DEFAULT_FI(T2** a, int b, int c, int d, int e) {
	return T(a[b][c]);
}
template<class T, class T3>
T3 DEFAULT_FO(T** a, int b, int c, int d, int e) {
	return a[b][c];
}
template<class T>
T* DEFAULT_FR(T* a, int M, int N) {
	auto x = new T[M];
	for (int i = 0; i < M; ++i)
		x[i] = a[i];
	return x;
}
template<class T>
T* DEFAULT_FC(T* a, int M, int N) {
	suto x = new T[M];
	for (int i = 0; i < M; ++i)
		x[i] = a[i];
	return x;
}

template<class T, class T2 = T, class T3 = T>
class trans2d {
public:
	int para;
	function<T*(T*, int, int)> fr, fc;
	function<T3(T**, int, int, int, int)> fi, fo;
	trans2d(int _para, function<T*(T*, int, int)> _fr = DEFAULT_FR<T>, function<T*(T*, int, int)> _fc = DEFAULT_FC<T>, \
		function<T3(T**, int, int, int, int)> _fo = DEFAULT_FO<T, T3>, function<T(T2**, int, int, int, int)> _fi = DEFAULT_FI<T, T2>)\
		:para(_para), fo(_fo), fi(_fi), fr(_fr), fc(_fc){}
	T3** operator()(T2** vec_, int mh, int mw) {
		racer rc(workers);
		bool reuse_mat = para & REUSE_MAT, reuse_row = para & REUSE_ROW, reuse_col = para & REUSE_COL, new_mat = para&NEW_MAT;
		bool skip_fi = (para & SKIP_FI) && reuse_mat, skip_row = para & SKIP_ROW, skip_col = para & SKIP_COL, skip_fo = (para&SKIP_FO) && (!new_mat);
		T** vec; T3** nmt;
		if (!reuse_mat)
			vec = aloc<T >(mh, mw);
		else
			vec = (T**)vec_;
		if (new_mat)
			nmt = aloc<T3>(mh, mw);
		else
			nmt = (T3**)vec;
		int mha = mh / workers + 1, mwa = mw / workers + 1;
		vector<thread> vot;
		for (int i = 0; i < workers; ++i) {
			vot.push_back(thread([=, &rc] {
				int begin = i*mha, end = mh>(i + 1)*mha ? (i + 1)*mha : mh;
				/// DO FI
				if (!skip_fi) {
					for (int j = begin; j < end; ++j)
						for (int k = 0; k < mw; ++k)
							vec[j][k] = fi(vec_, j, k, mh, mw);
				}
				rc.race();
				////DO ROW
				if (!skip_row) {
					auto barr = new T[mw];
					for (int j = begin; j < end; j++) {
						for (int k = 0; k < mw; ++k)
							barr[k] = vec[j][k];
						auto tret = fr(barr, mw, j);
						for (int k = 0; k < mw; ++k)
							vec[j][k] = tret[k];
						if (!reuse_row)
							delete[] tret;
					}
					delete[] barr;
				}
				rc.race();
				////DO COL
				begin = i*mwa, end = mw>(i + 1)*mwa ? (i + 1)*mwa : mw;
				if (!skip_col) {
					auto barr = new T[mh];
					for (int j = begin; j < end; j++) {
						for (int k = 0; k < mh; ++k)
							barr[k] = vec[k][j];
						auto tret = fc(barr, mh, j);
						for (int k = 0; k < mh; ++k)
							vec[k][j] = tret[k];
						if (!reuse_col)
							delete[] tret;
					}
					delete[] barr;
				}
				rc.race();
				///DO FO
				if (!skip_fo) {
					for (int k = 0; k < mh; ++k)
						for (int j = begin; j < end; ++j)
							nmt[k][j] = fo(vec, k, j, mh, mw);
				}
				return;
			}));
		}
		for (int i = 0; i < workers; ++i)
			vot[i].join();
		if (new_mat && (!reuse_mat))
			del(vec, mh);
		return nmt;
	}
};


trans2d<float> dct2d(REUSE_MAT | SKIP_FI | SKIP_FO, [](float* a1, int a2, int a3) {return simpledct(a1, a2); }, [](float* a1, int a2, int a3) {return simpledct(a1, a2); });
trans2d<complex<float> > fft2d(REUSE_MAT | SKIP_FI, [](complex<float>* a1, int a2, int a3) {return simplefft(a1, a2); }, \
	[](complex<float>* a1, int a2, int a3) {return simplefft(a1, a2); }, [](complex<float>** a1, int a2, int a3, int a4, int a5) {return a1[a2][a3] / sqrtf(a4*a5); });

trans2d<complex<float> > ifft2d(REUSE_MAT | NEW_MAT | SKIP_FI, [](complex<float>* a1, int a2, int a3) {return simplefft(a1, a2); }, \
	[](complex<float>* a1, int a2, int a3) {return simplefft(a1, a2); }, [](complex<float>** a1, int a2, int a3, int a4, int a5) {return a1[(a4 - a2) % a4][(a5 - a3) % a5] / sqrtf(a4*a5); });

