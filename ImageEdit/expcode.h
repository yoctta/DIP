#ifndef EXPCODE_H
#define EXPCODE_H

#define FRAME_HEIGHT 800
#define FRAME_WIDTH 800

int midptvalue(int** pixelmat, int mheight, int mwidth);
int maxvalue(int** pixelmat, int mheight, int mwidth);
int minvalue(int** pixelmat, int mheight, int mwidth);
float avgvalue(int** pixelmat, int mheight, int mwidth);
float varvalue(int** pixelmat, int mheight, int mwidth);
int* histogram(int** pixelmat, int mheight, int mwidth);
int** binaryimg(int** pixelmat, int mheight, int mwidth);
int** histogramequ(int** pixelmat, int mheight, int mwidth);
int** graystretch(int** pixelmat, int mheight, int mwidth);
int** centralize(int** framemat, int** pixelmat, int mheight, int mwidth);
int** rotation(int** framemat, int** pixelmat, int mheight, int mwidth);
int** moveimage(int** framemat, int** pixelmat, int mheight, int mwidth);
int** scaling(int** framemat, int** pixelmat, int mheight, int mwidth);
int** DFT(int** pixelmat, int mheight, int mwidth);
int** DCT(int** pixelmat, int mheight, int mwidth);
int** walsh(int** pixelmat, int mheight, int mwidth);
int** haar(int** pixelmat, int mheight, int mwidth);
int** medianfit(int** pixelmat, int mheight, int mwidth);
int** averagefit(int** pixelmat, int mheight, int mwidth);
int** lowpassfit(int** pixelmat, int mheight, int mwidth);
int** sobel(int** pixelmat, int mheight, int mwidth);
int** laplace(int** pixelmat, int mheight, int mwidth);
int** highpassfit(int** pixelmat, int mheight, int mwidth);
int** randomnoise(int** pixelmat, int mheight, int mwidth);
int** impulsenoise(int** pixelmat, int mheight, int mwidth);
int** inversefit(int** pixelmat, int mheight, int mwidth);
int** wienerfit(int** pixelmat, int mheight, int mwidth);
int** jpeg(int** pixelmat, int mheight, int mwidth);


#endif // EXPCODE_H
