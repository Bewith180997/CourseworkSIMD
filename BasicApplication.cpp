// BasicOpenCLApplication.cpp : Defines the entry point for the console application.
// : Results at the bottom of the document
// : Bewith180997



#include "stdafx.h"
#include "Chrono.h"
#include <math.h>
#include <thread>
#include <immintrin.h>
#include <mutex>

#define NB_THREADS 2

#define float8 __m256
#define int8 __m256i
#define mul8(a,b) _mm256_mul_ps(a,b)
#define add8(a,b) _mm256_add_ps(a,b)
#define div8(a,b) _mm256_div_ps(a,b)
#define set8(a) _mm256_set1_ps(a)
#define set8i(a) _mm256_set1_epi32(a)
#define sub8(a,b) _mm256_sub_ps(a,b)
#define min8(a,b) _mm256_min_ps(a,b)
#define cmp8(a,b) _mm256_cmp_ps(a, b, _CMP_LT_OS)
#define fToI(a) _mm256_cvtps_epi32(a)
#define iToF(a) _mm256_cvtepi32_ps(a)
#define compareGT(a,b) _mm256_cmpgt_epi32(a, b)
#define combine(a, b, mask) _mm256_blendv_ps(a, b, mask)


volatile int currentHeight = -1;
std::thread thr[NB_THREADS];
const int max_iterations_num = 255;

void SaveBMP(char* fname, unsigned char* image, int width, int height, int componentPerPixel = 1, int reverseColor = 0)
{
	FILE* destination;
	int i, j;
	int* pt;
	char name[512], hdr[0x36];
	unsigned char* imsource = new unsigned char[width * height * 3];
	//int al=(ImageSize*3)%4;

	if (componentPerPixel == 1)
		for (i = 0; i < width * height * 3; i++)
			imsource[i] = image[i / 3];
	else
		for (i = 0; i < width * height * 3; i++)
			imsource[i] = image[i];
	if (reverseColor)
		for (j = 0; j < height; j++)
			for (i = 0; i < width; i++)
			{
				unsigned char aux;
				aux = imsource[3 * (i + width * j)];
				imsource[3 * (i + width * j)] = imsource[3 * (i + width * j) + 2];
				imsource[3 * (i + width * j) + 2] = aux;
			}
	strcpy(name, fname);
	i = (int)strlen(name);
	if (!((i > 4) && (name[i - 4] == '.') && (name[i - 3] == 'b') && (name[i - 2] == 'm') && (name[i - 1] == 'p')))
	{
		name[i] = '.';
		name[i + 1] = 'b';
		name[i + 2] = 'm';
		name[i + 3] = 'p';
		name[i + 4] = 0;
	}
	if ((destination = fopen(name, "wb")) == NULL)
		perror("erreur de creation de fichier\n");
	hdr[0] = 'B';
	hdr[1] = 'M';
	pt = (int*)(hdr + 2);// file size
	*pt = 0x36 + width * height * 3;
	pt = (int*)(hdr + 6);//reserved
	*pt = 0x0;
	pt = (int*)(hdr + 10);// image address
	*pt = 0x36;
	pt = (int*)(hdr + 14);// size of [0E-35]
	*pt = 0x28;
	pt = (int*)(hdr + 0x12);// Image width
	*pt = width;
	pt = (int*)(hdr + 0x16);// Image heigth
	*pt = height;
	pt = (int*)(hdr + 0x1a);// color planes
	*pt = 1;
	pt = (int*)(hdr + 0x1c);// bit per pixel
	*pt = 24;
	for (i = 0x1E; i < 0x36; i++)
		hdr[i] = 0;
	fwrite(hdr, 0x36, 1, destination);
	fwrite(imsource, width * height * 3, 1, destination);
	fclose(destination);
	delete[] imsource;
}

typedef struct { float real; float im; } complex;
typedef struct { float real[8]; float im[8]; } complex8;

complex add(complex a, complex b)
{
	complex res;
	res.real = a.real + b.real;
	res.im = a.im + b.im;
	return res;
}
complex sub(complex a, complex b)
{
	complex res;
	res.real = a.real - b.real;
	res.im = a.im - b.im;
	return res;
}
complex mul(complex a, complex b)
{
	complex res;
	res.real = a.real * b.real - a.im * b.im;
	res.im = a.real * b.im + a.im * b.real;
	return res;
}

complex mulRef(float8* aReal, float8* bReal, float8* aIm, float8* bIm)
{
	complex res;
	res.real = (aReal->m256_f32)[0] * (bReal->m256_f32)[0] - (aIm->m256_f32)[0] * (bIm->m256_f32)[0];
	res.im = (aReal->m256_f32)[0] * (bIm->m256_f32)[0] + (aIm->m256_f32)[0] * (bReal->m256_f32)[0];
	return res;
}

void subC8(float8* aReal, float8* bReal, float8* aIm, float8* bIm)
{
	((float8*)bReal->m256_f32)[0] = ((float8*)sub8(*aReal, *bReal).m256_f32)[0];
	((float8*)bIm->m256_f32)[0] = ((float8*)sub8(*aIm, *bIm).m256_f32)[0];
}

void mulC8(float8* aReal, float8* bReal, float8* aIm, float8* bIm)
{
	float8* resReal = ((float8*)sub8(mul8(*aIm, *bIm), mul8(*aReal, *bReal)).m256_f32);
	float8* resIm = ((float8*)add8(mul8(*aReal, *bIm), mul8(*aIm, *bReal)).m256_f32);

	((float8*)bReal->m256_f32)[0] = ((float8*)resReal->m256_f32)[0];
	((float8*)bIm->m256_f32)[0] = ((float8*)resIm->m256_f32)[0];
}

float squaredNorm(complex c)
{
	return c.real * c.real + c.im * c.im;
}

float8 squaredNorm8(float8* cReal, float8* cIm)
{

	float8 squared = set8(0);
	squared = add8(mul8(*cReal, *cReal), mul8(*cIm, *cIm));
	return squared;
}


int Iterate(complex c)
{
	const int max_iterations = 255;
	complex z, a, l;
	a.real = 0.91;
	a.im = 0.;

	z = c;
	l.real = 4;
	l.im = 0;
	int i = 1;
	while (i < max_iterations)
	{

		z = mul(l, mul(z, sub(a, z)));
		if (squaredNorm(z) > 128)
			break;
		i += 2;
	}
	return (min(i, max_iterations));
}


int8 Iter8(float8* cPTRReal, float8* cPTRIm)
{
	
	const int8 max_iterations = set8i(max_iterations_num);
	complex8 z, a, l, b, d;

	//Instantiate all of the pointers
	float8* aPTRReal = (float8*)a.real;
	float8* aPTRIm = (float8*)a.im;
	float8* zPTRReal = (float8*)z.real;
	float8* zPTRIm = (float8*)z.im;
	float8* lPTRReal = (float8*)l.real;
	float8* lPTRIm = (float8*)l.im;
	float8* bPTRReal = (float8*)b.real;
	float8* bPTRIm = (float8*)b.im;
	float8* dPTRReal = (float8*)d.real;
	float8* dPTRIm = (float8*)d.im;

	//Set the values in each pointer
	((float8*)aPTRReal)[0] = set8(0.91f);
	((float8*)aPTRIm)[0] = set8(0.f);

	((float8*)zPTRReal) = cPTRReal;
	((float8*)zPTRIm) = cPTRIm;

	((float8*)lPTRReal)[0] = set8(4.f);
	((float8*)lPTRIm)[0] = set8(0.f);

	((float8*)bPTRReal) = bPTRReal;
	((float8*)bPTRIm) = bPTRIm;

	((float8*)dPTRReal) = dPTRReal;
	((float8*)dPTRIm) = dPTRIm;

	int8 i = set8i(1);
	int j = 1;

	float8 incrementArray = set8(0.f);
	float8 isAllSquared = set8(0.f);
	float8 addOn = set8(0.f);
	float8 leaveEarly = set8(0.f);
	float8 mask = set8(0.f);

	bool keepLooping = true;

	while (keepLooping)
	{
		//Reset to the appropriate values for the upcoming calculations
		((float8*)bPTRReal->m256_f32)[0] = ((float8*)zPTRReal->m256_f32)[0];
		((float8*)bPTRIm->m256_f32)[0] = ((float8*)zPTRIm->m256_f32)[0];
		((float8*)dPTRReal->m256_f32)[0] = ((float8*)zPTRReal->m256_f32)[0];
		((float8*)dPTRIm->m256_f32)[0] = ((float8*)zPTRIm->m256_f32)[0];

		//Perform the calculations, using pass by reference rather than pass by value
		subC8(aPTRReal, bPTRReal, aPTRIm, bPTRIm);
		mulC8(bPTRReal, dPTRReal, bPTRIm, dPTRIm);
		mulC8(lPTRReal, dPTRReal, lPTRIm, dPTRIm);

		((float8*)zPTRReal->m256_f32)[0] = ((float8*)dPTRReal->m256_f32)[0];
		((float8*)zPTRIm->m256_f32)[0] = ((float8*)dPTRIm->m256_f32)[0];

		isAllSquared = squaredNorm8(zPTRReal, zPTRIm);


		//Set the mask
		mask = cmp8(isAllSquared, set8(max_iterations_num), _CMP_LT_OS);
		incrementArray = combine(set8(0.f), set8(2.f), mask);
		leaveEarly = combine(set8(0.f), isAllSquared, mask);

		float* ptleaveEarly = (float*)(&leaveEarly);

		//Using the max() instruction could be better for this case
		if ((ptleaveEarly[0] + ptleaveEarly[1] + ptleaveEarly[2] + ptleaveEarly[3] + ptleaveEarly[4] + ptleaveEarly[5] + ptleaveEarly[6] + ptleaveEarly[7]) == 0)
		{
			break;
		}

		//Add 2 to the array if they need it
		((int8*)i.m256i_i32)[0] = (int8)fToI(add8(iToF(i), incrementArray));
		
		j += 2;
		if (j > 255)
		{
			keepLooping = false;
		}
	}
	return fToI(min8(iToF(i), iToF(max_iterations)));
}



void SimpleFractalDrawing(unsigned char* image, int dim[2], float range[2][2])
{
	//Chrono c;
	for (int j = 0; j < dim[1]; j++)
	{
		for (int i = 0; i < dim[0]; i++)
		{
			complex c;
			c.real = range[0][0] + (i + 0.5) * (range[0][1] - range[0][0]) / dim[0]; //Create x coordinates within the range [range[0][0] .. range[0][1]] 
			c.im = range[1][0] + (j + 0.5) * (range[1][1] - range[1][0]) / dim[1]; //Create x coordinates within the range [range[1][0] .. range[1][1]] 
			float f = 2 * Iterate(c);
			if (f > 255.)
			{
				f = 255.;
			}
			image[j * dim[0] + i] = f;
		}
	}
	//c.PrintElapsedTime("time CPU (Single) (s): ");
}

void SimpleFractalDrawingMT(unsigned char* image, int dim[2], float range[2][2], int threadNum, std::mutex* lock1)
{
	Chrono c;
	int iLength = dim[0];
	int jLength = dim[1];
	int thisHeight;
	int thisWidth;
	while (currentHeight < (jLength - 1)) 
	{

		//Get the row that this thread will next work on
		lock1->lock();
		currentHeight++;
		thisHeight = currentHeight;
		lock1->unlock();

		for (int i = 0; i < iLength; i++)
		{
			complex c;
			c.real = range[0][0] + (i + 0.5) * (range[0][1] - range[0][0]) / dim[0]; //Create x coordinates within the range [range[0][0] .. range[0][1]] 
			c.im = range[1][0] + (thisHeight + 0.5) * (range[1][1] - range[1][0]) / dim[1]; //Create x coordinates within the range [range[1][0] .. range[1][1]] 
			float f = 2 * Iterate(c);
			if (f > 255.)
				f = 255.;
			image[thisHeight * dim[0] + i] = f;
		}
	}
	c.PrintElapsedTime("time CPU (Multithreaded) Core (s): ");
}

void SimpleFractalDrawingSIMD(unsigned char* image, int dim[2], float range[2][2])
{
	//Chrono chr;
	int8 f1;
	float8 compareNum = set8(255);

	complex8 c;
	float8* cPTRReal = (float8*)c.real;
	float8* cPTRIm = (float8*)c.im;
	float8 index = { 0,1,2,3,4,5,6,7 };
	float8 set1;
	float8 set2;
	for (int j = 0; j < dim[1]; j++)
	{
		for (int i = 0; i < dim[0]; i = i + 8)
		{

			//Sets the Real value
			set1 = sub8(set8(range[0][1]), set8(range[0][0]));
			set1 = div8(set1, set8(dim[0]));
			set1 = mul8(set1, add8(add8(set8(i), index), set8(0.5f)));
			set1 = add8(set1, set8(range[0][0]));


			//Sets the Im value
			set2 = sub8(set8(range[1][1]), set8(range[1][0]));
			set2 = div8(set2, set8(dim[1]));
			set2 = mul8(set2, add8(set8(j), set8(0.5f)));
			set2 = add8(set2, set8(range[1][0]));


			((float8*)cPTRReal->m256_f32)[0] = set1;
			((float8*)cPTRIm->m256_f32)[0] = set2;


			((int8*)f1.m256i_i32)[0] = Iter8(cPTRReal, cPTRIm);
			f1 = fToI(mul8(set8(2.f), iToF(f1)));

			//Compares whether the stored values are greater than 255
			int8 all255 = set8i(255);
			int8 isLarger = compareGT(f1, all255);
			f1 = fToI(combine(iToF(f1), iToF(all255), iToF(isLarger)));


			image[j * dim[0] + (i + 0)] = f1.m256i_i32[0];
			image[j * dim[0] + (i + 1)] = f1.m256i_i32[1];
			image[j * dim[0] + (i + 2)] = f1.m256i_i32[2];
			image[j * dim[0] + (i + 3)] = f1.m256i_i32[3];
			image[j * dim[0] + (i + 4)] = f1.m256i_i32[4];
			image[j * dim[0] + (i + 5)] = f1.m256i_i32[5];
			image[j * dim[0] + (i + 6)] = f1.m256i_i32[6];
			image[j * dim[0] + (i + 7)] = f1.m256i_i32[7];

		}
	}
	//chr.PrintElapsedTime("time CPU SIMD (s): ");
}

void SimpleFractalDrawingSIMD_MT(unsigned char* image, int dim[2], float range[2][2], int threadNum, std::mutex* lock1)
{
	Chrono chr;
	float8 f;
	int8 f1;
	float8 compareNum;

	complex8 c;
	float8* cPTRReal = (float8*)c.real;
	float8* cPTRIm = (float8*)c.im;
	float8 index = { 0,1,2,3,4,5,6,7 };
	float8 set1;
	float8 set2;

	int iLength = dim[0];
	int jLength = dim[1];
	int thisHeight;
	int thisWidth;

	
	while (currentHeight < (jLength - 1))
	{
		//Get the row that this thread will next work on
		lock1->lock();
		currentHeight++;
		thisHeight = currentHeight;
		lock1->unlock();
		
		for (int i = 0; i < iLength; i = i + 8)
		{


			//Sets the Real value
			set1 = sub8(set8(range[0][1]), set8(range[0][0]));
			set1 = div8(set1, set8(iLength));
			set1 = mul8(set1, add8(add8(set8(i), index), set8(0.5f)));
			set1 = add8(set1, set8(range[0][0]));


			//Sets the Im value
			set2 = sub8(set8(range[1][1]), set8(range[1][0]));
			set2 = div8(set2, set8(jLength));
			set2 = mul8(set2, add8(set8(thisHeight), set8(0.5f)));
			set2 = add8(set2, set8(range[1][0]));


			((float8*)cPTRReal->m256_f32)[0] = set1;
			((float8*)cPTRIm->m256_f32)[0] = set2;


			((int8*)f1.m256i_i32)[0] = Iter8(cPTRReal, cPTRIm);
			f1 = fToI(mul8(set8(2.f), iToF(f1)));

			//Compares whether the stored values are greater than 255
			int8 all255 = set8i(255);
			int8 isLarger = compareGT(f1, all255);
			f1 = fToI(combine(iToF(f1), iToF(all255), iToF(isLarger)));

			image[thisHeight * jLength + (i + 0)] = f1.m256i_i32[0];
			image[thisHeight * jLength + (i + 1)] = f1.m256i_i32[1];
			image[thisHeight * jLength + (i + 2)] = f1.m256i_i32[2];
			image[thisHeight * jLength + (i + 3)] = f1.m256i_i32[3];
			image[thisHeight * jLength + (i + 4)] = f1.m256i_i32[4];
			image[thisHeight * jLength + (i + 5)] = f1.m256i_i32[5];
			image[thisHeight * jLength + (i + 6)] = f1.m256i_i32[6];
			image[thisHeight * jLength + (i + 7)] = f1.m256i_i32[7];


		}
	}
	chr.PrintElapsedTime("time CPU SIMD_MT Core: (s): ");
}


int main(int argc, char* argv[])
{

	std::mutex lock1;

	int dims[2] = { 1024,1024 };
	float range[2][2] = { {-0.003,0.008},{-0.0002,0.0005} };
	unsigned char* image = new unsigned char[dims[0] * dims[1]];

	//WORK ON SIMPLE, SINGLE CORE
	Chrono c;
	SimpleFractalDrawing(image, dims, range); //largest 64bit prime
	SaveBMP("fractal.bmp", image, dims[0], dims[1]);
	c.PrintElapsedTime("time CPU (Single Core) (s): ");
	for (int i = 0; i < dims[0] * dims[1]; i++)
		image[i] = 127; //resetting image to grey



	//WORK ON SIMPLE, MULTIPLE CORES
	c;
	for (int i = 0; i < NB_THREADS; i++)
	{
		thr[i] = std::thread(SimpleFractalDrawingMT, image, dims, range, i, &lock1);
	}

	for (int i = 0; i < NB_THREADS; i++)
	{
		thr[i].join();
	}
	c.PrintElapsedTime("time CPU (Multithreaded) Overall: (s): ");
	SaveBMP("fractalMT.bmp", image, dims[0], dims[1]);
	
	currentHeight = -1; //Set the currentHeight in preparation for the SIMD_MT

	for (int i = 0; i < dims[0] * dims[1]; i++)
		image[i] = 127; //resetting image to grey



	//WORK ON SIMD, SINGLE CORE
	c;
	SimpleFractalDrawingSIMD(image, dims, range);
	c.PrintElapsedTime("time CPU (SIMD Single Core) (s): ");
	SaveBMP("fractalSIMD.bmp", image, dims[0], dims[1]);

	for (int i = 0; i < dims[0] * dims[1]; i++)
		image[i] = 127;



	//WORK ON SIMD, MULTIPLE CORES
	c;
	for (int i = 0; i < NB_THREADS; i++)
	{
		thr[i] = std::thread(SimpleFractalDrawingSIMD_MT, image, dims, range, i, &lock1);
	}

	for (int i = 0; i < NB_THREADS; i++)
	{
		thr[i].join();
	}
	c.PrintElapsedTime("time CPU SIMD_MT Overall: (s): ");

	SaveBMP("fractalSIMD_MT.bmp", image, dims[0], dims[1]);
	delete[] image;
	return 0;
}
