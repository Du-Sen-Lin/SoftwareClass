#include <iostream>
#include <math.h>

#include <string.h>
#include <algorithm>
using namespace std;
#define PI 3.14
#define uchar unsigned char

#include "./gdal/gdal_priv.h"
#pragma comment(lib,"gdal_i.lib")
int main(int argc, const char * argv[])
{
	// insert code here...
	void RGBtoHIS(uchar r, uchar g, uchar b, float* HIS);
	void HIStoRGB(float H, float I, float S, float* m);
	void HistogramMatch(GDALDataset *HPan, GDALDataset *HI, const char* str);
	void stretch(GDALDataset *ss, const char* str);

	GDALAllRegister();         //先要进行注册
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");   //设置支持中文路径

														 //打开RGB波段影像与pan波段影像文件
	const char *str1 = "/Users/Star/Desktop/遥感定量理论/IHS/HISdata/landsat_hz_r.tif";
	GDALDataset *picR = (GDALDataset*)GDALOpen(str1, GA_ReadOnly);

	const char *str2 = "/Users/Star/Desktop/遥感定量理论/IHS/HISdata/landsat_hz_g.tif";
	GDALDataset *picG = (GDALDataset*)GDALOpen(str2, GA_ReadOnly);

	const char *str3 = "/Users/Star/Desktop/遥感定量理论/IHS/HISdata/landsat_hz_b.tif";
	GDALDataset *picB = (GDALDataset*)GDALOpen(str3, GA_ReadOnly);

	const char *str4 = "/Users/Star/Desktop/遥感定量理论/IHS/HISdata/landsat_hz_pan.tif";
	GDALDataset *picPan = (GDALDataset*)GDALOpen(str4, GA_ReadOnly);

	int width = picR->GetRasterXSize();
	int height = picR->GetRasterYSize();

	float *RDataBuff = new float[width];
	float *GDataBuff = new float[width];
	float *BDataBuff = new float[width];
	float *TotalBuff = new float[width];

	//创建存储文件的空文件
	GDALDriver *poDriver1 = (GDALDriver*)GDALGetDriverByName("GTiff");
	GDALDataset *picH = poDriver1->Create("/Users/Star/Desktop/遥感定量理论/IHS/workplace/H.tif", width, height, 1, GDT_Byte, NULL);
	GDALDriver *poDriver2 = (GDALDriver*)GDALGetDriverByName("GTiff");
	GDALDataset *picI = poDriver2->Create("/Users/Star/Desktop/遥感定量理论/IHS/workplace/I.tif", width, height, 1, GDT_Byte, NULL);
	GDALDriver *poDriver3 = (GDALDriver*)GDALGetDriverByName("GTiff");
	GDALDataset *picS = poDriver3->Create("/Users/Star/Desktop/遥感定量理论/IHS/workplace/S.tif", width, height, 1, GDT_Byte, NULL);
	GDALDriver *poDriver4 = (GDALDriver*)GDALGetDriverByName("GTiff");
	GDALDataset *picMix = poDriver4->Create("/Users/Star/Desktop/遥感定量理论/IHS/workplace/result.tif", width, height, 3, GDT_Byte, NULL);

	float *HDataBuff = new float[width];
	float *IDataBuff = new float[width];
	float *SDataBuff = new float[width];
	uchar m1, m2, m3;
	float m[3];

	for (int i = 0; i < height; i++)
	{
		//读取影像数据
		picR->GetRasterBand(1)->RasterIO(GF_Read, 0, i, width, 1, RDataBuff, width, 1, GDT_Float32, 0, 0);
		picG->GetRasterBand(1)->RasterIO(GF_Read, 0, i, width, 1, GDataBuff, width, 1, GDT_Float32, 0, 0);
		picB->GetRasterBand(1)->RasterIO(GF_Read, 0, i, width, 1, BDataBuff, width, 1, GDT_Float32, 0, 0);

		//IHS正变换
		for (int j = 0; j < width; j++)
		{
			m1 = RDataBuff[j];
			m2 = GDataBuff[j];
			m3 = BDataBuff[j];
			RGBtoHIS(m1, m2, m3, m);
			HDataBuff[j] = m[0];
			IDataBuff[j] = m[1];
			SDataBuff[j] = m[2];
		}

		//将IHS影像数据写入文件中
		picH->GetRasterBand(1)->RasterIO(GF_Write, 0, i, width, 1, HDataBuff, width, 1, GDT_Float32, 0, 0);
		picI->GetRasterBand(1)->RasterIO(GF_Write, 0, i, width, 1, IDataBuff, width, 1, GDT_Float32, 0, 0);
		picS->GetRasterBand(1)->RasterIO(GF_Write, 0, i, width, 1, SDataBuff, width, 1, GDT_Float32, 0, 0);
	}
	const char* path = "/Users/Star/Desktop/遥感定量理论/IHS/workplace/pan.tif";//全色图像

	//全色波段与I波段进行直方图匹配
	HistogramMatch(picPan, picI, path);


	float *RBuff = new float[width];
	float *GBuff = new float[width];
	float *BBuff = new float[width];

	float n[3];
	float n1, n2, n3;
	GDALDataset *change_I = (GDALDataset *)GDALOpen("/Users/Star/Desktop/遥感定量理论/IHS/workplace/pan.tif", GA_ReadOnly);

	//IHS空间进行拉伸
	path = "/Users/Star/Desktop/遥感定量理论/IHS/workplace/S_H.tif";
	stretch(picH, path);
	path = "/Users/Star/Desktop/遥感定量理论/IHS/workplace/S_I.tif";
	stretch(change_I, path);
	path = "/Users/Star/Desktop/遥感定量理论/IHS/workplace/S_S.tif";
	stretch(picS, path);

	//重新导入数据
	GDALDataset *trans_H = (GDALDataset *)GDALOpen("/Users/Star/Desktop/遥感定量理论/IHS/workplace/S_H.tif", GA_ReadOnly);
	GDALDataset *trans_I = (GDALDataset *)GDALOpen("/Users/Star/Desktop/遥感定量理论/IHS/workplace/S_I.tif", GA_ReadOnly);
	GDALDataset *trans_S = (GDALDataset *)GDALOpen("/Users/Star/Desktop/遥感定量理论/IHS/workplace/S_S.tif", GA_ReadOnly);

	for (int i = 0; i < height; i++)
	{
		picH->GetRasterBand(1)->RasterIO(GF_Read, 0, i, width, 1, HDataBuff, width, 1, GDT_Float32, 0, 0);
		trans_S->GetRasterBand(1)->RasterIO(GF_Read, 0, i, width, 1, SDataBuff, width, 1, GDT_Float32, 0, 0);
		change_I->GetRasterBand(1)->RasterIO(GF_Read, 0, i, width, 1, TotalBuff, width, 1, GDT_Float32, 0, 0);
		for (int j = 0; j < width; j++)
		{
			n1 = HDataBuff[j];
			n2 = TotalBuff[j];
			n3 = SDataBuff[j];
			HIStoRGB(n1, n2, n3, n);
			RBuff[j] = n[0];
			GBuff[j] = n[1];
			BBuff[j] = n[2];
		}

		picMix->GetRasterBand(1)->RasterIO(GF_Write, 0, i, width, 1, RBuff, width, 1, GDT_Float32, 0, 0);
		picMix->GetRasterBand(2)->RasterIO(GF_Write, 0, i, width, 1, GBuff, width, 1, GDT_Float32, 0, 0);
		picMix->GetRasterBand(3)->RasterIO(GF_Write, 0, i, width, 1, BBuff, width, 1, GDT_Float32, 0, 0);
	}

	//释放内存
	delete[] RBuff;
	delete[] GBuff;
	delete[] BBuff;
	delete[] RDataBuff;
	delete[] GDataBuff;
	delete[] BDataBuff;
	delete[] TotalBuff;
	delete[] HDataBuff;
	delete[] IDataBuff;
	delete[] SDataBuff;
	GDALClose((GDALDatasetH)trans_I);
	GDALClose((GDALDatasetH)picR);
	GDALClose((GDALDatasetH)picG);
	GDALClose((GDALDatasetH)picB);
	GDALClose((GDALDatasetH)picPan);
	GDALClose((GDALDatasetH)picS);
	GDALClose((GDALDatasetH)picI);
	GDALClose((GDALDatasetH)picH);
	GDALClose((GDALDatasetH)picMix);
	return 0;
}

//IHS正变换
void RGBtoHIS(uchar r, uchar g, uchar b, float* HIS)
{
	float H;
	float I;
	float S;

	float angle;
	angle = 0.0;   //初始化
	float R, G, B;

	//归一化
	R = ((float)r) / 255.0;
	G = ((float)g) / 255.0;
	B = ((float)b) / 255.0;

	//三角形变换
	I = (R + G + B) / 3.0;
	S = 1 - min(min(R, G), B) / I;
	angle = acos(0.5*(R - G + R - B) / (sqrt((R - G)*(R - G) + (R - B)*(G - B))));
	angle = angle*180.0 / PI;
	if (B>G)
	{
		H = 360 - angle;        //H量计算
	}
	else
	{
		H = angle;
	}

	H = H*255.0 / 360.0;
	I = I*255.0;
	S = S*255.0;

	HIS[0] = H;
	HIS[1] = I;
	HIS[2] = S;
}

//IHS反变换
void HIStoRGB(float H, float I, float S, float* RGB)
{
	float s = (S) / 255;
	float i = (I) / 255;
	float h = (H) * 360 / 255;
	if (h<0)
	{
		h = h + 360;
	}
	if (h<120 && h >= 0)
	{
		RGB[0] = i*(1.0 + ((s*cos(h*PI / 180)) / cos((60 - h)*PI / 180)));
		RGB[2] = i*(1.0 - s);
		RGB[1] = 3.0*i - (RGB[0] + RGB[2]);
	}

	if (h<240 && h >= 120)
	{
		h = h - 120;
		RGB[1] = i*(1.0 + ((s*cos(h*PI / 180)) / cos((60 - h)*PI / 180)));
		RGB[0] = i*(1.0 - s);
		RGB[2] = 3.0*i - (RGB[0] + RGB[1]);
	}
	if (h >= 240 && h <= 360)
	{
		h = h - 240;
		RGB[2] = i*(1.0 + ((s*cos(h*PI / 180)) / cos((60 - h)*PI / 180)));
		RGB[1] = i*(1.0 - s);
		RGB[0] = 3 * i - (RGB[1] + RGB[2]);
	}
	RGB[0] = (RGB[0] * 255);
	RGB[1] = (RGB[1] * 255);
	RGB[2] = (RGB[2] * 255);
}

//线性拉伸
void stretch(GDALDataset *ss, const char* str)
{
	int iWidth = ss->GetRasterXSize();
	int iHeight = ss->GetRasterYSize();
	int Size = iWidth * iHeight;

	//计算输入波段的累计直方图
	float fpPro[256];
	memset(fpPro, 0, sizeof(float) * 256);
	GDALRasterBand *poBand_I = ss->GetRasterBand(1);
	GUIntBig anHistogram_ss[256];
	poBand_I->GetHistogram(-0.5, 255.5, 256, anHistogram_ss, FALSE, FALSE, GDALDummyProgress, NULL);

	//直方图归一化
	for (int i = 0; i < 256; i++)
	{
		fpPro[i] = (float)anHistogram_ss[i] / Size;
	}

	//计算累计直方图
	float temp1[256];
	for (int i = 0; i < 256; i++)
	{
		if (i == 0)
			temp1[0] = fpPro[0];
		else
			temp1[i] = temp1[i - 1] + fpPro[i];
		fpPro[i] = temp1[i];
	}

	//前后各拉伸1%
	unsigned char Map[256];
	int flag1, flag2;
	for (int i = 0; i < 256; i++)
	{
		if ((fpPro[i] - 0.01) <= 0.001)
		{
			flag1 = i;
		}
		if ((fpPro[i] - 0.99) <= 0.001)
		{
			flag2 = i;
		}
	}
	//确定映射关系//
	for (int i = 0; i < 256; i++)
	{
		int c = 0;
		for (int j = 0; j < 256; j++)
		{
			if (i<flag1)
			{
				c = 0;
			}
			if (i>flag2)
			{
				c = 255;
			}
			if (i >= flag1 && i <= flag2)
			{
				c = (int)((float)(i - flag1) / (float)(flag2 - flag1) * 255);
			}
		}
		//建立灰度映射表//
		Map[i] = (uchar)c;
	}
	//一一进行映射
	unsigned char * poDataPan = new unsigned char[Size];
	unsigned char * poNewDataPan = new unsigned char[Size];
	poBand_I->RasterIO(GF_Read, 0, 0, iWidth, iHeight, poDataPan, iWidth, iHeight, GDT_Byte, 0, 0);
	int m = 0;
	for (int ny = 0; ny < iHeight; ny++)
	{
		for (int nx = 0; nx < iWidth; nx++)
		{
			unsigned char x = poDataPan[m];
			poNewDataPan[m] = (unsigned char)(Map[x] + 0.5f);
			m++;
		}
	}

	GDALDriver * poDriver = (GDALDriver *)GDALGetDriverByName("GTiff");
	GDALDataset * poDSResult = poDriver->Create(str, iWidth, iHeight, 1, GDT_Float32, NULL);
	poDSResult->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, iWidth, iHeight, poNewDataPan, iWidth, iHeight, GDT_Byte, 0, 0);
	GDALClose((GDALDatasetH)poDSResult);
}

//直方图匹配
void HistogramMatch(GDALDataset *HPan, GDALDataset *HI, const char* str)
{
	int iWidth = HI->GetRasterXSize();
	int iHeight = HI->GetRasterYSize();
	int Size = iWidth * iHeight;

	//计算I分量的累计直方图
	float fpPro[256];
	memset(fpPro, 0, sizeof(float) * 256);
	GDALRasterBand *poBand_I = HI->GetRasterBand(1);
	GUIntBig anHistogram_I[256];
	poBand_I->GetHistogram(-0.5, 255.5, 256,anHistogram_I, FALSE, FALSE, GDALDummyProgress, NULL);

	//直方图归一化
	for (int i = 0; i < 256; i++)
	{
		fpPro[i] = (float)anHistogram_I[i] / Size;
	}

	//计算累计直方图
	float temp1[256];
	for (int i = 0; i < 256; i++)
	{
		if (i == 0)
			temp1[0] = fpPro[0];
		else
			temp1[i] = temp1[i - 1] + fpPro[i];
		fpPro[i] = temp1[i];
	}

	//计算Pan波段的累计直方图
	float pPro[256];
	memset(pPro, 0, sizeof(float) * 256);
	GDALRasterBand *poBand_Pan = HPan->GetRasterBand(1);
	GUIntBig anHistogram_Pan[256];
	poBand_Pan->GetHistogram(-0.5, 255.5, 256, anHistogram_Pan, FALSE, FALSE, GDALDummyProgress, NULL);

	//直方图归一化
	for (int i = 0; i < 256; i++)
	{
		pPro[i] = (float)anHistogram_Pan[i] / Size;
	}

	//计算累计直方图
	float temp2[256];
	for (int i = 0; i < 256; i++)
	{
		if (i == 0)
			temp2[0] = pPro[0];
		else
			temp2[i] = temp2[i - 1] + pPro[i];
		pPro[i] = temp2[i];
	}

	//确定映射关系//
	unsigned char Map[256];
	for (int i = 0; i < 256; i++)
	{
		int c = 0;                           //最接近的规定直方图灰度级变量//
		float min_value = 1.0;              //最小差值变量//
											//枚举规定直方图各个灰度//
		for (int j = 0; j < 256; j++)
		{
			float now_value = 0.0;         //当前插值变量//
										   //差值计算//
			if (pPro[i] - fpPro[j] >= 0.0f)
				now_value = pPro[i] - fpPro[j];
			else
				now_value = fpPro[j] - pPro[i];
			//寻找最接近的规定直方图灰度级//
			if (now_value < min_value)
			{
				c = j;                      //最接近的灰度级//
				min_value = now_value;      //最小差值//
			}
		}
		//建立灰度映射表//
		Map[i] = (uchar)c;
	}
	//一一进行映射
	unsigned char * poDataPan = new unsigned char[Size];
	unsigned char * poNewDataPan = new unsigned char[Size];
	poBand_Pan->RasterIO(GF_Read, 0, 0, iWidth, iHeight, poDataPan, iWidth, iHeight, GDT_Byte, 0, 0);
	int m = 0;
	for (int ny = 0; ny < iHeight; ny++)
	{
		for (int nx = 0; nx < iWidth; nx++)
		{
			unsigned char x = poDataPan[m];
			poNewDataPan[m] = (unsigned char)(Map[x] + 0.5f);
			m++;
		}
	}

	GDALDriver * poDriver = (GDALDriver *)GDALGetDriverByName("GTiff");
	GDALDataset * poDSResult = poDriver->Create(str, iWidth, iHeight, 1, GDT_Float32, NULL);
	poDSResult->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, iWidth, iHeight, poNewDataPan, iWidth, iHeight, GDT_Byte, 0, 0);
	GDALClose((GDALDatasetH)poDSResult);
}
