#include<iostream>
using namespace std;
/***
��gdal&&lib�ļ�������ĿĿ¼��
��gdal18.dll������Ŀ��ReleaseĿ¼��
***/
#include "./gdal/gdal_priv.h"
#pragma comment(lib,"gdal_i.lib")
int main() {
	//���Ӵ���
	//����ͼ��
	GDALDataset* poSrcDS;
	//���ͼ��
	GDALDataset* poDstDS;
	//ͼ��Ŀ��Ⱥ͸߶�
	int imgXlen, imgYlen;
	//����ͼ��·��
	char* srcPath = "trees.jpg";
	//���ͼ��·��
	char* dstPath = "res.tif";
	//ͼ���ڴ�洢 byte 8λ 255 unsign char
	GByte* buffTmp;
	//ͼ�񲨶���
	int i, bandNum;
	//ע������
	GDALAllRegister();
	//��ͼ��
	poSrcDS = (GDALDataset*)GDALOpenShared(srcPath, GA_ReadOnly);

	//��ȡͼ����ȡ��߶ȺͲ�������
	imgXlen = poSrcDS->GetRasterXSize();
	imgYlen = poSrcDS->GetRasterYSize();
	bandNum = poSrcDS->GetRasterCount();
	//�����ȡ���
	cout << "Image X Length: " << imgXlen << endl;
	cout << "Image Y Length: " << imgYlen << endl;
	cout << "Band number :   " << bandNum << endl;
	//����ͼ��Ŀ��Ⱥ͸߶ȷ����ڴ�
	buffTmp = (GByte*)CPLMalloc(imgXlen*imgYlen * sizeof(GByte));
	//�������ͼ��
	poDstDS = GetGDALDriverManager()->GetDriverByName("GTiff")->Create(dstPath, imgXlen, imgYlen, bandNum, GDT_Byte, NULL);
	//һ�������ε����룬Ȼ��һ�������ε����
	for (i = 0;i < bandNum;i++) {
		poSrcDS->GetRasterBand(i + 1)->RasterIO(GF_Read,
			0, 0, imgXlen, imgYlen, buffTmp, imgXlen, imgYlen, GDT_Byte, 0, 0);
		poDstDS->GetRasterBand(i + 1)->RasterIO(GF_Write,
			0, 0, imgXlen, imgYlen, buffTmp, imgXlen, imgYlen, GDT_Byte, 0, 0);
		printf("... ... band %d processing ... ...\n", i);
	}
	//����ڴ�
	CPLFree(buffTmp);
	//�ر�dataset
	GDALClose(poDstDS);
	GDALClose(poSrcDS);

	system("PAUSE");
	return 0;
}