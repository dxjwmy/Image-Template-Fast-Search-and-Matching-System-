#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include<fstream>
#include <cstdio>
#include <string.h>   
#include <math.h> 
#include <cmath>
#include <stdlib.h>     
#include <malloc.h>   
#include <vector>
#include <string>
#include <iomanip>
#include <chrono>
using namespace std;
using namespace std::chrono;
#define PI 3.1415926
#define   WIDTHBYTES(bits) (((bits)+31)/32*4)//����ʹͼ������ռ�ֽ���Ϊ4byte�ı���  
typedef unsigned char  BYTE;
typedef unsigned short WORD;
typedef unsigned long  DWORD; 
typedef int LONG;

struct true_data
{
    int start_x;
    int start_y;
    double height;
    double width;
};

typedef struct tagBITMAPFILEHEADER {
    DWORD  bfSize;          //�ļ���С  
    WORD   bfReserved1;     //�����֣�������  
    WORD   bfReserved2;     //�����֣�ͬ��  
    DWORD  bfOffBits;       //ʵ��λͼ���ݵ�ƫ���ֽ�������ǰ�������ֳ���֮��  
} BITMAPFILEHEADER;
typedef struct tagBITMAPINFOHEADER {
    //public:  
    DWORD   biSize;             //ָ���˽ṹ��ĳ��ȣ�Ϊ40  
    LONG    biWidth;            //λͼ��  
    LONG    biHeight;           //λͼ��  
    WORD    biPlanes;           //ƽ������Ϊ1  
    WORD    biBitCount;         //������ɫλ����������1��2��4��8��16��24���µĿ�����32  
    DWORD   biCompression;      //ѹ����ʽ��������0��1��2������0��ʾ��ѹ��  
    DWORD   biSizeImage;        //ʵ��λͼ����ռ�õ��ֽ���  
    LONG    biXPelsPerMeter;    //X����ֱ���  
    LONG    biYPelsPerMeter;    //Y����ֱ���  
    DWORD   biClrUsed;          //ʹ�õ���ɫ�������Ϊ0�����ʾĬ��ֵ(2^��ɫλ��)  
    DWORD   biClrImportant;     //��Ҫ��ɫ�������Ϊ0�����ʾ������ɫ������Ҫ��  
} BITMAPINFOHEADER;

class BMP
{
private:
    BITMAPFILEHEADER bitHead;
    BITMAPINFOHEADER  bitInfoHead;
    int width;
    int height;
    WORD fileType;
public:
    BYTE* pColorData;
    BMP();
    BMP(int width, int height);
    ~BMP();
    void setBitHead(BITMAPFILEHEADER bitHead);
    void setBitInfoHead(BITMAPINFOHEADER  bitInfoHead);
    void setPColorData(BYTE* pColorData);
    BITMAPFILEHEADER getBitHead();
    BITMAPINFOHEADER getBitInfoHead();
    BYTE* getPColorData();
    WORD getBitCount();
    void setfileType(WORD fileType);
    WORD getfileType();
    int rows();
    int cols();
    vector<vector<double>>rgb2gray;
    double mu = 0;
    double sigma = 0;
};
BMP::BMP()
{
}
BMP::BMP(int width, int height) :width(width), height(height)    //�в����Ĺ��캯����ͨ����ʼ���б�ֵ
{
    this->bitHead.bfSize = width * height * 3 + 54;
    this->bitHead.bfReserved1 = 0;     //�����֣�������  
    this->bitHead.bfReserved2 = 0;     //�����֣�ͬ��  
    this->bitHead.bfOffBits = 54;
    this->fileType = 19778;
    this->bitInfoHead.biSize = 40;             //ָ���˽ṹ��ĳ��ȣ�Ϊ40  
    this->bitInfoHead.biWidth = width;            //λͼ��  
    this->bitInfoHead.biHeight = height;           //λͼ��  
    this->bitInfoHead.biPlanes = 1;           //ƽ������Ϊ1  
    this->bitInfoHead.biBitCount = 24;         //������ɫλ����������1��2��4��8��16��24���µĿ�����32  
    this->bitInfoHead.biCompression = 0;      //ѹ����ʽ��������0��1��2������0��ʾ��ѹ��  
    this->bitInfoHead.biSizeImage = width * height * 3;        //ʵ��λͼ����ռ�õ��ֽ���  
    this->bitInfoHead.biXPelsPerMeter = 0;    //X����ֱ���  
    this->bitInfoHead.biYPelsPerMeter = 0;    //Y����ֱ���  
    this->bitInfoHead.biClrUsed = 0;          //ʹ�õ���ɫ�������Ϊ0�����ʾĬ��ֵ(2^��ɫλ��)  
    this->bitInfoHead.biClrImportant = 0;
    //дλͼ������
    int write_width = (((width * 24) + 31) / 32 * 4);//����дλͼ��ʵ�ʿ�ȣ��ֽڣ���ȷ����Ϊ4byte�ı���   
    this->pColorData = (BYTE*)malloc(write_width * height);//�����ڴ�ռ�洢ͼ����֮������
    this->rgb2gray.resize(height, vector<double>(width));
    memset(this->pColorData, 0, write_width * height);
}
BMP::~BMP()
{
}
void BMP::setfileType(WORD fileType)
{
    this->fileType = fileType;
}
WORD BMP::getfileType()
{
    return this->fileType;
}
BITMAPFILEHEADER BMP::getBitHead()
{
    return this->bitHead;
}

BITMAPINFOHEADER BMP::getBitInfoHead()
{
    return this->bitInfoHead;
}
BYTE* BMP::getPColorData()
{
    return this->pColorData;
}
void BMP::setBitHead(BITMAPFILEHEADER bitHead)
{
    this->bitHead = bitHead;
}

void BMP::setBitInfoHead(BITMAPINFOHEADER  bitInfoHead)
{
    this->bitInfoHead = bitInfoHead;
    this->width = bitInfoHead.biWidth;
    this->height = bitInfoHead.biHeight;
}
void BMP::setPColorData(BYTE* pColorData)
{
    this->pColorData = pColorData;
}

int BMP::rows()
{
    return this->height;
}
int BMP::cols()
{
    return this->width;
}
WORD BMP::getBitCount()
{
    return this->bitInfoHead.biBitCount;
}

int imread(const char strFile[50], BMP& bmp)
{
    FILE* pfile = fopen(strFile, "rb");//�ļ���ͼ��
    BITMAPFILEHEADER  bitHead;
    BITMAPINFOHEADER  bitInfoHead;

    WORD fileType;
    fread(&fileType, 1, sizeof(WORD), pfile);
    bmp.setfileType(fileType);
    if (fileType != 0x4d42)
    {
        std::cout << "file is not .bmp file!" << std::endl;
        return 0;

    }
    fread(&bitHead, 1, sizeof(tagBITMAPFILEHEADER), pfile);
    fread(&bitInfoHead, 1, sizeof(BITMAPINFOHEADER), pfile);
    int width = bitInfoHead.biWidth;
    int height = bitInfoHead.biHeight;
    bmp.rgb2gray.resize(height, vector<double>(width));
    //�����ڴ�ռ��Դͼ�����ڴ�     
    int l_width = WIDTHBYTES(width * bitInfoHead.biBitCount);//����λͼ��ʵ�ʿ�Ȳ�ȷ����Ϊ4byte�ı���  
    BYTE* pColorData = (BYTE*)malloc(height * l_width);//�����ڴ�ռ�洢ͼ������  
    memset(pColorData, 0, height * l_width);
    long nData = height * l_width; //��λͼ������Ϣ����������     
    fread(pColorData, 1, nData, pfile);//ͼ�����ͨ�������ⲿ�����ݼ���ʵ��  
    bmp.setBitHead(bitHead);
    bmp.setBitInfoHead(bitInfoHead);
    bmp.setPColorData(pColorData);
    fclose(pfile);

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            int data_r = static_cast<int>(pColorData[(y)*width * 3 + (x) * 3 + 2]);
            int data_g = static_cast<int>(pColorData[(y)*width * 3 + (x) * 3 + 1]);
            int data_b = static_cast<int>(pColorData[(y)*width * 3 + (x) * 3 + 0]);
            bmp.rgb2gray[y][x] = (double)(data_r * 0.299 + data_g * 0.587 + data_b * 0.114);
        }
    }
    if (height < 20 && width < 20)
    {
        for (int y = 0; y < height; y++)
        {
            for (int x = 0; x < width; x++)
            {

                bmp.mu += bmp.rgb2gray[y][x];
            }
        }
        double size = (height * width);
        bmp.mu /= size;

        for (int y = 0; y < height; y++)
        {
            for (int x = 0; x < width; x++)
            {

                bmp.sigma += ((bmp.rgb2gray[y][x] - bmp.mu), 2);
            }
        }
        bmp.sigma = pow(bmp.sigma / size, 0.5);
    }
}
void imwrite(const char strFile[50], BMP& bmp)
{
    FILE* wfile = fopen(strFile, "wb");//���ļ�Ϊ�洢�޸ĺ�ͼ����׼��  
    WORD fileType = bmp.getfileType();
    fwrite(&fileType, 1, sizeof(WORD), wfile);
    BITMAPFILEHEADER writebitHead = bmp.getBitHead();
    BITMAPINFOHEADER  writebitInfoHead = bmp.getBitInfoHead();
    int mywritewidth = WIDTHBYTES(writebitInfoHead.biWidth * writebitInfoHead.biBitCount);//BMPͼ��ʵ��λͼ�������Ŀ��Ϊ4byte�ı������ڴ˼���ʵ�����������  

    long write_nData = mywritewidth * writebitInfoHead.biHeight;//��ȡ��λͼ���������ȶ���  
    fwrite(&writebitHead, 1, sizeof(tagBITMAPFILEHEADER), wfile);//д��λͼ�ļ�ͷ��Ϣ������ļ�  
    fwrite(&writebitInfoHead, 1, sizeof(BITMAPINFOHEADER), wfile);//д��λͼ��Ϣͷ��Ϣ������ļ�
    BYTE* PColorData = bmp.getPColorData();
    fwrite(PColorData, 1, write_nData, wfile);   //��������ͼ��������д���ļ� 
    fclose(wfile);

}

BMP p_change(BMP bmp, const char change_file[50], double height_change, double width_change)
{
    int width = bmp.cols();
    int height = bmp.rows();
    int DRAW_HEIGHT = height * height_change;
    int DRAW_WIDTH = width * width_change;
    BMP change_bmp(DRAW_WIDTH, DRAW_HEIGHT);
    int l_width = WIDTHBYTES(width * 24);//����λͼ��ʵ�ʿ�Ȳ�ȷ����Ϊ4byte�ı���  
    //дλͼ������
    int write_width = WIDTHBYTES(DRAW_WIDTH * 24);//����дλͼ��ʵ�ʿ�ȣ��ֽڣ���ȷ����Ϊ4byte�ı���   
    /*******************ͼ������******************/
    for (int hnum = 0; hnum < DRAW_HEIGHT; hnum++)
    {
        for (int wnum = 0; wnum < DRAW_WIDTH; wnum++)
        {
            double d_original_img_hnum = hnum * height / (double)DRAW_HEIGHT;
            double d_original_img_wnum = wnum * width / (double)DRAW_WIDTH;
            int i_original_img_hnum = d_original_img_hnum;
            int i_original_img_wnum = d_original_img_wnum;
            double distance_to_a_x = d_original_img_wnum - i_original_img_wnum;//��ԭͼ������a���ˮƽ����  
            double distance_to_a_y = d_original_img_hnum - i_original_img_hnum;//��ԭͼ������a��Ĵ�ֱ����  
            int original_point_a = i_original_img_hnum * l_width + i_original_img_wnum * 3;//����λ��ƫ��������Ӧ��ͼ��ĸ����ص�RGB�����,�൱�ڵ�A    
            int original_point_b = i_original_img_hnum * l_width + (i_original_img_wnum + 1) * 3;//����λ��ƫ��������Ӧ��ͼ��ĸ����ص�RGB�����,�൱�ڵ�B  
            int original_point_c = (i_original_img_hnum + 1) * l_width + i_original_img_wnum * 3;//����λ��ƫ��������Ӧ��ͼ��ĸ����ص�RGB�����,�൱�ڵ�C   
            int original_point_d = (i_original_img_hnum + 1) * l_width + (i_original_img_wnum + 1) * 3;//����λ��ƫ��������Ӧ��ͼ��ĸ����ص�RGB�����,�൱�ڵ�D   
            if (i_original_img_hnum + 1 == DRAW_HEIGHT - 1)
            {
                original_point_c = original_point_a;
                original_point_d = original_point_b;
            }
            if (i_original_img_wnum + 1 == DRAW_WIDTH - 1)
            {
                original_point_b = original_point_a;
                original_point_d = original_point_c;
            }

            int pixel_point = hnum * write_width + wnum * 3;//ӳ��߶ȱ任ͼ������λ��ƫ����  
            change_bmp.pColorData[pixel_point] =
                bmp.pColorData[original_point_a] * (1 - distance_to_a_x) * (1 - distance_to_a_y) +
                bmp.pColorData[original_point_b] * distance_to_a_x * (1 - distance_to_a_y) +
                bmp.pColorData[original_point_c] * distance_to_a_y * (1 - distance_to_a_x) +
                bmp.pColorData[original_point_d] * distance_to_a_y * distance_to_a_x;
            change_bmp.pColorData[pixel_point + 1] =
                bmp.pColorData[original_point_a + 1] * (1 - distance_to_a_x) * (1 - distance_to_a_y) +
                bmp.pColorData[original_point_b + 1] * distance_to_a_x * (1 - distance_to_a_y) +
                bmp.pColorData[original_point_c + 1] * distance_to_a_y * (1 - distance_to_a_x) +
                bmp.pColorData[original_point_d + 1] * distance_to_a_y * distance_to_a_x;
            change_bmp.pColorData[pixel_point + 2] =
                bmp.pColorData[original_point_a + 2] * (1 - distance_to_a_x) * (1 - distance_to_a_y) +
                bmp.pColorData[original_point_b + 2] * distance_to_a_x * (1 - distance_to_a_y) +
                bmp.pColorData[original_point_c + 2] * distance_to_a_y * (1 - distance_to_a_x) +
                bmp.pColorData[original_point_d + 2] * distance_to_a_y * distance_to_a_x;

        }
    }
    /*******************ͼ������******************/

    //imwrite(change_file, change_bmp);

    for (int y = 0; y < DRAW_HEIGHT; y++)
    {
        for (int x = 0; x < DRAW_WIDTH; x++)
        {
            int data_r = static_cast<int>(change_bmp.pColorData[(y)*DRAW_WIDTH * 3 + (x) * 3 + 2]);
            int data_g = static_cast<int>(change_bmp.pColorData[(y)*DRAW_WIDTH * 3 + (x) * 3 + 1]);
            int data_b = static_cast<int>(change_bmp.pColorData[(y)*DRAW_WIDTH * 3 + (x) * 3 + 0]);
            change_bmp.rgb2gray[y][x] = (data_r * 0.299 + data_g * 0.587 + data_b * 0.114);
        }
    }

    //cout << "Done!" << endl;
    return change_bmp;
}

double get_NCC(BMP test, BMP target, int start_y, int start_x)
{
    double target_mu = 0;
    double target_sigma = 0;
    double ans = 0;
    int height = target.rows();
    int width = target.cols();
    int test_height = test.rows();
    int test_width = test.cols();
    double size = test_width * test_height;
    for (int y = 0; y + start_y < height && y < test_height; y++)
    {
        for (int x = 0; x + start_x < width && x < test_width; x++)
        {
            target_mu += target.rgb2gray[y + start_y][x + start_x];
        }
    }
    target_mu /= size;

    for (int y = 0; y + start_y < height && y < test_height; y++)
    {
        for (int x = 0; x + start_x < width && x < test_width; x++)
        {
            target_sigma += pow((target.rgb2gray[y + start_y][x + start_x] - target_mu), 2);
        }
    }
    target_sigma = pow(target_sigma / size, 0.5);
    for (int y = 0; y + start_y < height && y < test_height; y++)
    {
        for (int x = 0; x + start_x < width && x < test_width; x++)
        {
            ans += (target.rgb2gray[y + start_y][x + start_x] - target_mu) * (test.rgb2gray[y][x] - test.mu);
        }
    }
    double delt = size * test.sigma * target_sigma;
    ans = ans / delt;
    return ans;
}

void gaussianBlur(BMP& bmp, int radius, double sigma) {
    int width = bmp.cols();
    int height = bmp.rows();
    int channels = bmp.getBitCount() / 8;

    // ������˹��
    std::vector<std::vector<double>> kernel(2 * radius + 1, std::vector<double>(2 * radius + 1));
    double sum = 0.0;
    for (int i = -radius; i <= radius; i++) {
        for (int j = -radius; j <= radius; j++) {
            kernel[i + radius][j + radius] = exp(-(i * i + j * j) / (2 * sigma * sigma)) / (2 * PI * sigma * sigma);
            sum += kernel[i + radius][j + radius];
        }
    }

    // ��һ����˹��
    for (int i = 0; i < 2 * radius + 1; i++) {
        for (int j = 0; j < 2 * radius + 1; j++) {
            kernel[i][j] /= sum;
        }
    }

    // ������ʱͼ������
    BYTE* tempData = (BYTE*)malloc(width * height * channels);

    // ��ÿ�����ؽ��и�˹ģ������
    for (int y = radius; y < height - radius; y++) {
        for (int x = radius; x < width - radius; x++) {
            for (int c = 0; c < channels; c++) {
                double newValue = 0.0;
                for (int i = -radius; i <= radius; i++) {
                    for (int j = -radius; j <= radius; j++) {
                        newValue += bmp.pColorData[(y + i) * width * channels + (x + j) * channels + c] * kernel[i + radius][j + radius];
                    }
                }
                tempData[y * width * channels + x * channels + c] = static_cast<BYTE>(newValue);
            }
        }
    }

    // �����������ݿ�����ԭʼ����
    memcpy(bmp.pColorData, tempData, width * height * channels);

    // �ͷ���ʱ����
    free(tempData);
}

void match(BMP target_bmp, BMP test_bmp, string p_target, string p_test, double d_height, double d_width, int& start_y, int& start_x)
{
    BMP change_target_bmp;
    string change_targets = "change_" + p_target;

    const char* target = p_target.c_str(), * test = p_test.c_str(), * change_target = change_targets.c_str();

    // ��ȡλͼ������
    BYTE* pColorData = target_bmp.getPColorData();//���ԭͼ��������
    // ��˹ģ��
    gaussianBlur(target_bmp, 2, 3);


    // ԭͼ���� ����С��
    change_target_bmp = p_change(target_bmp, change_target, d_height, d_width);
    int height = change_target_bmp.rows();
    int width = change_target_bmp.cols();
    int test_hight = test_bmp.rows();
    int test_width = test_bmp.cols();
    int target_height = target_bmp.rows();
    int target_width = target_bmp.cols();

    /*ͼ��ƥ�����*/
    double max = -1.0;

    for (int i = 0; i < height - test_hight; i++)
    {

        for (int j = 0; j < width - test_width; j++)
        {

            double NCC = get_NCC(test_bmp, change_target_bmp, i, j);

            if (NCC > max)
            {
                max = NCC;
                start_y = i;
                start_x = j;
            }
        }
    }

    start_y = start_y / d_height;
    start_x = start_x / d_width;

}

void true_message(true_data& tru)
{
    tru.start_x = 533; tru.start_y = 721; tru.height = 159; tru.width = 105;
}

int main()
{
    ofstream outfile("output.txt");
    auto start = high_resolution_clock::now();
    double d_height, d_width;
    string input1 = "input1.bmp";    string input2 = "input2.bmp";
    BMP input1_bmp, input2_bmp;
    imread(input1.c_str(), input1_bmp);    imread(input2.c_str(), input2_bmp);

    true_data tru;  true_message(tru);
    d_height = ((double)16 / tru.height)-0.02;
    d_width = ((double)16 / tru.width)-0.02;

    double size = tru.width * tru.height;
    int start_x = 0, start_y = 0;
    double size_match = 0;
    int final_x = 0, final_y = 0;
    double max = 0;
    outfile << input1 << ":  ";
    double Precision = 0;
    for (int j = 0; j < 3; j++) {
        bool judge = 0;
        for (int k = 0; k < 3; k++) {
            double d = 0.02;
            cout << "��" << j * 3 + k + 1 << "��ƥ��    " ;
            match(input1_bmp, input2_bmp, input1, input2, d_height + d * j, d_width + d * k, start_y, start_x);
            if (abs(start_x - tru.start_x) > tru.width || abs(start_y - tru.start_y) > tru.height)
                Precision = 0;
            else
            {
                size_match = (tru.width - abs(start_x - tru.start_x)) * (tru.height - abs(start_y - tru.start_y));
                Precision = size_match / size;
            }
            cout << Precision ;
            if (Precision > max)
            {
                max = Precision;
                final_x = start_x;
                final_y = start_y;
            }
                     
        }
        cout << endl;
        if (judge)
            break;
        
    }
    if (max == 0)
    {
        cout << "ƥ��ʧ��" << endl;;
        outfile << "ƥ��ʧ��" << endl;
    }
    else
    {
        cout << max << endl;
        cout << "�ɹ�" << endl;
        outfile << "��" << final_x << "��" << final_y << ")     " << max << "    " << max / (2 - max) << endl;
    }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    cout << "����������ʱ��: "
        << duration.count() << " ��" << endl;
    outfile << "����������ʱ��: "
        << duration.count() << " ��" << endl;
    system("pause");

    return 0;
} 