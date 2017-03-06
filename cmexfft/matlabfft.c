
#include "mex.h"
#include<math.h>
#include<windows.h>

#define SIZE 4096
#define pi 3.141592654

int databit;     
int bandbit;   
int tempbit;

typedef struct {
	 LONG64 real;               
	 LONG64 image;          
} Bcomplex;

 int addr0[SIZE];  
 int addr1[SIZE];  
 int addr2[SIZE];  
 int addr3[SIZE];  
 int addr4[SIZE];  
 int addr5[SIZE];  
 int addr6[SIZE];  

Bcomplex DataRom;
Bcomplex DataRom0[16];
Bcomplex DataRom1[64];
Bcomplex DataRom2[256];
Bcomplex DataRom3[1024];
Bcomplex DataRom4[4096];

int Addr0Rom[SIZE];
int Addr1Rom[SIZE];
int Addr2Rom[SIZE];
int Addr3Rom[SIZE];
int Addr4Rom[SIZE];

Bcomplex Idata[SIZE];
Bcomplex Odata[SIZE];

void LevelFftAddr();
void RomAddrGenerate();
void RomDataFromTxt();
void BufflyFFT();
void ComputeData(Bcomplex *x1, Bcomplex *x2, Bcomplex *x3, Bcomplex *x4,
	                           Bcomplex r1, Bcomplex r2, Bcomplex r3, Bcomplex r4,int inbit);


void LevelFftAddr() {

	for ( int i = 0;i < SIZE;++i) {
		addr0[i] = ((int)(i / 4) + ((int)(i % 4) * 1024)) & 0x0fff;
	}

	for ( int i = 0;i < SIZE;++i) {
		int temp;
		int temp1;
		temp = i & 0x0c00;  
		temp1 = i & 0x03ff; 
		addr1[i] = (temp + temp1 / 4 + (temp1 % 4) * 256) & 0x0fff;
	}

	for (int i = 0;i < SIZE;++i) {
		int temp;
		int temp1;
		temp = i & 0x0f00;  
		temp1 = i & 0x00ff;  
		addr2[i] = (temp + temp1 / 4 + (temp1 % 4) * 64) & 0x0fff;
	}

	for (int i = 0;i < SIZE;++i) {
		int temp;
		int temp1;
		temp = i & 0x0fc0;   
		temp1 = i & 0x003f;
		addr3[i] = (temp + temp1 / 4 + (temp1 % 4) * 16) & 0x0fff;
	}

	for (int i = 0;i < SIZE;++i) {
		int temp;
		int temp1;
		temp = i & 0x0ff0;
		temp1 = i & 0x000f;
		addr4[i] = (temp + temp1 / 4 + (temp1 % 4) * 4) & 0x0fff;
	}

	for (int i = 0;i < SIZE;++i) {
		addr5[i] = i & 0x0fff;
	}

	for (int i = 0;i < SIZE;++i) {
		int temp;
		int temp1;
		int temp2;
		int temp3;
		int temp4;
		int temp5;
		temp = i % 4;         
		temp1 = (i / 4) % 4;  
		temp2 = (i / 16) % 4;  
		temp3 = (i / 64) % 4; 
		temp4 = (i / 256) % 4; 
		temp5 = (i / 1024) % 4;
		addr6[i] = (temp * 1024 + temp1 * 256 + temp2 * 64 + temp3 * 16 + temp4 * 4 + temp5) & 0x0fff;
	}
}


void RomDataFromTxt()
{
	int real;
	int imag;
	int k;
	int A = (int)pow(2, bandbit-1) - 1;

	DataRom.real = A;
	DataRom.image = 0;

	k = 0;
	for (int i = 0;i < 4;++i)
	{
		for (int j = 0;j < 4;++j)
		{
			real = (int)(A*cos(-2 * 256 * i*j*pi / SIZE));
			imag = (int)(A*sin(-2 * 256 * i*j*pi / SIZE));
			DataRom0[k].real = real;
			DataRom0[k].image = imag;
			k = k + 1;
		}
	}
	k = 0;
	int num;
	for (int L = 0;L < 4;++L)
	{
		for (int i = 0;i < 4;++i)
		{
			num = L + i * 4;
			for (int j = 0;j < 4;++j)
			{
				real = (int)(A*cos(-2 * 64 * num*j*pi / SIZE));
				imag = (int)(A*sin(-2 * 64 * num*j*pi / SIZE));
				DataRom1[k].real = real;
				DataRom1[k].image = imag;
				k = k + 1;
			}
		}
	}

	k = 0;
	num = 0;
	for (int L1 = 0;L1 < 4;++L1)
	{
		for (int L = 0;L < 4;++L)
		{
			for (int i = 0;i < 4;++i)
			{
				num = L1 + (L + i * 4) * 4;
				for (int j = 0;j < 4;++j)
				{
					real = (int)(A*cos(-2 * 16 * num * j * pi / SIZE));
					imag = (int)(A*sin(-2 * 16 * num*j*pi / SIZE));
					DataRom2[k].real = real;
					DataRom2[k].image = imag;
					k = k + 1;
				}
			}
		}
	}

	k = 0;
	num = 0;
	for (int L2 = 0;L2 < 4;++L2)
	{
		for (int L1 = 0;L1 < 4;++L1)
		{
			for (int L = 0;L < 4;++L)
			{
				for (int i = 0;i < 4;++i)
				{
					num = L2 + (L1 + (L + i * 4) * 4) * 4;
					for (int j = 0;j < 4;++j)
					{
						real = (int)(A*cos(-2 * 4 * num * j * pi / SIZE));
						imag = (int)(A*sin(-2 * 4 * num * j * pi / SIZE));
						DataRom3[k].real = real;
						DataRom3[k].image = imag;
						k = k + 1;
					}
				}
			}
		}
	}

	k = 0;
	num = 0;
	for (int L3 = 0;L3 < 4;++L3)
	{
		for (int L2 = 0;L2 < 4;++L2)
		{
			for (int L1 = 0;L1 < 4;++L1)
			{
				for (int L = 0;L < 4;++L)
				{
					for (int i = 0;i < 4;++i)
					{
						num = L3 + (L2 + (L1 + (L + i * 4) * 4) * 4) * 4;
						for (int j = 0;j < 4;++j)
						{
							real = (int)(A*cos(-2 * 1 * num * j * pi / SIZE));
							imag = (int)(A*sin(-2 * 1 * num * j * pi / SIZE));
							DataRom4[k].real = real;
							DataRom4[k].image = imag;
							k = k + 1;
						}
					}
				}
			}
		}
	}
}


void RomAddrGenerate() {
	for (int i = 0;i < SIZE;++i)
	{
		 int temp;
		 int temp1;
		temp = i / 1024 & 0x00000fff;
		temp1 = i % 4 & 0x00000fff;
		Addr0Rom[i] = temp * 4 + temp1; 
	}

	for (int i = 0;i < SIZE;++i)
	{
		 int temp;
		 int temp1;
		temp = i / 256;
		temp1 = i % 4;
		Addr1Rom[i] = temp * 4 + temp1;
	}

	for (int i = 0;i < SIZE;++i)
	{
		int temp;
		int temp1;
		temp = i / 64;
		temp1 = i % 4;
		Addr2Rom[i] = temp * 4 + temp1; 
	}

	for (int i = 0;i < SIZE;++i) {

		int temp;
		int temp1;
		temp = i / 16;
		temp1 = i % 4;                     
		Addr3Rom[i] = temp * 4 + temp1;                
	}

	for (int i = 0;i < SIZE;++i)                  
	{
		//ad_reg6;
		Addr4Rom[i] = i;
	}
}


void ComputeData(Bcomplex *x1, Bcomplex *x2, Bcomplex *x3, Bcomplex *x4, Bcomplex r1,
	Bcomplex r2, Bcomplex r3, Bcomplex r4,int inbit)
{

	LONG64 brw1pr, biw1pi,
		 crw2pr, ciw2pi,
		 drw3pr, diw3pi,
		 brw1pi, biw1pr,
		 crw2pi, ciw2pr,
		 drw3pi, diw3pr;

	LONG64  mer, mei,
   		 mfr, mfi,
		 mgr, mgi,
		 mhr, mhi;

	LONG64
		er, ei,
		fr, fi,
		gr, gi,
		hr, hi;

	brw1pr = x2->real*r2.real; biw1pi = x2->image*r2.image;
	crw2pr = x3->real*r3.real; ciw2pi = x3->image*r3.image;
	drw3pr = x4->real*r4.real; diw3pi = x4->image*r4.image;
	brw1pi = x2->real*r2.image; biw1pr = x2->image*r2.real;
	crw2pi = x3->real*r3.image; ciw2pr = x3->image*r3.real;
	drw3pi = x4->real*r4.image; diw3pr = x4->image*r4.real;


	mer =   x1->real * pow(2, bandbit - 1)     + brw1pr  - biw1pi   + crw2pr - ciw2pi + drw3pr - diw3pi;
	mei =   x1->image * pow(2, bandbit - 1) + brw1pi  + biw1pr  + crw2pi + ciw2pr + drw3pi + diw3pr;
	mfr =    x1->real * pow(2, bandbit - 1)     + brw1pi  + biw1pr   - crw2pr + ciw2pi - drw3pi - diw3pr; 
	mfi =    x1->image *pow(2, bandbit - 1)   - brw1pr  + biw1pi   - crw2pi - ciw2pr + drw3pr - diw3pi;
	mgr =  x1->real * pow(2, bandbit - 1)       - brw1pr  + biw1pi   + crw2pr - ciw2pi - drw3pr + diw3pi;
	mgi =  x1->image * pow(2, bandbit - 1)   - brw1pi  - biw1pr    + crw2pi + ciw2pr - drw3pi - diw3pr;
	mhr =  x1->real * pow(2, bandbit - 1)       - brw1pi  - biw1pr     - crw2pr + ciw2pi + drw3pi + diw3pr;
	mhi =  x1->image *pow(2, bandbit - 1)   + brw1pr - biw1pi      - crw2pi - ciw2pr - drw3pr + diw3pi;

	LONG64 temp;
	LONG64 temp1;
	LONG64 temp2;
	LONG64 temp3;
	LONG64 temp4;


	temp4 = (int)pow(2, bandbit - 1); 
	temp3 = (int)pow(2, inbit - 1);



	temp2 = mer;
	temp = temp2 / temp4;
	temp1 = temp % 2;
	temp = (temp2 >> bandbit);

	if (temp >= temp3 - 1)
		er = temp3 - 1;
	else if (temp <= -1*temp3)
		er = -1*temp3;
	else
	{
		er = temp;
		if (temp1 == 1) {
			er = er + 1;
		}
	}


	temp2 = mei; 
	temp = mei / temp4;
	temp1 = temp % 2;
	temp = (temp2 >> bandbit);
	if (temp >= temp3 - 1)
		ei = temp3 - 1;
	else if (temp <= -temp3)
		ei = -temp3;
	else
	{
		ei = temp;
		if (temp1 == 1) {
			ei = ei + 1;
		}
	}

	temp2 = mfr;
	temp = temp2 / temp4;
	temp1 = temp % 2;
	temp = (temp2 >> bandbit);
	if (temp >= temp3 - 1)
		fr = temp3 - 1;
	else if (temp <= -temp3)
		fr = -temp3;
	else
	{
		fr = temp;
		if (temp1 == 1) {
			fr = fr + 1;
		}
	}

	temp2 = mfi; 
	temp = temp2 / temp4;
	temp1 = temp % 2;
	temp = (temp2 >> bandbit);
	if (temp >= temp3 - 1)
		fi = temp3 - 1;
	else if (temp <= -temp3)
		fi = -temp3;
	else
	{
		fi = temp;
		if (temp1 == 1) {
			fi = fi + 1;
		}
	}


	temp2 = mgr; 
	temp = temp2 / temp4;
	temp1 = temp % 2;
	temp = (temp2 >> bandbit);
	if (temp >= temp3 - 1)
		gr = temp3 - 1;
	else if (temp <= -temp3)
		gr = -temp3;
	else
	{
		gr = temp;
		if (temp1 == 1) {
			gr = gr + 1;
		}
	}
	

	temp2 = mgi; 
	temp = temp2 / temp4;
	temp1 = temp % 2;
	temp = (temp2 >> bandbit);

	if (temp >= temp3 - 1)
		gi = temp3 - 1;
	else if (temp <= -temp3)
		gi = -temp3;
	else
	{
		gi = temp;
		if (temp1 == 1) {
			gi = gi + 1;
		}
	}

	temp2 = mhr; 
	temp = temp2 / temp4;
	temp1 = temp % 2;
	temp = (temp2 >> bandbit);
	if (temp >= temp3 - 1)
		hr = temp3 - 1;
	else if (temp <= -temp3)
		hr = -temp3;
	else
	{
		hr = temp;
		if (temp1 == 1) {
			hr = hr + 1;
		}
	}


	temp2 = mhi; 
	temp = temp2 / temp4;
	temp1 = temp % 2;
	temp = (temp2 >> bandbit);

	if (temp >= temp3 - 1)
		hi = temp3 - 1;
	else if (temp <= -temp3)
		hi = -temp3;
	else
	{
		hi = temp;
		if (temp1 == 1) {
			hi = hi + 1;
		}
	}
	x1->real = er;
	x1->image = ei;

	x2->real = fr;
	x2->image = fi;

	x3->real = gr;
	x3->image = gi;

	x4->real = hr;
	x4->image = hi;
}

void BufflyFFT() {

	for (int i = 0;i < SIZE;i = i + 4) {
		Bcomplex x1, x2, x3, x4;
		Bcomplex r1, r2, r3, r4;
		x1 = Idata[addr0[i + 0]];
		x2 = Idata[addr0[i + 1]];
		x3 = Idata[addr0[i + 2]];
		x4 = Idata[addr0[i + 3]];

		r1 = DataRom;
		r2 = DataRom;
		r3 = DataRom;
		r4 = DataRom;

		ComputeData(&x1, &x2, &x3, &x4, r1, r2, r3, r4,tempbit);
		Idata[addr0[i + 0]] = x1;
		Idata[addr0[i + 1]] = x2;
		Idata[addr0[i + 2]] = x3;
		Idata[addr0[i + 3]] = x4;
	}

	for (int i = 0;i < SIZE;i = i + 4) {
		Bcomplex x1, x2, x3, x4;
		Bcomplex r1, r2, r3, r4;
		x1 = Idata[addr1[i + 0]];
		x2 = Idata[addr1[i + 1]];
		x3 = Idata[addr1[i + 2]];
		x4 = Idata[addr1[i + 3]];

		r1 = DataRom0[Addr0Rom[i + 0]];
		r2 = DataRom0[Addr0Rom[i + 1]];
		r3 = DataRom0[Addr0Rom[i + 2]];
		r4 = DataRom0[Addr0Rom[i + 3]];

		ComputeData(&x1, &x2, &x3, &x4, r1, r2, r3, r4,tempbit);
		Idata[addr1[i + 0]] = x1;
		Idata[addr1[i + 1]] = x2;
		Idata[addr1[i + 2]] = x3;
		Idata[addr1[i + 3]] = x4;
	}

	for (int i = 0;i < SIZE;i = i + 4) {
		Bcomplex x1, x2, x3, x4;
		Bcomplex r1, r2, r3, r4;
		x1 = Idata[addr2[i + 0]];
		x2 = Idata[addr2[i + 1]];
		x3 = Idata[addr2[i + 2]];
		x4 = Idata[addr2[i + 3]];

		r1 = DataRom1[Addr1Rom[i + 0]];
		r2 = DataRom1[Addr1Rom[i + 1]];
		r3 = DataRom1[Addr1Rom[i + 2]];
		r4 = DataRom1[Addr1Rom[i + 3]];

		ComputeData(&x1, &x2, &x3, &x4, r1, r2, r3, r4,tempbit);
		Idata[addr2[i + 0]] = x1;
		Idata[addr2[i + 1]] = x2;
		Idata[addr2[i + 2]] = x3;
		Idata[addr2[i + 3]] = x4;
	}

	for (int i = 0;i < SIZE;i = i + 4) {
		Bcomplex x1, x2, x3, x4;
		Bcomplex r1, r2, r3, r4;
		x1 = Idata[addr3[i + 0]];
		x2 = Idata[addr3[i + 1]];
		x3 = Idata[addr3[i + 2]];
		x4 = Idata[addr3[i + 3]];

		r1 = DataRom2[Addr2Rom[i + 0]];
		r2 = DataRom2[Addr2Rom[i + 1]];
		r3 = DataRom2[Addr2Rom[i + 2]];
		r4 = DataRom2[Addr2Rom[i + 3]];

		ComputeData(&x1, &x2, &x3, &x4, r1, r2, r3, r4,tempbit);
		Idata[addr3[i + 0]] = x1;
		Idata[addr3[i + 1]] = x2;
		Idata[addr3[i + 2]] = x3;
		Idata[addr3[i + 3]] = x4;
	}

	for (int i = 0;i < SIZE;i = i + 4) {
		Bcomplex x1, x2, x3, x4;
		Bcomplex r1, r2, r3, r4;
		x1 = Idata[addr4[i + 0]];
		x2 = Idata[addr4[i + 1]];
		x3 = Idata[addr4[i + 2]];
		x4 = Idata[addr4[i + 3]];

		r1 = DataRom3[Addr3Rom[i + 0]];
		r2 = DataRom3[Addr3Rom[i + 1]];
		r3 = DataRom3[Addr3Rom[i + 2]];
		r4 = DataRom3[Addr3Rom[i + 3]];

		ComputeData(&x1, &x2, &x3, &x4, r1, r2, r3, r4,tempbit);
		Idata[addr4[i + 0]] = x1;
		Idata[addr4[i + 1]] = x2;
		Idata[addr4[i + 2]] = x3;
		Idata[addr4[i + 3]] = x4;
	}
	for (int i = 0;i < SIZE;i = i + 4) {
		Bcomplex x1, x2, x3, x4;
		Bcomplex r1, r2, r3, r4;

		x1 = Idata[addr5[i + 0]];
		x2 = Idata[addr5[i + 1]];
		x3 = Idata[addr5[i + 2]];
		x4 = Idata[addr5[i + 3]];

		r1 = DataRom4[Addr4Rom[i + 0]];
		r2 = DataRom4[Addr4Rom[i + 1]];
		r3 = DataRom4[Addr4Rom[i + 2]];
		r4 = DataRom4[Addr4Rom[i + 3]];

		ComputeData(&x1, &x2, &x3, &x4, r1, r2, r3, r4,databit);
		Idata[addr5[i + 0]] = x1;
		Idata[addr5[i + 1]] = x2;
		Idata[addr5[i + 2]] = x3;
		Idata[addr5[i + 3]] = x4;
	}

	for (int i = 0;i < SIZE;i = i + 4) {
		Odata[i + 0] = Idata[addr6[i + 0]];
		Odata[i + 1] = Idata[addr6[i + 1]];
		Odata[i + 2] = Idata[addr6[i + 2]];
		Odata[i + 3] = Idata[addr6[i + 3]];
	}
}


void fftport(double *data_in,double *data_out)
{
    double temp[8192];
    bandbit = (int)data_in[0];
    databit =  (int)data_in[1];
    tempbit =(int)data_in[2];
    for (int i = 0; i<SIZE; ++i)
    {
        Idata[i].real = (long long)data_in[3+i];
        Idata[i].image = (long long)data_in[3+i+4096];
    }
    LevelFftAddr();
    RomDataFromTxt();
    RomAddrGenerate();
    BufflyFFT();
    for (int i = 0; i < SIZE; ++i) {
        temp[i] = (double)Odata[i].real;
        temp[i+4096] = (double)Odata[i].image;
    }
    
    for (int i = 0; i < SIZE; ++i) {
        data_out[i] = temp[i];
        data_out[i+4096] = temp[i+4096];
    }
}

void mexFunction( int nlhs, mxArray *plhs[],int nrhs,  mxArray *prhs[])
{
    double *dataCursor;
    double *dataoutCursor;
    plhs[0] = mxCreateDoubleMatrix(8192,1,mxREAL);
    
    dataCursor = mxGetPr(prhs[0]);
    dataoutCursor = mxGetPr(plhs[0]);
    
    fftport(dataCursor,dataoutCursor);
    
}





