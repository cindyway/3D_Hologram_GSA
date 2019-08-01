
//Reconstruction of holograms using Fourier Transform

//Read object from bmp image

#include "stdafx.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <conio.h>
#include <math.h>
#include <time.h>
#include "fft.h"

#define pi 3.14159265

extern void twiddle();
extern void FFT();
extern void IFFT();
extern void FFT_SHIFT();

//Size of the image
int ROW,COL;

//Size of the hologram
int H_ROW=1024;
int H_COL=1024;
int index;

double sdata[2048][2048],cdata[2048][2048];
double ht_r[2048][2048], ht_i[2048][2048], h_r[2048][2048], h_i[2048][2048], ph[2048][2048], s[2048][2048], s_r[2048][2048], s_i[2048][2048], new_r1[2048][2048];
double dx,dy, max_p, min_p;

unsigned char temp[10000],temp1[10000];

char outfile[30], ifile[30];	
								
char hfile[30],ofile[30];
unsigned char header[256];

double lamda,dam_pix_no;
double laser_angle=0.28;

FILE *in, *out, *head, *hptr, *rout, *src;


typedef struct pt
	{
		double x;
		double y;
		double z;
	} pt;

typedef struct cp
	{
		double real;
		double img;
	} cp;

struct cp fn[2048][2048], Q[2048][2048], t_image[2048][2048], B[2048];

//Variables for fft 
struct cmpx P[5000],w[5000];
int N;

void init();
void downsample(int DS);
void ref_gen();
void Transform_image();
void fn_gen(struct pt mp);
void holo_gen();
void phase_only();
void holo_recon();
void normalize();
void save_image();


void init()
{
	int x,y,z;

	//Load output filename
	printf("\nInput output image file ");
	scanf("%s",outfile);
	strcat(outfile,".bmp");
	
	//Load input image
	printf("\nInput image file ");
	scanf("%s",ifile);
	strcat(ifile,".bmp");
	if((src=fopen(ifile,"rb"))==NULL) 
	{
		printf("\nInput image file not found");
		exit(1);
	}
	fread(header, sizeof(char),54,src);
	ROW=header[23];ROW=ROW*256+header[22];
	COL=header[19];COL=COL*256+header[18];

	z=COL*3;
	for(x=ROW-1;x>=0;x--)
	{
		fread(temp,sizeof(char),z,src);
		for(y=0;y<z;y+=3)
		{
			sdata[x][y/3]=temp[y];	//Normal reading of hologram pixel - no change
		}
	}
	fclose(src);


	//Clear hologram buffer
	for(x=0;x<H_ROW;x++)
		for(y=0;y<H_COL;y++)
			h_r[x][y]=h_i[x][y]=0;
	//clear new buffer
	

}


void downsample(int DS)
{
	int i,j;
	
	
}

void ref_gen()
{
	int i;
	double y,wn,stp,d;
	
	wn=2*pi/lamda;
	laser_angle=laser_angle*pi/180;

	d=2*pi; 
	stp=d/8; 
	for(i=0;i<=H_ROW;i++)
	{
		y=(double)i*dy;					
		
		B[i].real=cos(wn*y*sin(laser_angle));		
		B[i].img=sin(wn*y*sin(laser_angle));
	}

}

void Transform_image()
{
	int i,j;
	int h_row,h_col,i_row,i_col;
	
	h_row=H_ROW/2;h_col=H_COL/2;
	i_row=ROW/2;i_col=COL/2;

	for(i=0;i<H_ROW;i++)
		for(j=0;j<H_COL;j++)
			t_image[i][j].real =t_image[i][j].img =0;
	
	for(i=h_row-i_row;i<h_row+i_row;i++)
		for(j=h_col-i_col;j<h_col+i_col;j++)
		{
			t_image[i][j].real =cdata[i-h_row+i_row][j-h_col+i_col];
			t_image[i][j].img =0;
		}
			

	N=H_COL; twiddle();
	for(i=0;i<H_ROW;i++)
	{
		for(j=0;j<H_COL;j++)
		{
			P[j].real = t_image[i][j].real;
			P[j].imag = t_image[i][j].img;
		}
		FFT();
		for(j=0;j<H_COL;j++)
		{
			Q[i][j].real = P[j].real;
			Q[i][j].img = P[j].imag;
		}
	}
	N=H_ROW; twiddle();
	for(j=0;j<H_COL;j++)
	{
		for(i=0;i<H_ROW;i++)
		{
			P[i].real = Q[i][j].real;
			P[i].imag = Q[i][j].img;
		}
		FFT();
		for(i=0;i<H_ROW;i++)
		{
			t_image[i][j].real = P[i].real;
			t_image[i][j].img = P[i].imag;
		}
	}
}
void Transform_image2()
{
	int i, j;
	int h_row, h_col, i_row, i_col;

	h_row = H_ROW / 2; h_col = H_COL / 2;
	i_row = ROW / 2; i_col = COL / 2;

	for (i = 0; i < H_ROW; i++)
		for (j = 0; j < H_COL; j++)
			t_image[i][j].real = t_image[i][j].img = 0;

	for (i = h_row - i_row; i < h_row + i_row; i++)
		for (j = h_col - i_col; j < h_col + i_col; j++)
		{
			t_image[i][j].real = s_r[i][j];
			t_image[i][j].img =  s_i[i][j];
		}


	N = H_COL; twiddle();
	for (i = 0; i < H_ROW; i++)
	{
		for (j = 0; j < H_COL; j++)
		{
			P[j].real = t_image[i][j].real;
			P[j].imag = t_image[i][j].img;
		}
		FFT();
		for (j = 0; j < H_COL; j++)
		{
			Q[i][j].real = P[j].real;
			Q[i][j].img = P[j].imag;
		}
	}
	N = H_ROW; twiddle();
	for (j = 0; j < H_COL; j++)
	{
		for (i = 0; i < H_ROW; i++)
		{
			P[i].real = Q[i][j].real;
			P[i].imag = Q[i][j].img;
		}
		FFT();
		for (i = 0; i < H_ROW; i++)
		{
			t_image[i][j].real = P[i].real;
			t_image[i][j].img = P[i].imag;
		}
	}
}

//Generate a Fresnel zone plate corresponding to a point source at x=y=0
void fn_gen(struct pt mp)	//only mp.z is used 
{
	int i,j;
	double d1,m,n,z,wn,h_row,h_col;

	wn=2*pi/lamda;
	z=mp.z*mp.z;

	for(i=0;i<H_ROW;i++)
		for(j=0;j<H_COL;j++)
			fn[i][j].real=fn[i][j].img=0;

	h_row=H_ROW/2;h_col=H_COL/2;
	for(i=0;i<H_ROW;i++)
	{
		m=i-h_row;m*=dx;	//Define orgin to be the center of the image array
		for(j=0;j<H_COL;j++)
		{
			n=j-h_col;n*=dy;	//Define orgin to be the center of the image array
			d1=sqrt(m*m+n*n+z);
			fn[i][j].real=cos(d1*wn)/d1;
			fn[i][j].img=sin(d1*wn)/d1;
		}
	}
	
	
	N=H_COL; twiddle();
	for(i=0;i<H_ROW;i++)
	{
		for(j=0;j<H_COL;j++)
		{
			P[j].real = fn[i][j].real;
			P[j].imag = fn[i][j].img;
		}
		FFT();
		for(j=0;j<H_COL;j++)
		{
			Q[i][j].real = P[j].real;
			Q[i][j].img = P[j].imag;
		}
	}
	N=H_ROW; twiddle();
	for(j=0;j<H_COL;j++)
	{
		for(i=0;i<H_ROW;i++)
		{
			P[i].real = Q[i][j].real;
			P[i].imag = Q[i][j].img;
		}
		FFT();
		for(i=0;i<H_COL;i++)
		{
			fn[i][j].real = P[i].real;
			fn[i][j].img = P[i].imag;
		}
	}
}

void holo_gen()
{
	int i,j;
	
	for(i=0;i<H_ROW;i++)
	{
		for(j=0;j<H_COL;j++)
		{
			ht_r[i][j]=fn[i][j].real *t_image[i][j].real -fn[i][j].img *t_image[i][j].img;
			ht_i[i][j]=fn[i][j].real *t_image[i][j].img  +fn[i][j].img *t_image[i][j].real;
		}
	}
	N=H_COL; twiddle();
	for(i=0;i<H_ROW;i++)
	{
		for(j=0;j<H_COL;j++)
		{
			P[j].real = ht_r[i][j];
			P[j].imag = ht_i[i][j];
		}
		IFFT();FFT_SHIFT();
		for(j=0;j<H_COL;j++)
		{
			Q[i][j].real = P[j].real;
			Q[i][j].img = P[j].imag;
		}
	}
	N=H_ROW; twiddle();
	for(j=0;j<H_COL;j++)
	{
		for(i=0;i<H_ROW;i++)
		{
			P[i].real = Q[i][j].real;
			P[i].imag = Q[i][j].img;
		}
		IFFT();FFT_SHIFT();
		for(i=0;i<H_COL;i++)
		{
			ht_r[i][j] = P[i].real;
			ht_i[i][j] = P[i].imag;
		}
	}

	for(i=0;i<H_ROW;i++)
	{
		for(j=0;j<H_COL;j++)
		{
			h_r[i][j]+=ht_r[i][j];			//Accumulate h_r[][] with sub-hologram at depth=mp.z
			h_i[i][j]+=ht_i[i][j];
		}
	}
}


void phase_only()
{
	int i, j;
	double r1, q1, q2, q3, q4, mag;

	q1 = pi / 2;
	q2 = pi;
	q3 = pi / 2 + pi;
	q4 = 2 * pi;


	for (i = 0; i < H_ROW; i++)
	{
		for (j = 0; j < H_COL; j++)
		{
			mag = h_r[i][j] * h_r[i][j] + h_i[i][j] * h_i[i][j];
			if (mag > 0) mag = sqrt(mag);
			if (h_r[i][j] == 0)
			{
				if (h_i[i][j] > 0) r1 = q1;
				if (h_i[i][j] < 0) r1 = q3;
			}
			else
			{
				r1 = h_i[i][j] / h_r[i][j];
				if (r1 < 0)r1 = -r1;
				r1 = atan(r1);
				if ((h_r[i][j] < 0) && (h_i[i][j] >= 0)) r1 = pi - r1;
				if ((h_r[i][j] < 0) && (h_i[i][j] < 0)) r1 += pi;
				if ((h_r[i][j] > 0) && (h_i[i][j] < 0)) r1 = q4 - r1;

			}
			h_r[i][j] = cos(r1);
			h_i[i][j] = sin(r1);
			ph[i][j] = r1 / q4;	//Normalize to [0,1] and store phase of hologram in array ph[][]
		}
	}
}
void holo_recon()
{
	int i,j;
		
	N=H_COL; twiddle();
	for(i=0;i<H_ROW;i++)
	{
		for(j=0;j<H_COL;j++)
		{
			P[j].real = h_r[i][j];
			P[j].imag = h_i[i][j];
		}

		FFT();
		for(j=0;j<H_COL;j++)
		{
			Q[i][j].real = P[j].real;
			Q[i][j].img = P[j].imag;
		}
	}
	N=H_ROW; twiddle();
	for(j=0;j<H_COL;j++)
	{
		for(i=0;i<H_ROW;i++)
		{
			P[i].real = Q[i][j].real;
			P[i].imag = Q[i][j].img;
		}
		FFT();
		for(i=0;i<H_ROW;i++)
		{
			h_r[i][j] = P[i].real;
			h_i[i][j] = P[i].imag;
		}
	}

	for(i=0;i<H_ROW;i++)
	{
		for(j=0;j<H_COL;j++)
		{
			ht_r[i][j]=fn[i][j].real *h_r[i][j] +fn[i][j].img *h_i[i][j];
			ht_i[i][j]=fn[i][j].real *h_i[i][j] -fn[i][j].img *h_r[i][j];
		}
	}

	

	N=H_COL; twiddle();
	for(i=0;i<H_ROW;i++)
	{
		for(j=0;j<H_COL;j++)
		{
			P[j].real = ht_r[i][j];
			P[j].imag = ht_i[i][j];
		}
		IFFT();FFT_SHIFT();
		for(j=0;j<H_COL;j++)
		{
			Q[i][j].real = P[j].real;
			Q[i][j].img = P[j].imag;
		}
	}
	N=H_ROW; twiddle();
	for(j=0;j<H_COL;j++)
	{
		for(i=0;i<H_ROW;i++)
		{
			P[i].real = Q[i][j].real;
			P[i].imag = Q[i][j].img;
		}
		IFFT();FFT_SHIFT();
		for(i=0;i<H_COL;i++)
		{
			ht_r[i][j] = P[i].real;
			ht_i[i][j] = P[i].imag;
		}
	}

	for(i=0;i<H_ROW;i++)
	{
		for(j=0;j<H_COL;j++)
		{
			h_r[i][j]=ht_r[i][j]*ht_r[i][j]+ht_i[i][j]*ht_i[i][j];			//Accumulate h_r[][] with sub-hologram at depth=mp.z
			if(h_r[i][j]>0) h_r[i][j]=sqrt(h_r[i][j]);
		}
	}
}







//Hologram is bipolar, add an offset to make all values positive, then normalize to
//the range [0,255]
void normalize()
{
	int i,j;
	double range, mean;


	//Normalize Real part
	max_p=-1e9;
	min_p=1e9;
	mean=0;

	for(i=0;i<H_ROW;i++)
		for(j=0;j<H_COL;j++)
		{
			if(h_r[i][j]>max_p) max_p=h_r[i][j];
			if(h_r[i][j]<min_p) min_p=h_r[i][j];
		}
	
	range=max_p-min_p;

	for(i=0;i<H_ROW;i++)
		for(j=0;j<H_COL;j++)
			h_r[i][j]=255*(h_r[i][j]-min_p)/range;
	printf("\nMax=%g, Min=%g",max_p,min_p);			

}



//This procedure save the hologram image[][] to a bmp file for displaying the source image
void save_image()
{
	int row,col,j,k;

	//Create output file with size = Hologram dimension
	header[22]=H_ROW%256;header[23]=H_ROW/256;
	header[18]=H_COL%256;header[19]=H_COL/256;
	
	if((out=fopen(outfile,"wb"))==NULL) exit(1);
	fwrite(header,sizeof(char),54,out);

	
	printf("\nHologram size: Row=%d, Column=%d",H_ROW,H_COL);
	k=COL*3;

	normalize();	//normalize h_r[i][j] to range [0,255]
	
	//Convert output image to bmp format and stored in outfile
	for(row=H_ROW-1;row>=0;row--)
	{
		for(col=0;col<H_COL;col++)
		{
			j=col*3;
			temp[j]=temp[j+1]=temp[j+2]=h_r[row][col];
		}
		fwrite(temp,sizeof(char),H_COL*3,out);
	}
	fclose(out);
	printf("\nSource image successfully loaded and saved in bmp format");

}


int main(int argc, char* argv[])
{
	struct pt mp;
	int i, j, dist_L, dist_R;
	double q4;
	q4 = 2 * pi;

	dx = dy = 8.1e-6;	//Pixel size of the hologram
	lamda = 633e-9;	//Set wavelength of reference beam

	printf("\nPixel size = %g", dx);
	printf("\nWavelength = %g\n", lamda);

	//Load object intensity image
	init();

	downsample(12);	//Downsample source image

	printf("\nImage size:    ROW=%d, COL=%d", ROW, COL);
	printf("\nHologram size: H_ROW=%d, H_COL=%d", H_ROW, H_COL);

	printf("\nInput distance between left side of image and hologram (mm) (suggest 500)");
	scanf("%d", &dist_L);
	printf("\nInput distance between right side of image and hologram (mm) (suggest 550)");
	scanf("%d", &dist_R);



	//Generate Double depth image hologram (left and right side)

		//Hologram of the left side of the image
	mp.z = (double)dist_L / 1000;
	fn_gen(mp);
	for (i = 0; i < ROW; i++)
		for (j = 0; j < COL; j++)
			if (j < COL / 2) cdata[i][j] = sdata[i][j];
			else cdata[i][j] = 0;
	Transform_image();
	holo_gen();

	//Hologram of the right side of the image
	mp.z = (double)dist_R / 1000;
	fn_gen(mp);
	for (i = 0; i < ROW; i++)
		for (j = 0; j < COL; j++)
			if (j >= COL / 2) cdata[i][j] = sdata[i][j];
			else cdata[i][j] = 0;
	Transform_image();
	holo_gen();


	//Retain only the phase component
	phase_only();

	//After generating the hologram, reconstruct the image at a given focal plane

	//Input the distance of the focal plane
	printf("\nInput reconstruction distance (mm) ");
	scanf("%d", &dist_L);
	mp.z = (double)dist_L / 1000;
	fn_gen(mp);

	holo_recon();	//Reconstruct the hologram at the focal plane 
	

	
	for (int k = 0; k<5; k++) {

		int h_row, h_col, i_row, i_col;

		h_row = H_ROW / 2; h_col = H_COL / 2;
		i_row = ROW / 2; i_col = COL / 2;
		for (i = h_row - i_row; i < h_row + i_row; i++)
			for (j = h_col - i_col; j < h_col + i_col; j++)
			{
				s[i][j] = sdata[i - h_row + i_row][j - h_col + i_col];
				new_r1[i][j] = ht_i[i][j] / ht_r[i][j];
				if (new_r1[i][j] < 0) new_r1[i][j] = -new_r1[i][j];
				new_r1[i][j] = atan(new_r1[i][j]);
				s_r[i][j] = s[i][j] * cos(new_r1[i][j]);
				s_i[i][j] = s[i][j] * sin(new_r1[i][j]);
			}

		//Hologram  of the image
		mp.z = 0.525;
		fn_gen(mp);

		Transform_image2();
		holo_gen();


		//Retain only the phase component
		phase_only();


		holo_recon();

	}
	int h_row, h_col, i_row, i_col;

	h_row = H_ROW / 2; h_col = H_COL / 2;
	i_row = ROW / 2; i_col = COL / 2;
	for (i = h_row - i_row; i < h_row + i_row; i++)
		for (j = h_col - i_col; j < h_col + i_col; j++)
		{
			s[i][j] = sdata[i - h_row + i_row][j - h_col + i_col];
			new_r1[i][j] = ht_i[i][j] / ht_r[i][j];
			if (new_r1[i][j] < 0) new_r1[i][j] = -new_r1[i][j];
			new_r1[i][j] = atan(new_r1[i][j]);
			s_r[i][j] = s[i][j] * cos(new_r1[i][j]);
			s_i[i][j] = s[i][j] * sin(new_r1[i][j]);
		}
	
	mp.z = 0.525;
	fn_gen(mp);
	Transform_image2();
	holo_gen();
	phase_only();


	printf("\nImage reconstructed from POH \n");

	save_image();	//Save reconstructed image
	
	printf("\n\n------ Process completed, press any key to continue --------\n");
	getche();
	return 0;
}






