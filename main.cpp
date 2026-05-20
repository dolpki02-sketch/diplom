#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h>

#include "include/fda_poliphase_short.h"

#include <chrono>
#include <thread>
using namespace std::chrono;
using namespace std;

const float pi = 3.1415926;
const int N_Q = 80000;
int N_Q_corr = 2 * N_Q;
int Nfft_Q = 1048576;
int i_nois;
//int ppt;
float nois;
float corr_pick;
//float corr_pick_;
float frequency_position[33];
float corr_position[33];
unsigned int rad;
float f_car_est_0;
int cor_Max_search;
//
float Fs, Fd;
//
int N_samples;
//
short in_poli_i[BL >> 1];
short in_poli_q[BL >> 1];
//
short I_buf_s[8]; // resampler buffer I
short Q_buf_s[8]; // resampler buffer Q
// Receiver CIC filter I
int64_t f1_inti1, f1_inti2, f1_inti3, f1_inti4, f1_combi1, f1_combi2, f1_combi3, f1_combi4, f1_last_combi1, f1_last_combi2,
f1_last_combi3, f1_last_combi4;
// Receiver CIC filter Q
int64_t f1_intq1, f1_intq2, f1_intq3, f1_intq4, f1_combq1, f1_combq2, f1_combq3, f1_combq4, f1_last_combq1, f1_last_combq2,
f1_last_combq3, f1_last_combq4;
short* s_short = (short*)malloc(4480000 * sizeof(short));
float* corr_Q = (float*)malloc(N_Q_corr * sizeof(float));
float* fft_src_q = (float*)malloc(2 * Nfft_Q * sizeof(float));
float* fft_corr_q = (float*)malloc(2 * Nfft_Q * sizeof(float));
float* corr_q = (float*)malloc(Nfft_Q * sizeof(float));
float* w_f = (float*)malloc(2 * Nfft_Q * sizeof(float));
float* fft_ref_q = (float*)malloc(2 * Nfft_Q * sizeof(float));
float* Q_corr = (float*)malloc(N_Q * sizeof(float));
float* Q_corr_2 = (float*)malloc(2 * N_Q * sizeof(float));
short* Q_corr_2_s = (short*)malloc(2 * N_Q * sizeof(short));
bool* Q_corr_2_ = (bool*)malloc(2 * N_Q * sizeof(bool));
float* fft_res_ref_q = (float*)malloc(2 * Nfft_Q * sizeof(float));
float* ref0_q = (float*)malloc(2 * Nfft_Q * sizeof(float));
float* fft_in_q = (float*)malloc(2 * Nfft_Q * sizeof(float));
float* fft_res_q = (float*)malloc(2 * Nfft_Q * sizeof(float));
float* tmp_x2 = (float*)malloc(4 * N_Q_corr * sizeof(float));
float* tmp_y2 = (float*)malloc(4 * N_Q_corr * sizeof(float));
float* fft_src = (float*)malloc(2 * Nfft_Q * sizeof(float));
float* res_koef = (float*)malloc(5000 * sizeof(float));
//
const int L_degree = 8;
float Lagr_[8];
float Lagr[8][1024];
short Lagr_RX_short[8][1024];
//
uint32_t n_kp_0; // int resampler coef
//
short re_f_c[1024];
short im_f_c[1024];
//
unsigned int getRadix(unsigned long);
void tw_gen(float*, int);
void DSPF_sp_fftSPxSP_cn(int, float*, float*, float*,/* unsigned char *brev*/int, int, int);
void DSPF_sp_ifftSPxSP_cn(int, float*, float*, float*, /*unsigned char *brev,*/ int, int, int);
float samples_cic_fs_s(FILE*, float*, int, int, float, float, float, float);
int file_to_buff_converter(FILE*, short*, float);
int get_samples_(short*, int, float, float, float, float*, float*);
int get_samples_1(short*, int, float, float, float, float*, float*);
int search_short_fft_Q(float*, float*, int, int);

typedef struct {
	int corr_Maximorum;
	float max_corr;
} SearchRes;

void part_search(bool k, SearchRes* result);

//
void calculate_lagrange_coeff(const float x, float* lagrange_coeff, const unsigned int lagrange_degree)
{
	int k, m;
	const int w = lagrange_degree / 2;

	const int offset = (lagrange_degree & 1) ^ 1;

	int n = 0;

	for (k = -w + offset; k <= w; k++)
	{
		float lambda = 1.0;

		for (m = -w + offset; m <= w; m++)
		{
			if (k != m)
			{
				const float numerator = (x - m);
				const float denominator = (k - m);

				lambda *= numerator;
				lambda /= denominator;
			}
		}

		lagrange_coeff[n++] = lambda;
	}
}
//

int main()
{
	srand(time(NULL));
	FILE* fPtr;
	FILE* f_float;
	FILE* f_short_I, * f_short_Q;
	char str1[5];
	// float amp_norm[91];
	 //float frequency_position[33];
	 //float corr_position[33];
	int k, i;
	short N_sym, z_s, i_freq;
	int corr_Max;
	float car_0;
	float x;
	float car_shift, car_shift_Hz;
	int N_;
	for (k = 0; k < 1024; k++)
	{
		x = (float)k / (float)1024;
		calculate_lagrange_coeff(x, Lagr_, L_degree);
		for (i = 0; i < L_degree; i++)
		{
			Lagr[i][k] = Lagr_[L_degree - 1 - i];
			Lagr_RX_short[i][k] = (short)(Lagr[i][k] * 32767);
		}
	}
	fPtr = fopen("/home/root/projects/gska_search_parralel/81840.txt", "r");
	if (fPtr == NULL)
	{
		printf("File %s Open failed\n", "81840.txt");
		exit(2);
	}
	for (i = 0; i < N_Q; i++)
	{
		fgets(str1, 7, fPtr);
		if (str1[0] == 116) x = 1;
		if (str1[0] == 102) x = -1;
		Q_corr[i] = x;
		for (k = 0; k < 2; k++)
		{
			Q_corr_2[i * 2 + k] = x;
			Q_corr_2_s[i * 2 + k] = (short)x;
			if (x == 1)  Q_corr_2_[i * 2 + k] = 1; else Q_corr_2_[i * 2 + k] = 0;
		}
	}
	fclose(fPtr);
	for (i = 0; i < Nfft_Q; i++)
	{
		if (i < N_Q_corr)
		{
			fft_ref_q[2 * i] = Q_corr_2[i];
			fft_ref_q[2 * i + 1] = Q_corr_2[i];
		}
		else
		{
			fft_ref_q[2 * i] = 0;
			fft_ref_q[2 * i + 1] = 0;
		}
	}
	Fd = 280.0;// 280 MHz
	Fs = 20.0; // 20 Mchip
	n_kp_0 = (uint32_t)((4 * Fs / Fd) * pow(2, 32));// 2*Fs/Fd/2
	rad = getRadix(Nfft_Q);
	tw_gen(w_f, Nfft_Q);
	DSPF_sp_fftSPxSP_cn(Nfft_Q, fft_ref_q, w_f, fft_res_ref_q, rad, 0, Nfft_Q);
	//Reference sequence - complex conjugate of ref samples FFT
	for (i = 0; i < Nfft_Q; i++)
	{
		ref0_q[2 * i] = -fft_res_ref_q[2 * i];
		ref0_q[2 * i + 1] = fft_res_ref_q[2 * i + 1];
	}
	for (i = 0; i < 1024; i++)
	{
		re_f_c[i] = (short)(32767 * sin(2 * pi * (float)i / 1024.0));//*sin(_2pi*(float)i/1024.0));
		im_f_c[i] = (short)(32767 * cos(2 * pi * (float)i / 1024.0));
	}
	i_nois = 90; // set nois
	if (i_nois < 0 || i_nois>90)
	{
		printf("%s is out of range\n", "i");
		exit(2);
	}
	nois = (float)i_nois - 40.0; // start -40 dB
	printf("%f\n", nois);
	N_sym = 4; //symbols for search
	if (N_sym < 4)
	{
		printf("%s value is less than 4. It's out of range!!\n", "N_sym");
		exit(2);
	}
	car_0 = ((float)rand() - (float)(RAND_MAX >> 1)) / ((float)(RAND_MAX)) * 20.0e+6;// car_0 in Hz +/- 20 MHz
	//car_0 = 20.0e+06; // 10 MHz
	//car_0 = 0.0e+06;
	car_shift_Hz = ((float)rand() - (float)(RAND_MAX >> 1)) / ((float)(RAND_MAX)) * 2000.0;// car_shift in Hz +/- 2 kHz
	car_shift_Hz = 100;
	car_0 /= 1.0e+6; // car_0 in MHzcar_shift_Hz = ((float)rand()-(float)(RAND_MAX>>1))/((float)(RAND_MAX))*2000.0;// car_shift in Hz +/- 2 kHz
	car_shift = car_shift_Hz / 1.0e+06; // car_shift in MHz
	//ppt = 0;
	N_ = 0;
	f_float = fopen("f_float.bin", "wb");
	f_short_I = fopen("f_short_I.bin,", "wb");
	f_short_Q = fopen("f_short_Q.bin,", "wb");
	x = samples_cic_fs_s(f_float, Q_corr_2, N_Q_corr, N_sym, car_0 + car_shift, Fs, Fd, nois);
	fclose(f_float);
	f_float = fopen("f_float.bin", "rb");
	N_samples = file_to_buff_converter(f_float, s_short, x);
	fclose(f_float);
	i = get_samples_1(s_short, N_samples, car_0, Fs, Fd, tmp_x2, tmp_y2);
	for (i = 0; i < Nfft_Q; i++)
	{
		if (i < 4 * N_Q_corr) // 4 symbols
		{
			fft_in_q[2 * i] = tmp_x2[i + N_];
			fft_in_q[2 * i + 1] = tmp_y2[i + N_];
			z_s = (short)(511 * tmp_x2[i + N_]);
			fwrite(&z_s, sizeof(short), 1, f_short_I);
			z_s = (short)(511 * tmp_y2[i + N_]);
			fwrite(&z_s, sizeof(short), 1, f_short_Q);
		}
		else
		{
			fft_in_q[2 * i] = 0;
			fft_in_q[2 * i + 1] = 0;
		}
	}
	DSPF_sp_fftSPxSP_cn(Nfft_Q, fft_in_q, w_f, fft_res_q, rad, 0, Nfft_Q);
	/*
	//i = 0;
	for( i_freq=-32;i_freq<2;i_freq+=2 )
	{
	 corr_Max = search_short_fft_Q(ref0_q,fft_res_q,i_freq,Nfft_Q);
	 frequency_position[(i_freq>>1)+16] = corr_pick;
	 corr_position[(i_freq>>1)+16]  = corr_Max;
	 //i++;
	}
	max_corr0 = 0;
	corr_Maximorum0 = 0;
	for( i=0;i<17;i++)
	{
	  if( frequency_position[i]>max_corr0 )
	  {
		max_corr0 = frequency_position[i];
		corr_Maximorum0 = i;
	  }
	}
   // i = 0;
	for( i_freq=2;i_freq<34;i_freq+=2 )
	{
	 corr_Max = search_short_fft_Q(ref0_q,fft_res_q,i_freq,Nfft_Q);
	 frequency_position[(i_freq>>1)+16] = corr_pick;
	 corr_position[(i_freq>>1)+16]  = corr_Max;
	// i++;
	}
	max_corr1 = 0;
	corr_Maximorum1 = 0;
	for( i=17;i<33;i++)
	{
	  if( frequency_position[i]>max_corr1 )
	  {
		max_corr1 = frequency_position[i];
		corr_Maximorum1 = i;
	  }
	}
	*/

	printf("calculation start\n");
	auto time_start = high_resolution_clock::now();
	
	SearchRes sr[2];
	//part_search(0, &sr[0]);
	thread th1(part_search, 0, &sr[0]);
	part_search(1, &sr[1]);
	th1.join();

	int max_search_ind = (sr[1].max_corr >= sr[0].max_corr) ? sr[1].corr_Maximorum : sr[0].corr_Maximorum;
	f_car_est_0 = 2 * (max_search_ind - 16) * (2 * Fs * 1.0e+06 / Nfft_Q);// frequency in Hz
	x = (f_car_est_0 - car_shift_Hz) / 250.0; //

	auto time_stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(time_stop - time_start);
	
	printf("%f\n", f_car_est_0);
	printf("%i\n", corr_position[max_search_ind]);
	printf("calculation end (%f sec)", duration.count() / 1000000.0);
	return 0;
}

void part_search(bool k, SearchRes* result)
{
	int i_freq, i, start;
	int corr_Max;

	if (k == 0) 
		start = -32; 
	else 
		start = 0;

	for (i_freq = start; i_freq < (start + 34); i_freq += 2)
	{
		corr_Max = search_short_fft_Q(ref0_q, fft_res_q, i_freq, Nfft_Q);
		frequency_position[(i_freq >> 1) + 16] = corr_pick;
		corr_position[(i_freq >> 1) + 16] = corr_Max;
		//i++;
	}

	if (k == 0) 
		start = 0; 
	else 
		start = 17;

	result->corr_Maximorum = 0;
	result->max_corr = 0;
	
	for (i = start; i < (start + 16); i++)
	{
		if (frequency_position[i] > result->max_corr)
		{
			result->max_corr = frequency_position[i];
			result->corr_Maximorum = i;
		}
	}

	printf("search %d finished (%d, %f)\n", k, result->corr_Maximorum, result->max_corr);
}


int search_short_fft_Q(float* ref_0, float* fft_res, int k_0, int n_fft)
{
	int i, N_corr;
	float x;
	if (k_0 >= 0)
	{
		for (i = 2 * k_0; i < 2 * n_fft; i += 2)
		{
			fft_src_q[i - 2 * k_0] = fft_res[i]; //x_i[i-2*k_0] = fft_res[i];
			fft_src_q[i - 2 * k_0 + 1] = fft_res[i + 1]; //x_i[i-2*k_0+1] = fft_res[i+1];
		}
		for (i = 0; i < 2 * k_0; i += 2)
		{
			fft_src_q[2 * n_fft + i - 2 * k_0] = fft_res[i]; //x_i[2*N_fft+i-2*k_0] = fft_res[i];
			fft_src_q[2 * n_fft + i - 2 * k_0 + 1] = fft_res[i + 1]; // x_i[2*N_fft+i-2*k_0+1] = fft_res[i+1]
		}
	}
	if (k_0 < 0)
	{
		for (i = 0; i < 2 * n_fft + 2 * k_0; i += 2)
		{
			fft_src_q[i - 2 * k_0] = fft_res[i]; // x_i[i-2*k_0] = fft_res[i];
			fft_src_q[i - 2 * k_0 + 1] = fft_res[i + 1]; // x_i[i-2*k_0+1] = fft_res[i+1];
		}
		for (i = 0; i < -2 * k_0; i += 2)
		{
			fft_src_q[i] = fft_res[2 * n_fft + i + 2 * k_0]; // x_i[i] = fft_res[2*N_fft+i+2*k_0];
			fft_src_q[i + 1] = fft_res[2 * n_fft + i + 2 * k_0 + 1]; // x_i[i+1] = fft_res[2*N_fft+i+2*k_0+1];
		}
	}
	// FFT by reference FFT multiplication
	for (i = 0; i < n_fft; i++)
	{
		fft_corr_q[2 * i] = fft_src_q[2 * i + 1] * ref_0[2 * i + 1] - fft_src_q[2 * i] * ref_0[2 * i];
		fft_corr_q[2 * i + 1] = fft_src_q[2 * i + 1] * ref_0[2 * i] + fft_src_q[2 * i] * ref_0[2 * i + 1];
	}
	// Inverce FFT
	DSPF_sp_ifftSPxSP_cn(n_fft, fft_corr_q, w_f, fft_src, rad, 0, n_fft);
	for (i = 0; i < n_fft; i++)
	{
		corr_q[i] = sqrtf(fft_src[2 * i] * fft_src[2 * i] + fft_src[2 * i + 1] * fft_src[2 * i + 1]);
	}
	for (i = 0; i < N_Q_corr; i++)
	{
		corr_Q[i] = corr_q[i];
		corr_Q[i] += corr_q[n_fft + i - N_Q_corr];
		corr_Q[i] += corr_q[n_fft + i - 2 * N_Q_corr];
		corr_Q[i] += corr_q[n_fft + i - 3 * N_Q_corr];
		corr_Q[i] += corr_q[n_fft + i - 4 * N_Q_corr];
	}
	x = 0;
	N_corr = 0;
	for (i = 0; i < N_Q_corr; i++)
	{
		if (corr_Q[i] > x)
		{
			x = corr_Q[i];
			N_corr = i;
		}
	}
	corr_pick = x;
	return N_corr;
}
