#include <stdio.h>
#include <stdint.h>
#include <math.h>

#define pi 3.1415926
#define _2pi 2*3.1415926

float noise()
{
    short v;
    float noi;
    v=rand();
    noi=((float)(v)-(float)(16383.5))/(float)(32768);
    return (noi);
}



float samples_cic_fs_s(FILE *f_x,float *Spread,int N_spread,int N_0,float f_car,float Fs,float Fd,float S_N)
{
    float z_max;
    float x;
    float wx;
    int jj,nn,kk,mm,i,n,ii;
    int i_cic,i_cic_shift;
    short k_random;
    float ns,mu; // resampler accumulator
    float kp; //
    float x_res;
    extern float Lagr[8][1024];
    float I_buf[8];
    float d0; // carrier frequency is part of symbol rate
    float defi;
    float k_res;
    // CIC filter
    int64_t f0_inti1,f0_inti2,f0_inti3,f0_inti4;
    int64_t f0_combi1,f0_combi2,f0_combi3,f0_combi4;
    int64_t f0_last_combi1,f0_last_combi2,f0_last_combi3,f0_last_combi4;
    f0_inti1 = 0;
    f0_inti2 = 0;
    f0_inti3 = 0;
    f0_inti4 = 0;
    f0_combi1 = 0;
    f0_combi2 = 0;
    f0_combi3 = 0;
    f0_combi4 = 0;
    f0_last_combi1 = 0;
    f0_last_combi2 = 0;
    f0_last_combi3 = 0;
    f0_last_combi4 = 0;
    d0 =2*pi*(f_car+70.0)/Fd;
    i = 1;
    while (Fs*powf(2,i+1)<(Fd/2))
    {
        x = Fs*powf(2,i+1);
        i++;
    }
    i--;
    i_cic = 1<<i; // interpolation coeffissient
    switch(i_cic)
    {
        case 2:
        i_cic_shift = 4;
            break;
        case 4:
        i_cic_shift = 8;
            break;
        case 8:
        i_cic_shift = 12;
            break;
        case 16:
        i_cic_shift = 16;
            break;
        case 32:
        i_cic_shift = 20;
            break;
        case 64:
        i_cic_shift = 24;
           break;
        case 128:
        i_cic_shift = 28;
            break;
        case 256:
        i_cic_shift = 32;
        break;
        case 512:
        i_cic_shift = 36;
        break;
        case 1024:
        i_cic_shift = 40;
        break;
    }
    z_max = 0;
    x = Fd/(float)i_cic;
    k_res = 2*Fs/x;
    kp = k_res;  // resampling coef. The Fs rises in 1/k_res
    for( i =0;i<8;i++)
    {
           I_buf[i] = 0;
    }
    float z=log10f(sqrt(2*k_res)); // C/N = ((Es/No)*/2)*k_res for fd = 4Fs*k_res (2 - is half Fd)=> z = lg(sqrt(2*k_res))
    z-=S_N/20;     // S_N - signal to nois ratio
    float gs=pow(10,z);
    kk = 0;
    ns = 0;
    k_random=rand();
    defi=_2pi*((float)(k_random)-(float)(16383.5))/(float)(32768);// rand float from -2pi to 2pi
    n = 0;
    ii = 0;
for(jj=0;jj<N_0;jj++)
{
nn=rand();
nn&=1;
if(nn==0) {wx=0.707;}
if(nn==1) {wx=-0.707;}
    kk = 0;
    while (kk<N_spread)
    {
        if(ns>1)
        {
           I_buf[n] = wx*Spread[kk];
            n++;
            if( n>7) n -= 8;
             ns-=1;
             kk++;
        }
         else
        {
           mu = ns;
           int l=(int)(mu*1023);
           x_res = 0;
           for(int i=0;i<8;i++)
           {
             int i0 = i+n;
             if(i0>7) i0-=8;
             x_res+=Lagr[7-i][l]*I_buf[i0];
           }
           ns+=kp;
           //Comb
           f0_combi1 = (int64_t) (x_res*1023)-f0_last_combi1;
           f0_combi2 = f0_combi1 - f0_last_combi2;
           f0_combi3 = f0_combi2 - f0_last_combi3;
           f0_combi4 = f0_combi3 - f0_last_combi4;
           // update
           f0_last_combi1 = (int64_t) (x_res*1023);
           f0_last_combi2 = f0_combi1;
           f0_last_combi3 = f0_combi2;
           f0_last_combi4 = f0_combi3;
           // interpolation by i_cic
           for(mm=0;mm<i_cic;mm++)
           {
            //Integrated
             f0_inti1 += (f0_combi4);
             f0_inti2 += f0_inti1;
             f0_inti3 += f0_inti2;
             f0_inti4 += f0_inti3;
             defi+=d0;
             if(defi>_2pi) defi-=_2pi;
             if(defi<0) defi+=_2pi;
             z=0;
             for(nn=0;nn<12;nn++) z+=gs*noise();
             x = ((float)(f0_inti4>>i_cic_shift)/1023.0)*cosf(defi)+((float)(f0_inti4>>i_cic_shift)/1023.0)*sinf(defi)+z;
             //x = ((float)(f0_inti4>>i_cic_shift)/1023.0)*cosf(defi)+((float)(f0_inti4>>i_cic_shift)/1023.0)*sinf(defi);
             if(fabsf(x)>=z_max) z_max = fabsf(x);
             fwrite(&x,4,1,f_x);
             if(ii<5000)
             {
              z = (float)(f0_inti4>>i_cic_shift)/(float)1023;
             }
             ii++;
            }
          }
         }
    }
return z_max;
}
//}




