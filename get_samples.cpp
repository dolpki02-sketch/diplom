#include <stdint.h>
#include <math.h>

#define pi 3.1415926
#define _2pi 2*3.1415926

#include "include/fda_poliphase_short.h"
#include "include/fda_poliphase.h"

int get_samples_1(short *f_x,int N,float f_car,float Fs,float Fd,float *i_tmp,float *q_tmp)
{
   int i,jj,n,l,m,nn;
   int i_cic,i_cic_,i_cic_shift,k_decim;
   float k_res,ns;
   float x_res,y_res;
   short x_s;
   //extern float tmp_x_in[655360];
   short wx_s,wy_s,wx_x_s,wy_y_s,wx_x_s_,wy_y_s_;
   int z_x,z_y;
   short in_poli_i[BL>>1];
   short in_poli_q[BL>>1];
   float s_i[BL_p];
   short s_q[BL_p];
   extern float Lagr[8][1024];
   float defi_dem,d_dem,d0,defi;
   float I_buf[8]; // resampler buffer I
   float Q_buf[8]; // resampler buffer Q
   // CIC filter I
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
   // CIC filter Q
   int64_t f0_intq1,f0_intq2,f0_intq3,f0_intq4;
   int64_t f0_combq1,f0_combq2,f0_combq3,f0_combq4;
   int64_t f0_last_combq1,f0_last_combq2,f0_last_combq3,f0_last_combq4;
   f0_intq1 = 0;
   f0_intq2 = 0;
   f0_intq3 = 0;
   f0_intq4 = 0;
   f0_combq1 = 0;
   f0_combq2 = 0;
   f0_combq3 = 0;
   f0_combq4 = 0;
   f0_last_combq1 = 0;
   f0_last_combq2 = 0;
   f0_last_combq3 = 0;
   f0_last_combq4 = 0;
   for( i=0;i<BL>>1;i++ )
   {
       in_poli_i[i] = 0;
       in_poli_q[i] = 0;
   }
   for( i=0;i<BL_p;i++ )
   {
       s_i[i] = 0;
       s_q[i] = 0;
   }
   for( i=0;i<8;i++ )
   {
     I_buf[i] = 0;
     Q_buf[i] = 0;
   }
   i = 1;
   k_res = 2*Fs/140;
   while (k_res<0.25)
   {
       k_res  *= 2.0;
       i++;
   }
   i_cic = i;
   i_cic_= 1;
   if(i_cic>1)
   {
       i_cic_=1<<(i-1);// decimation coefficient
       switch(i_cic_)
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
   }
   d_dem = 2*pi*f_car/(Fd/2);
   defi_dem = 0;
   defi = 0;
   d0 = pi/2.0f;
   n = 0;
   l = 0;
   m = 0;
   nn = 0;
   k_decim = 0;
   ns = 0;
   for(jj=0;jj<N;jj++)
   {
       //fread(&x_s,2,1,f_in); // f_in is short
       x_s = *f_x++;
      // if(jj<500000) tmp_x_in[jj] = (float)x_s;
       if( jj%2==0 )
       {
         z_x = 0;
         for( i=0;i<BL>>1;i++ )
         {
          z_x +=  (int)(in_poli_i[i]*B[BL-2-2*i]);
         }
         z_y = 0;
         for( i=0;i<BL>>1;i++ )
         {
          z_y +=  (int)(in_poli_q[i]*B[BL-1-2*i]);
         }
         for( i=(BL>>1)-1;i>0;i-- )
         {
          in_poli_i[i] = in_poli_i[i-1];
         }
         wx_s = (short)(z_x>>15);
         wy_s = (short)(z_y>>15);
         if( nn%2==0 )
         {
          in_poli_i[0] = x_s;
         }
         else
         {
           in_poli_i[0] = -x_s;
         }
         nn++;
       }
       else
       {
         z_x = 0;
         for( i=0;i<(BL>>1);i++ )
         {
          z_x += (int)(in_poli_i[i]*B[BL-1-2*i]);
         }
         z_y = 0;
         for( i=0;i<(BL>>1);i++ )
         {
          z_y +=  (int)(in_poli_q[i]*B[BL-2-2*i]);
         }
         for( i=(BL>>1)-1;i>0;i-- )
         {
          in_poli_q[i] = in_poli_q[i-1];
         }
         wx_s = (short)(z_x>>15);
         wy_s = (short)(z_y>>15);
         if( l%2==0 )
         {
          in_poli_q[0] = x_s;
         }
         else
         {
           in_poli_q[0] = -x_s;
         }
         l++;
       }
       if( jj%2==0 )
       {
         defi_dem+=d_dem;
         if(defi_dem>_2pi) defi_dem-=_2pi;
         if(defi_dem<0) defi_dem+=_2pi;
         //
         wx_x_s = (short)(((int)(2*wy_s*(short)(32767*sin(defi_dem))) - (int)(2*wx_s*(short)(32767*cos(defi_dem))))>>15);
         wy_y_s = (short)(((int)(2*wx_s*(short)(32767*sin(defi_dem))) + (int)(2*wy_s*(short)(32767*cos(defi_dem))))>>15);
         if(i_cic>1)
         {
             // I chanel
             //Integrated
             //f0_inti1 += (int64_t )(wx_x*1023);
             f0_inti1 += (int64_t )(wx_x_s);
             f0_inti2 += f0_inti1;
             f0_inti3 += f0_inti2;
             f0_inti4 += f0_inti3;
             // Q chanel
             //Integrated
             //f0_intq1 += (int64_t )(wy_y*1023);
             f0_intq1 += (int64_t )(wy_y_s);
             f0_intq2 += f0_intq1;
             f0_intq3 += f0_intq2;
             f0_intq4 += f0_intq3;
             // decimation by dec_cic
             //Comb
             if( k_decim%i_cic_==0 )
             {
                 //I chanel
                 f0_combi1 = f0_inti4-f0_last_combi1;
                 f0_combi2 = f0_combi1 - f0_last_combi2;
                 f0_combi3 = f0_combi2 - f0_last_combi3;
                 f0_combi4 = f0_combi3 - f0_last_combi4;
                 // update
                 f0_last_combi1 = f0_inti4;
                 f0_last_combi2 = f0_combi1;
                 f0_last_combi3 = f0_combi2;
                 f0_last_combi4 = f0_combi3;
                 wx_x_s_ =(short) ((f0_combi4>>i_cic_shift));
                 //tmp_x2[k_decim>>i_cic] = wx_x_;
                 //Q chanel
                 f0_combq1 = f0_intq4-f0_last_combq1;
                 f0_combq2 = f0_combq1 - f0_last_combq2;
                 f0_combq3 = f0_combq2 - f0_last_combq3;
                 f0_combq4 = f0_combq3 - f0_last_combq4;
                 // update
                 f0_last_combq1 = f0_intq4;
                 f0_last_combq2 = f0_combq1;
                 f0_last_combq3 = f0_combq2;
                 f0_last_combq4 = f0_combq3;
                 wy_y_s_ =(short)((f0_combq4>>i_cic_shift));
                 //tmp_y2[k_decim>>i_cic] = wy_y_;
             }
         }
         else
         {
            // no CIC filters
             wx_x_s_ = wx_x_s;
             wy_y_s_ = wy_y_s;
         }
         if( k_decim%i_cic_==0  )
         {
             if(ns<1)
             {
                 I_buf[n] = wx_x_s_/1024.0;
                 Q_buf[n] = wy_y_s_/1024.0;
                 n++;
                 if( n>7) n -= 8;
                 ns+=k_res;
             }
            else
             {
                 ns-=1;
                // mu = ns;
                 int l=(int)(ns*1023);
                 x_res = 0;
                 y_res = 0;
                 for(i=0;i<8;i++)
                 {
                     int i0 = i+n;
                     if(i0>7) i0-=8;
                     x_res+=Lagr[i][l]*I_buf[i0];
                     y_res+=Lagr[i][l]*Q_buf[i0];
                 }
                 ns+=k_res;
                 i_tmp[m] = x_res;
                 q_tmp[m] = y_res;
                 m++;
              }
             // resampler
         }
         k_decim++;
       }
   }
   //f_x-=N; // return pointer
   return m;
}

