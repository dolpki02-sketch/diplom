#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <inttypes.h>


int file_to_buff_converter(FILE *f_in,short *s_out,float s_ou_max)
{
    int i,m;
    short s_ou_16;
    float x;
    // Read file size
    int64_t fileSize;
    if(f_in  == NULL)
    {
        printf("File %s Open failed\n", "f_in");
        exit (2);
    }
    fileSize = 0;
    while(getc(f_in)!=EOF) fileSize ++;
    fseek(f_in,0,SEEK_SET);
    m = fileSize>>2;
    for( i=0;i<m;i++)
    {
      fread(&x,4,1,f_in);
      x/=s_ou_max;
      s_ou_16 = (short)(x*1023); // 12 bit, to MSB for overflow
      s_out[i]=s_ou_16;
    }
    return m;
}

