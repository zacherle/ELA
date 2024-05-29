//#include <stdio.h>
#include <stdlib.h>

int compare (const void * a, const void * b)
{
  float fa = *(const float*) a;
  float fb = *(const float*) b;
  return (fa > fb) - (fa < fb);
}


// #define COMPARE(a, b) (((a) > (b)) - ((a) < (b)))


float maxgap(int naz, float* az) {

float gap,daz;
int i;

    qsort (az, naz, sizeof(float), compare);

    gap = az[0]+360.0-az[naz-1];

    for (i=1; i<naz; i++) {
        daz=az[i]-az[i-1];
        if (daz > gap) {
            gap=daz;
            }
        }
    return gap;
}
