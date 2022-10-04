#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include "image.h"

float get_pixel(image im, int x, int y, int c)
{
    // make sure we get the correct channel number
    assert(c >= 0 && c < im.c);
    // deal with boundary problems
	x = x < 0 ? 0 : x;
	x = x >= im.w ? im.w-1 : x;
	y = y < 0 ? 0 : y;
	y = y >= im.h ? im.h-1 : y;
	// return
	return im.data[x + im.w*y + im.w*im.h*c];
}

void set_pixel(image im, int x, int y, int c, float v)
{
    // TODO Fill this in
    // deal with the boundary problem
    assert(x >= 0 && x < im.w && y >= 0 && y < im.h && c >= 0 && c < im.c);
	// set the pixel
	im.data[x + im.w*y + im.w*im.h*c] = v;
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    int i,j,k;
    // go over the whole image
    for(k = 0; k < im.c; ++k){
        for(j = 0; j < im.h; ++j){
            for(i = 0; i < im.w; ++i){
                copy.data[i + im.w*j + im.w*im.h*k] = im.data[i + im.w*j + im.w*im.h*k];
            }
        }
    }
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    int i,j;
    for(j = 0; j < im.h; ++j){
        for(i = 0; i < im.w; ++i){
            gray.data[i + im.w*j] =0.299*im.data[i + im.w*j]+ 0.587*im.data[i + im.w*j + + im.w*im.h] + 0.114* im.data[i + im.w*j + im.w*im.h*2];
        }
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    // TODO Fill this in
    assert(c>=0 && c<im.c); // needs to be a valid channel
     int i,j;
     for(j = 0; j < im.h; ++j){
        for(i = 0; i < im.w; ++i){
            im.data[i + im.w*j +  im.w*im.h*c ] += v ;
           }
     }
}

void clamp_image(image im)
{
    // TODO Fill this in
    int i,j,k;
    for(k = 0; k < im.c; ++k){
        for(j = 0; j < im.h; ++j){
            for(i = 0; i < im.w; ++i){
               im.data[i + im.w*j + im.w*im.h*k] = im.data[i + im.w*j + im.w*im.h*k] < 0 ? 0 : im.data[i + im.w*j + im.w*im.h*k];
               im.data[i + im.w*j + im.w*im.h*k] = im.data[i + im.w*j + im.w*im.h*k] > 1 ? 1 : im.data[i + im.w*j + im.w*im.h*k];
            }
        }
    }

}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

float H_fun(float R, float G, float B, float V, float C) {
    // first case
	if (C == 0)
		return 0;
	float H;
	// find the case / mean // with in 0,1
	if (V == B)
		H = (R - G) / C + 4;
	if (V == R)
		H = (G - B) / C;
	if (V == G)
		H = (B - R) / C + 2;
	return H < 0 ? H / 6 + 1 : H / 6;
}
void rgb_to_hsv(image im)
{
    // TODO Fill this in
    int i,j;
    for(i = 0; i < im.w; ++i){
        for(j = 0; j < im.h; ++j){
            float R,G,B,V,C,S,H;
            R = im.data[i + im.w*j];
            G = im.data[i + im.w*j + im.w*im.h];
            B = im.data[i + im.w*j + im.w*im.h*2];
            V = three_way_max(R,G,B);
            C = V- three_way_min(R,G,B);
            S = V == 0? 0 : C/V;
            H = H_fun(R,G,B,V,C);
            im.data[i + im.w*j] = H;
            im.data[i + im.w*j + im.w*im.h] = S;
            im.data[i + im.w*j + im.w*im.h*2] = V;
        }

        }
}

void hsv_to_rgb(image im)
{
    // TODO Fill this in
  assert(im.c==3);
  for (int i = 0; i < im.w; i++) {
      for (int j = 0; j < im.h; j++) {
          float H = im.data[i + im.w*j];
          float S = im.data[i + im.w*j + im.w*im.h];
          float V = im.data[i + im.w*j + im.w*im.h*2];
          float C = V * S;
          float X = C * (1 - fabs(fmod(6 * H, 2.0) - 1));
          float m = V - C;
          float temp = 1.0 / 6.0;
          if (H >= 0 && H < temp) {
              im.data[i + im.w*j ] = C + m;
              im.data[i + im.w*j + im.w*im.h] = X + m;
              im.data[i + im.w*j + im.w*im.h*2] = m;
          }
          if (H >= temp && H < 2 * temp) {
              im.data[i + im.w*j ] = X + m;
              im.data[i + im.w*j + im.w*im.h] = C + m;
              im.data[i + im.w*j + im.w*im.h*2] = m;
          }
          if (H >= 2 * temp && H < 3 * temp) {
              im.data[i + im.w*j ] = m;
              im.data[i + im.w*j + im.w*im.h] = C + m;
              im.data[i + im.w*j + im.w*im.h*2] = X + m;
          }
          if (H >= 3 * temp && H < 4 * temp) {
              im.data[i + im.w*j ] = m;
              im.data[i + im.w*j + im.w*im.h] = X + m;
              im.data[i + im.w*j + im.w*im.h * 2] = C + m;
          }
          if (H >= 4 * temp && H < 5 * temp) {
              im.data[i + im.w*j ] = X + m;
              im.data[i + im.w*j + im.w*im.h] = m;
              im.data[i + im.w*j + im.w*im.h*2] = C + m;
          }
          if (H >= 5 * temp && H < 1) {
              im.data[i + im.w*j] = C + m;
              im.data[i + im.w*j + im.w*im.h] = m;
              im.data[i + im.w*j + im.w*im.h*2] = X + m;
          }
      }
  }
}

void scale_image(image im, int c, float v)
{
    // give the channel assert:

    assert(c>=0 && c < im.c);
     int i,j;
     for(j = 0; j < im.h; ++j){
        for(i = 0; i < im.w; ++i){
            im.data[i + im.w*j +  im.w*im.h*c ] *= v ;
           }
        }
}

void rgb_to_HCL(image im)
{
    // assert channel
    assert(im.c == 3);
    int i,j;
     for(j = 0; j < im.h; ++j){
        for(i = 0; i < im.w; ++i){
            float gama = 2.2;
            float R = pow(im.data[i + im.w*j], gama);
            float G = pow(im.data[i + im.w*j +  im.w*im.h], gama);
            float B = pow(im.data[i + im.w*j +  im.w*im.h*2], gama);

            //SET THE xyz, XYZ
            float X = 0.5767309*R + 0.1855540*G + 0.1881852*B;
			float Y = 0.2973769*R + 0.6273491*G + 0.0752741*B;
			float Z = 0.0270343*R + 0.0706872*G + 0.9911085*B;
			float x = 0.9504;
			float y = 1.0000;
			float z = 1.0888;


			//set u,v
			float uu = 4 * X / (X + 15 * Y + 3 * Z);
			float vv = 9 * Y / (X + 15 * Y + 3 * Z);
			float ur = 4 * x / (x + 15 * y + 3 * z);
			float vr = 9 * y / (x + 15 * y + 3 * z);


			// find the L,C,H
			float L= Y/y > 0.008856 ? 116 * pow(Y/y, 1.0 / 3.0) - 16 : 903.3*Y/y;
			float u = 13 * L*(uu - ur);
			float v = 13 * L*(vv - vr);
			float C=sqrt(pow(u, 2) + pow(v, 2));
			float H=atan2(13 * L*(vv - vr),13 * L*(uu - ur));
			im.data[i + im.w*j] = L;
            im.data[i + im.w*j + im.w*im.h] = C;
            im.data[i + im.w*j + im.w*im.h*2] = H;
           }
        }


}
