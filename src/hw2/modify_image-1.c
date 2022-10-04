#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

/******************************** Resizing *****************************
  To resize we'll need some interpolation methods and a function to create
  a new image and fill it in with our interpolation methods.
************************************************************************/
// add new method
float clamp_pixel(image im, int x, int y, int c)
    {
    assert(c<im.c && c>=0);
    if(x>=im.w)x=im.w-1;
    if(y>=im.h)y=im.h-1;
    if(x<0)x=0;
    if(y<0)y=0;
    return im.data[x+ im.w*y+im.w*im.h*c];
    }

float nn_interpolate(image im, float x, float y, int c)
{
    // TODO
    /***********************************************************************
      This function performs nearest-neighbor interpolation on image "im"
      given a floating column value "x", row value "y" and integer channel "c",
      and returns the interpolated value.
    ************************************************************************/
    int near_x = roundf(x);
	int near_y = roundf(y);

	if (near_x>=im.w) near_x=im.w-1;
    if (near_y>=im.h) near_y=im.h-1;
    if (near_x<0) near_x=0;
    if (near_y<0) near_y=0;

	return im.data[near_x+im.w*near_y+im.w*im.h*c];
}

image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix the return line)
    /***********************************************************************
      This function uses nearest-neighbor interpolation on image "im" to a new
      image of size "w x h"
    ************************************************************************/
    image new_im = make_image(w,h,im.c);

    float scale_w = (im.w * 1.0) / w;
    float scale_h = (im.h * 1.0) / h;
    for (int c = 0; c < im.c; c++) {
        for (int i = 0; i < w; i++) {
		  for (int j = 0; j < h; j++) {
			  new_im.data[i+ w*j+w*h*c] = nn_interpolate(im,scale_w*(i+0.5)-0.5, scale_h*(j+0.5)-0.5, c);
		  }
	    }
	  }
    return new_im;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // TODO
    /***********************************************************************
      This function performs bilinear interpolation on image "im" given
      a floating column value "x", row value "y" and integer channel "c".
      It interpolates and returns the interpolated value.
    ************************************************************************/
    int x1 = floor(x);
	int x2 = ceil(x);
	int y1 = floor(y);
	int y2 = ceil(y);

	//v1
	float v1 = clamp_pixel(im, x1, y1, c);
	//v2
	float v2 = clamp_pixel(im, x2, y1, c);
	//v3
	float v3 = clamp_pixel(im, x1, y2, c);
	//v4
    float v4 = clamp_pixel(im, x2, y2, c);
    return v1*(x2-x)*(y2-y)+v2*(x-x1)*(y2-y)+v3*(x2-x)*(y-y1)+v4*(x-x1)*(y-y1);
  }

image bilinear_resize(image im, int w, int h)
{
    // TODO
    /***********************************************************************
      This function uses bilinear interpolation on image "im" to a new image
      of size "w x h". Algorithm is same as nearest-neighbor interpolation.
    ************************************************************************/
    image new_im = make_image(w,h,im.c);

    float scale_w = (im.w * 1.0) / w;
    float scale_h = (im.h * 1.0) / h;
    for (int c = 0; c < im.c; c++) {
        for (int i = 0; i < w; i++) {
		  for (int j = 0; j < h; j++) {
			  new_im.data[i+ w*j+w*h*c] = bilinear_interpolate(im,scale_w*(i+0.5)-0.5, scale_h*(j+0.5)-0.5, c);
		  }
	    }
	  }


    return new_im;
}


/********************** Filtering: Box filter ***************************
  We want to create a box filter. We will only use square box filters.
************************************************************************/

void l1_normalize(image im)
{
    // TODO
    /***********************************************************************
      This function divides each value in image "im" by the sum of all the
      values in the image and modifies the image in place.
    ************************************************************************/
    int i, j,m,n, c;
    float sum;
    for (c = 0; c < im.c; c++){
        sum = 0;
        for ( i = 0; i < im.w; i++) {
			for (j = 0; j < im.h; j++) {
				sum += im.data[i + j*im.w+ c*im.w*im.h];
			}
    }
    for (m=0; m < im.w; m++) {
			for (n = 0; n < im.h; n++) {
				im.data[m + n*im.w+ c*im.w*im.h] /= sum;
			}
		}
    }
}

image make_box_filter(int w)
{
    // TODO
    /***********************************************************************
      This function makes a square filter of size "w x w". Make an image of
      width = height = w and number of channels = 1, with all entries equal
      to 1. Then use "l1_normalize" to normalize your filter.
    ************************************************************************/
    assert(w%2); // w needs to be odd
    image im = make_image(w,w,1);
    for (int i = 0; i < w; i++) {
	  for (int j = 0; j < w; j++) {
		  im.data[i+ j*w]= 1.0 / (w*w);
	  }
  }

    l1_normalize(im);
  return im;
}

image convolve_image(image im, image filter, int preserve)
{
    // TODO
    /***********************************************************************
      This function convolves the image "im" with the "filter". The value
      of preserve is 1 if the number of input image channels need to be 
      preserved. Check the detailed algorithm given in the README.  
    ************************************************************************/
  assert(filter.c ==1);
  int mid = filter.w / 2;
  image new_im;

  if (preserve==1) {
	  new_im = make_image(im.w, im.h, im.c);
	  for (int c = 0; c < im.c; c++) {
		  for (int i = 0; i < im.w; i++) {
			  for (int j = 0; j < im.h; j++) {
				  float q = 0;
				  for (int m = 0; m < filter.w; m++) {
					  for (int n = 0; n < filter.w; n++) {
						  q += filter.data[m+n*filter.w] *clamp_pixel(im,i+m-mid, j+n-mid, c);
					  }
				  }
				  new_im.data[i + j*new_im.w + c*new_im.w*new_im.h] = q;
			  }
		  }
	  }
  }
  else {
	  new_im = make_image(im.w, im.h, 1);
	  for (int i = 0; i < im.w; i++) {
		  for (int j = 0; j < im.h; j++) {
			  float q = 0;
			  for (int c = 0; c < im.c; c++) {
				  for (int m = 0; m < filter.w; m++) {
					  for (int n = 0; n < filter.w; n++) {
						  q += filter.data[m + n*filter.w]*clamp_pixel(im, i+m-mid, j+n-mid, c);
					  }
				  }
			  }
			  new_im.data[i + j*new_im.w] = q;
		  }
	  }
  }
  return new_im;
  }


image make_highpass_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with highpass filter values using image.data[]
    ************************************************************************/
    image new_im = make_box_filter(3);
    new_im.data[0] = new_im.data[2] = new_im.data[2*3] = new_im.data[2+2*3] = 0;
	new_im.data[1] = new_im.data[1*3] = new_im.data[2+1*3] = new_im.data[1+ 2*3]= -1;
	new_im.data[1+ 3] = 4;

    return new_im;
}

image make_sharpen_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with sharpen filter values using image.data[]
    ************************************************************************/
    image new_im = make_box_filter(3);
    new_im.data[0] = new_im.data[2] = new_im.data[2*3] = new_im.data[2+2*3] = 0;
	new_im.data[1] = new_im.data[1*3] = new_im.data[2+1*3] = new_im.data[1+ 2*3]= -1;
	new_im.data[1+ 3] = 5;

    return new_im;
}

image make_emboss_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with emboss filter values using image.data[]
    ************************************************************************/
    image new_im = make_box_filter(3);
    new_im.data[0] = -2;
    new_im.data[2] = new_im.data[2*3] = 0;
    new_im.data[2+2*3] = 2;
	new_im.data[1] = new_im.data[1*3] = -1;
	new_im.data[2+1*3] = new_im.data[1+ 2*3]= new_im.data[1+ 3] = 1;

    return new_im;
}

// Question 2.3.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: the emboss filter and sharpen filter need preserve. The high pass filter don't need preserve. This is because we need
// to keep the color and distinguish the channels after emboss and sharpen filter in our case
// but high pass only need to keep the edges and curves and minor difference between channels.

// Question 2.3.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: we need use "clamp" as the post processing after filters. Need to keep value range from 0 to 1 by using the clamp.


float compute_gaussian(float sigma, int x, int y) {
	return expf(-(pow(x, 2) + pow(y, 2)) / (2 * pow(sigma, 2))) / (TWOPI * pow(sigma, 2));
}

image make_gaussian_filter(float sigma)
{
    // TODO
    /***********************************************************************
      sigma: a float number for the Gaussian.
      Create a Gaussian filter with the given sigma. Note that the kernel size 
      is the next highest odd integer from 6 x sigma. Return the Gaussian filter.
    ************************************************************************/
    int w = ceil(6 * sigma);
	w = w % 2 == 0 ? w + 1 : w;
	image new_im = make_box_filter(w);
	int offset = w / 2;
	for (int i = 0; i < w; i++) {
		for (int j = 0; j < w; j++) {
			new_im.data[i + j*new_im.w] = compute_gaussian(sigma, (i-offset), (j-offset));
		}
	}
	l1_normalize(new_im);

  return new_im;

}

image add_image(image a, image b)
{
    // TODO
    /***********************************************************************
      The input images a and image b have the same height, width, and channels.
      Sum the given two images and return the result, which should also have
      the same height, width, and channels as the inputs. Do necessary checks.
    ************************************************************************/
  assert(a.w==b.w && a.h==b.h && a.c==b.c); // assure images are the same size
  int w = a.w;
  int h = a.h;
  int c = a.c;
  image im = make_image(w, h, c);
  for (int k = 0; k < c; k++) {
	  for (int i = 0; i < w; i++) {
		  for (int j = 0; j < h; j++) {
			  im.data[i + j*im.w + k*im.h*im.w] = a.data[i + j*im.w + k*im.h*im.w] + b.data[i + j*im.w + k*im.h*im.w];
		  }
	  }
  }

  return im;

}

image sub_image(image a, image b)
{
  assert(a.w==b.w && a.h==b.h && a.c==b.c); // assure images are the same size
  int w = b.w;
  int h = b.h;
  int c = b.c;
  image im = make_image(w, h, c);
  // TODO: Implement subtraction
  for (int k = 0; k < c; k++) {
	  for (int i = 0; i < w; i++) {
		  for (int j = 0; j < h; j++) {
			 im.data[i + j*im.w + k*im.h*im.w] = a.data[i + j*im.w + k*im.h*im.w] - b.data[i + j*im.w + k*im.h*im.w];
		  }
	  }
  }

  return im;
}

image make_gx_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 Sobel Gx filter and return it
    ************************************************************************/
    image new_im = make_box_filter(3);
    new_im.data[0] = new_im.data[2*3] = -1;
    new_im.data[1] = new_im.data[1+1*3] = new_im.data[1+ 2*3] = 0;
	new_im.data[2] = new_im.data[2 + 2*3] = 1;
	new_im.data[3] = -2;
	new_im.data[2 + 1*3] = 2;


    return new_im;
}

image make_gy_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 Sobel Gy filter and return it
    ************************************************************************/
    image new_im = make_box_filter(3);
    new_im.data[0] = new_im.data[2] = -1;
    new_im.data[3] = new_im.data[1+1*3] = new_im.data[2+ 1*3] = 0;
	new_im.data[2*3] = new_im.data[2 + 2*3] = 1;
	new_im.data[1] = -2;
	new_im.data[1 + 2*3] = 2;


    return new_im;
}

void feature_normalize(image im)
{
    // TODO
    /***********************************************************************
      Calculate minimum and maximum pixel values. Normalize the image by
      subtracting the minimum and dividing by the max-min difference.
    ************************************************************************/
      float l = 50000;
      float r = -50000;
      for (int c = 0; c < im.c; c++) {
          for (int i = 0; i < im.w; i++) {
              for (int j = 0; j < im.h; j++) {
                  l = (l < im.data[i + j*im.w + c * im.w *im.h]) ? l : im.data[i + j*im.w + c * im.w *im.h] ;
                  r = (r > im.data[i + j*im.w + c * im.w *im.h]) ? r : im.data[i + j*im.w + c * im.w *im.h] ;
              }
          }
      }
      float range = r - l;
      for (int c = 0; c < im.c; c++) {
          for (int i = 0; i < im.w; i++) {
              for (int j = 0; j < im.h; j++) {
                  im.data[i + j*im.w + c * im.w *im.h] = range == 0 ? 0 : (im.data[i + j*im.w + c * im.w *im.h] - l) / range;
              }
          }
      }

  }

image *sobel_image(image im)
{
    // TODO
    /***********************************************************************
      Apply Sobel filter to the input image "im", get the magnitude as sobelimg[0]
      and gradient as sobelimg[1], and return the result.
    ************************************************************************/
    image *sobelimg = calloc(2, sizeof(image));

    image m_x = convolve_image(im, make_gx_filter(), 0);
	image m_y = convolve_image(im, make_gy_filter(), 0);
	image mag = make_image(im.w, im.h,1);
	image grad = make_image(im.w, im.h,1);
	for (int i = 0; i < im.w; i++) {
		for (int j = 0; j < im.h; j++) {
			mag.data[i + j*mag.w] = sqrtf(powf(m_x.data[i + j*m_x.w], 2) + powf(m_y.data[i + j*m_y.w], 2));
			grad.data[i + j*grad.w] = atan2(m_y.data[i + j*m_y.w], m_x.data[i + j*m_x.w]);
		}
	}
	sobelimg[0] = mag;
	sobelimg[1] = grad;

    return sobelimg;
}

image colorize_sobel(image im)
{
  // TODO
  /***********************************************************************
    Create a colorized version of the edges in image "im" using the 
    algorithm described in the README.
  ************************************************************************/
    image *sobel_res = calloc(2, sizeof(image));
    image filter = make_gaussian_filter(4);
    image new_im = convolve_image(im, filter, 1);
    clamp_image(new_im);

    sobel_res = sobel_image(new_im);
	image mag = sobel_res[0];
	image grad = sobel_res[1];
	image hsv_img = make_image(mag.w, mag.h, 3);

	//Nomalize the magnitude:
	feature_normalize(mag);
	feature_normalize(grad);
	//Nomalize the angle and get the hsv image
	assert(grad.w == mag.w);
	assert(grad.h == mag.h);
	for (int i = 0; i < grad.w; i++) {
		for (int j = 0; j < grad.h; j++) {
			//grad.data[i * j*grad.w] = grad.data[i * j*grad.w] / TWOPI + 0.5;
			hsv_img.data[i + j*hsv_img.w] = grad.data[i + j*grad.w];
			hsv_img.data[i + j*hsv_img.w + hsv_img.w * hsv_img.h] = mag.data[i + j * mag.w];
			hsv_img.data[i + j*hsv_img.w + hsv_img.w * hsv_img.h*2] = mag.data[i + j * mag.w];
		}
	}

	//Gaussian smooth the image
	hsv_to_rgb(hsv_img);

  return hsv_img;

}

// EXTRA CREDIT: Median filter

/*
image apply_median_filter(image im, int kernel_size)
{
  return make_image(1,1,1);
}
*/

// SUPER EXTRA CREDIT: Bilateral filter
float compute_gaussian1(float sigma, float error) {
	return expf(-(powf(error, 2)) / (2 * pow(sigma, 2))) / (2 * M_PI * pow(sigma, 2));
}

image apply_bilateral_filter(image im, float sigma1, float sigma2)
{
  int w = ceil(6 * sigma1);
	w = w % 2 == 0 ? w + 1 : w;
	int k = w / 2;
	image bf = make_box_filter(w);
	image new_im = make_image(im.w, im.h, im.c);
	image spatial = make_gaussian_filter(sigma1);

	for (int c = 0; c < im.c; c++) {
		for (int x = 0; x < im.w; x++) {
			for (int y = 0; y < im.h; y++) {
				float N = 0;
				for (int i = 0; i < w; i++) {
					for (int j = 0; j < w; j++) {
						bf.data[i + j*bf.w] = spatial.data[i + j * spatial.w] * compute_gaussian1(sigma2, (clamp_pixel(im, x, y, c) - clamp_pixel(im, x + i - k, y + j - k, c)));
						N += bf.data[i + j*bf.w];
					}
				}
				l1_normalize(bf);
				//Convolution
				float sum = 0;
				for (int m = 0; m < bf.w; m++) {
					for (int n = 0; n < bf.w; n++) {
						sum += bf.data[m +  bf.w * n] * clamp_pixel(im, x + m - k, y + n - k, c);
					}
				}
				new_im.data[x + new_im.w*y + new_im.w * new_im.h * c] = sum;
			 }
		}
	}

  return new_im;
}
