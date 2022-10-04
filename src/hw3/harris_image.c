#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"
#include <time.h>
#define TWOPI 6.2831853

// Frees an array of descriptors.
// descriptor *d: the array.
// int n: number of elements in array.

void free_descriptors(descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        free(d[i].data);
    }
    free(d);
}

// Create a feature descriptor for an index in an image.
// image im: source image.
// int i: index in image for the pixel we want to describe.
// returns: descriptor for that index.
descriptor describe_index(image im, int i)
{
    int w = 5;
    descriptor d;
    d.p.x = i%im.w;
    d.p.y = i/im.w;
    d.data = calloc(w*w*im.c, sizeof(float));
    d.n = w*w*im.c;
    int c, dx, dy;
    int count = 0;
    // If you want you can experiment with other descriptors
    // This subtracts the central value from neighbors
    // to compensate some for exposure/lighting changes.
    for(c = 0; c < im.c; ++c){
        float cval = im.data[c*im.w*im.h + i];
        for(dx = -w/2; dx < (w+1)/2; ++dx){
            for(dy = -w/2; dy < (w+1)/2; ++dy){
                float val = get_pixel(im, i%im.w+dx, i/im.w+dy, c);
                d.data[count++] = cval - val;
            }
        }
    }
    return d;
}

// Marks the spot of a point in an image.
// image im: image to mark.
// ponit p: spot to mark in the image.
void mark_spot(image im, point p)
{
    int x = p.x;
    int y = p.y;
    int i;
    for(i = -9; i < 10; ++i){
        set_pixel(im, x+i, y, 0, 1);
        set_pixel(im, x, y+i, 0, 1);
        set_pixel(im, x+i, y, 1, 0);
        set_pixel(im, x, y+i, 1, 0);
        set_pixel(im, x+i, y, 2, 1);
        set_pixel(im, x, y+i, 2, 1);
    }
}

// Marks corners denoted by an array of descriptors.
// image im: image to mark.
// descriptor *d: corners in the image.
// int n: number of descriptors to mark.
void mark_corners(image im, descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        mark_spot(im, d[i].p);
    }
}

// Creates a 1d Gaussian filter.
// float sigma: standard deviation of Gaussian.
// returns: single row image of the filter.

image make_1d_gaussian(float sigma)
{
    // TODO: generate 1 d gaussian
    int w = (int)roundf(sigma * 6) + 1;
    if (w%2 == 0){
        w++;
    }
    image g_filter = make_image(w, 1, 1);

    int mid = (int)w / 2;

    for (int i=0; i<w; i++){
        float a1 = 1./(TWOPI*(pow(sigma,2)));
        float a2 = -((i-mid)*(i-mid))/(pow(sigma,2)*2);
          g_filter.data[i] = a1 * exp(a2);
    }

    l1_normalize(g_filter);

    return g_filter;
}



// Smooths an image using separable Gaussian filter.
// image im: image to smooth.
// float sigma: std dev. for Gaussian.
// returns: smoothed image.
image smooth_image(image im, float sigma)
{
    // TODO: use two convolutions with 1d gaussian filter.
    image filter1 = make_1d_gaussian(sigma);
    image filter2 = make_image(1,filter1.w,1);
    for(int i=0; i<filter1.w; i++){
        filter2.data[i] = filter1.data[i];
    }
    image temp = convolve_image(im, filter1, 1);
    image img = convolve_image(temp, filter2, 1);
    return img;
}

// Calculate the structure matrix of an image.
// image im: the input image.
// float sigma: std dev. to use for weighted sum.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          third channel is IxIy.
image structure_matrix(image im, float sigma)
{
    image S = make_image(im.w, im.h, 3);
    assert((im.c==1 || im.c==3) && "only grayscale or rgb supported");
    // convert to grayscale if necessary
//    if(im.c==1)im=im;
//    else im=rgb_to_grayscale(im);

    // TODO: calculate structure matrix for im.
    image sobel_x = make_gx_filter();
    image sobel_y = make_gy_filter();
    image Ix = convolve_image(im, sobel_x, 0);
    image Iy = convolve_image(im, sobel_y, 0);
    for (int i = 0; i < im.w; i++) {
       for (int j = 0; j < im.h; j++) {
          S.data[i + j *S.w ] = pow(Ix.data[i + j * Ix.w],2);
          S.data[i + j *S.w + S.w * S.h]= pow(Iy.data[i + j * Iy.w],2);
          S.data[i + j *S.w + S.w * S.h*2] = Ix.data[i + j * Ix.w]* Iy.data[i + j * Iy.w];
        }
      }
      S = smooth_image(S, sigma);
      return S;
}


// Estimate the cornerness of each pixel given a structure matrix S.
// image S: structure matrix for an image.
// returns: a response map of cornerness calculations.
image cornerness_response(image S)
{
    image R = make_image(S.w, S.h, 1);
    // TODO: fill in R, "cornerness" for each pixel using the structure matrix.
    // We'll use formulation det(S) - alpha * trace(S)^2, alpha = .06.

    for (int i = 0; i < R.w; i++) {
        for (int j = 0; j < R.h; j++) {
          float det = S.data[i+j*S.w] *  S.data[i+j*S.w+S.w*S.h] - pow(S.data[i+j*S.w+S.w*S.h*2], 2);
          float trace = S.data[i+j*S.w] +  S.data[i+j*S.w+S.w*S.h] ;
          R.data[i+j*S.w] = det - 0.06 * trace * trace;
        }
      }

    return R;
}


// Perform non-max supression on an image of feature responses.
// image im: 1-channel image of feature responses.
// int w: distance to look for larger responses.
// returns: image with only local-maxima responses within w pixels.
image nms_image(image im, int w)
{
    image r = copy_image(im);
    // perform NMS on the response map.
    // for every pixel in the image:
    //     for neighbors within w:
    //         if neighbor response greater than pixel response:
    //             set response to be very low (I use -999999 [why not 0??])
    for(int x=w; x<r.w-w; x++){
        for(int y=w; y<r.h-w; y++){
            int idc = 0;
            for(int m=x-w; m<=x+w; m++) {
                for(int n=y-w; n<=y+w; n++){
                    if(get_pixel(r, m, n, 0) > get_pixel(r, x, y, 0)){
                        set_pixel(r, x, y, 0, -999999);
                        idc = 1;
                        break;
                    }
                }
                if(idc == 1)
                    break;
            }
        }
    }
    return r;
}





// Perform harris corner detection and extract features from the corners.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
// int *n: pointer to number of corners detected, should fill in.
// returns: array of descriptors of the corners in the image.
descriptor *harris_corner_detector(image im, float sigma, float thresh, int nms, int *n)
{
    // Calculate structure matrix
    image S = structure_matrix(im, sigma);

    // Estimate cornerness
    image R = cornerness_response(S);

    // Run NMS on the responses
    image Rnms = nms_image(R, nms);


    //TODO: count number of responses over threshold
    int count = 0; // change this
    for(int i =0; i < Rnms.w * Rnms.h*Rnms.c;i++){
        if(Rnms.data[i] > thresh){
            count +=1;
        }
    }
    *n = count; // <- set *n equal to number of corners in image.
    descriptor *d = calloc(count, sizeof(descriptor));
    descriptor *points = d;
    //TODO: fill in array *d with descriptors of corners, use describe_index.
    for(int i=0; i<Rnms.w*Rnms.h*Rnms.c; i++){
        if(Rnms.data[i] > thresh){
            *points = describe_index(im, i);
            points += 1;
        }
    }

    free_image(S);
    free_image(R);
    free_image(Rnms);
    return d;
}

// Find and draw corners on an image.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
void detect_and_draw_corners(image im, float sigma, float thresh, int nms)
{
    int n = 0;
    descriptor *d = harris_corner_detector(im, sigma, thresh, nms, &n);
    mark_corners(im, d, n);
}



