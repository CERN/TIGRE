/*-------------------------------------------------------------------------
 * CUDA function for optimized proton CT radiographies
 * The full method is described in Kaser et al.: Integration of proton imaging into the TIGRE toolbox (submitted to ZMP)
 * and based on the method of Collins-Fekete (https://doi.org/10.1088/0031-9155/61/23/8232)
 */
 
/*--------------------------------------------------------------------------
 This file is part of the TIGRE Toolbox
 
 Copyright (c) 2015, University of Bath and 
                     CERN-European Organization for Nuclear Research
                     All rights reserved.

 License:            Open Source under BSD. 
                     See the full license at
                     https://github.com/CERN/TIGRE/blob/master/LICENSE

 Contact:            tigre.toolbox@gmail.com
 Codes:              https://github.com/CERN/TIGRE/
 Coded by:           Stefanie Kaser, Benjamin Kirchmayer 
--------------------------------------------------------------------------*/

#include <cuda.h>
#include "mex.h"
#include <cuda_runtime_api.h>
#include "improvedForwardProjections.hpp"
#include <algorithm>
#include <math.h>

#define cudaCheckErrors(msg) \
do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
                mexPrintf("%s \n",msg);\
                mexErrMsgIdAndTxt("ImprovedForwardProj:",cudaGetErrorString(__err));\
        } \
} while (0)


__device__ int SolvePolynomial(float*x, float a, float b, float c){
    // Calculates real roots of a third-order polynomial function using Vieta's method and Cardano's method
    // We obtain a polynomial of the form x³ + ax² + bx + c = 0 and reduce it to z³+pz+q = 0 
    // Herefore, we have to make a substitution: x = z - a/3
    float p = b - a*a / 3.0;
    float q = 2*a*a*a/27.0 - a*b / 3.0 + c;
    float disc = q*q/4.0 + p*p*p/27.0;
    if(disc > 0){
        float u = cbrt(-0.5*q + sqrt(disc)); 
        float v = cbrt(-0.5*q - sqrt(disc)); 
        x[0] = u + v - a/3.0; // don't forget to substitute back z --> x
        return 1;
    }
    else if(disc == 0 && p == 0){
        x[0] = -a/3.0; // don't forget to substitute back z --> x
        return 1;
    }
    else if(disc == 0 && p != 0){
        x[0] = 3.0*q/p - a/3.0; // don't forget to substitute back z --> x
        x[1] = -3.0*q/(2.0*p) - a/3.0; 
        return 2;
    }
    else{
        x[0] = -sqrt(-4.0 * p / 3.0) * cos(1./3. * acos(-0.5*q*sqrt(-27./(p*p*p))) + pi/3.0) - a/3.0; // don't forget to substitute back z --> x
        x[1] = sqrt(-4.0 * p / 3.0) * cos(1./3. * acos(-0.5*q*sqrt(-27./(p*p*p)))) - a/3.0;
        x[2] = -sqrt(-4.0 * p / 3.0) * cos(1./3. * acos(-0.5*q*sqrt(-27./(p*p*p))) - pi/3.0) - a/3.0;
        return 3;
    }
}

__device__ float cspline(float t, float a, float b, float c, float d){

    return a*(t*t*t) + b*(t*t) + c*t +d;

}

__device__ void SimpleSort(float* arr, int size_arr){
    // Insertion sorting method
    float curr_elem;
    int j;
    
    for (int i=1; i<size_arr; i++){
    
        curr_elem = arr[i];
        j = i-1;    // minimum is zero

        while(j>=0 && curr_elem<arr[j]){
            arr[j+1] = arr[j];
            j = j-1;
        }//j
        arr[j+1] = curr_elem;
    }//i 
  }


__device__ int hullEntryExit(float* HullIntercept, float* position, float* direction, int in_or_out, float* hullparams, float detOff){
  float a = hullparams[0];
  float b = hullparams[1];
  float alpha = hullparams[2];
  float h = hullparams[3];
  float kx = direction[0];
  float dx = position[0] - kx*detOff;
  float pref_z2 = b*b*kx*kx*cos(alpha)*cos(alpha) - 2.0 * b*b*kx*cos(alpha)*sin(alpha) + b*b*sin(alpha)*sin(alpha) \
          + a*a*kx*kx*sin(alpha)*sin(alpha) + 2.0 * a*a*kx*cos(alpha)*sin(alpha) + a*a*cos(alpha)*cos(alpha);

  float pref_z = b*b*2.0*kx*dx*cos(alpha)*cos(alpha) - 2.0*b*b*dx*cos(alpha)*sin(alpha) + \
           a*a*2.0*kx*dx*sin(alpha)*sin(alpha) + 2.0*a*a*dx*cos(alpha)*sin(alpha);

  float pref = b*b*dx*dx*cos(alpha)*cos(alpha) + a*a*dx*dx*sin(alpha)*sin(alpha) - a*a*b*b;

  float p = pref_z/pref_z2;
  float q = pref/pref_z2;
  float disc = (p/2.0) * (p/2.0) - q;
  
  if(disc>0){

    float z_1 = -p/2.0 + sqrt(disc);
    float z_2 = -p/2.0 - sqrt(disc);
    float z_solve;

    if(in_or_out == 1){
      z_solve = min(z_1, z_2);
    }
    else {
      z_solve = max(z_1, z_2);
    }

  	float x_solve = kx*z_solve + dx;

    float ky = direction[1];
    float dy = position[1] - ky*detOff;
    float y_solve = ky*z_solve + dy;

    if(-h/2 <= y_solve && y_solve <= h/2){

    HullIntercept[0] = x_solve;
    HullIntercept[1] = y_solve;
    HullIntercept[2] = z_solve;

    return 0;
    }
    else{
    float z1_h = (1.0/ky) * (0.5*h-dy); 
    float z2_h = (1.0/ky) * (-0.5*h-dy);  

    if(in_or_out == 1){
      z_solve = min(z1_h, z2_h);
      if(dy > 0){y_solve = -h*0.5;}
      else{y_solve = h*0.5;}
      x_solve = kx*z_solve + dx;
    }
    else {
      z_solve = max(z1_h, z2_h);
      if(dy < 0){y_solve = -h*0.5;}
      else{y_solve = h*0.5;}
      x_solve = kx*z_solve + dx;
    }
    
    if(min(z_1, z_2) <= z_solve && z_solve <= max(z_1, z_2)){

    HullIntercept[0] = x_solve;
    HullIntercept[1] = y_solve;
    HullIntercept[2] = z_solve;

    return 0;
    }

    else{return 1;}}
  }
else{return 1;}
}


__device__ int MinMax(float* solutions, float a, float b, float c){
    float p = 2*b/(3*a);
    float q = c / (3*a);
    float disc = 0.25*p*p - q;
    if (disc > 0){
        solutions[0] = -0.5*p + sqrt(disc);
        solutions[1] = -0.5*p - sqrt(disc);
        return 0;
    }
    solutions[0] = -1;
    solutions[1] = -1;
    return 1;
}


__device__ int calcInterceptsLinear(float* LinInterceptsVec, float* start, float* stop, float* direction, float* pix, int maxIntercep, bool* protFlag){
  float boundary;
  int counter = 0;
  int nx, ny;
  nx = int(abs(stop[0] - start[0])/pix[0]);
  ny = int(abs(stop[1] - start[1])/pix[1]);
    if(nx+ny>=maxIntercep){
        *protFlag = false;
        return 1;}
  
  if (int(stop[0]/pix[0]) == int(start[0]/pix[0]) && int(stop[1]/pix[1]) == int(start[1]/pix[1])) {
  *protFlag = true;
  return 0;
  }
          
  if (int(stop[0]/pix[0]) != int(start[0]/pix[0])) {
    float k = direction[0];
    float d = start[0] - k*start[2];
    boundary = trunc( ((stop[0] > start[0]) ? stop[0]:start[0])/pix[0])*pix[0];

    for (int ix=0; ix<nx; ix++){
        if(ix != 0){
          boundary = boundary - pix[0];
        }
        float intercept = (boundary - d) / k;

        if(intercept > start[2] && intercept < stop[2]){
          LinInterceptsVec[ix] = intercept; 
          counter++;
          if (counter >= maxIntercep){
              *protFlag = false;
              return counter;}
        }
    }
  }

  if (int(stop[1]/pix[1]) != int(start[1]/pix[1])) {
    float k = direction[1];
    float d = start[1] - k*start[2];
    boundary = trunc( ((stop[1] > start[1]) ? stop[1]:start[1])/pix[1])*pix[1];
    for (int iy=nx; iy<nx+ny; iy++){
        if(iy != nx){
          boundary = boundary - pix[1];
        }
    float intercept = (boundary - d) / k;
    if(intercept > start[2] && intercept < stop[2]){
      LinInterceptsVec[iy] = intercept; 
      counter++;
      if(counter >= maxIntercep){
          *protFlag = false;
          return counter;}
    }
  }
  }
  int diff = maxIntercep - counter;
  for(int j = 0; j<diff; j++){
    LinInterceptsVec[counter+j] = 2*abs(stop[2]-start[2]); //Just ensure that array Element is larger than total distance                      
  }
  SimpleSort(LinInterceptsVec, maxIntercep);
  for(int j = 0; j<diff; j++){
    LinInterceptsVec[counter+j] = 0; // Set value back to zero (just for safety...)                     
  } 
  *protFlag = true;
  return counter;
}

        


__device__ int calcIntercepts(float* InterceptsVec ,float* a, float* b, \
                      float* c, float* d, float* pos1, float* pixelSize, bool* protFlag, int maxIntercep){
                          
            /*Calculates channel Intercepts and the lengths the proton (ion) has spent in the
              corresponding channel.
              Returns 1 if proton is accepted and 0 if it is rejected due to too many Intercepts
            */
                    
      float oneX, oneY, zeroX, zeroY;
	  zeroX = d[0];
	  oneX = pos1[0];
	  zeroY = d[1];
	  oneY = pos1[1];


          int status, nx, ny;
          float IntercepX[3];
          float IntercepY[3];
          float solutions[2];
          float boundary;
          // counter has to be implemented despite the initial discrimination because one can not state beforehand if
          // the cubic spline has more than one Intercept with the channel boundary
          int counter=0;
        

          //Check how many Intercepts will occur approximately
          int test = MinMax(solutions, a[0], b[0], c[0]);
           if (test == 0){
           if (solutions[0] < 1 && solutions[0] > 0){
               float cand = a[0] * solutions[0]*solutions[0]*solutions[0] + b[0] * solutions[0]*solutions[0] + c[0] * solutions[0] + d[0];
               if (cand > d[0] && cand > pos1[0]){
               (oneX > zeroX) ? oneX:zeroX=cand;
               }
               else if(cand < d[0] && cand < pos1[0]){
                (oneX < zeroX) ? oneX:zeroX=cand;
               }
           }

           if (solutions[1] < 1 && solutions[1] > 0){
               float cand = a[0] * solutions[1]*solutions[1]*solutions[1] + b[0] * solutions[1]*solutions[1] + c[0] * solutions[1] + d[0];
               if (cand > oneX && cand > zeroX){
                (oneX > zeroX) ? oneX:zeroX=cand;
               }
               else if(cand < oneX && cand < zeroX){
                (oneX < zeroX) ? oneX:zeroX=cand;
               }
           }
           }

          
           test = MinMax(solutions, a[1], b[1], c[1]);
           if (test == 0){
           if (solutions[0] < 1 && solutions[0] > 0){
               float cand = a[1] * solutions[0]*solutions[0]*solutions[0] + b[1] * solutions[0]*solutions[0] + c[1] * solutions[0] + d[1];
               if (cand > d[1] && cand > pos1[1]){
               (oneY > zeroY) ? oneY:zeroY=cand;
               }
               else if(cand < d[1] && cand < pos1[1]){
                (oneY < zeroY) ? oneY:zeroY=cand;
               }
           }

           if (solutions[1] < 1 && solutions[1] > 0){
               float cand = a[1] * solutions[1]*solutions[1]*solutions[1] + b[1] * solutions[1]*solutions[1] + c[1] * solutions[1] + d[1];
               if (cand > oneY && cand > zeroY){
                (oneY > zeroY) ? oneY:zeroY=cand;
               }
               else if(cand < oneY && cand < zeroY){
                (oneY < zeroY) ? oneY:zeroY=cand;
               }
           }
           } 

          nx = int(abs(oneX - zeroX) / pixelSize[0]);
          ny = int(abs(oneY - zeroY) / pixelSize[1]);
          if (nx + ny == 0) {
          *protFlag = true;
          return 0;
         }

          if ((nx + ny) <= maxIntercep){ 
          
              if (int(oneX/pixelSize[0]) != int(zeroX/pixelSize[0])) {
                boundary = trunc( ((oneX > zeroX) ? oneX:zeroX)/pixelSize[0])*pixelSize[0];
                for (int ix=0; ix<nx; ix++){
                  if(ix != 0){
                    boundary = boundary - pixelSize[0];
                  }
                  //Start from the largest pixel boundary and propagate to the smallest
                  status = SolvePolynomial(IntercepX, b[0]/a[0], c[0]/a[0], d[0]/a[0] - boundary/a[0]);
                  for (int kx=0; kx < status; kx++ ){
                    if(IntercepX[kx]< 1. && IntercepX[kx] > 0. ){
                      if (counter >=maxIntercep){break;}
                      InterceptsVec[counter] = IntercepX[kx];
                      counter++;
                    }
                  }//kx
                 if (counter >=maxIntercep){break;}     
                }
              }

                if ( int(oneY/pixelSize[1]) != int(zeroY/pixelSize[1])) {
                  boundary = trunc( ((oneY > zeroY) ? oneY:zeroY)/pixelSize[1])*pixelSize[1];
                  for (int iy=0; iy<ny; iy++){ 
                    if(iy != 0){
                        boundary = boundary - pixelSize[1];
                    }
                    //Start from the largest pixel boundary and propagate to the smallest
                    status = SolvePolynomial(IntercepY, b[1]/a[1], c[1]/a[1], d[1]/a[1] - boundary/a[1]);
                    for (int ky=0; ky < status; ky++ ){
                      if ((IntercepY[ky]< 1.) &&  (IntercepY[ky] > 0.) ){
                        if (counter >=maxIntercep){break;}
                        InterceptsVec[counter] = IntercepY[ky];
                        counter++;
                      }
                     }//ky
                    if (counter >=maxIntercep){break;}
                    }
                  }

                  if (counter >= maxIntercep){ // || counter == 0){ 
                    *protFlag = false;
                    return counter;
                  }else{
                      

                    int diff = maxIntercep - counter;
                    for(int j = 0; j<diff; j++){
                        InterceptsVec[counter+j] = 2. + (float)j; //Just ensure that array Element is larger than 1                        
                      }     
                   
                    SimpleSort(InterceptsVec, maxIntercep);
                    *protFlag = true;
                    return counter;
                  }

          }else{
          // Too many channel Intercepts - Proton neglected 
          // Discrimination is implemented to neglect protons with large entry angles
          // and to reduce the size of the array that has to be allocated for each thread
          *protFlag = false;
          return counter;
          }
        }


__global__ void ParticleKernel(float* dhist1, float* dhist2, float* devicePosIn, float* devicePosOut, float* devicedirIn, \
                               float* devicedirOut ,float* p_wepl,int* numOfEntries, int* detectSizeX, int* detectSizeY, \
                               float* pix, float* detectDistIn, float* detectDistOut, float *ein, float *hull, float *reject){
    /*Calculate Spline Parameters
    c = deviceDirIn / d = devicePosIn (pos0)
    */
    
    // int customsize = int(50/(*pixelSize));
    /*float *tInterceptsVec;  ---> this is too slow! 7 s instead of 1.5 s
    tInterceptsVec = new float[customsize]; 
    delete[] tInterceptsVec;*/
    /*float *ptr; ---> this is too slow! 7.3s instead of 1.5 s
    ptr = (float*) malloc(customsize * sizeof(float));
    free(ptr);*/
            
    unsigned int protonIndex = blockIdx.x*blockDim.x  + threadIdx.x;
    float dimX, dimY, lk, lenX, lenY;
    float lenZ = abs(*detectDistIn) + abs(*detectDistOut);
    dimX = (float) *detectSizeX;
    dimY = (float) *detectSizeY;

    //Dereference input parameters
    int entries, dSizeX, dSizeY;
    // float pix;
    
    entries = *numOfEntries;
    dSizeX = *detectSizeX;
    dSizeY = *detectSizeY;
    // pix = *pixelSize;
            
            
    if(hull[3] == 0){
    lenX = sqrt((devicePosOut[protonIndex] - devicePosIn[protonIndex]) * (devicePosOut[protonIndex] - devicePosIn[protonIndex]) \
            + lenZ*lenZ); 
    lenY = sqrt((devicePosOut[protonIndex + entries] - devicePosIn[protonIndex + entries]) * (devicePosOut[protonIndex + entries] - devicePosIn[protonIndex + entries]) \
            + lenZ*lenZ);
   
    float lambda0, lambda1, ref_wepl;
    ref_wepl = 10 * 0.00244 * powf(*ein, 1.75);
    lambda0 = 1.01 + 0.43 * (p_wepl[protonIndex]/ref_wepl) * (p_wepl[protonIndex]/ref_wepl);
    lambda1 = 0.99 - 0.46 * (p_wepl[protonIndex]/ref_wepl) * (p_wepl[protonIndex]/ref_wepl);

    float a[2], b[2], c[2], d[2], pos1[2];
    
    //Allocate memory for all pointers
    // Calculate optimized xdir_in
    devicedirIn[protonIndex] = devicedirIn[protonIndex] \
            / sqrt(devicedirIn[protonIndex]*devicedirIn[protonIndex] + 1.0);    //  ... dz = 1!
    devicedirIn[protonIndex] = devicedirIn[protonIndex] * lenX * lambda0;
    
    // Calculate optimized ydir_in
    devicedirIn[protonIndex + entries] = devicedirIn[protonIndex + entries] \
            / sqrt(devicedirIn[protonIndex + entries]*devicedirIn[protonIndex + entries] + 1.0);  // ... dz = 1!
    devicedirIn[protonIndex + entries] = devicedirIn[protonIndex + entries] * lenY * lambda0;
    
    // Calculate optimized xdir_out
    devicedirOut[protonIndex] = devicedirOut[protonIndex] \
            / sqrt(devicedirOut[protonIndex]*devicedirOut[protonIndex] + 1.0); //  ... dz = 1!
    devicedirOut[protonIndex] = devicedirOut[protonIndex] * lenX * lambda1;
    
    // Calculate optimized ydir_out
    devicedirOut[protonIndex + entries] = devicedirOut[protonIndex + entries] \
            / sqrt(devicedirOut[protonIndex + entries]*devicedirOut[protonIndex + entries] + 1.0); // ... dz = 1!
    devicedirOut[protonIndex + entries] = devicedirOut[protonIndex + entries] * lenY * lambda1;
            
    // Calculate spline parameters
    a[0] = devicePosIn[protonIndex]*2. + devicedirIn[protonIndex] - 2.*devicePosOut[protonIndex] + devicedirOut[protonIndex];
    a[1] = devicePosIn[protonIndex + entries]*2. + devicedirIn[protonIndex + entries] - \
    2.*devicePosOut[protonIndex + entries] +  devicedirOut[protonIndex + entries];

    b[0] = -3.*devicePosIn[protonIndex] -2.*devicedirIn[protonIndex] + 3.*devicePosOut[protonIndex] - devicedirOut[protonIndex];
    b[1] = -3.*devicePosIn[protonIndex + entries] -2.* devicedirIn[protonIndex + entries] \
    + 3.*devicePosOut[protonIndex + entries] - devicedirOut[protonIndex + entries];

    c[0] = devicedirIn[protonIndex];
    c[1] = devicedirIn[protonIndex + entries];

    d[0] = devicePosIn[protonIndex];
    d[1] = devicePosIn[protonIndex + entries];

    pos1[0] = devicePosOut[protonIndex];
    pos1[1] = devicePosOut[protonIndex + entries];
    
    /* --------------------------------------------------------------------------------- */
    /* ------------------------ Start without Hull (CS only)  -------------------------- */
    /* --------------------------------------------------------------------------------- */ 
    int count;
    bool status = false;
    float InterceptsVec[vecSizeCS] = {0}; 
    
    count = calcIntercepts(InterceptsVec, a, b, c, d, pos1, pix, &status, vecSizeCS);
       
    if (status) { 
        int indX, indY, linInd;
        float tOld = 0.0;
         if (count==0){ 
           indX = int(pos1[0]/pix[0]+dimX/2.); // REPLACE: pos1 by pos0
           indY = int(pos1[1]/pix[1]+dimY/2.);

           if ((0 <= indX) && (indX < dSizeX) && (0 <= indY) && (indY < dSizeY)){ 
               linInd = indY + indX*(dSizeY);  
               atomicAdd(&dhist1[linInd], p_wepl[protonIndex]);
               atomicAdd(&dhist2[linInd], 1.0f);
           }

         } 
         else{
            for(int i= 0; i<=count; i++){
              lk = (InterceptsVec[i]- tOld)*lenZ;
              if(tOld == 0){
                indX = int(d[0]/pix[0] +dimX/2);
                indY = int(d[1]/pix[1] +dimY/2);
                linInd = indY + indX*(dSizeY); 

                if ((0 <= indX) && (indX < dSizeX) && (0 <= indY) && (indY < dSizeY)){
                    linInd = indY + indX*(dSizeY);
                    atomicAdd(&dhist1[linInd], (lk/lenZ)*(lk/lenZ)*p_wepl[protonIndex]);
                    atomicAdd(&dhist2[linInd], (lk/lenZ)*(lk/lenZ));
                }
                tOld = InterceptsVec[i];

              }else if(i == count){
                lk = lenZ - InterceptsVec[i-1]*lenZ;
                indX = int(pos1[0]/pix[0] +dimX/2);
                indY = int(pos1[1]/pix[1] +dimY/2);

                if ((0 <= indX) && (indX < dSizeX) && (0 <= indY) && (indY < dSizeY)){
                    linInd = indY + indX*(dSizeY); 
                    atomicAdd(&dhist1[linInd], (lk/lenZ)*(lk/lenZ)*p_wepl[protonIndex]);
                    atomicAdd(&dhist2[linInd], (lk/lenZ)*(lk/lenZ));
                }

              }else{
                indX = int(cspline(InterceptsVec[i] - eps, a[0], b[0], c[0], d[0])/pix[0] +dimX/2);
                indY = int(cspline(InterceptsVec[i] - eps, a[1], b[1], c[1], d[1])/pix[1] +dimY/2);

                if ((0 <= indX) && (indX < dSizeX) && (0 <= indY) && (indY < dSizeY)){
                    linInd = indY + indX*(dSizeY); 
                    atomicAdd(&dhist1[linInd], (lk/lenZ)*(lk/lenZ)*p_wepl[protonIndex]);
                    atomicAdd(&dhist2[linInd], (lk/lenZ)*(lk/lenZ));
                }
                tOld = InterceptsVec[i];
              }

            }//i
         }//if - Intercepts
     }
    else{
        atomicAdd(reject, 1.0);
    }
/* ------------------------ End no Hull calculation (CS only)  -------------------------- */
    }

else{
    // WEIGHTING FACTORS FOR CHANNELS I 
    float weight_air_in = 0.00479; 
    float weight_air_out = 0.00479; 

    float HullIn[3], HullOut[3], initpos[3], exitpos[3];  
    float initdir[2], exitdir[2]; 
            
    initpos[0] = devicePosIn[protonIndex];
    initpos[1] = devicePosIn[protonIndex + entries];
    initpos[2] = *detectDistIn;

    exitpos[0] = devicePosOut[protonIndex];
    exitpos[1] = devicePosOut[protonIndex + entries];
    exitpos[2] = *detectDistOut;

    initdir[0] = devicedirIn[protonIndex];
    initdir[1] = devicedirIn[protonIndex + entries];

    exitdir[0] = devicedirOut[protonIndex];
    exitdir[1] = devicedirOut[protonIndex + entries];

    int check = hullEntryExit(HullIn, initpos, initdir, 1, hull, *detectDistIn);

    if(check == 0){
        check = hullEntryExit(HullOut, exitpos, exitdir, 0, hull, *detectDistOut);
    }

    if(check == 0 && HullOut[2] > HullIn[2]){            
        /* --------------------------------------------------------------------------------- */
        /* ------------------------ Start with Hull + SL outside  -------------------------- */
        /* --------------------------------------------------------------------------------- */
        const int hullIntercep = int(vecSizeCS);  
        const int airIntercepIn = int(vecSizeIn);   
        const int airIntercepOut = int(vecSizeOut);   
        bool status1 = false;
        bool status2 = false; 
        bool status3 = false;
        
        int countIn, countHull, countOut;
        float InterceptsVecOut[airIntercepOut] = {0}; 
        float InterceptsVecIn[airIntercepIn] = {0};
        float InterceptsVecHull[hullIntercep] = {0}; 
        lenX = sqrt((HullOut[0] - HullIn[0])*(HullOut[0] - HullIn[0]) + (HullOut[2] - HullIn[2])*(HullOut[2] - HullIn[2])); 
        lenY = sqrt((HullOut[1] - HullIn[1])*(HullOut[1] - HullIn[1]) + (HullOut[2] - HullIn[2])*(HullOut[2] - HullIn[2]));

        countIn = calcInterceptsLinear(InterceptsVecIn, initpos, HullIn, initdir, pix, airIntercepIn, &status1);
        countOut = calcInterceptsLinear(InterceptsVecOut, HullOut, exitpos, exitdir, pix, airIntercepOut, &status2);

        /* ------------ CUBIC SPLINE PREPARATIONS ---------------- */
        float lambda0, lambda1, ref_wepl;
        ref_wepl = 10 * 0.00244 * powf(*ein, 1.75);
        lambda0 = 1.01 + 0.43 * (p_wepl[protonIndex]/ref_wepl)*(p_wepl[protonIndex]/ref_wepl);
        lambda1 = 0.99 - 0.46 * (p_wepl[protonIndex]/ref_wepl)*(p_wepl[protonIndex]/ref_wepl);

        float a[2], b[2], c[2], d[2], pos1[2];

        //Allocate memory for all pointers
        // Calculate optimized xdir_in
	devicedirIn[protonIndex] = devicedirIn[protonIndex] \
                / sqrt(devicedirIn[protonIndex]*devicedirIn[protonIndex] + 1.0);    // ... dz = 1! 
        devicedirIn[protonIndex] = devicedirIn[protonIndex] * lenX * lambda0;

        // Calculate optimized ydir_in
	devicedirIn[protonIndex + entries] = devicedirIn[protonIndex + entries] \
                / sqrt(devicedirIn[protonIndex + entries]*devicedirIn[protonIndex + entries] + 1.0);   // ... dz = 1! 
        devicedirIn[protonIndex + entries] = devicedirIn[protonIndex + entries] * lenY * lambda0;

        // Calculate optimized xdir_out
	devicedirOut[protonIndex] = devicedirOut[protonIndex] \
                / sqrt(devicedirOut[protonIndex]*devicedirOut[protonIndex] + 1.0); // ... dz = 1!
        devicedirOut[protonIndex] = devicedirOut[protonIndex] * lenX * lambda1;

        // Calculate optimized ydir_out
	devicedirOut[protonIndex + entries] = devicedirOut[protonIndex + entries] \
                / sqrt(devicedirOut[protonIndex + entries]*devicedirOut[protonIndex + entries] + 1.0); // ... dz = 1!
        devicedirOut[protonIndex + entries] = devicedirOut[protonIndex + entries] * lenY * lambda1;

        // Calculate spline parameters
        a[0] = HullIn[0]*2. + devicedirIn[protonIndex] - 2.*HullOut[0] + devicedirOut[protonIndex];
        a[1] = HullIn[1]*2. + devicedirIn[protonIndex + entries] - \
        2.*HullOut[1] +  devicedirOut[protonIndex + entries];

        b[0] = -3.*HullIn[0] -2.*devicedirIn[protonIndex] + 3.*HullOut[0] - devicedirOut[protonIndex];
        b[1] = -3.*HullIn[1] -2.* devicedirIn[protonIndex + entries] \
        + 3.*HullOut[1] - devicedirOut[protonIndex + entries];

        c[0] = devicedirIn[protonIndex];
        c[1] = devicedirIn[protonIndex + entries];

        d[0] = HullIn[0];
        d[1] = HullIn[1];

        pos1[0] = HullOut[0];
        pos1[1] = HullOut[1];

        countHull = calcIntercepts(InterceptsVecHull, a, b, c, d, pos1, pix, &status3, hullIntercep);
        /* -------------------- End CS Preparations! -------------- */

        if(status1 && status2 && status3){
        float tOld = initpos[2];
        int indX, indY, linInd;

        // WEIGHTING FACTORS FOR CHANNELS II
        float weight_water = 1;  // p_wepl[protonIndex]/(len_b*weight_air_in);

        // ---------------------------------------- Start with SL from detector to hull
        if (countIn == 0){
        indX = int(initpos[0]/pix[0] + dimX/2.);
        indY = int(initpos[1]/pix[1] + dimY/2.);
        lk = HullIn[2] - initpos[2];
        if ((0 <= indX) && (indX < dSizeX) && (0 <= indY) && (indY < dSizeY)){ 
           linInd = indY + indX*(dSizeY);  
           atomicAdd(&dhist1[linInd], weight_air_in*(lk/lenZ)*(lk/lenZ)*p_wepl[protonIndex]);
           atomicAdd(&dhist2[linInd], weight_air_in*(lk/lenZ)*(lk/lenZ));
            }
        }

        else{
        for(int i= 0; i<=countIn; i++){
           lk = InterceptsVecIn[i] - tOld;
           if(i == 0){
             indX = int(initpos[0]/pix[0] + dimX/2.);
             indY = int(initpos[1]/pix[1] + dimY/2.);
             if ((0 <= indX) && (indX < dSizeX) && (0 <= indY) && (indY < dSizeY) && (0 < lk) && (lk < (HullIn[2]-initpos[2]))){
             linInd = indY + indX*(dSizeY);
             atomicAdd(&dhist1[linInd], weight_air_in*(lk/lenZ)*(lk/lenZ)*p_wepl[protonIndex]);
             atomicAdd(&dhist2[linInd], weight_air_in*(lk/lenZ)*(lk/lenZ));
             tOld = InterceptsVecIn[i];
             }   
           }
           else if(i == countIn){
             lk = HullIn[2] - InterceptsVecIn[i-1];
             indX = int(HullIn[0]/pix[0] + dimX/2.);
             indY = int(HullIn[1]/pix[1] + dimY/2.);
             if ((0 <= indX) && (indX < dSizeX) && (0 <= indY) && (indY < dSizeY) && (0 < lk) && (lk < (HullIn[2]-initpos[2]))){
             linInd = indY + indX*(dSizeY);
             atomicAdd(&dhist1[linInd], weight_air_in*(lk/lenZ)*(lk/lenZ)*p_wepl[protonIndex]);
             atomicAdd(&dhist2[linInd], weight_air_in*(lk/lenZ)*(lk/lenZ));
             }
           }

           else{
             indX = int(((initdir[0]*(InterceptsVecIn[i]-eps) + (initpos[0] - initdir[0] * initpos[2])))/pix[0] + dimX/2.);
             indY = int(((initdir[1]*(InterceptsVecIn[i]-eps) + (initpos[1] - initdir[1] * initpos[2])))/pix[1] + dimY/2.);
             if ((0 <= indX) && (indX < dSizeX) && (0 <= indY) && (indY < dSizeY) && (0 < lk) && (lk < (HullIn[2]-initpos[2]))){
             linInd = indY + indX*(dSizeY);
             atomicAdd(&dhist1[linInd], weight_air_in*(lk/lenZ)*(lk/lenZ)*p_wepl[protonIndex]);
             atomicAdd(&dhist2[linInd], weight_air_in*(lk/lenZ)*(lk/lenZ));
             tOld = InterceptsVecIn[i];
             }
            }
           }
          }   // end else
        // --------------------------- CS within hull
        
             tOld = 0.0;
             if (countHull==0){ 
               indX = int(HullIn[0]/pix[0] + dimX/2.); 
               indY = int(HullIn[1]/pix[1] + dimY/2.);
               lk = HullOut[2] - HullIn[2];
               if ((0 <= indX) && (indX < dSizeX) && (0 <= indY) && (indY < dSizeY)){ 
                   linInd = indY + indX*(dSizeY);  
                   atomicAdd(&dhist1[linInd], weight_water*(lk/lenZ)*(lk/lenZ)*p_wepl[protonIndex]);
                   atomicAdd(&dhist2[linInd], weight_water*(lk/lenZ)*(lk/lenZ));
               }

             } else{
                for(int i= 0; i<=countHull; i++){
                  lk = (InterceptsVecHull[i] - tOld)*(HullOut[2] - HullIn[2]);
                  if(tOld == 0){
                    indX = int(d[0]/pix[0] + dimX/2.);
                    indY = int(d[1]/pix[1] + dimY/2.);
                    linInd = indY + indX*(dSizeY); 

                    if ((0 <= indX) && (indX < dSizeX) && (0 <= indY) && (indY < dSizeY) && (0 < lk) && (lk < (HullOut[2]-HullIn[2]))){
                        linInd = indY + indX*(dSizeY);
                        atomicAdd(&dhist1[linInd], weight_water*(lk/lenZ)*(lk/lenZ)*p_wepl[protonIndex]);
                        atomicAdd(&dhist2[linInd], weight_water*(lk/lenZ)*(lk/lenZ));
                    }
                    tOld = InterceptsVecHull[i]; 

                  }else if(i == countHull){
                    lk = (HullOut[2] - HullIn[2]) - InterceptsVecHull[i-1]*(HullOut[2] - HullIn[2]);
                    indX = int(pos1[0]/pix[0] + dimX/2.);
                    indY = int(pos1[1]/pix[1] + dimY/2.);

                    if ((0 <= indX) && (indX < dSizeX) && (0 <= indY) && (indY < dSizeY) && (0 < lk) && (lk < (HullOut[2]-HullIn[2]))){
                        linInd = indY + indX*(dSizeY); 
                        atomicAdd(&dhist1[linInd], weight_water*(lk/lenZ)*(lk/lenZ)*p_wepl[protonIndex]);
                        atomicAdd(&dhist2[linInd], weight_water*(lk/lenZ)*(lk/lenZ));
                    }

                  }else{
                    indX = int(cspline(InterceptsVecHull[i] -eps, a[0], b[0], c[0], d[0])/pix[0] + dimX/2.);
                    indY = int(cspline(InterceptsVecHull[i] -eps, a[1], b[1], c[1], d[1])/pix[1] + dimY/2.);

                    if ((0 <= indX) && (indX < dSizeX) && (0 <= indY) && (indY < dSizeY) && (0 < lk) && (lk < (HullOut[2]-HullIn[2]))){
                        linInd = indY + indX*(dSizeY); 
                        atomicAdd(&dhist1[linInd], weight_water*(lk/lenZ)*(lk/lenZ)*p_wepl[protonIndex]);
                        atomicAdd(&dhist2[linInd], weight_water*(lk/lenZ)*(lk/lenZ));
                    }
                    tOld = InterceptsVecHull[i];
                  }

             }//i
         }

        // --------------------------- SL from hull to detector
        tOld = HullOut[2];
        if (countOut == 0){
        indX = int(exitpos[0]/pix[0] + dimX/2.);
        indY = int(exitpos[1]/pix[1] + dimY/2.);
        lk = exitpos[2] - HullOut[2];
        if ((0 <= indX) && (indX < dSizeX) && (0 <= indY) && (indY < dSizeY)){ 
           linInd = indY + indX*(dSizeY);  
           atomicAdd(&dhist1[linInd], weight_air_out*(lk/lenZ)*(lk/lenZ)*p_wepl[protonIndex]);
           atomicAdd(&dhist2[linInd], weight_air_out*(lk/lenZ)*(lk/lenZ));
            }
        }

        else{
        for(int i= 0; i<=countOut; i++){
           lk = abs(InterceptsVecOut[i] - tOld);
           if(i == 0){
             indX = int(HullOut[0]/pix[0] + dimX/2.);
             indY = int(HullOut[1]/pix[1] + dimY/2.);
             if ((0 <= indX) && (indX < dSizeX) && (0 <= indY) && (indY < dSizeY) && (0 < lk) && (lk < (exitpos[2]-HullOut[2]))){
             linInd = indY + indX*(dSizeY);  
             atomicAdd(&dhist1[linInd], weight_air_out*(lk/lenZ)*(lk/lenZ)*p_wepl[protonIndex]);
             atomicAdd(&dhist2[linInd], weight_air_out*(lk/lenZ)*(lk/lenZ));
             tOld = InterceptsVecOut[i];
             }   
           }
           else if(i == countOut){
             lk = exitpos[2] - InterceptsVecOut[i-1];
             indX = int(exitpos[0]/pix[0] + dimX/2.);
             indY = int(exitpos[1]/pix[1] + dimY/2.);
             if ((0 <= indX) && (indX < dSizeX) && (0 <= indY) && (indY < dSizeY) && (0 < lk) && (lk < (exitpos[2]-HullOut[2]))){
             linInd = indY + indX*(dSizeY);
             atomicAdd(&dhist1[linInd], weight_air_out*(lk/lenZ)*(lk/lenZ)*p_wepl[protonIndex]);
             atomicAdd(&dhist2[linInd], weight_air_out*(lk/lenZ)*(lk/lenZ));
             }
           }

           else{
             indX = int(((exitdir[0]*(InterceptsVecOut[i]-eps) + (HullOut[0] - exitdir[0] * HullOut[2])))/pix[0] + dimX/2.);
             indY = int(((exitdir[1]*(InterceptsVecOut[i]-eps) + (HullOut[1] - exitdir[1] * HullOut[2])))/pix[1] + dimY/2.);
             if ((0 <= indX) && (indX < dSizeX) && (0 <= indY) && (indY < dSizeY) && (0 < lk) && (lk < (exitpos[2]-HullOut[2]))){
             linInd = indY + indX*(dSizeY);
             atomicAdd(&dhist1[linInd], weight_air_out*(lk/lenZ)*(lk/lenZ)*p_wepl[protonIndex]);
             atomicAdd(&dhist2[linInd], weight_air_out*(lk/lenZ)*(lk/lenZ));
             tOld = InterceptsVecOut[i];
             }
            }
           }
          }   // end else
        }
        else{
        atomicAdd(reject, 1.0);
    }

        /* --------------------------- End Hull + SL outside ------------------------------- */
        
        }  

    else{   
    
    /* --------------------------------------------------------------------------------- */
    /* ----------------------------- Start with SL only!  ------------------------------ */
    /* --------------------------------------------------------------------------------- */ 
    int count;
    bool status = false;
    float InterceptsVec[vecSizeCS] = {0}; 
    
    float initpos[3], exitpos[3]; 
    float mydir[2]; 
    initpos[0] = devicePosIn[protonIndex];
    initpos[1] = devicePosIn[protonIndex + entries];
    initpos[2] = *detectDistIn;
    exitpos[0] = devicePosOut[protonIndex];
    exitpos[1] = devicePosOut[protonIndex + entries];
    exitpos[2] = *detectDistOut;

    mydir[0] = (exitpos[0] - initpos[0])/lenZ;
    mydir[1] = (exitpos[1] - initpos[1])/lenZ;  // dz = 1
    count = calcInterceptsLinear(InterceptsVec, initpos, exitpos, mydir, pix, vecSizeCS, &status);
            
       
    if (status) { 
        int indX, indY, linInd;
        float tOld = initpos[2];
         if (count==0){ 
           indX = int(initpos[0]/pix[0] + dimX/2.); 
           indY = int(initpos[1]/pix[1] + dimY/2.);

           if ((0 <= indX) && (indX < dSizeX) && (0 <= indY) && (indY < dSizeY)){
               linInd = indY + indX*(dSizeY);  
               atomicAdd(&dhist1[linInd], weight_air_out*p_wepl[protonIndex]);
               atomicAdd(&dhist2[linInd], weight_air_out*1.0f);
           }

         } else{
            for(int i= 0; i<=count; i++){
              lk = InterceptsVec[i] - tOld;
              if(tOld == initpos[2]){
                indX = int(initpos[0]/pix[0] + dimX/2.);
                indY = int(initpos[1]/pix[1] + dimY/2.);
                linInd = indY + indX*(dSizeY); 

                if ((0 <= indX) && (indX < dSizeX) && (0 <= indY) && (indY < dSizeY) && (0 < lk) && (lk < lenZ)){
                    linInd = indY + indX*(dSizeY);
                    atomicAdd(&dhist1[linInd], weight_air_out*(lk/lenZ)*(lk/lenZ)*p_wepl[protonIndex]);
                    atomicAdd(&dhist2[linInd], weight_air_out*(lk/lenZ)*(lk/lenZ));
                }
                tOld = InterceptsVec[i];

              }else if(i == count){
                lk = exitpos[2] - InterceptsVec[i-1];
                indX = int(exitpos[0]/pix[0] + dimX/2.);
                indY = int(exitpos[1]/pix[1] + dimY/2.);

                if ((0 <= indX) && (indX < dSizeX) && (0 <= indY) && (indY < dSizeY) && (0 < lk) && (lk < lenZ)){
                    linInd = indY + indX*(dSizeY); 
                    atomicAdd(&dhist1[linInd], weight_air_out*(lk/lenZ)*(lk/lenZ)*p_wepl[protonIndex]);
                    atomicAdd(&dhist2[linInd], weight_air_out*(lk/lenZ)*(lk/lenZ));
                }

              }else{
                indX = int(((mydir[0]*(InterceptsVec[i]-eps) + (initpos[0] - mydir[0] * (initpos[2]))))/pix[0] + dimX/2.);
                indY = int(((mydir[1]*(InterceptsVec[i]-eps) + (initpos[1] - mydir[1] * (initpos[2]))))/pix[1] + dimY/2.);

                if ((0 <= indX) && (indX < dSizeX) && (0 <= indY) && (indY < dSizeY) && (0 < lk) && (lk < lenZ)){
                    linInd = indY + indX*(dSizeY); 
                    atomicAdd(&dhist1[linInd], weight_air_out*(lk/lenZ)*(lk/lenZ)*p_wepl[protonIndex]);
                    atomicAdd(&dhist2[linInd], weight_air_out*(lk/lenZ)*(lk/lenZ));
                }
                tOld = InterceptsVec[i];
              }

            } //i
         }//if - Intercepts
     }
    else{
        // *reject += 1;
        atomicAdd(reject, 1.0);
    }
    /* ------------------------------ End SL only! ------ -------------------------- */
    }   
   }
}

__global__ void sumHist(float* hist, float* histNorm){
    
    unsigned int index = blockIdx.x*blockDim.x  + threadIdx.x;
    hist[index] = hist[index]/histNorm[index];
}

__host__ void ParticleProjections(float * outProjection, float* posIn, float* posOut, float* dirIn, float* dirOut, \
                                  float* p_wepl, int numOfEntries, int detectSizeX, int detectSizeY, float* pixelSize, \
                                  float detectDistIn, float detectDistOut, float ein, float* ch_param){

    /*
    Detect Size = 400x400
    Prepare Input for GPU*/

    const int sizeInputs = 2*numOfEntries*sizeof(float);
    const int detectorMem = detectSizeX*detectSizeY*sizeof(float);
    float reject = 0.0;

    float *dPosIn, *dPosOut, *ddirIn, *ddirOut, *dhist1, *dhist2, *d_wepl, *dHull;
    int *dnumEntries, *ddetectorX, *ddetectorY;
    float *dpixelSize, *dDetectDistIn, *dDetectDistOut, *dEin, *dReject;

    float *hist1, *hist2;
    hist1 = new float[detectSizeX*detectSizeY];
    hist2 = new float[detectSizeX*detectSizeY];
    for(int i = 0; i<detectSizeX*detectSizeY; i++){
        hist1[i] = 0.f;
        hist2[i]= 0.f;
    
    }

    //Allocate Memory on GPU
    cudaMalloc( (void**) &dPosIn, sizeInputs );
    cudaMalloc( (void**) &dPosOut, sizeInputs );
    cudaMalloc( (void**) &ddirIn, sizeInputs );
    cudaMalloc( (void**) &ddirOut, sizeInputs );
    cudaMalloc( (void**) &d_wepl, numOfEntries*sizeof(float));
    cudaMalloc( (void**) &dhist1, detectorMem );
    cudaMalloc( (void**) &dhist2, detectorMem );
    cudaMalloc( (void**) &dnumEntries, sizeof(int));
    cudaMalloc( (void**) &ddetectorX, sizeof(int));
    cudaMalloc( (void**) &ddetectorY, sizeof(int));
    cudaMalloc( (void**) &dpixelSize, 2*sizeof(float));
    cudaMalloc( (void**) &dDetectDistIn, sizeof(float));
    cudaMalloc( (void**) &dDetectDistOut, sizeof(float));
    cudaMalloc( (void**) &dEin, sizeof(float));
    cudaMalloc( (void**) &dReject, sizeof(float));
    cudaMalloc( (void**) &dHull, 5*sizeof(float));
    cudaError_t _err_alloc = cudaGetLastError();
    mexPrintf("%s \n", cudaGetErrorString(_err_alloc));
    cudaCheckErrors("GPU Allocation failed!");

    //Copy Arrays to GPU
    cudaMemcpy(dPosIn, posIn,sizeInputs ,cudaMemcpyHostToDevice);
    cudaMemcpy(dPosOut, posOut,sizeInputs,cudaMemcpyHostToDevice);
    cudaMemcpy(ddirIn, dirIn,sizeInputs,cudaMemcpyHostToDevice);
    cudaMemcpy(ddirOut, dirOut,sizeInputs,cudaMemcpyHostToDevice);
    cudaMemcpy(d_wepl, p_wepl, numOfEntries*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dnumEntries, &numOfEntries,sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(ddetectorX, &detectSizeX, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(ddetectorY, &detectSizeY, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dpixelSize, pixelSize, 2*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dDetectDistIn, &detectDistIn, sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dDetectDistOut, &detectDistOut, sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dEin, &ein, sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dReject, &reject, sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dHull, ch_param, 5*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dhist1, hist1, detectorMem, cudaMemcpyHostToDevice);
    cudaMemcpy(dhist2, hist2, detectorMem, cudaMemcpyHostToDevice);
    cudaCheckErrors("Host to device transport failed!");



    dim3 grid(floor(numOfEntries/maxthreads),1,1);
    dim3 block(maxthreads,1,1);

    
    ParticleKernel<<<grid, block>>>(dhist1, dhist2, dPosIn, dPosOut, ddirIn, ddirOut, d_wepl, dnumEntries, ddetectorX, ddetectorY, \
            dpixelSize, dDetectDistIn, dDetectDistOut, dEin, dHull, dReject);
    cudaError_t _err = cudaGetLastError();
    mexPrintf("%s \n", cudaGetErrorString(_err));
    cudaCheckErrors("Kernel fail!");
    
    //dim3 grid_sum((int)floor(detectSizeX*detectSizeY/64),1,1);
    //dim3 block_sum(64,1,1);
    //sumHist<<<grid_sum, block_sum>>>(dhist1, dhist2);
        
    //Copy result from device to host
    //cudaMemcpy(outProjection, dhist1,detectorMem ,cudaMemcpyDeviceToHost);
    cudaMemcpy(hist1, dhist1,detectorMem ,cudaMemcpyDeviceToHost);
    cudaMemcpy(hist2, dhist2,detectorMem ,cudaMemcpyDeviceToHost);
    cudaMemcpy(&reject, dReject,sizeof(float) ,cudaMemcpyDeviceToHost);
    //cudaError_t _errcp = cudaGetLastError();
    //mexPrintf("%s \n", cudaGetErrorString(_errcp));
    cudaCheckErrors("Device to host transport failed!");
    
    for(int j = 0; j<detectSizeX*detectSizeY; j++){
        outProjection[j] = hist1[j]/hist2[j]; 
    }

    std::cout << "Particles rejected [%]: " << 100*reject/numOfEntries << std::endl;

    cudaFree(dPosIn);
    cudaFree(dPosOut);
    cudaFree(ddirIn);
    cudaFree(ddirOut);
    cudaFree(dhist1);
    cudaFree(dhist2);
    cudaFree(d_wepl);
    cudaFree(dnumEntries);
    cudaFree(ddetectorX);
    cudaFree(ddetectorY);
    cudaFree(dpixelSize);
    cudaFree(dDetectDistIn);
    cudaFree(dDetectDistOut);
    cudaFree(dEin);
    cudaFree(dReject);
    cudaFree(dHull);

    delete(hist1);
    delete(hist2);
    // delete(&reject);


}
