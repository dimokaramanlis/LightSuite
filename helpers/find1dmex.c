/*************************************************************************
 * MEX ROUTINE FIND1DMEX.C
 * B = find1dmex(A, count)
 *
 * A is 3-d arrays
 * FIND along the second dimension
 * COUNT is scacar, >0: find most first COUNT
 *                  <0: find most last -COUNT
 *
 * Method: there is no copy of array
 *
 * Compilation:
 *  >> mex -O -v find1dmex.c % add -largeArrayDims on 64-bit computer
 *
 * Author Bruno Luong <brunoluong@yahoo.com>
 * Last update: 27/June/2009 
 * 19-May-2010: change count to mwSignedIndex 
 ************************************************************************/

#include "mex.h"
#include "matrix.h"

/* Uncomment this on older Matlab version where size_t has not been
 * defined */
/*
 * #define mwSize int
 * #define size_t int
 */

/* Define correct type depending on platform 
  You might have to modify here depending on your compiler */
#if defined(_MSC_VER) || defined(__BORLANDC__)
typedef __int64 int64;
typedef __int32 int32;
typedef __int16 int16;
typedef __int8 int08;
#else /* LINUX + LCC, CAUTION: not tested by the author */
typedef long long int int64;
typedef long int int32;
typedef short int16;
typedef char int08;
#endif

/* Macro use for find engine */
/* pA is a pointer to a specific type of data, specified by pAtype */
#define FIRST_ENGINE(pA, pAtype) \
    for (j=0; j<n; j++) \
    { \
        pA = (pAtype*)PrA + j*km; \
        pB = PrB + j*kc; \
        for (i=0; i<k; i++) { \
            c = 0; \
            for (p=0; p<m; p++) \
                if (pA[p*k] != 0) { \
                    pB[c*k] = (double)(p+1); \
                    if ((++c)==count) break; \
                } \
            pA++; \
            pB++; \
        } \
    } \
    break;

/* LAST engine, the only difference is the reverse for-loop on p */
#define LAST_ENGINE(pA, pAtype) \
    for (j=0; j<n; j++) \
    { \
        pA = (pAtype*)PrA + j*km; \
        pB = PrB + j*kc; \
        for (i=0; i<k; i++) { \
            c = 0; \
            for (p=m; p--;) \
                if (pA[p*k] != 0) { \
                    pB[c*k] = (double)(p+1); \
                    if ((++c)==count) break; \
                } \
            pA++; \
            pB++; \
        } \
    } \
    break;    

/* Engine selection depending of FIRST or LAST option */
#define FIND_ENGINE(pA, pAtype, islast) \
if (islast) \
    {LAST_ENGINE(pA, pAtype)} \
else \
    {FIRST_ENGINE(pA, pAtype)}

/* First engine for complex data, chech real and imaginary parts (pr/pi) */
#define CMPLX_FIRST_ENGINE(pr, pi, pAtype) \
    for (j=0; j<n; j++) \
    { \
        pr = (pAtype*)PrA + j*km; \
        pi = (pAtype*)PiA + j*km; \
        pB = PrB + j*kc; \
        for (i=0; i<k; i++) { \
            c = 0; \
            for (p=0; p<m; p++) \
                if (pr[p*k] != 0 || pi[p*k] != 0) { \
                    pB[c*k] = (double)(p+1); \
                    if ((++c)==count) break; \
                } \
            pr++; pi++; \
            pB++; \
        } \
    } \
    break;

/* LAST engine for complex data, the only difference is the reverse
 * for-loop on p */    
#define CMPLX_LAST_ENGINE(pr, pi, pAtype) \
    for (j=0; j<n; j++) \
    { \
        pr = (pAtype*)PrA + j*km; \
        pi = (pAtype*)PiA + j*km; \
        pB = PrB + j*kc; \
        for (i=0; i<k; i++) { \
            c = 0; \
            for (p=m; p--;) \
                if (pr[p*k] != 0 || pi[p*k] != 0) { \
                    pB[c*k] = (double)(p+1); \
                    if ((++c)==count) break; \
                } \
            pr++; pi++; \
            pB++; \
        } \
    } \
    break;

/* Complex engine selection depending of FIRST or LAST option */    
#define CMPLX_FIND_ENGINE(pr, pi, pAtype, islast) \
if (islast) \
    {CMPLX_LAST_ENGINE(pr, pi, pAtype)} \
else \
    {CMPLX_FIRST_ENGINE(pr, pi, pAtype)}    

/* Define the name for Input/Output ARGUMENTS */
#define A prhs[0]
#define COUNT prhs[1]
#define B plhs[0]

/* Gateway of find1dmex */
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {

    mxClassID ClassID;
    mwSize m, n, k;
    const mwSize* dimA;
    mwSize dimB[3];
    mwSize i, j, p, c;
    mwSignedIndex count;
    mwIndex km, kc;
    double *prAdouble, *piAdouble, *PrA, *PiA;
    double *pB, *PrB;
    float *prAsingle, *piAsingle;
    int64 *pA64;
    int32 *pA32;
    int16 *pA16;
    int08 *pA08;
    int islast;

    /* Check number of arguments */
    if (nrhs!=2)
        mexErrMsgTxt("FIND1DMEX: two arguments are required.");
    
    /* Comment these 8 lines if you are sure about your MEX installation
     * and want to remove unecessary overhead checking */
    if (sizeof(int08) != 1)
        mexErrMsgTxt("FIND1DMEX: incorrect int08 definition (modify MEX file is required)");
    if (sizeof(int16) != 2)
        mexErrMsgTxt("FIND1DMEX: incorrect int16 definition (modify MEX file is required)");
    if (sizeof(int32) != 4)
        mexErrMsgTxt("FIND1DMEX: incorrect int32 definition (modify MEX file is required)");
    if (sizeof(int64) != 8)
        mexErrMsgTxt("FIND1DMEX: incorrect int64 definition (modify MEX file is required)");

    /* Get class of input matrix */
    ClassID = mxGetClassID(A);

    /* Get the size, MUST BE two or three, no check */
    dimA = mxGetDimensions(A);
    k = dimA[0];
    m = dimA[1];
    if (mxGetNumberOfDimensions(A)<3)
        n = 1; /* third dimension is singleton */
    else
        n = dimA[2];
    
    /* Get data pointers */
    PrA = mxGetPr(A);
    PiA = mxGetPi(A);
    km = k*m;
    
    count = (mwSignedIndex)(*(mxGetPr(COUNT)));
    if (count<0) /* by convention count<0 triggers find with 'LAST' option */
    {
        islast = 1;
        count = -count;
    }
    else islast = 0;
      
    /* Generate output array B, prefilled with zeros */
    dimB[0] = k; dimB[1] = count; dimB[2] = n;
    B = mxCreateNumericArray(3, dimB, mxDOUBLE_CLASS, mxREAL);
    if (B==NULL)
        mexErrMsgTxt("FIND1DMEX: out of memory.");
    PrB = mxGetPr(B);
    kc = k*count;
    
    /* Limit count to m */
    if (count>m) count = m;
    
    /* Nothing to do */
    if (count==0)
        return;
    
    /* Call the engine */
    switch (ClassID) {
        case mxDOUBLE_CLASS:
            if (PiA==NULL)
                {FIND_ENGINE(prAdouble, double, islast);}
            else
                {CMPLX_FIND_ENGINE(prAdouble, piAdouble, double, islast);}                
        case mxSINGLE_CLASS:
            if (PiA==NULL)
                {FIND_ENGINE(prAsingle, float, islast);}
            else
                {CMPLX_FIND_ENGINE(prAsingle, piAsingle, float, islast);}
        case mxINT64_CLASS:
        case mxUINT64_CLASS:
            FIND_ENGINE(pA64, int64, islast);
        case mxINT32_CLASS:
        case mxUINT32_CLASS:
            FIND_ENGINE(pA32, int32, islast);
        case mxCHAR_CLASS:  
        case mxINT16_CLASS:
        case mxUINT16_CLASS:
            FIND_ENGINE(pA16, int16, islast);
        case mxLOGICAL_CLASS:
        case mxINT8_CLASS:
        case mxUINT8_CLASS:
            FIND_ENGINE(pA08, int08, islast);
        default:
            mexErrMsgTxt("FIND1DMEX: Class not supported.");
    }            
    
    return;

} /* FIND1DMEX */
