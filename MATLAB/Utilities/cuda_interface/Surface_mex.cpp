#include "mex.h"
#include "../../../Common/CUDA/Siddon_surface_projection.hpp"
#include <cmath>
#include <limits>

// Internal packed interface; use Ax_surface / Atb_surface for public layouts.
static mwSize product(mwSize a, mwSize b) {
    if (b && a > std::numeric_limits<mwSize>::max() / b)
        mexErrMsgIdAndTxt("TIGRE:surface:size", "Array dimensions overflow addressable size.");
    return a * b;
}

void mexFunction(int nlhs, mxArray** out, int nrhs, const mxArray** in) {
    if (nrhs != 9 || nlhs != 1)
        mexErrMsgIdAndTxt("TIGRE:surface:arguments", "Expected nine inputs and one output.");
    for (int k = 0; k < 9; ++k) {
        bool validType = k == 7 ? mxIsUint8(in[k]) : (k == 2 ? mxIsDouble(in[k]) : mxIsSingle(in[k]));
        if (!validType || mxIsComplex(in[k]) || mxIsSparse(in[k]))
            mexErrMsgIdAndTxt("TIGRE:surface:type", "Expected full real single arrays, dimensions double, mask uint8.");
    }
    if (mxGetNumberOfElements(in[0]) != 1 || mxGetNumberOfElements(in[8]) != 1)
        mexErrMsgIdAndTxt("TIGRE:surface:scalar", "Mode and device must be scalar.");
    for (int k = 2; k <= 4; ++k)
        if (mxGetNumberOfElements(in[k]) != 3)
            mexErrMsgIdAndTxt("TIGRE:surface:geometry", "Volume vectors need three entries.");
    auto ptr = [&](int k) { return static_cast<const float*>(mxGetData(in[k])); };
    const double* n = mxGetPr(in[2]);
    const float *d = ptr(3), *o = ptr(4);
    for (int k = 0; k < 3; ++k)
        if (!std::isfinite(n[k]) || n[k] < 1 || n[k] > std::numeric_limits<int>::max() || n[k] != std::floor(n[k]) ||
            !std::isfinite(d[k]) || d[k] <= 0 || !std::isfinite(o[k]))
            mexErrMsgIdAndTxt("TIGRE:surface:geometry", "Invalid volume geometry.");
    float mode = ptr(0)[0], device = ptr(8)[0];
    if ((mode != 0 && mode != 1) || !std::isfinite(device) || device < 0 ||
            double(device) > std::numeric_limits<int>::max() || device != std::floor(device))
        mexErrMsgIdAndTxt("TIGRE:surface:mode", "Invalid mode or zero-based device.");
    bool back = mode == 1;
    const mwSize* sz = mxGetDimensions(in[6]);
    mwSize nd = mxGetNumberOfDimensions(in[6]);
    if (nd > 4 || sz[0] != 3)
        mexErrMsgIdAndTxt("TIGRE:surface:pixels", "Packed pixels must be 3 x columns x rows x views.");
    mwSize cols = sz[1], rows = nd > 2 ? sz[2] : 1, pv = nd > 3 ? sz[3] : 1;
    mwSize views = mxGetN(in[5]);
    if (mxGetNumberOfDimensions(in[5]) != 2 || mxGetM(in[5]) != 3 || !views || !rows || !cols ||
            (pv != 1 && pv != views) || rows > std::numeric_limits<int>::max() ||
            cols > std::numeric_limits<int>::max() || views > std::numeric_limits<int>::max())
        mexErrMsgIdAndTxt("TIGRE:surface:views", "Invalid source/pixel view dimensions.");
    mwSize plane = product(rows, cols), count = mxGetNumberOfElements(in[7]);
    mwSize mv = count / plane;
    if (count != product(mv, plane) || (mv != 1 && mv != views))
        mexErrMsgIdAndTxt("TIGRE:surface:mask", "Invalid packed mask size.");
    const unsigned char* mask = static_cast<const unsigned char*>(mxGetData(in[7]));
    for (mwSize j = 0; j < count; ++j)
        if (mask[j] > 1) mexErrMsgIdAndTxt("TIGRE:surface:mask", "Mask must be binary.");
    mwSize vox = product(product(mwSize(n[0]), mwSize(n[1])), mwSize(n[2]));
    mwSize samples = product(plane, views);
    product(back ? vox : samples, sizeof(float));
    if (mxGetNumberOfElements(in[1]) != (back ? samples : vox))
        mexErrMsgIdAndTxt("TIGRE:surface:data", "Input data size mismatch.");
    const int finiteInputs[] = {1, 5, 6};
    for (int k : finiteInputs)
        for (mwSize j = 0; j < mxGetNumberOfElements(in[k]); ++j)
            if (!std::isfinite(ptr(k)[j]))
                mexErrMsgIdAndTxt("TIGRE:surface:finite", "Packed data and coordinates must be finite.");
    mwSize dims[3] = {back ? mwSize(n[0]) : cols, back ? mwSize(n[1]) : rows, back ? mwSize(n[2]) : views};
    out[0] = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);
    auto fn = back ? siddon_surface_backprojection : siddon_surface_projection;
    int status = fn(ptr(1), int(n[0]), int(n[1]), int(n[2]), d[0], d[1], d[2], o[0], o[1], o[2],
        int(views), int(rows), int(cols), ptr(5), ptr(6), int(pv), mask, int(mv),
        static_cast<float*>(mxGetData(out[0])), int(device));
    if (status)
        mexErrMsgIdAndTxt("TIGRE:surface:CUDA", "Surface CUDA operation failed; see CUDA diagnostic above.");
}
