
// maxgap: real(C_FLOAT) function maxgap(naz, az) BIND(C)
extern "C" float maxgap(int naz, const float* az)
{
    if(naz <= 0) return 0.0f;
    std::vector<float> arr(az, az + naz);
    std::sort(arr.begin(), arr.end());
    float max_gap = 0.0f;
    for (int i = 1; i < naz; i++) {
         float gap = arr[i] - arr[i-1];
         if(gap > max_gap) max_gap = gap;
    }
    float circular_gap = 360.0f - arr.back() + arr.front();
    if(circular_gap > max_gap) max_gap = circular_gap;
    return max_gap;
}
