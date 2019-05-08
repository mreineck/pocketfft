#include <complex>
#include <vector>
#include "pocketfft_hdronly.h"

using namespace std;

using pocketfft_private::shape_t;
using pocketfft_private::stride_t;

int main()
  {
  shape_t shape{128,128,320};
  stride_t stride(shape.size()), rstride(shape.size());
  size_t tmp=sizeof(complex<double>);
  size_t rtmp=sizeof(double);
  for (int i=shape.size()-1; i>=0; --i)
    {
    stride[i]=tmp;
    tmp*=shape[i];
    rstride[i]=rtmp;
    rtmp*=shape[i];
    }
  size_t ndata=1;
  for (size_t i=0; i<shape.size(); ++i)
    ndata*=shape[i];
  vector<complex<double>> data(ndata), data2(ndata);
  vector<double> rdata(ndata), rdata2(ndata);
  shape_t axes;
  for (size_t i=0; i<shape.size(); ++i)
    axes.push_back(i);
  for (size_t rep=0; rep<10; ++rep)
    pocketfft_c2c<double>(shape, stride, stride, axes, true,
      data.data(), data.data(), 1.);
  for (size_t rep=0; rep<10; ++rep)
    pocketfft_r2r_fftpack<double>(shape, stride, stride, 0, true,
      data.data(), data.data(), 1.);
  for (size_t rep=0; rep<10; ++rep)
    pocketfft_r2r_hartley<double>(shape, stride, stride, axes,
      data.data(), data.data(), 1.);
  }
