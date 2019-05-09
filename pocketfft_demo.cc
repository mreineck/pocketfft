#include <complex>
#include <vector>
#include "pocketfft_hdronly.h"

using namespace std;
using namespace pocketfft;

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
    c2c<double>(shape, stride, stride, axes, true,
      data.data(), data.data(), 1.);
  for (size_t rep=0; rep<10; ++rep)
    r2r_fftpack<double>(shape, rstride, rstride, 0, true,
      rdata.data(), rdata.data(), 1.);
  for (size_t rep=0; rep<10; ++rep)
    r2r_hartley<double>(shape, rstride, rstride, axes,
      rdata.data(), rdata.data(), 1.);
  }
