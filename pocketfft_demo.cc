#include <complex>
#include <vector>
#include "pocketfft.h"

using namespace std;

int main()
  {
  vector<size_t> shape{2048,2048};
  vector<ptrdiff_t> stride(shape.size());
  size_t tmp=sizeof(complex<double>);
  for (int i=shape.size()-1; i>=0; --i)
    {
    stride[i]=tmp;
    tmp*=shape[i];
    }
  size_t ndata=1;
  for (size_t i=0; i<shape.size(); ++i)
    ndata*=shape[i];
  vector<complex<double>> data(ndata);
  vector<size_t> axes;
  for (size_t i=0; i<shape.size(); ++i)
    axes.push_back(i);
  for (size_t rep=0; rep<10; ++rep)
    pocketfft_complex(shape.size(), shape.data(), stride.data(), stride.data(),
      axes.size(), axes.data(), true, data.data(), data.data(), 1., true);

  }
