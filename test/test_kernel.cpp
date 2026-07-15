// Standalone kernel unit test: verifies \int W d^Dx = 1 and \int r dW/dr d^Dx = -D for every
// (kernel, D). Build (from repo root):
//   x86_64-conda-linux-gnu-c++ -std=c++17 -O2 -fopenmp -I include \
//     -I /u/kpala/opt/cabana-0.8.0/include -I $CONDA/ctop10/include test/test_kernel.cpp src/kernel.cpp \
//     -o /tmp/test_kernel -L $CONDA/ctop10/lib -lkokkoscore -lkokkossimd && /tmp/test_kernel

#include "kernel.h"
#include <cstdio>
#include <cmath>
using namespace ccake;
template<unsigned int D>
void check(KernelType kt, const char* name){
  double h=1.0, s=SPHkernel<D>::support(kt); int N=800000; double dr=s/N;
  double IW=0, Igrad=0;
  for(int i=0;i<N;i++){ double r=(i+0.5)*dr;
    double W=SPHkernel<D>::kernel(r,h,kt);
    double rel[3]={r,0,0}, g[3]={0,0,0};
    SPHkernel<D>::gradKernel(rel,r,h,g,kt);   // g[0] = dW/dr
    double w=(D==1)?2.0:(D==2)?2*M_PI*r:4*M_PI*r*r;
    IW += W*w*dr;  Igrad += r*g[0]*w*dr; }
  printf("  %-13s D=%u  intW=%.6f (1)   int r dW/dr=%.6f (%d)\n",name,D,IW,Igrad,-(int)D);
}
int main(){
  const char* kn[]={"cubic","quintic","wendland_c2","wendland_c4","wendland_c6"};
  KernelType ks[]={KernelType::CUBIC,KernelType::QUINTIC,KernelType::WENDLAND_C2,KernelType::WENDLAND_C4,KernelType::WENDLAND_C6};
  for(int j=0;j<5;j++){ check<1>(ks[j],kn[j]); check<2>(ks[j],kn[j]); check<3>(ks[j],kn[j]); }
  return 0;
}
