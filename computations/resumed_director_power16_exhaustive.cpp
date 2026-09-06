// Exact exhaustive balanced-profile test of F_16 P_7 F_16.
// Compile: g++ -O3 -std=c++17 this_file.cpp -o /home/math/quadra/tmp/power16
// This is a finite lower test, NOT an upper bound for regularized signing norms.
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>
int mul(int a,int b){int c=0;for(;b;b>>=1){if(b&1)c^=a;a<<=1;if(a&16)a^=19;}return c;}
int pw(int a,int e){int c=1;for(;e;e>>=1,a=mul(a,a))if(e&1)c=mul(c,a);return c;}
int tr(int a){return a^pw(a,2)^pw(a,4)^pw(a,8);}
int main(){
 const auto start=std::chrono::steady_clock::now();
 int W[16][16],T[16][16];
 for(int i=0;i<16;i++)for(int j=0;j<16;j++){
  assert(tr(mul(i,j))==0||tr(mul(i,j))==1);
  W[i][j]=tr(mul(i,j))?-1:1;
 }
 for(int i=0;i<16;i++)for(int j=0;j<16;j++){
  T[i][j]=0;for(int a=0;a<16;a++)T[i][j]+=W[i][pw(a,7)]*W[a][j];
 }
 for(int i=0;i<16;i++)for(int j=0;j<16;j++){
  int c=0;for(int a=0;a<16;a++)c+=T[a][i]*T[a][j];
  assert(c==(i==j?256:0));
 }
 for(int i=0;i<16;i++){int s=0;for(int j=0;j<16;j++)s+=T[i][j];assert(s==16);}
 std::vector<unsigned> masks;
 std::vector<std::array<int,16>> values;
 for(unsigned mask=0;mask<65536;mask++)if(__builtin_popcount(mask)==8){
  masks.push_back(mask);std::array<int,16> a{};
  for(int i=0;i<16;i++)for(int j=0;j<16;j++)a[i]+=T[i][j]*((mask>>j&1)?1:-1);
  values.push_back(a);
 }
 int best=-1;unsigned bx=0,by=0;std::uint64_t pairs=0;
 // Simultaneous reversal of both inputs preserves the objective; fix first bit of first.
 for(size_t x=0;x<values.size();x++)if(masks[x]&1){
  for(size_t y=0;y<values.size();y++){
   int v=0;for(int j=0;j<16;j++)v+=std::abs(-values[x][j]+2*values[y][j])+std::abs(2*values[x][j]+4*values[y][j]);
   ++pairs;if(v>best){best=v;bx=masks[x];by=masks[y];}
  }
 }
 const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 std::cout<<"{\n  \"field\": \"GF(16), x^4+x+1\",\n  \"power\": 7,\n  \"balanced_profiles\": "<<masks.size()<<",\n  \"pairs_after_global_reversal\": "<<pairs<<",\n  \"exact_numerator\": "<<best<<",\n  \"exact_denominator\": 512,\n  \"value\": "<<double(best)/512<<",\n  \"witness_masks\": ["<<bx<<", "<<by<<"],\n  \"elapsed_seconds\": "<<seconds<<"\n}\n";
}
