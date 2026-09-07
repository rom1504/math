#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

// Exact projective Gray scan; command arguments: order and upper triangle bits.
// Bit '1' means +1, '0' means -1, pairs in lexicographic i<j order.
int main(int argc,char**argv){
  if(argc!=3)return 2;
  const int n=std::atoi(argv[1]); const std::string bits=argv[2];
  if(n<2||n>30||bits.size()!=size_t(n*(n-1)/2))return 3;
  std::vector<std::vector<int>> a(n,std::vector<int>(n));
  int t=0,q=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
    a[i][j]=a[j][i]=bits[t++]=='1'?1:-1;q+=a[i][j];
  }
  std::vector<int>x(n,1),h(n);for(int i=0;i<n;i++)for(int j=0;j<n;j++)h[i]+=a[i][j];
  int lo=q,hi=q;uint64_t ilo=0,ihi=0;
  const uint64_t end=uint64_t(1)<<(n-1);
  for(uint64_t k=1;k<end;k++){
    const int i=__builtin_ctzll(k)+1,s=x[i];q-=2*s*h[i];
    for(int j=0;j<n;j++)h[j]-=2*s*a[j][i];x[i]=-s;
    if(q>hi){hi=q;ihi=k^(k>>1);}if(q<lo){lo=q;ilo=k^(k>>1);}
  }
  std::cout<<"{\"n\":"<<n<<",\"min\":"<<lo<<",\"max\":"<<hi
    <<",\"cap\":"<<std::max(hi,-lo)<<",\"min_gray_spin_code\":"<<ilo
    <<",\"max_gray_spin_code\":"<<ihi<<",\"projective_states\":"<<end<<"}\n";
}
