// Exact Gray-code cap enumeration; compile outputs into the workspace tmp directory.
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

int main(int argc, char **argv) {
    int n = argc > 1 ? std::atoi(argv[1]) : 6;
    int k = argc > 2 ? std::atoi(argv[2]) : 4;
    if ((n != 2 && n != 5 && n != 6) || (k != 1 && k != 4)) return 2;
    const int seed5[5][5] = {{0,-1,1,-1,1},{-1,0,-1,1,1},
                            {1,-1,0,1,-1},{-1,1,1,0,-1},{1,1,-1,-1,0}};
    std::vector<int> A(n*n, 0);
    if (n == 2) {
        A = {1,1,1,-1};
    } else if (n == 5) {
        for (int i=0;i<n;++i) for(int j=0;j<n;++j) A[i*n+j]=seed5[i][j];
    } else {
        for (int i=0;i<n;++i) for(int j=0;j<n;++j) A[i*n+j]=(i!=j);
        int bit=0, code=220;
        for(int i=1;i<n;++i) for(int j=i+1;j<n;++j) {
            A[i*n+j]=A[j*n+i]=1-2*((code>>bit)&1); ++bit;
        }
    }
    const int N=n*k;
    std::vector<int> B(N*N), x(N,1), h(N,0), xmin, xmax;
    for(int a=0;a<k;++a) for(int b=0;b<k;++b)
        for(int i=0;i<n;++i) for(int j=0;j<n;++j)
            B[(a*n+i)*N+b*n+j]=(k==1 ? 1 : 1-2*(a==b))*A[i*n+j];
    int E=0;
    for(int i=0;i<N;++i) for(int j=0;j<N;++j) { h[i]+=B[i*N+j]; E+=B[i*N+j]; }
    E/=2;
    int lo=E, hi=E; xmin=x; xmax=x;
    const uint64_t end=uint64_t(1)<<(N-1);
    for(uint64_t t=1;t<end;++t) {
        int i=1+__builtin_ctzll(t), old=x[i];
        E+=-2*old*h[i]+2*B[i*N+i];
        for(int j=0;j<N;++j) h[j]-=2*old*B[j*N+i];
        x[i]=-old;
        if(E<lo) {lo=E; xmin=x;}
        if(E>hi) {hi=E; xmax=x;}
    }
    std::cout << "n="<<n<<" k="<<k<<" min="<<lo<<" max="<<hi<<" Q="<<std::max(-lo,hi)<<"\n";
    std::cout << "minimizer:"; for(int v:xmin) std::cout<<" "<<v; std::cout<<"\n";
    std::cout << "maximizer:"; for(int v:xmax) std::cout<<" "<<v; std::cout<<"\n";
}
