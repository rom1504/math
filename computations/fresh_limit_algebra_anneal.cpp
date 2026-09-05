// Lower witnesses only, with fixed mt19937 seed 912734. No optimality claim.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

int main(int argc,char**argv) {
    const int n=argc>1?std::atoi(argv[1]):5, d=argc>2?std::atoi(argv[2]):2;
    const int runs=argc>3?std::atoi(argv[3]):200;
    if((n!=2 && n!=5 && n!=6)||d<0||d>4)return 2;
    int k=1;for(int i=0;i<d;++i)k*=4;
    const int N=n*k;
    const int seed5[5][5]={{0,-1,1,-1,1},{-1,0,-1,1,1},{1,-1,0,1,-1},{-1,1,1,0,-1},{1,1,-1,-1,0}};
    std::vector<int>A(n*n,0),B(N*N),x(N),h(N),bestx;
    if(n==2){A={1,1,1,-1};}
    else if(n==5){for(int i=0;i<n;++i)for(int j=0;j<n;++j)A[i*n+j]=seed5[i][j];}
    else{for(int i=0;i<n;++i)for(int j=0;j<n;++j)A[i*n+j]=(i!=j);int bit=0;
         for(int i=1;i<n;++i)for(int j=i+1;j<n;++j){A[i*n+j]=A[j*n+i]=1-2*((220>>bit)&1);++bit;}}
    for(int a=0;a<k;++a)for(int b=0;b<k;++b){int v=1,aa=a,bb=b;
       for(int t=0;t<d;++t){v*=1-2*(aa%4==bb%4);aa/=4;bb/=4;}
       for(int i=0;i<n;++i)for(int j=0;j<n;++j)B[(a*n+i)*N+b*n+j]=v*A[i*n+j];}
    std::mt19937 gen(912734);std::uniform_real_distribution<double>u(0,1);
    int best=-1;
    for(int run=0;run<runs;++run){
        for(int&i:x)i=2*(gen()%2)-1;
        int E=0;for(int i=0;i<N;++i){h[i]=0;for(int j=0;j<N;++j)h[i]+=B[i*N+j]*x[j];E+=x[i]*h[i];}E/=2;
        for(int sweep=0;sweep<150;++sweep){double temp=2*std::sqrt(N)*std::pow(.005,sweep/149.0);
            for(int p=0;p<N;++p){int i=gen()%N,old=x[i],delta=-2*old*h[i]+2*B[i*N+i];
                if(delta>=0||u(gen)<std::exp(delta/temp)){E+=delta;x[i]=-old;for(int j=0;j<N;++j)h[j]-=2*old*B[j*N+i];}}
        }
        bool changed=true;while(changed){changed=false;for(int i=0;i<N;++i)if(x[i]*h[i]-B[i*N+i]<0){int old=x[i];E+=-2*old*h[i]+2*B[i*N+i];x[i]=-old;for(int j=0;j<N;++j)h[j]-=2*old*B[j*N+i];changed=true;}}
        if(E>best){best=E;bestx=x;}
    }
    std::cout<<"n="<<n<<" k="<<k<<" runs="<<runs<<" lower_Q="<<best<<" normalized_seed="<<best/std::pow(k,1.5)<<"\n";
    std::cout<<"witness:";for(int v:bestx)std::cout<<" "<<v;std::cout<<"\n";
}
