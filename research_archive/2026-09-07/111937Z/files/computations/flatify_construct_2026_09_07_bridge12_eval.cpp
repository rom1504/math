#include <cstdint>
#include <cstdlib>
extern "C" uint64_t evaluate(const int8_t* k,const int8_t* e,int size,int p,int q,int polarity,int cutoff){
  int best=-1; uint32_t count=0;
  for(int i=0;i<size;i++){
    const int a=e[i^p]; const int8_t* row=k+i*size;
    for(int j=0;j<size;j++){
      int value=std::abs(a+polarity*int(e[j^q]))+int(row[j]);
      if(value>cutoff)return uint64_t(value)<<32;
      if(value>best){best=value;count=1;}else if(value==best)count++;
    }
  }
  return (uint64_t(best)<<32)|count;
}
extern "C" uint64_t worst_pair(const int8_t* k,const int8_t* e,int size,int p,int q,int polarity){
  int best=-1;uint32_t arg=0;
  for(int i=0;i<size;i++)for(int j=0;j<size;j++){
    int value=std::abs(int(e[i^p])+polarity*int(e[j^q]))+int(k[i*size+j]);
    if(value>best){best=value;arg=i*size+j;}
  }
  return(uint64_t(best)<<32)|arg;
}
