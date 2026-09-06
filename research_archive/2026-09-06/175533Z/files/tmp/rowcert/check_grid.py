exec(open('/home/math/quadra/tmp/rowcert/check_c.py').read().split('for N in [6')[0])
for N in range(4,10):
 for beta in [.25,.5,1,2,4]:
  m,n,L,T=pressure(N,beta); C=cmat(T);rho=max(np.linalg.eigvalsh(C)) if m>1 else 0
  passes=[la for la in [.25,.5,1,2,4] if la*rho<4]
  print(N,beta,m,n,rho,passes)
