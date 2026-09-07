"""Fast alternating witnesses for tensor averaging-reflection tests (not upper bounds)."""
import numpy as np

rng=np.random.default_rng(20260906)

def transform(X,p,t):
    Y=X.copy()
    for axis in range(t):
        block=Y.reshape(len(Y),2**axis,2,2**(t-axis-1),2)
        mean=p*block[:,:,0,:,:]+(1-p)*block[:,:,1,:,:]
        block[:,:,:,:,:]=2*mean[:,:,None,:,:]-block
    return Y

for t in [3,4,5,6,8]:
    for p in [.1,.2,.3,.363726,.4,.45]:
        mu=np.array([1.])
        for _ in range(t):mu=np.kron(mu,[p,1-p])
        for r in [1.,2.]:
            B=np.array([[-1,r],[r,r*r]])
            X=rng.choice([-1.,1.],size=(96,2**t,2))
            X[:4,:,:]=np.array([[-1,-1],[-1,1],[1,-1],[1,1]])[:,None,:]
            for _ in range(35):
                Y=np.where(transform(X,p,t)@B>=0,1.,-1.)
                X=np.where(transform(Y,p,t)@B>=0,1.,-1.)
            field=transform(X,p,t)@B
            vals=(np.abs(field).sum(axis=2)*mu).sum(axis=1)/2
            best=vals.max()
            if (r==2 and best>3.50000001) or (r==1 and best>1.2873):
                print('improvement',t,p,r,best,flush=True)
        print('done',t,p,flush=True)
