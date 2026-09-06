"""Standalone integer certificate for the resumed campaign; no solver needed."""
import numpy as np

sq = {i*i % 11 for i in range(1,11)}
F = np.ones((12,12),dtype=np.int64)
F[1:,0] = -1
for i in range(11):
    for j in range(11):
        F[i+1,j+1] = 1 if i == j or (j-i)%11 in sq else -1
colors = np.array([
 [0,1,2,1,1,2,2,2,1,0,0,0],
 [0,2,0,2,0,1,1,2,2,1,0,1],
 [2,2,0,1,0,1,0,1,2,2,1,0],
 [2,2,1,0,0,1,0,0,2,2,1,1],
 [2,0,0,1,1,2,2,0,0,1,1,2],
 [0,1,2,1,0,2,2,2,0,1,0,1],
 [2,0,0,0,2,2,1,1,1,2,1,0],
 [1,0,1,2,2,1,0,0,1,0,2,2],
 [1,0,2,0,2,0,1,1,0,2,2,1],
 [1,1,1,2,2,0,1,0,0,0,2,2],
 [1,1,2,0,1,0,0,2,2,1,2,0],
 [0,2,1,2,1,0,2,1,1,0,0,2]],dtype=np.int64)
assert np.array_equal(F.T@F,12*np.eye(12,dtype=np.int64))
for c in range(3):
    Z=F*(colors==c)
    assert np.array_equal(Z.T@F+F.T@Z,8*np.eye(12,dtype=np.int64))
    assert np.count_nonzero(colors==c)==48
H=np.empty((144,144),dtype=np.int64)
for i in range(12):
    for j in range(12):
        for l in range(12):
            for k in range(12):
                H[12*i+j,12*l+k]=F[i,j]*F[i,k]*F[l,k]*F[l,j]
assert np.array_equal(H,H.T)
assert np.array_equal(H@H,144*np.eye(144,dtype=np.int64))
assert np.array_equal(H@np.ones(144,dtype=np.int64),12*np.ones(144,dtype=np.int64))
for c in range(3):
    P=(colors.ravel()==c).astype(np.int64)
    assert np.array_equal(H@P,8*np.ones(144,dtype=np.int64)-12*P)
B=np.ones((3,3),dtype=np.int64)-2*np.eye(3,dtype=np.int64)
X=B[colors.ravel(),:]
energy=int(np.sum(X*(H@X@B)))
assert energy == 9792
print('PASS: quotient 2E3-I3; quadratic energy=9792; q=4896; R(B)=17/6')
B2=np.ones((3,3),dtype=np.int64)
B2[0,0]=B2[1,1]=-1
F0=np.array([[1,1,1],[1,-1,-1],[-1,1,-1]],dtype=np.int64)
left=F0[colors.ravel(),:]
right=-F0[np.array([0,2,1])[colors.ravel()],:]
cross=int(np.sum(right*(H@left@B2)))
assert cross == 9792
# R4=(J4-2I4), z=(left,right,right,left); verify the 8*cross identity
# using the four blocks, without constructing the 1728-by-1728 matrix.
blocks=[left,right,right,left]
lift_energy=sum((1 if a!=b else -1)*int(np.sum(blocks[a]*(H@blocks[b]@B2)))
                for a in range(4) for b in range(4))
assert lift_energy == 8*cross == 78336
print('PASS: second seed bilinear=9792; R4-lift quadratic=78336; R(B2)=17/6')
