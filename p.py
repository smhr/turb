import numpy as np
n = 3
#B = np.random.rand(n, n, n)
B = np.arange(1,28).astype('float64')
A = np.reshape(B, (n, n, n), order="F") # smhr: It seems order is not important, what is important is to how to write in file.
print(A[:,:,:])
print(A.sum())
A.T.tofile('data.bin')
#np.save('data.npy', A)
