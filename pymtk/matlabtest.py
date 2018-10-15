import numpy as np
from matrizestencil import banded

k = 4
s = np.arange(1,5)

for j in range(k//2, k-1):
    print("j", j, "- k//2", k//2, "+ 1")
    a = np.zeros(int(j - k//2 + 1))
    print("k", k, "- j", j, "-1")
    b = np.zeros(k - j - 1)
    a = np.concatenate((a, s, b))
    print(a)

