import numpy as np
import math
import matplotlib.pyplot as plt

EbN0dB = 2.00
CODE_LENGTH = 1024
Res = [[] for i in range(CODE_LENGTH)]
for phi in range(128,CODE_LENGTH):
    with open("./calculatedLLR/error/SN{0:.2f}/phi{1:04d}.txt".format(EbN0dB, phi)) as f:
        for s_line in f:
            Res[phi].append(math.fabs(float(s_line)))

fig, ax = plt.subplots()
bp = ax.boxplot(Res)
plt.title('Error: EbN0 [dB] =2.00')
plt.xlabel('Channel index')
plt.ylabel('calculatedLLR')
# Y軸のメモリのrange
plt.ylim([0, 1000])
plt.grid()
# 描画
plt.show()
