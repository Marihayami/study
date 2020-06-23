import matplotlib.pyplot as plt
import numpy as np
import math
import statistics

EbN0dB = 2.00
CODE_LENGTH = 1024
axis = [i for i in range(CODE_LENGTH)]
NE = [[] for i in range(CODE_LENGTH)]
NN = [[] for i in range(CODE_LENGTH)]
NEave = []
NNave = []
cond = [False for i in range(CODE_LENGTH)]

# 情報ビット
info = [False for i in range(CODE_LENGTH)]
isFrozen = []
with open("./isFrozen/SN{0:.2f}Len{1:04d}.txt".format(EbN0dB, CODE_LENGTH)) as f:
    for s_line in f:
        isFrozen.append(int(s_line))
for i in range(CODE_LENGTH):
    if isFrozen[i] == 0:
        info[i] = True

# NE
for phi in range(0, CODE_LENGTH):
    with open("./calculatedLLR/preNonError/thisError/SN{0:.2f}/phi{1:04d}.txt".format(EbN0dB, phi)) as f:
        for s_line in f:
            NE[phi].append(math.fabs(float(s_line)))
    if len(NE[phi]) == 0:
        NEave.append(None)
    else:
        NEave.append(statistics.mean(NE[phi]))
        cond[phi] = True

# NN
for phi in range(0, CODE_LENGTH):
    with open("./calculatedLLR/preNonError/thisNonError/SN{0:.2f}/phi{1:04d}.txt".format(EbN0dB, phi)) as f:
        for s_line in f:
            NN[phi].append(math.fabs(float(s_line)))
    if len(NN[phi]) == 0:
        NNave.append(None)
    else:
        NNave.append(statistics.mean(NN[phi]))

info = np.array(info)
cond = np.array(cond)
axis = np.array(axis)
NEave = np.array(NEave)
NNave = np.array(NNave)
fig = plt.figure()

option = np.logical_and(info, cond)

plt.plot(axis[option], NEave[option], marker="o", color="blue",
         linestyle="", label="preNonError_thisError")
plt.plot(axis[option], NNave[option], marker="o", color="purple",
         linestyle="", label="preNonError_thisNonError")

plt.title("Information Bit NE vs NN")
plt.xlabel("Channel Index")
plt.ylabel("average absLLR in SCDecoder")
plt.legend()
plt.grid(True)

fig.savefig("img.png")
