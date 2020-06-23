import matplotlib.pyplot as plt
import numpy as np
import math
import statistics

EbN0dB = 2.00
CODE_LENGTH = 1024
axis = [i for i in range(CODE_LENGTH)]
EE = [[] for i in range(CODE_LENGTH)]
NE = [[] for i in range(CODE_LENGTH)]
EEave = []
NEave = []

# 情報ビット
info = [False for i in range(CODE_LENGTH)]
isFrozen = []
with open("./isFrozen/SN{0:.2f}Len{1:04d}.txt".format(EbN0dB, CODE_LENGTH)) as f:
    for s_line in f:
        isFrozen.append(int(s_line))
for i in range(CODE_LENGTH):
    if isFrozen[i] == 0:
        info[i] = True

# EE
for phi in range(0, CODE_LENGTH):
    with open("./calculatedLLR/preError/thisError/SN{0:.2f}/phi{1:04d}.txt".format(EbN0dB, phi)) as f:
        for s_line in f:
            EE[phi].append(math.fabs(float(s_line)))
    if len(EE[phi]) == 0:
        EEave.append(None)
    else:
        EEave.append(statistics.mean(EE[phi]))

# NE
for phi in range(0, CODE_LENGTH):
    with open("./calculatedLLR/preNonError/thisError/SN{0:.2f}/phi{1:04d}.txt".format(EbN0dB, phi)) as f:
        for s_line in f:
            NE[phi].append(math.fabs(float(s_line)))
    if len(NE[phi]) == 0:
        NEave.append(None)
    else:
        NEave.append(statistics.mean(NE[phi]))

info = np.array(info)
axis = np.array(axis)
EEave = np.array(EEave)
NEave = np.array(NEave)
fig = plt.figure()

option = np.logical_not(info)

plt.plot(axis[option], EEave[option], marker="o", color="blue",
         linestyle="", label="preError_thisError")
plt.plot(axis[option], NEave[option], marker="o", color="purple",
         linestyle="", label="preNonError_thisError")

plt.title("Frozen Bit EE vs NE")
plt.xlabel("Channel Index")
plt.ylabel("average absLLR in SCDecoder")
plt.legend()
plt.grid(True)

fig.savefig("img.png")
