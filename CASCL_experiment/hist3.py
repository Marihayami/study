import matplotlib.pyplot as plt
import numpy as np
import math
import statistics

EbN0dB = 2.00
CODE_LENGTH = 1024


# 情報ビット
info = [False for i in range(CODE_LENGTH)]
isFrozen = []
with open("./isFrozen/SN{0:.2f}Len{1:04d}.txt".format(EbN0dB, CODE_LENGTH)) as f:
    for s_line in f:
        isFrozen.append(int(s_line))
for i in range(CODE_LENGTH):
    if isFrozen[i] == 0:
        info[i] = True

for phi in range(0, CODE_LENGTH):
    if info[phi] == True:
        continue
    EN = []
    NN = []
    with open("./calculatedLLR/preError/thisNonError/SN{0:.2f}/phi{1:04d}.txt".format(EbN0dB, phi)) as f:
        for s_line in f:
            EN.append(math.fabs(float(s_line)))
        if len(EN) == 0:
            continue

    with open("./calculatedLLR/preNonError/thisNonError/SN{0:.2f}/phi{1:04d}.txt".format(EbN0dB, phi)) as f:
        for s_line in f:
            NN.append(math.fabs(float(s_line)))

    EN = np.array(EN)
    NN = np.array(NN)
    """
    fig = plt.figure()
    labels = ["preError_thisNonError", "preNonError_thisNonError"]
    plt.hist([EN, NN], label=labels, log=True)
    plt.title("Frozen Bit EN vs NN:{0:04d}".format(phi))
    plt.xlabel("LLR")
    plt.ylabel("count")
    plt.legend()
    plt.grid(True)
    fig.savefig(
        "./images/FrozenENvsNN/SN{0:.2f}phi{1:04d}.png".format(EbN0dB, phi))
    """
