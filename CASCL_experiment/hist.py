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

# NE
for phi in range(0, CODE_LENGTH):
    if info[phi] == False:
        continue
    EE = []
    EN = []
    with open("./calculatedLLR/preNonError/thisError/SN{0:.2f}/phi{1:04d}.txt".format(EbN0dB, phi)) as f:
        for s_line in f:
            EE.append(math.fabs(float(s_line)))
    if len(EE) == 0:
        continue

    with open("./calculatedLLR/preNonError/thisNonError/SN{0:.2f}/phi{1:04d}.txt".format(EbN0dB, phi)) as f:
        for s_line in f:
            EN.append(math.fabs(float(s_line)))
    EE = np.array(EE)
    EN = np.array(EN)
    fig = plt.figure()
    labels = ["preError_thisError", "preNonError_thisError"]
    plt.hist([EE,EN], label=labels, log=True)
    plt.title("Information Bit EE vs EN:{0:04d}".format(phi))
    plt.xlabel("LLR")
    plt.ylabel("count")
    plt.legend()
    plt.grid(True)
    fig.savefig(
        "./images/InformationEEvsEN/SN{0:.2f}phi{1:04d}.png".format(EbN0dB, phi))
