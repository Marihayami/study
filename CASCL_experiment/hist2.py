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
    NE = []
    NN = []
    with open("./calculatedLLR/preNonError/thisError/SN{0:.2f}/phi{1:04d}.txt".format(EbN0dB, phi)) as f:
        for s_line in f:
            NE.append(math.fabs(float(s_line)))
    if len(NE) == 0:
        continue

    with open("./calculatedLLR/preNonError/thisNonError/SN{0:.2f}/phi{1:04d}.txt".format(EbN0dB, phi)) as f:
        for s_line in f:
            NN.append(math.fabs(float(s_line)))

    with open("./calculatedLLR/preError/thisError/SN{0:.2f}/phi{1:04d}.txt".format(EbN0dB, phi)) as f:
        for s_line in f:
            EE.append(math.fabs(float(s_line)))

    with open("./calculatedLLR/preError/thisNonError/SN{0:.2f}/phi{1:04d}.txt".format(EbN0dB, phi)) as f:
        for s_line in f:
            EN.append(math.fabs(float(s_line)))

    EE = np.array(EE)
    EN = np.array(EN)
    NE = np.array(NE)
    NN = np.array(NN)
    fig = plt.figure()
    labels = ["preError_thisError", "preError_thisNonError",
              "preNonError_thisError", "preNonError_thisNonError"]
    plt.hist([EE, EN, NE, NN], label=labels, log=True)
    plt.title("Information Bit EE vs EN vs NE vs NN:{0:04d}".format(phi))
    plt.xlabel("LLR")
    plt.ylabel("count")
    plt.legend()
    plt.grid(True)
    fig.savefig(
        "./images/InformationEEvsENvsNEvsNN/SN{0:.2f}phi{1:04d}.png".format(EbN0dB, phi))
