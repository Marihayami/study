import matplotlib.pyplot as plt
import math
import statistics

EbN0dB = 2.00
CODE_LENGTH = 1024
axis = [i for i in range(CODE_LENGTH)]
EE = [[] for i in range(CODE_LENGTH)]
EN = [[] for i in range(CODE_LENGTH)]
NE = [[] for i in range(CODE_LENGTH)]
NN = [[] for i in range(CODE_LENGTH)]
EEave = [None for i in range(CODE_LENGTH)]
ENave = [None for i in range(CODE_LENGTH)]
NEave = [None for i in range(CODE_LENGTH)]
NNave = [None for i in range(CODE_LENGTH)]

for phi in range(0, CODE_LENGTH):
    with open("./calculatedLLR/preError/thisError/SN{0:.2f}/phi{1:04d}.txt".format(EbN0dB, phi)) as f:
        for s_line in f:
            EE[phi].append(math.fabs(float(s_line)))
    if len(EE[phi]) == 0:
        EEave[phi]=None
    else:
        EEave[phi] = statistics.mean(EE[phi])


for phi in range(0, CODE_LENGTH):
    with open("./calculatedLLR/preError/thisNonError/SN{0:.2f}/phi{1:04d}.txt".format(EbN0dB, phi)) as f:
        for s_line in f:
            EN[phi].append(math.fabs(float(s_line)))
    if len(EN[phi]) == 0:
        ENave[phi]=None
    else:
        ENave[phi] = statistics.mean(EN[phi])


for phi in range(0, CODE_LENGTH):
    with open("./calculatedLLR/preNonError/thisError/SN{0:.2f}/phi{1:04d}.txt".format(EbN0dB, phi)) as f:
        for s_line in f:
            NE[phi].append(math.fabs(float(s_line)))
    if len(NE[phi]) == 0:
            NEave[phi]=None
    else:
        NEave[phi] = statistics.mean(NE[phi])
        
for phi in range(0, CODE_LENGTH):
    with open("./calculatedLLR/preNonError/thisNonError/SN{0:.2f}/phi{1:04d}.txt".format(EbN0dB, phi)) as f:
        for s_line in f:
            NN[phi].append(math.fabs(float(s_line)))
    if len(NN[phi]) == 0:
            NNave[phi]=None
    else:
        NNave[phi] = statistics.mean(NN[phi])

fig = plt.figure()

plt.plot(axis, EEave, marker="o", color="red", linestyle="",label="preError_thisError")
plt.plot(axis, ENave, marker="o", color="green", linestyle="",label="preError_thisNonError")
plt.plot(axis, NEave, marker="o", color="blue", linestyle="",label="preNonError_thisError")
plt.plot(axis, NNave, marker="o", color="purple", linestyle="",label="preNonError_thisNonError")
plt.title("All Bit EE vs EN vs NE vs NN")
plt.xlabel("Channel Index")
plt.ylabel("average absLLR in SCDecoder")
plt.legend()
plt.grid(True)

fig.savefig("img.png")
