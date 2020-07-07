#include "bits/stdc++.h"
using namespace std;
pair<double, double> boxMuller(double);
pair<double, double> complexSum(pair<double, double>, pair<double, double>);
pair<double, double> complexDifference(pair<double, double>, pair<double, double>);
pair<double, double> complexProduct(pair<double, double>, pair<double, double>);
pair<double, double> complexQuotient(pair<double, double>, pair<double, double>);
pair<double, double> complexConjugate(pair<double, double>);
pair<double, double> scalarProduct(pair<double, double>, double);
pair<double, double> scalarQuotient(pair<double, double>, double);
void awgn(int, double, vector<pair<double, double>>, vector<pair<double, double>> &);
void bpskModulation(int, double, vector<bool>, vector<pair<double, double>> &);
void bpskDemodulation(int, double, vector<pair<double, double>>, vector<double> &);
void crcEncode(vector<bool> &, vector<int>, vector<int>, vector<int>);
void errorCount(int, double &, double &, vector<bool>, vector<bool>);
void generateBitInverse(int, vector<int> &);
void generateChanellLLR(int, double, vector<double> &);
void generateFrozenArray(int, int, vector<bool> &, vector<double>);
void generateInformationArray(vector<bool> &, vector<bool> &, vector<int>, vector<int>);
void generateInformationCount(int, vector<bool>, vector<int> &);
void generateInputArray(int, vector<bool>, vector<bool> &, vector<bool>, vector<int>);
void initializeCRC(vector<int>, vector<int> &);
void polarEncode(int, vector<bool>, vector<bool> &, vector<int>);
bool boolSum(bool, bool);
bool checkCRC(int, int, int, vector<bool>);
bool errorCount(int, vector<bool>, vector<bool>);
int randBin();
double bpskLLR(double, pair<double, double>);
double complexNorm(pair<double, double>);
class CASCL_Decoder
{
    int SEGMENT_LENGTH;
    int crcNum;
    int infoNum;

    int STAGE_LENGTH;
    int CODE_LENGTH;
    int CRC_LENGTH;
    int INFO_LENGTH;
    int LIST_LENGTH;

    vector<bool> CRCInformation;
    vector<bool> isFrozen;
    vector<int> BitInverse;

    vector<int> CRC_LENGTHS;
    vector<int> crcPolynomial;
    vector<int> crcPosition;
    vector<int> INFO_LENGTHS;

    vector<double> Boundary;

    vector<vector<vector<bool>>> Codeword;

    stack<int> inactivePathIndices;
    vector<bool> activePath;
    vector<vector<vector<vector<double>>>> arrayPointer_P;
    vector<vector<vector<vector<bool>>>> arrayPointer_C;
    vector<stack<int>> inactiveArrayIndices_P;
    vector<vector<int>> pathIndexToArrayIndex_P;
    vector<vector<int>> arrayReferenceCount_P;
    vector<stack<int>> inactiveArrayIndices_C;
    vector<vector<int>> pathIndexToArrayIndex_C;
    vector<vector<int>> arrayReferenceCount_C;

public:
    CASCL_Decoder(int CODE_LENGTH_COPY, int CRC_LENGTH_COPY, int INFO_LENGTH_COPY, int LIST_LENGTH_COPY, int SEGMENT_LENGTH_COPY, vector<bool> CRCInformation_COPY, vector<bool> &EstimatedCRCInformation, vector<bool> isFrozen_COPY, vector<int> BitInverse_COPY, vector<int> CRC_LENGTHS_COPY, vector<int> crcPolynomial_COPY, vector<int> crcPosition_COPY, vector<int> INFO_LENGTHS_COPY, vector<double> Boundary_COPY, vector<double> ReceivedLLR)
    {
        crcNum = 0;
        infoNum = 0;

        /*Initialization*/
        InitializeDataStructures(CODE_LENGTH_COPY, CRC_LENGTH_COPY, INFO_LENGTH_COPY, LIST_LENGTH_COPY, SEGMENT_LENGTH_COPY, CRCInformation_COPY, isFrozen_COPY, BitInverse_COPY, CRC_LENGTHS_COPY, crcPolynomial_COPY, crcPosition_COPY, INFO_LENGTHS_COPY, Boundary_COPY);
        int l = assignInitialPath();
        vector<vector<double>> &P_0 = getArrayPointer_P(0, l);
        for (int beta = 0; beta < CODE_LENGTH; beta++)
        {
            P_0[beta][1] = 1.0 / (1.0 + exp(ReceivedLLR[beta]));
            P_0[beta][0] = 1 - P_0[beta][1];
        }
        /*Initialization-end*/

        /*Main loop*/
        for (int phi = 0; phi < CODE_LENGTH; phi++)
        {
            bool phiMod = phi % 2;
            recursivelyCalcP(STAGE_LENGTH - 1, phi);
            if (isFrozen[phi])
            {
                for (int l = 0; l < LIST_LENGTH; l++)
                {
                    if (activePath[l] == false)
                        continue;
                    vector<vector<bool>> &C_m = getArrayPointer_C(STAGE_LENGTH - 1, l);
                    C_m[0][phiMod] = 0;
                }
            }
            else
            {
                continuePaths_UnfrozenBit(phi);
                infoNum++;
            }

            if (phiMod == 1)
                recursivelyUpdateC(STAGE_LENGTH - 1, phi);

            if (phi == crcPosition[crcNum])
            {
                getCodeword(EstimatedCRCInformation);
                crcNum++;
                infoNum = 0;
            }
        }
        /*Main loop-end*/
    }
    void InitializeDataStructures(int CODE_LENGTH_COPY, int CRC_LENGTH_COPY, int INFO_LENGTH_COPY, int LIST_LENGTH_COPY, int SEGMENT_LENGTH_COPY, vector<bool> CRCInformation_COPY, vector<bool> isFrozen_COPY, vector<int> BitInverse_COPY, vector<int> CRC_LENGTHS_COPY, vector<int> crcPolynomial_COPY, vector<int> crcPosition_COPY, vector<int> INFO_LENGTHS_COPY, vector<double> Boundary_COPY)
    {
        SEGMENT_LENGTH = SEGMENT_LENGTH_COPY;

        CODE_LENGTH = CODE_LENGTH_COPY;
        CRC_LENGTH = CRC_LENGTH_COPY;
        INFO_LENGTH = INFO_LENGTH_COPY;
        LIST_LENGTH = LIST_LENGTH_COPY;

        CRC_LENGTHS.resize(SEGMENT_LENGTH);
        crcPolynomial.resize(SEGMENT_LENGTH);
        crcPosition.resize(SEGMENT_LENGTH);
        INFO_LENGTHS.resize(SEGMENT_LENGTH);

        Codeword.resize(SEGMENT_LENGTH);
        for (int i = 0; i < SEGMENT_LENGTH; i++)
        {
            CRC_LENGTHS[i] = CRC_LENGTHS_COPY[i];
            crcPolynomial[i] = crcPolynomial_COPY[i];
            crcPosition[i] = crcPosition_COPY[i];
            INFO_LENGTHS[i] = INFO_LENGTHS_COPY[i];
            Codeword.at(i).resize(LIST_LENGTH);
        }
        for (int i = 0; i < SEGMENT_LENGTH; i++)
            for (int l = 0; l < LIST_LENGTH; l++)
                Codeword.at(i).at(l).resize(CRC_LENGTHS.at(i) + INFO_LENGTHS.at(i));

        isFrozen.resize(CODE_LENGTH);
        BitInverse.resize(CODE_LENGTH);
        Boundary.resize(CODE_LENGTH);
        for (int i = 0; i < CODE_LENGTH; i++)
        {
            isFrozen[i] = isFrozen_COPY[i];
            BitInverse[i] = BitInverse_COPY[i];
            Boundary[i] = Boundary_COPY[i];
        }
        for (int i = 0; i < CRCInformation_COPY.size(); i++)
        {
            CRCInformation.push_back(CRCInformation_COPY.at(i));
        }

        STAGE_LENGTH = log2(CODE_LENGTH) + 1;
        activePath.resize(LIST_LENGTH, false);
        inactiveArrayIndices_P.resize(STAGE_LENGTH);
        pathIndexToArrayIndex_P.resize(STAGE_LENGTH);
        arrayReferenceCount_P.resize(STAGE_LENGTH);
        inactiveArrayIndices_C.resize(STAGE_LENGTH);
        pathIndexToArrayIndex_C.resize(STAGE_LENGTH);
        arrayReferenceCount_C.resize(STAGE_LENGTH);
        arrayPointer_P.resize(STAGE_LENGTH);
        arrayPointer_C.resize(STAGE_LENGTH);

        for (int lambda = 0; lambda < STAGE_LENGTH; lambda++)
        {
            pathIndexToArrayIndex_P[lambda].resize(LIST_LENGTH);
            arrayReferenceCount_P[lambda].resize(LIST_LENGTH);
            pathIndexToArrayIndex_C[lambda].resize(LIST_LENGTH);
            arrayReferenceCount_C[lambda].resize(LIST_LENGTH);
            arrayPointer_P[lambda].resize(LIST_LENGTH);
            arrayPointer_C[lambda].resize(LIST_LENGTH);
        }

        for (int lambda = 0; lambda < STAGE_LENGTH; lambda++)
        {
            for (int s = 0; s < LIST_LENGTH; s++)
            {
                inactiveArrayIndices_P[lambda].push(s);
                inactiveArrayIndices_C[lambda].push(s);
                arrayPointer_P[lambda][s].resize(1 << (STAGE_LENGTH - 1 - lambda));
                arrayPointer_C[lambda][s].resize(1 << (STAGE_LENGTH - 1 - lambda));
                for (int beta = 0; beta < (1 << (STAGE_LENGTH - 1 - lambda)); beta++)
                {
                    arrayPointer_P[lambda][s][beta].resize(2);
                    arrayPointer_C[lambda][s][beta].resize(2);
                }
            }
        }

        for (int l = 0; l < LIST_LENGTH; l++)
            inactivePathIndices.push(l);
    }
    int assignInitialPath(void)
    {
        int l, s;
        l = inactivePathIndices.top();
        inactivePathIndices.pop();
        activePath[l] = true;
        for (int lambda = 0; lambda < STAGE_LENGTH; lambda++)
        {
            s = inactiveArrayIndices_P[lambda].top();
            inactiveArrayIndices_P[lambda].pop();
            pathIndexToArrayIndex_P[lambda][l] = s;
            arrayReferenceCount_P[lambda][s] = 1;

            s = inactiveArrayIndices_C[lambda].top();
            inactiveArrayIndices_C[lambda].pop();
            pathIndexToArrayIndex_C[lambda][l] = s;
            arrayReferenceCount_C[lambda][s] = 1;
        }
        return l;
    }
    int clonePath(int l)
    {
        int l_, s;
        l_ = inactivePathIndices.top();
        inactivePathIndices.pop();
        activePath[l_] = true;
        for (int lambda = 0; lambda < STAGE_LENGTH; lambda++)
        {
            s = pathIndexToArrayIndex_P[lambda][l];
            pathIndexToArrayIndex_P[lambda][l_] = s;
            arrayReferenceCount_P[lambda][s]++;

            s = pathIndexToArrayIndex_C[lambda][l];
            pathIndexToArrayIndex_C[lambda][l_] = s;
            arrayReferenceCount_C[lambda][s]++;
        }
        return l_;
    }
    void killPath(int l)
    {
        int s;
        activePath[l] = false;
        inactivePathIndices.push(l);
        for (int lambda = 0; lambda < STAGE_LENGTH; lambda++)
        {
            s = pathIndexToArrayIndex_P[lambda][l];
            arrayReferenceCount_P[lambda][s]--;
            if (arrayReferenceCount_P[lambda][s] == 0)
                inactiveArrayIndices_P[lambda].push(s);

            s = pathIndexToArrayIndex_C[lambda][l];
            arrayReferenceCount_C[lambda][s]--;
            if (arrayReferenceCount_C[lambda][s] == 0)
                inactiveArrayIndices_C[lambda].push(s);
        }
    }
    vector<vector<double>> &getArrayPointer_P(int lambda, int l)
    {
        int s, s_;
        s = pathIndexToArrayIndex_P[lambda][l];
        if (arrayReferenceCount_P[lambda][s] == 1)
            s_ = s;
        else
        {
            s_ = inactiveArrayIndices_P[lambda].top();
            inactiveArrayIndices_P[lambda].pop();
            for (int i = 0; i < (1 << (STAGE_LENGTH - 1 - lambda)); i++)
            {
                for (int j = 0; j < 2; j++)
                    arrayPointer_P[lambda][s_][i][j] = arrayPointer_P[lambda][s][i][j];
            }
            arrayReferenceCount_P[lambda][s]--;
            arrayReferenceCount_P[lambda][s_] = 1;
            pathIndexToArrayIndex_P[lambda][l] = s_;
        }
        return arrayPointer_P[lambda][s_];
    }
    vector<vector<bool>> &getArrayPointer_C(int lambda, int l)
    {
        int s, s_;
        s = pathIndexToArrayIndex_C[lambda][l];
        if (arrayReferenceCount_C[lambda][s] == 1)
            s_ = s;
        else
        {
            s_ = inactiveArrayIndices_C[lambda].top();
            inactiveArrayIndices_C[lambda].pop();
            for (int i = 0; i < (1 << (STAGE_LENGTH - 1 - lambda)); i++)
            {
                for (int j = 0; j < 2; j++)
                    arrayPointer_C[lambda][s_][i][j] = arrayPointer_C[lambda][s][i][j];
            }
            arrayReferenceCount_C[lambda][s]--;
            arrayReferenceCount_C[lambda][s_] = 1;
            pathIndexToArrayIndex_C[lambda][l] = s_;
        }
        return arrayPointer_C[lambda][s_];
    }
    int recursivelyCalcP(int lambda, int phi)
    {
        if (lambda == 0)
            return 0;
        int psi = phi / 2;
        bool phiMod = phi % 2;
        if (phiMod == 0)
            recursivelyCalcP(lambda - 1, psi);
        double sigma = 0;
        for (int l = 0; l < LIST_LENGTH; l++)
        {
            if (activePath[l] == false)
                continue;
            vector<vector<double>> &P_lambda = getArrayPointer_P(lambda, l);
            vector<vector<double>> &P_lambda_1 = getArrayPointer_P(lambda - 1, l);
            vector<vector<bool>> &C_lambda = getArrayPointer_C(lambda, l);
            for (int beta = 0; beta < (1 << (STAGE_LENGTH - 1 - lambda)); beta++)
            {
                if (phiMod == 0)
                {
                    for (int u_dash = 0; u_dash < 2; u_dash++)
                    {
                        double cal = 0;
                        for (int u_dash_dash = 0; u_dash_dash < 2; u_dash_dash++)
                            cal += (P_lambda_1[2 * beta][(u_dash ^ u_dash_dash)] / 2) * P_lambda_1[2 * beta + 1][u_dash_dash];
                        P_lambda[beta][u_dash] = cal;
                        sigma = max(sigma, P_lambda[beta][u_dash]);
                    }
                }
                else
                {
                    int u_dash = C_lambda[beta][0];
                    for (int u_dash_dash = 0; u_dash_dash < 2; u_dash_dash++)
                    {
                        P_lambda[beta][u_dash_dash] = (P_lambda_1[2 * beta][(u_dash ^ u_dash_dash)] / 2) * P_lambda_1[2 * beta + 1][u_dash_dash];
                        sigma = max(sigma, P_lambda[beta][u_dash_dash]);
                    }
                }
            }
        }
        for (int l = 0; l < LIST_LENGTH; l++)
        {
            if (activePath[l] == false)
                continue;
            vector<vector<double>> &P_lambda = getArrayPointer_P(lambda, l);
            for (int beta = 0; beta < (1 << (STAGE_LENGTH - 1 - lambda)); beta++)
            {
                for (int u = 0; u < 2; u++)
                    P_lambda[beta][u] /= sigma;
            }
        }
        return 0;
    }
    void recursivelyUpdateC(int lambda, int phi)
    {
        int psi = phi / 2;
        bool psiMod = psi % 2;
        for (int l = 0; l < LIST_LENGTH; l++)
        {
            if (activePath[l] == false)
                continue;
            vector<vector<bool>> &C_lambda = getArrayPointer_C(lambda, l);
            vector<vector<bool>> &C_lambda_1 = getArrayPointer_C(lambda - 1, l);
            for (int beta = 0; beta < (1 << (STAGE_LENGTH - 1 - lambda)); beta++)
            {
                C_lambda_1[2 * beta][psiMod] = (C_lambda[beta][0] ^ C_lambda[beta][1]);
                C_lambda_1[2 * beta + 1][psiMod] = C_lambda[beta][1];
            }
        }
        if (psiMod == 1)
            recursivelyUpdateC(lambda - 1, psi);
    }
    void continuePaths_UnfrozenBit(int phi)
    {
        bool phiMod = phi % 2;
        set<int> S;
        for (int l = 0; l < LIST_LENGTH; l++)
        {
            if (activePath.at(l) == false)
                continue;
            if (S.find(l) != S.end())
                continue;
            vector<vector<bool>> &C_m = getArrayPointer_C(STAGE_LENGTH - 1, l);
            vector<vector<double>> &P_m = getArrayPointer_P(STAGE_LENGTH - 1, l);
            double LLR = log(P_m.at(0).at(0) / P_m.at(0).at(1));
            int sizeStack = inactivePathIndices.size();

            bool bit;
            if (P_m.at(0).at(0) >= P_m.at(0).at(1))
                bit = 0;
            else
                bit = 1;

            if ((fabs(LLR) <= Boundary.at(phi)) && (0 < sizeStack))
            {

                C_m.at(0).at(phiMod) = bit;
                Codeword.at(crcNum).at(l).at(infoNum) = (bit ^ 1);
                int l_ = clonePath(l);
                C_m = getArrayPointer_C(STAGE_LENGTH - 1, l_);
                C_m.at(0).at(phiMod) = (bit ^ 1);
                for (int i = 0; i < infoNum; i++)
                    Codeword.at(crcNum).at(l_).at(i) = Codeword.at(crcNum).at(l).at(i);
                Codeword.at(crcNum).at(l_).at(infoNum) = (bit);
                S.insert(l_);
            }
            else
            {
                C_m.at(0).at(phiMod) = bit;
                Codeword.at(crcNum).at(l).at(infoNum) = bit;
            }
        }
    }
    /*void continuePaths_UnfrozenBit(int phi)
    {
        bool phiMod = phi % 2;
        vector<vector<double>> probForks(LIST_LENGTH, vector<double>(2));
        vector<vector<bool>> contForks(LIST_LENGTH, vector<bool>(2));
        vector<pair<double, pair<int, int>>> list;
        list.resize(2 * LIST_LENGTH);
        int i = 0;
        for (int l = 0; l < LIST_LENGTH; l++)
        {
            if (activePath[l] == true)
            {
                vector<vector<double>> &P_m = getArrayPointer_P(STAGE_LENGTH - 1, l);
                probForks[l][0] = P_m[0][0];
                probForks[l][1] = P_m[0][1];
                i++;
            }
            else
            {
                probForks[l][0] = -1.0;
                probForks[l][1] = -1.0;
            }
        }
        int rho = min(2 * i, LIST_LENGTH);
        int count = 0;
        for (int l = 0; l < LIST_LENGTH; l++)
        {
            for (int b = 0; b < 2; b++)
            {
                list[2 * l + b].first = probForks[l][b];
                list[2 * l + b].second.first = l;
                list[2 * l + b].second.second = b;
            }
        }
        sort(list.rbegin(), list.rend());
        for (i = 0; i < (2 * LIST_LENGTH); i++)
        {
            if (count < rho)
            {
                contForks[list[i].second.first][list[i].second.second] = true;
                count++;
            }
            else
                contForks[list[i].second.first][list[i].second.second] = false;
        }
        for (int l = 0; l < LIST_LENGTH; l++)
        {
            if (activePath[l] == false)
                continue;
            if (contForks[l][0] == false && contForks[l][1] == false)
                killPath(l);
        }
        for (int l = 0; l < LIST_LENGTH; l++)
        {
            if (contForks[l][0] == false && contForks[l][1] == false)
                continue;
            vector<vector<bool>> &C_m = getArrayPointer_C(STAGE_LENGTH - 1, l);
            if (contForks[l][0] == true && contForks[l][1] == true)
            {
                C_m[0][phiMod] = 0;

                Codeword.at(crcNum).at(l).at(infoNum) = 1;

                int l_ = clonePath(l);
                C_m = getArrayPointer_C(STAGE_LENGTH - 1, l_);
                C_m[0][phiMod] = 1;

                for (int i = 0; i < infoNum; i++)
                    Codeword.at(crcNum).at(l_).at(i) = Codeword.at(crcNum).at(l).at(i);
                Codeword.at(crcNum).at(l_).at(infoNum) = 0;
            }
            else
            {
                if (contForks[l][0] == true)
                {
                    C_m[0][phiMod] = 0;

                    Codeword.at(crcNum).at(l).at(infoNum) = 0;
                }

                else
                {
                    C_m[0][phiMod] = 1;

                    Codeword.at(crcNum).at(l).at(infoNum) = 1;
                }
            }
        }
    }*/
    void getCodeword(vector<bool> &EstimatedCRCInformation)
    {
        vector<pair<double, int>> list;
        list.resize(LIST_LENGTH);
        for (int l = 0; l < LIST_LENGTH; l++)
        {
            if (activePath.at(l) == false)
                continue;
            vector<vector<bool>> &C_m = getArrayPointer_C(STAGE_LENGTH - 1, l);
            vector<vector<double>> &P_m = getArrayPointer_P(STAGE_LENGTH - 1, l);
            list.at(l).first = P_m.at(0).at(C_m.at(0).at(1));
            list.at(l).second = l;
        }
        sort(list.rbegin(), list.rend());
        bool flag = true;
        if (crcPolynomial.at(crcNum) != 0)
        {
            for (int l = 0; l < LIST_LENGTH; l++)
            {
                if (activePath.at(list.at(l).second) == false)
                    continue;
                if (checkCRC(CRC_LENGTHS.at(crcNum), INFO_LENGTHS.at(crcNum), crcPolynomial.at(crcNum), Codeword.at(crcNum).at(list.at(l).second)) == true)
                {
                    flag = false;
                    for (int i = 0; i < INFO_LENGTHS.at(crcNum) + CRC_LENGTHS.at(crcNum); i++)
                        EstimatedCRCInformation.push_back(Codeword.at(crcNum).at(list.at(l).second).at(i));
                    for (int i = 0; i < LIST_LENGTH; i++)
                    {
                        if (activePath.at(i) == false)
                            continue;
                        if (i == list.at(l).second)
                            continue;
                        killPath(i);
                    }
                    break;
                }
            }
        }
        if (flag)
        {
            for (int i = 0; i < (INFO_LENGTHS.at(crcNum) + CRC_LENGTHS.at(crcNum)); i++)
                EstimatedCRCInformation.push_back(Codeword.at(crcNum).at(list.at(0).second).at(i));
            for (int i = 0; i < LIST_LENGTH; i++)
            {
                if (activePath.at(i) == false)
                    continue;
                if (i == list.at(0).second)
                    continue;
                killPath(i);
            }
        }
    }
    /*void getCodeword(vector<bool> &EstimatedCRCInfo)
    {
        vector<pair<double, pair<int, int>>> list;
        list.resize(LIST_LENGTH);
        for (int l = 0; l < LIST_LENGTH; l++)
        {
            if (activePath[l] == false)
                continue;
            vector<vector<bool>> &C_m = getArrayPointer_C(STAGE_LENGTH - 1, l);
            vector<vector<double>> &P_m = getArrayPointer_P(STAGE_LENGTH - 1, l);
            list[l].first = P_m[0][C_m[0][1]];
            list[l].second.first = l;
        }
        sort(list.rbegin(), list.rend());
        bool flag = true;
        if (crcPolynomial.at(crcNum) != 0)
        {
            for (int l = 0; l < LIST_LENGTH; l++)
            {
                //C_0更新
                vector<vector<bool>> &C_0 = getArrayPointer_C(0, list[l].second.first);
                vector<bool> C, U, damy;
                C.resize(CODE_LENGTH);
                U.resize(CODE_LENGTH);
                for (int beta = 0; beta < CODE_LENGTH; beta++)
                    C[beta] = C_0[beta][0];
                polarEncode(CODE_LENGTH, C, U, BitInverse);
                int infoCount = 0;
                for (int beta = 0; beta < CODE_LENGTH; beta++)
                {
                    if (isFrozen[beta] == 0)
                    {

                        damy.push_back(U[beta]);
                    }
                }
                if (checkCRC(CRC_LENGTHS.at(crcNum), INFO_LENGTHS.at(crcNum), crcPolynomial.at(crcNum), damy) == true)
                {
                    cout << endl;
                    cout << "list: " << list.at(l).second.first << endl;
                    cout << endl;

                    for (int i = 0; i < damy.size(); i++)
                    {
                        EstimatedCRCInfo.push_back(damy.at(i));

                        cout << i << " d:" << damy.at(i) << " c:" << Codeword.at(crcNum).at(list.at(l).second.first).at(i);
                        if (damy.at(i) != Codeword.at(crcNum).at(list.at(l).second.first).at(i))
                            cout << "Error! ";
                        cout << endl;
                    }
                    flag = false;
                    break;
                }
            }
        }
        if (flag)
        {
            cout << "None" << endl;

            vector<vector<bool>> &C_0 = getArrayPointer_C(0, list[0].second.first);
            vector<bool> C, U;
            C.resize(CODE_LENGTH);
            U.resize(CODE_LENGTH);
            for (int beta = 0; beta < CODE_LENGTH; beta++)
                C[beta] = C_0[beta][0];
            polarEncode(CODE_LENGTH, C, U, BitInverse);
            int infoCount = 0;
            for (int beta = 0; beta < CODE_LENGTH; beta++)
            {
                if (isFrozen[beta] == 0)
                {
                    EstimatedCRCInfo.push_back(U[beta]);
                }
            }
        }
    }*/
};