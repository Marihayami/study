#include "SCFlip.hpp"
pair<double, double> boxMuller(double standartDeviation)
{
    pair<double, double> res;
    double u1 = double(random()) / RAND_MAX;
    double u2 = double(random()) / RAND_MAX;
    double v1 = sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
    double v2 = sqrt(-2 * log(u1)) * sin(2 * M_PI * u2);
    double ni = standartDeviation * v1;
    double nq = standartDeviation * v2;
    res = make_pair(ni, nq);
    return res;
}
pair<double, double> complexSum(pair<double, double> s, pair<double, double> t)
{
    pair<double, double> res;
    res.first = s.first + t.first;
    res.second = s.second + t.second;
    return res;
}
pair<double, double> complexDifference(pair<double, double> s, pair<double, double> t)
{
    pair<double, double> res;
    res.first = s.first - t.first;
    res.second = s.second - t.second;
    return res;
}
pair<double, double> complexProduct(pair<double, double> s, pair<double, double> t)
{
    pair<double, double> res;
    res.first = s.first * t.first - s.second * t.second;
    res.second = s.second * t.first + s.first * t.second;
    return res;
}
pair<double, double> complexQuotient(pair<double, double> s, pair<double, double> t)
{
    pair<double, double> res;
    double square_norm = pow(t.first, 2) + pow(t.second, 2);
    pair<double, double> t_conjugate = complexConjugate(t);
    s = scalarQuotient(s, square_norm);
    res = complexProduct(s, t_conjugate);
    return res;
}
pair<double, double> complexConjugate(pair<double, double> s)
{
    pair<double, double> res;
    res.first = s.first;
    res.second = -s.second;
    return res;
}
pair<double, double> scalarProduct(pair<double, double> s, double t)
{
    pair<double, double> res;
    res.first = s.first * t;
    res.second = s.second * t;
    return res;
}
pair<double, double> scalarQuotient(pair<double, double> s, double t)
{
    pair<double, double> res;
    res.first = s.first / t;
    res.second = s.second / t;
    return res;
}
void awgn(int CODE_LENGTH, double standartDeviation, vector<pair<double, double>> SentSymbol, vector<pair<double, double>> &ReceivedSymbol)
{
    for (int j = 0; j < CODE_LENGTH; j++)
        ReceivedSymbol[j] = complexSum(SentSymbol[j], boxMuller(standartDeviation));
}
void bpskModulation(int CODE_LENGTH, double A, vector<bool> EncodedArray, vector<pair<double, double>> &SentSymbol)
{
    for (int j = 0; j < CODE_LENGTH; j++)
    {
        if (EncodedArray[j] == 0)
            SentSymbol[j].first = -A;
        else
            SentSymbol[j].first = A;
    }
}
void bpskDemodulation(int CODE_LENGTH, double standartDeviation, vector<pair<double, double>> ReceivedSymbol, vector<double> &ReceivedLLR)
{
    for (int j = 0; j < CODE_LENGTH; j++)
        ReceivedLLR[j] = bpskLLR(standartDeviation, ReceivedSymbol[j]);
}
void crcEncode(int CRC_LENGTH, int INFO_LENGTH, int polynominal, vector<bool> &CRCInformation)
{
    int i;
    long int dividend = 0;
    for (i = 0; i < CRC_LENGTH; i++)
    {
        dividend <<= 1;
        dividend += CRCInformation[i];
    }
    for (i = CRC_LENGTH; i < CRC_LENGTH + INFO_LENGTH; i++)
    {
        dividend <<= 1;
        dividend += CRCInformation[i];
        if ((dividend >> CRC_LENGTH) == 0)
            continue;
        dividend ^= polynominal;
    }
    for (i = 0; i < CRC_LENGTH; i++)
        CRCInformation[INFO_LENGTH + i] = (dividend >> (CRC_LENGTH - 1 - i)) & 1;
}
void errorCount(int INFO_LENGTH, double &bitError, double &flameError, vector<bool> Information, vector<bool> EstimatedInformation)
{
    bool flag = true;
    for (int j = 0; j < INFO_LENGTH; j++)
    {
        if (EstimatedInformation[j] != Information[j])
        {
            bitError++;
            if (flag)
            {
                flameError++;
                flag = false;
            }
        }
    }
}
void generateBitInverse(int CODE_LENGTH, vector<int> &BitInverse)
{
    int i, j, p1, p2, ni;
    int n = (int)(log(CODE_LENGTH) / log(2.0) + 0.5);
    BitInverse[0] = 0;
    p1 = p2 = 1;
    for (i = 1; i <= n; i++)
    {
        p2 *= 2;
        ni = CODE_LENGTH / p2;
        for (j = p1; j < p2; j++)
        {
            BitInverse[j] = BitInverse[j - p1] + ni;
        }
        p1 = p2;
    }
}
void generateChanellLLR(int CODE_LENGTH, double EbN0dB, vector<double> &ChannelLLR)
{
    FILE *fp;
    char filename[256];
    sprintf(filename, "./ChannelLLR/%.2lfLen%d.txt", EbN0dB, CODE_LENGTH);
    if ((fp = fopen(filename, "r")) == NULL)
    {
        printf("Error occured6\n");
        exit(1);
    }
    for (int i = 0; i < CODE_LENGTH; i++)
        fscanf(fp, "%lf", &ChannelLLR[i]);
    fclose(fp);
}
void generateFrozenArray(int CODE_LENGTH, int frozenNum, vector<bool> &isFrozen, vector<double> ChannelLLR)
{
    vector<pair<double, int>> v;
    for (int i = 0; i < CODE_LENGTH; i++)
        v.push_back(make_pair(ChannelLLR[i], i));
    sort(v.begin(), v.end());
    for (int i = 0; i < frozenNum; i++)
        isFrozen[int(v[i].second)] = 1;
}
void generateInformationArray(int CRC_LENGTH, int INFO_LENGTH, vector<bool> &Information, vector<bool> &CRCInformation)
{
    for (int i = 0; i < INFO_LENGTH; i++)
    {
        // Information[i] = randBin();
        Information[i] = 0;
        CRCInformation[i] = Information[i];
    }
    for (int i = INFO_LENGTH; i < INFO_LENGTH + CRC_LENGTH; i++)
    {
        CRCInformation[i] = 0;
    }
}
void generateInformationCount(int CODE_LENGTH, vector<bool> isFrozen, vector<int> &InformationCount)
{
    for (int i = 0, index = 0; i < CODE_LENGTH; i++)
    {
        if (isFrozen[i])
            continue;
        InformationCount[i] = index++;
    }
}
void generateInputArray(int CODE_LENGTH, vector<bool> CRCInformation, vector<bool> &InputArray, vector<bool> isFrozen, vector<int> InformationCount)
{
    for (int i = 0; i < CODE_LENGTH; i++)
    {
        if (isFrozen[i])
            InputArray[i] = 0;
        else
            InputArray[i] = CRCInformation[InformationCount[i]];
    }
}
void initializeCRC(int CRC_LENGTH, int &crcPolynomial)
{
    if (CRC_LENGTH == 4)
    {
        crcPolynomial = 0x9;
        crcPolynomial = (crcPolynomial << 1) + 1;
    }
    else if (CRC_LENGTH == 8)
    {
        crcPolynomial = 0xA6;
        crcPolynomial = (crcPolynomial << 1) + 1;
    }
    else if (CRC_LENGTH == 10)
    {
        crcPolynomial = 0x327;
        crcPolynomial = (crcPolynomial << 1) + 1;
    }
    else if (CRC_LENGTH == 16)
    {
        crcPolynomial = 0x8810;
        crcPolynomial = (crcPolynomial << 1) + 1;
    }
    else if (CRC_LENGTH == 0)
    {
        crcPolynomial = 0;
    }
    else
    {
        printf("Error Occured\n");
        exit(1);
    }
}
void polarEncode(int CODE_LENGTH, vector<bool> InputArray, vector<bool> &EncodedArray, vector<int> BitInverse)
{
    set<int> Check;
    int left, right;
    bool cl, cr;
    for (int i = 0; i < CODE_LENGTH; i++)
        EncodedArray[i] = InputArray[i];
    int stageNum = log2(CODE_LENGTH);
    for (int i = 0; i < stageNum; i++)
    {
        int jumpNum = 1 << i;
        for (int j = 0; j < jumpNum; j++)
        {
            int slideNum = CODE_LENGTH / (2 * jumpNum);
            for (int k = 0; k < slideNum; k++)
            {
                if (k == 0)
                {
                    left = 0 + j * slideNum * 2;
                    right = slideNum + j * slideNum * 2;
                }
                cl = EncodedArray[left];
                cr = EncodedArray[right];
                EncodedArray[left] = boolSum(cl, cr);
                left++;
                right++;
            }
        }
    }
    for (int i = 0; i < CODE_LENGTH; i++) //ビット反転処理
    {
        if (Check.find(i) == Check.end())
        {
            bool tem;
            tem = EncodedArray[i];
            int j = BitInverse[i];
            EncodedArray[i] = EncodedArray[j];
            EncodedArray[j] = tem;
            Check.insert(i);
            Check.insert(j);
        }
    }
}
bool boolSum(bool s, bool t)
{
    bool res;
    int cal = (s + t) % 2;
    res = cal;
    return res;
}
bool checkCRC(int CRC_LENGTH, int INFO_LENGTH, int crcPolynomial, vector<bool> EstimatedCRCInformation)
{

    int i;
    long int dividend = 0;

    for (i = 0; i < CRC_LENGTH; i++)
    {
        dividend <<= 1;
        dividend += EstimatedCRCInformation[i];
    }

    for (i = CRC_LENGTH; i < CRC_LENGTH + INFO_LENGTH; i++)
    {
        dividend <<= 1;
        dividend += EstimatedCRCInformation[i];

        if ((dividend >> CRC_LENGTH) == 0)
            continue;
        dividend ^= crcPolynomial;
    }
    if (dividend == 0)
        return true;
    else
        return false;
}
bool errorCount(int INFO_LENGTH, vector<bool> Information, vector<bool> EstimatedInformation)
{
    for (int j = 0; j < INFO_LENGTH; j++)
    {
        if (EstimatedInformation[j] != Information[j])
            return false;
    }
    return true;
}
int randBin()
{
    return (random() & 1);
}
double bpskLLR(double standartDeviation, pair<double, double> ReceivedSymbol)
{
    double res = -2.0 * ReceivedSymbol.first / standartDeviation / standartDeviation;
    return res;
}
double complexNorm(pair<double, double> s)
{
    double res;
    res = sqrt(s.first * s.first + s.second * s.second);
    return res;
}
