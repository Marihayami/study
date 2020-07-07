#include "SCFlip.hpp"
signed main()
{
    srandom(time(0));

    /*パラメータ*/
    double R = 0.5; //コードレート

    int LIMIT_REPEAT_NUM = 1000;

    int CODE_LENGTH = 1024;
    int CRC_LENGTH = 16;
    int CRCINFO_LENGTH = CODE_LENGTH * R;
    int FROZEN_LENGTH = CODE_LENGTH - CRCINFO_LENGTH;
    int INFO_LENGTH = CRCINFO_LENGTH - CRC_LENGTH;
    int LIST_LENGTH = 32;

    double Rate = double(INFO_LENGTH) / CODE_LENGTH; //実行ビットレート

    int segNum = 6;
    vector<int> CRC_LENGTHS = {5, 4, 4, 1, 1, 1};
    vector<int> crcPolynomial(segNum);
    vector<int> crcPosition = {511, 767, 895, 959, 991, 1023};
    /*int segNum = 1;
    vector<int> CRC_LENGTHS = {16};
    vector<int> crcPolynomial(segNum);
    vector<int> crcPosition = {1023};*/

    /*パラメータ-end*/

    /*宣言*/
    // int crcPolynomial;
    int loopNum;
    double EbN0;
    double standartDeviation;
    double bitError;
    double flameError;
    /*宣言-end*/

    /*定義*/

    int M = 2;    //変調多値数
    int Ns = 1;   //1シンボルあたりのサンプル数
    double A = 1; //振幅
    double SN_WIDTH = 0.25;
    /*定義-end*/

    /*SN範囲を指定*/
    double SN_MIN;
    double SN_MAX;
    printf("---------------------------------------------------------------------\n");
    printf("please tell me the value of SN_MIN:");
    scanf("%lf", &SN_MIN);
    printf("please tell me the value of SN_MAX:");
    scanf("%lf", &SN_MAX);
    printf("\n");
    /*SN範囲を指定ーend*/

    /*CRC初期設定*/
    initializeCRC(CRC_LENGTHS, crcPolynomial);
    /*CRC初期設定-end*/

    /*パラメータ出力*/
    printf("CODE_LENGTH: %d\n", CODE_LENGTH);
    printf("CRC_LENGTH: %d\n", CRC_LENGTH);
    printf("CRCINFO_LENGTH: %d\n", CRCINFO_LENGTH);
    printf("LIST_LENGTH: %d\n", LIST_LENGTH);
    printf("FROZEN_LENGTH: %d\n", FROZEN_LENGTH);
    printf("INFO_LENGTH: %d\n", INFO_LENGTH);
    printf("Rate: %f\n", Rate);
    printf("LIMIT_REPEAT_NUM: %d\n", LIMIT_REPEAT_NUM);
    printf("SN_MIN: %.2f\n", SN_MIN);
    printf("SN_MAX: %.2f\n", SN_MAX);
    printf("SN_WIDTH: %.2f\n", SN_WIDTH);
    printf("\n");
    /*パラメータ出力-end*/

    /*逆ビット列生成*/
    vector<int> BitInverse(CODE_LENGTH);
    generateBitInverse(CODE_LENGTH, BitInverse);
    /*逆ビット列生成-end*/

    /*Simulation*/
    for (double EbN0dB = SN_MIN; EbN0dB <= SN_MAX; EbN0dB += SN_WIDTH)
    {

        /*各ループ初期設定*/
        loopNum = 0;
        bitError = 0;
        flameError = 0;
        EbN0 = pow(10.0, EbN0dB / 10.0);
        standartDeviation = sqrt(A * A * Ns / (2 * (Rate)*log2(M) * EbN0));
        /*各ループ初期設定-end*/

        /*各チャネル特性を取得*/
        vector<double> ChannelLLR(CODE_LENGTH);
        generateChanellLLR(CODE_LENGTH, EbN0dB, ChannelLLR);
        /*各チャネル特性を取得-end*/

        /*凍結ビット決定*/
        vector<bool> isFrozen(CODE_LENGTH);
        generateFrozenArray(CODE_LENGTH, FROZEN_LENGTH, isFrozen, ChannelLLR);
        /*凍結ビット決定-end*/

        /*InformationCount作成*/
        vector<int> InformationCount(CODE_LENGTH, -1);
        generateInformationCount(CODE_LENGTH, isFrozen, InformationCount);
        /*InformationCount作成-end*/

        vector<int> INFO_LENGTHS(segNum, 0);
        int index = 0;
        for (int i = 0; i < CODE_LENGTH; i++)
        {
            if (isFrozen[i])
                continue;
            INFO_LENGTHS[index]++;
            if (i == crcPosition[index])
            {
                INFO_LENGTHS[index] -= CRC_LENGTHS[index];
                index++;
            }
        }

        vector<double> Boundary(CODE_LENGTH);
        FILE *fp;
        char filename[256];
        sprintf(filename, "./Boundary/SN%.2lf.txt", EbN0dB);
        if ((fp = fopen(filename, "r")) == NULL)
        {
            printf("Error occured6\n");
            exit(1);
        }
        for (int i = 0; i < CODE_LENGTH; i++)
            fscanf(fp, "%lf", &Boundary[i]);
        fclose(fp);

        /*Loop*/
        for (int i = 0; i < LIMIT_REPEAT_NUM; i++)
        {
            /*配列宣言*/
            vector<bool> Information;
            vector<bool> CRCInformation;
            vector<bool> InputArray(CODE_LENGTH);
            vector<bool> EncodedArray(CODE_LENGTH);
            vector<bool> EstimatedCRCInformation;
            vector<bool> EstimatedInformation;
            vector<double> ReceivedLLR(CODE_LENGTH);
            vector<pair<double, double>> SentSymbol(CODE_LENGTH);
            vector<pair<double, double>> ReceivedSymbol(CODE_LENGTH);
            /*配列宣言-end*/

            /*情報ビット列生成*/
            generateInformationArray(Information, CRCInformation, CRC_LENGTHS, INFO_LENGTHS);
            /*情報ビット列生成-end*/

            /*CRC付与*/
            crcEncode(CRCInformation, CRC_LENGTHS, crcPolynomial, INFO_LENGTHS);
            /*CRC付与-end*/

            /*入力ビット列生成*/
            generateInputArray(CODE_LENGTH, CRCInformation, InputArray, isFrozen, InformationCount);
            /*入力ビット列生成-end*/

            /*Polar Encode*/
            polarEncode(CODE_LENGTH, InputArray, EncodedArray, BitInverse);
            /*Polar Encode-end*/

            /*Modulation*/
            bpskModulation(CODE_LENGTH, A, EncodedArray, SentSymbol);
            /*Modulation-end*/

            /*AWGN*/
            awgn(CODE_LENGTH, standartDeviation, SentSymbol, ReceivedSymbol);
            /*AWGN-end*/

            /*Demodulation*/
            bpskDemodulation(CODE_LENGTH, standartDeviation, ReceivedSymbol, ReceivedLLR);
            /*Demodulation-end*/

            /*Decoder*/
            CASCL_Decoder *Decoder;
            Decoder = new CASCL_Decoder(CODE_LENGTH, CRC_LENGTH, INFO_LENGTH, LIST_LENGTH, segNum, CRCInformation, EstimatedCRCInformation, isFrozen, BitInverse, CRC_LENGTHS, crcPolynomial, crcPosition, INFO_LENGTHS, Boundary, ReceivedLLR);
            delete Decoder;

            /*Decoder-end*/

            int cal = 0;
            for (int i = 0; i < segNum; i++)
            {
                for (int j = 0; j < INFO_LENGTHS[i]; j++)
                {
                    EstimatedInformation.push_back(EstimatedCRCInformation.at(cal + j));
                }
                cal += INFO_LENGTHS[i] + CRC_LENGTHS[i];
            }

            /*Count*/
            errorCount(INFO_LENGTH, bitError, flameError, Information, EstimatedInformation);
            /*Count-end*/

            loopNum++;
            if (loopNum >= 1000 && flameError >= 100)
                break;
        }
        /*Loop-end*/

        /*SN出力*/
        printf("%lf,%lf\n", EbN0dB, flameError / loopNum);
        /*SN出力-end*/
    }
    /*Simulation-end*/
    printf("---------------------------------------------------------------------\n");
    return 0;
}
