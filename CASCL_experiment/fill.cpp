#include "bits/stdc++.h"
using namespace std;
int main()
{
    double EbN0dB = 2.00;
    int CODE_LENGTH = 1024;
    for (int phi = 0; phi < CODE_LENGTH; phi++)
    {
        FILE *fp;
        char filename[256];
        sprintf(filename, "./calculatedLLR/preError/thisError/SN%.2lf/phi%04d.txt", EbN0dB, phi);
        if ((fp = fopen(filename, "a")) == NULL)
        {
            printf("Error occured1\n");
            exit(1);
        }
        fclose(fp);
    }
    for (int phi = 0; phi < CODE_LENGTH; phi++)
    {
        FILE *fp;
        char filename[256];
        sprintf(filename, "./calculatedLLR/preError/thisNonError/SN%.2lf/phi%04d.txt", EbN0dB, phi);
        if ((fp = fopen(filename, "a")) == NULL)
        {
            printf("Error occured2\n");
            exit(1);
        }
        fclose(fp);
    }
    for (int phi = 0; phi < CODE_LENGTH; phi++)
    {
        FILE *fp;
        char filename[256];
        sprintf(filename, "./calculatedLLR/preNonError/thisError/SN%.2lf/phi%04d.txt", EbN0dB, phi);
        if ((fp = fopen(filename, "a")) == NULL)
        {
            printf("Error occured3\n");
            exit(1);
        }
        fclose(fp);
    }

    for (int phi = 0; phi < CODE_LENGTH; phi++)
    {
        FILE *fp;
        char filename[256];
        sprintf(filename, "./calculatedLLR/preNonError/thisNonError/SN%.2lf/phi%04d.txt", EbN0dB, phi);
        if ((fp = fopen(filename, "a")) == NULL)
        {
            printf("Error occured4\n");
            exit(1);
        }
        fclose(fp);
    }
}
