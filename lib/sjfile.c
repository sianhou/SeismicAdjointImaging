//
// Created by hsa on 09/12/16.
//

#include "sjfile.h"

//! Show information
void sjbasicinformation() {
    printf("\n");
    printf("Author:  Hou, Sian\n");
    printf("E-mail:  sianhou1987@outlook.com\n");
    printf("Website: https://github.com/housian1987/gfdxyz\n");
}

//! Process standard input
int sgpsisep(char *TempString, char *TempValue1, char *TempValue2) {
    int i, j;
    int nchar, num;

    nchar = strlen(TempString);

    i = 0;
    num = 0;
    while (TempString[i] != '=') {
        TempValue1[i] = TempString[i];
        i++;
        num++;
    }
    for (i++, j = 0; i < nchar; i++, j++) {
        TempValue2[j] = TempString[i];
    }
    return num;
}

int sjgetint(int argc, char **argv, char *option, int *value) {
    int i, j;
    int flag, numchar = 0, numopt = 0;
    char *TempValue1;
    char *TempValue2;

    numopt = strlen(option);
    for (i = 1; i < argc; i++) {
        flag = 1;
        TempValue1 = (char *) calloc(2048, sizeof(char));
        TempValue2 = (char *) calloc(2048, sizeof(char));

        numchar = sgpsisep(argv[i], TempValue1, TempValue2);
        if (numopt == numchar) {
            j = 0;
            while (j < numchar) {
                if (TempValue1[j] != option[j]) flag = 0;
                j++;
            }
            if (flag == 1) {
                *value = atoi(TempValue2);

                free((void *) TempValue1);
                free((void *) TempValue2);
                return 1;
            }
        }
        free((void *) TempValue1);
        free((void *) TempValue2);
    }
    return 0;
}

int sjgetfloat(int argc, char **argv, char *option, float *value) {
    int i, j;
    int flag, numchar = 0, numopt = 0;
    char *TempValue1;
    char *TempValue2;

    numopt = strlen(option);
    for (i = 1; i < argc; i++) {
        flag = 1;
        TempValue1 = (char *) calloc(2048, sizeof(char));
        TempValue2 = (char *) calloc(2048, sizeof(char));

        numchar = sgpsisep(argv[i], TempValue1, TempValue2);
        if (numopt == numchar) {
            j = 0;
            while (j < numchar) {
                if (TempValue1[j] != option[j]) flag = 0;
                j++;
            }
            if (flag == 1) {
                *value = atof(TempValue2);

                free((void *) TempValue1);
                free((void *) TempValue2);
                return 1;
            }
        }
        free((void *) TempValue1);
        free((void *) TempValue2);
    }
    return 0;
}

int sjgetchar(int argc, char **argv, char *option, char *value) {
    int i, j, k;
    int flag = 1, numchar = 0, numopt = 0;
    char *TempValue1;
    char *TempValue2;

    numopt = strlen(option);
    for (i = 1; i < argc; i++) {
        flag = 1;
        TempValue1 = (char *) calloc(2048, sizeof(char));
        TempValue2 = (char *) calloc(2048, sizeof(char));

        numchar = sgpsisep(argv[i], TempValue1, TempValue2);
        if (numopt == numchar) {
            j = 0;
            while (j < numchar) {
                if (TempValue1[j] != option[j]) flag = 0;
                j++;
            }
            if (flag == 1) {
                *value = TempValue2[0];

                free((void *) TempValue1);
                free((void *) TempValue2);
                return 1;
            }
        }
        free((void *) TempValue1);
        free((void *) TempValue2);
    }
    return 0;
}

int sjgetstring(int argc, char **argv, char *option, char **value) {
    int i, j;
    int flag, numchar = 0, numopt = 0;
    char *TempValue1;
    char *TempValue2;

    numopt = strlen(option);
    for (i = 1; i < argc; i++) {
        flag = 1;
        TempValue1 = (char *) calloc(2048, sizeof(char));
        TempValue2 = (char *) calloc(2048, sizeof(char));

        numchar = sgpsisep(argv[i], TempValue1, TempValue2);
        if (numopt == numchar) {
            j = 0;
            while (j < numchar) {
                if (TempValue1[j] != option[j]) flag = 0;
                j++;
            }
            if (flag == 1) {
                (*value) = (char *) malloc((strlen(TempValue2) + 1) * sizeof(char));
                for (j = 0; j < strlen(TempValue2); j++) {
                    (*value)[j] = TempValue2[j];
                }
                (*value)[j] = '\0';
                free((void *) TempValue1);
                free((void *) TempValue2);
                return 1;
            }
        }
        free((void *) TempValue1);
        free((void *) TempValue2);
    }
    return 0;
}

//! Read and write GFDXYZ file
long int sggetfilesize(void *filename) {
    FILE *fp;
    if ((fp = fopen(filename, "rb")) != NULL) {
        long int size;
        fseek(fp, 0, SEEK_END);
        size = ftell(fp);
        fclose(fp);
        return size;
    } else {
        printf("Error in sggetfilesize() - can't open the file.\n");
        exit(0);
    }
}

//! Write and read SU file
int sjgetsun1(size_t size, char *inputname) {
    segy tracehead = {0};
    FILE *fpsegy;
    if ((fpsegy = fopen(inputname, "rb")) == NULL) {
        printf("ERROR: Cannot open file in function sjreadsu() !\n");
        exit(0);
    }

    //! Read SEGY header
    fread(&tracehead, sizeof(segy), 1, fpsegy);

    //! Fclose
    fclose(fpsegy);

    return tracehead.ns;
}

int sjgetsun2(size_t size, char *inputname) {
    segy tracehead = {0};
    FILE *fpsegy;
    if ((fpsegy = fopen(inputname, "rb")) == NULL) {
        printf("ERROR: Cannot open file in function sjreadsu() !\n");
        exit(0);
    }

    //! Read SEGY header
    fread(&tracehead, sizeof(segy), 1, fpsegy);

    //! Calculate n2
    int n2 = sggetfilesize(inputname) / (tracehead.ns * size + sizeof(segy));

    //! Fclose
    fclose(fpsegy);

    return n2;
}

int sjwritesu(void *ptr, size_t n2, size_t n1, size_t size, float d1, int ifappend, char *outputname) {
    int ii;
    segy tracehead = {0};

    if (ifappend == 0) { //! Creative a new file or overwrite a old file
        //! Open file
        FILE *fpsegy;
        if ((fpsegy = fopen(outputname, "wb")) == NULL) {
            printf("ERROR: Cannot open file in function sjwritesu() !\n");
            exit(0);
        }

        //! Write n1 and d1
        tracehead.ns = (unsigned short) n1;
        tracehead.dt = (unsigned short) (d1 * 1e6);

        //! Write su file
        for (ii = 0; ii < n2; ii++) {
            tracehead.tracl = tracehead.tracr = ii + 1;
            tracehead.trid = 1;
            //! Write su head
            fwrite(&tracehead, sizeof(segy), 1, fpsegy);
            //! Write data
            fwrite(ptr + ii * n1 * size, size, n1, fpsegy);
        }
        fclose(fpsegy);
        return 1;
    } else { //! append data to an existed file
        //! Open file
        FILE *fpsegy;
        if ((fpsegy = fopen(outputname, "ab")) == NULL) {
            printf("ERROR: Cannot open file in function sjwritesu() !\n");
            exit(0);
        }

        //! Write n1 and d1
        tracehead.ns = (unsigned short) n1;
        tracehead.dt = (unsigned short) (d1 * 1e6);

        //! Calculate tracl and tracr
        int tmp = sjgetsun2(size, outputname);

        //! Write su file
        fseek(fpsegy, 0, SEEK_END);
        for (ii = 0; ii < n2; ++ii) {
            tracehead.tracl = tracehead.tracr = tmp + ii + 1;
            tracehead.trid = 1;
            //! Write su head
            fwrite(&tracehead, sizeof(segy), 1, fpsegy);
            //! Write data
            fwrite(ptr + ii * n1 * size, size, n1, fpsegy);
        }
        fclose(fpsegy);
        return 1;
    }
}

int sjreadsu(void *ptr, size_t n2, size_t n1, size_t size, size_t offsetn2, size_t offsetn1, char *inputname) {
    int ii, ns;
    FILE *fpsegy;
    if ((fpsegy = fopen(inputname, "rb")) == NULL) {
        printf("ERROR: Cannot open file in function sjreadsu() !\n");
        exit(0);
    }

    //! Get ns
    ns = sjgetsun1(size, inputname);

    //! Read SU data
    fseek(fpsegy, offsetn2 * (ns * size + 240L), SEEK_SET); //! Offset in n2 direction
    for (ii = 0; ii < n2; ++ii) {
        fseek(fpsegy, 240L, 1); //! Fseek SU header
        fseek(fpsegy, offsetn1 * size, 1); //! Offset in n1 direction
        fread(ptr + ii * n1 * size, size, n1, fpsegy); //! Read data
    }
    fclose(fpsegy);
    return 1;
}

int sjwritesuall(float *ptr, int n2, int n1, float dt, char *outputname) {
    int ii;
    segy tracehead = {0};
    //! Open file
    FILE *fpsegy;
    if ((fpsegy = fopen(outputname, "wb")) == NULL) {
        printf("ERROR: Cannot open file in function sjwritesuall() !\n");
        exit(0);
    }
    //! Write nsample and dt
    tracehead.ns = (unsigned short) n1;
    tracehead.dt = (unsigned short) (dt * 1e6);

    //! Write SU file
    for (ii = 0; ii < n2; ii++) {
        tracehead.tracl = tracehead.tracr = ii + 1;
        tracehead.trid = 1;
        //! Write SU head
        fwrite(&tracehead, sizeof(segy), 1, fpsegy);
        //! Write data
        fwrite(ptr + ii * n1, sizeof(float), n1, fpsegy);
    }
    fclose(fpsegy);
    return 1;
}

int sjreadsuall(float *ptr, int n2, int n1, char *inputname) {
    int ii;
    FILE *fpsegy;
    if ((fpsegy = fopen(inputname, "rb")) == NULL) {
        printf("ERROR: Cannot open file in function readsegy1() !\n");
        exit(0);
    }

    //! Read SU data
    for (ii = 0; ii < n2; ++ii) {
        //! Fseek SU header
        fseek(fpsegy, 240L, 1);
        //! Read data
        fread(ptr + ii * n1, sizeof(float), n1, fpsegy);
    }
    fclose(fpsegy);
    return 1;
}

//! Survey
int sjgetsvynum() {
    return sizeof(sury) / sizeof(int);
}

int sjgetsvyns(char *inputname) {
    return sjgetsun2(sizeof(int), inputname);
}

int sjgetsvynr(char *inputname) {
    return (sjgetsun1(sizeof(int), inputname) - sjgetsvynum()) / 3;
}

int sjreadsurvey(int is, sury *svy, int *ry, int *rx, int *rz, char *inputname) {
    //! Get the ns
    int nsvy = sjgetsvynum();
    int ns = sjgetsvyns(inputname);
    int nr = sjgetsvynr(inputname);
    if (is < ns) {
        //! Get suy
        sjreadsu(svy, 1, nsvy, sizeof(int), is, 0, inputname);
        //! get ry rx rz
        sjreadsu(ry, 1, nr, sizeof(int), is, nsvy + 0 * nr, inputname);
        sjreadsu(rx, 1, nr, sizeof(int), is, nsvy + 1 * nr, inputname);
        sjreadsu(rz, 1, nr, sizeof(int), is, nsvy + 2 * nr, inputname);
    } else {
        printf("ERROR: Is exceed ns in function sgbin2su()!");
        return 0;
    }
}
