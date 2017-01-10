//
// Created by hsa on 09/12/16.
//

#include "sjfile.h"

//! Show information
void sjbasicinformation() {
    printf("\n");
    printf("Author:     Hou, Sian\n");
    printf("E-mail:     sianhou1987@outlook.com\n");
    printf("Website:    https://github.com/housian1987/SeismicAdjointModeling\n");
    exit(0);
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

//! Read and write file
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
int sjssurvey_init(sjssurvey *ptr) {
    ptr->sy = -1;
    ptr->sx = -1;
    ptr->sz = -1;
    ptr->y0 = -1;
    ptr->x0 = -1;
    ptr->z0 = -1;
    ptr->ny = -1;
    ptr->nx = -1;
    ptr->nz = -1;
    ptr->gny = -1;
    ptr->gnx = -1;
    ptr->gnz = -1;
    ptr->nr = -1;
    ptr->tr = -1;
    ptr->ns = -1;
    ptr->maxnr = -1;
    ptr->nparas = 14;
    ptr->ry = NULL;
    ptr->rx = NULL;
    ptr->rz = NULL;
    ptr->surveyfile = NULL;
    return 1;
}

int sjssurvey_display(sjssurvey *ptr) {
    printf("Display survey information:\n");
    printf("  surveyfile:  %s\n", ptr->surveyfile);
    return 1;
}

int sjssurvey_readis(sjssurvey *ptr, int is) {

    int nr = ptr->maxnr;
    int np = ptr->nparas;

    if (is < ptr->ns) {
        //! Allocate memory
        sjcheckfree1d((void *) ptr->ry);
        sjcheckfree1d((void *) ptr->rx);
        sjcheckfree1d((void *) ptr->rz);
        ptr->ry = (int *) sjalloc1d(nr, sizeof(int));
        ptr->rx = (int *) sjalloc1d(nr, sizeof(int));
        ptr->rz = (int *) sjalloc1d(nr, sizeof(int));

        //! Get parameters
        sjreadsu(ptr, 1, np, sizeof(int), is, 0, ptr->surveyfile);

        //! get ry rx rz
        sjreadsu(ptr->ry, 1, nr, sizeof(int), is, np + 0 * nr, ptr->surveyfile);
        sjreadsu(ptr->rx, 1, nr, sizeof(int), is, np + 1 * nr, ptr->surveyfile);
        sjreadsu(ptr->rz, 1, nr, sizeof(int), is, np + 2 * nr, ptr->surveyfile);
        return 1;
    } else {
        printf("ERROR: Is exceed ns in function sjssurvey_readis()!");
        return 0;
    }
}

int sjssurvey_write(sjssurvey *ptr, int ifappend) {
    int nr = ptr->maxnr;
    int np = ptr->nparas;
    int *tmp = (int *) sjalloc1d(np + 3 * nr, sizeof(int));
    memcpy(tmp + 0 * np + 0 * nr, ptr, np * sizeof(int));
    memcpy(tmp + 1 * np + 0 * nr, ptr->ry, nr * sizeof(int));
    memcpy(tmp + 1 * np + 1 * nr, ptr->rx, nr * sizeof(int));
    memcpy(tmp + 1 * np + 2 * nr, ptr->rz, nr * sizeof(int));
    sjwritesu((void *) tmp, 1, np + 3 * nr, sizeof(int), 1.0f, ifappend, ptr->surveyfile);
    return 1;
}

int sjssurvey_getparas(sjssurvey *ptr, int argc, char **argv) {

    if (argc == 1) {
        printf("* survey:      Input filename of seismic survey.\n");
        return 0;
    } else {
        if (!sjmgets("survey", ptr->surveyfile)) {
            printf("ERROR: Should set survey file!\n");
            return 0;
        };

        //! Get ns
        ptr->ns = sjgetsun2(sizeof(int), ptr->surveyfile);

        //! Get maxnr
        ptr->maxnr = (sjgetsun1(sizeof(int), ptr->surveyfile) - ptr->nparas) / 3;

        //! Get parameters
        sjreadsu(ptr, 1, ptr->nparas, sizeof(int), 0, 0, ptr->surveyfile);

        return 1;
    }
}

//! Geology
int sjsgeo_init(sjsgeology *ptr) {
    ptr->gvp2d = NULL;
    ptr->gvs2d = NULL;
    ptr->vp2d = NULL;
    ptr->vp2d = NULL;
    ptr->ipp2d = NULL;
    ptr->nipp2d = NULL;
    ptr->gipp2d = NULL;
    ptr->vpfile = NULL;
    ptr->vsfile = NULL;
    ptr->ippfile = NULL;
    ptr->lsippfile = NULL;
    return 1;
}

int sjsgeo_display(sjsgeology *ptr) {
    printf("Display geo information:\n");
    printf("  vpfile:      %s\n", ptr->vpfile);
    printf("  vsfile:      %s\n", ptr->vsfile);
    printf("  ippfile:     %s\n", ptr->ippfile);
    printf("  lsippfile:   %s\n", ptr->lsippfile);
    return 1;
}

int sjsgeo_getparas2d(sjsgeology *ptr, int argc, char **argv, char *info) {
    if (strcmp(info, "vp") == 0) {
        if (argc == 1) {
            printf("* vp:          2D seismic P-velocity.\n");
            return 0;
        } else {
            if (!sjmgets("vp", ptr->vpfile)) {
                printf("ERROR: Should set vp file!\n");
                exit(0);
            }
            return 1;
        }
    }
    if (strcmp(info, "vs") == 0) {
        if (argc == 1) {
            printf("* vs:          2D seismic S-velocity.\n");
            return 0;
        } else {
            if (!sjmgets("vs", ptr->vsfile)) {
                printf("ERROR: Should set vs file!\n");
                exit(0);
            }
            return 1;
        }
    }
    if (strcmp(info, "ipp") == 0) {
        if (argc == 1) {
            printf("* ipp:         Image file of P-P wave.\n");
            return 0;
        } else {
            if (!sjmgets("ipp", ptr->ippfile)) {
                printf("ERROR: Should set ipp file!\n");
                exit(0);
            }
            return 1;
        }
    }
    if (strcmp(info, "lsipp") == 0) {
        if (argc == 1) {
            printf("* lsipp:       Least square image file of P-P wave.\n");
            return 0;
        } else {
            if (!sjmgets("lsipp", ptr->lsippfile)) {
                printf("ERROR: Should set lsipp file!\n");
                exit(0);
            }
            return 1;
        }
    }
}

//! Wavefield
int sjswave_init(sjswave *ptr) {
    ptr->recy = NULL;
    ptr->recx = NULL;
    ptr->recz = NULL;
    ptr->snapy2d = NULL;
    ptr->snapx2d = NULL;
    ptr->snapz2d = NULL;
    ptr->recyfile = NULL;
    ptr->recxfile = NULL;
    ptr->reczfile = NULL;
    return 1;
}

int sjswave_display(sjswave *ptr) {
    printf("Display wave information:\n");
    printf("  recyfile:    %s\n", ptr->recyfile);
    printf("  recxfile:    %s\n", ptr->recxfile);
    printf("  reczfile:    %s\n", ptr->reczfile);
    return 1;
}

int sjswave_getparas(sjswave *ptr, int argc, char **argv, char *info) {
    if (strcmp(info, "recy") == 0) {
        if (argc == 1) {
            printf("* recy:        Seicmic record in y (cxline) - direction.\n");
            return 0;
        } else {
            if (!sjmgets("recy", ptr->recyfile)) {
                printf("ERROR: Should set recy!\n");
                exit(0);
            }
            return 1;
        }
    }
    if (strcmp(info, "recx") == 0) {
        if (argc == 1) {
            printf("* recx:        Seicmic record in x (inline) - direction.\n");

            return 0;
        } else {
            if (!sjmgets("recx", ptr->recxfile)) {
                printf("ERROR: Should set recx!\n");
                exit(0);
            }
            return 1;
        }
    }
    if (strcmp(info, "recz") == 0) {
        if (argc == 1) {
            printf("* recz:        Seicmic record in z (depth) - direction.\n");
            return 0;
        } else {
            if (!sjmgets("recz", ptr->reczfile)) {
                printf("ERROR: Should set recz!\n");
                exit(0);
            }
            return 1;
        }
    }

    return 0;
}

//! Option
int sjsoption_init(sjsoption *ptr) {
    ptr->nt = -1;
    ptr->k1 = -1;
    ptr->jsnap = -1;
    ptr->nsnap = -1;
    ptr->srcrange = -1;
    ptr->srctrunc = -1;
    ptr->ystacksrc = -1;
    ptr->dt = -1.0f;
    ptr->fp = -1.0f;
    ptr->amp = -1.0f;
    ptr->srcdecay = -1.0f;
    ptr->nb = -1;
    ptr->ds = -1.0f;
    ptr->ycutdirect = -1;
    ptr->ycalscatter = -1;
    return 1;
}

int sjsoption_display(sjsoption *ptr) {
    printf("Display option information:\n");
    printf("  nt:          %d\n", ptr->nt);
    printf("  k1:          %d\n", ptr->k1);
    printf("  jsnap:       %d\n", ptr->jsnap);
    printf("  nsnap:       %d\n", ptr->nsnap);
    printf("  srcrange:    %d\n", ptr->srcrange);
    printf("  srctrunc:    %d\n", ptr->srctrunc);
    printf("  ystacksrc:   %d\n", ptr->ystacksrc);
    printf("  dt:          %f\n", ptr->dt);
    printf("  fp:          %f\n", ptr->fp);
    printf("  amp:         %f\n", ptr->amp);
    printf("  srcdecay:    %f\n", ptr->srcdecay);
    printf("  nb:          %d\n", ptr->nb);
    printf("  ds:          %f\n", ptr->ds);
    printf("  ycutdirect:  %d\n", ptr->ycutdirect);
    printf("  ycalscatter: %d\n", ptr->ycalscatter);
    return 1;
}

int sjsoption_getparas(sjsoption *ptr, int argc, char **argv) {
    if (argc == 1) {
        printf("* nt:          Total calculate time.\n");
        printf("  k1:          Peak position of wavelet, default = 50.\n");
        printf("  jsnap:       Interval of snap, default = 2.\n");
        printf("  srcrange:    Total range of source, default = 10.\n");
        printf("  srctrunc:    Total time of source, default = 301.\n");
        printf("  ystacksrc:   Process source, 0: boundary condition;\n");
        printf("                               1: stacking (default);\n");
        printf("  dt:          Interval of time(s), default = 0.001.\n");
        printf("  fp:          Peak frequency of wavelet, default = 20.\n");
        printf("  amp:         Peak amplitude of wavelet, default = 1.0.\n");
        printf("  srcdecay:    Decay of source, default = 0.4.\n");
        printf("  nb:          Range of ABC, default = 15.\n");
        printf("  ds:          Interval of space(m), default = 10.0.\n");
        printf("  ycutdirect:  Cut direct wave, 0: don't cut;\n");
        printf("                                1: cut direct (default).\n");
        printf("  ycalscatter: Calculate scatter, 0: don't calculate (default);\n");
        printf("                                  1: calculate.\n");
        return 0;
    } else {
        if (!sjmgeti("nt", ptr->nt)) {
            printf("ERROR: Should set nt!\n");
            exit(0);
        }
        if (!sjmgeti("k1", ptr->k1)) ptr->k1 = 50;
        if (!sjmgeti("jsnap", ptr->jsnap)) ptr->jsnap = 2;
        ptr->nsnap = (ptr->nt - 1) / ptr->jsnap + 1;
        if (!sjmgeti("srcrange", ptr->srcrange)) ptr->srcrange = 10;
        if (!sjmgeti("srctrunc", ptr->srctrunc)) ptr->srctrunc = 301;
        if (!sjmgeti("ystacksrc", ptr->ystacksrc)) ptr->ystacksrc = 1;
        if (!sjmgetf("dt", ptr->dt)) ptr->dt = 0.001;
        if (!sjmgetf("fp", ptr->fp)) ptr->fp = 20.0;
        if (!sjmgetf("amp", ptr->amp)) ptr->amp = 1.0;
        if (!sjmgetf("srcdecay", ptr->srcdecay)) ptr->srcdecay = 0.4;
        if (!sjmgeti("nb", ptr->nb)) ptr->nb = 15;
        if (!sjmgetf("ds", ptr->ds)) ptr->ds = 10.0;
        if (!sjmgeti("ycutdirect", ptr->ycutdirect)) ptr->ycutdirect = 1;
        if (!sjmgeti("ycalscatter", ptr->ycalscatter)) ptr->ycalscatter = 0;
        return 1;
    }
}
