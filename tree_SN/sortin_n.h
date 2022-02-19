
#include <stddef.h>
#define NDIM 2
#define NPOINT 12
#define NETSIZE 16
typedef struct knode knode;
typedef struct kpoint kpoint;
#define TYPE double

struct kpoint {
    TYPE coord[NDIM];
};

struct knode {
    kpoint split;
    knode *left;
    knode *right;
    unsigned short int axis;
};

#define SWAP(A, B)               \
    {                            \
        kpoint t = A;            \
        A.coord[0] = B.coord[0]; \
        A.coord[1] = B.coord[1]; \
        B.coord[0] = t.coord[0]; \
        B.coord[1] = t.coord[1]; \
    }

#define SWAP_N(A, B, AXIS)                   \
    {                                        \
        if (A.coord[AXIS] > B.coord[AXIS]) { \
            kpoint t = A;                    \
            A.coord[0] = B.coord[0];         \
            A.coord[1] = B.coord[1];         \
            B.coord[0] = t.coord[0];         \
            B.coord[1] = t.coord[1];         \
        }                                    \
    }
void sortin_n_3(kpoint *data, int axis);
void sortin_n_4(kpoint *data, int axis);
void sortin_n_5(kpoint *data, int axis);
void sortin_n_6(kpoint *data, int axis);
void sortin_n_7(kpoint *data, int axis);
void sortin_n_8(kpoint *data, int axis);
void sortin_n_9(kpoint *data, int axis);
void sortin_n_10(kpoint *data, int axis);
void sortin_n_11(kpoint *data, int axis);
void sortin_n_12(kpoint *data, int axis);
void sortin_n_13(kpoint *data, int axis);
void sortin_n_14(kpoint *data, int axis);
void sortin_n_15(kpoint *data, int axis);
void sortin_n_16(kpoint *data, int axis);
void sortin_n_17(kpoint *data, int axis);
void sortin_n_18(kpoint *data, int axis);
void sortin_n_19(kpoint *data, int axis);
void sortin_n_20(kpoint *data, int axis);
void sortin_n_21(kpoint *data, int axis);
void sortin_n_22(kpoint *data, int axis);
void sortin_n_23(kpoint *data, int axis);
void sortin_n_24(kpoint *data, int axis);
void sortin_n_25(kpoint *data, int axis);
void sortin_n_26(kpoint *data, int axis);
void sortin_n_27(kpoint *data, int axis);
void sortin_n_28(kpoint *data, int axis);
void sortin_n_29(kpoint *data, int axis);
void sortin_n_30(kpoint *data, int axis);
void sortin_n_31(kpoint *data, int axis);
void sortin_n_32(kpoint *data, int axis);

extern void (*sortin_network[32])(kpoint *x, int y);