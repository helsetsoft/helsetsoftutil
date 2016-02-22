VI 0 0 3 0 0 2
MODULE NR,0 0 0
FILE 0,../../nr.f90
GMODPROC SPRSAX: SPRSAX_DP
GMODPROC SPRSAX: SPRSAX_SP
GMODPROC SPRSTX: SPRSTX_DP
GMODPROC SPRSTX: SPRSTX_SP
GMODPROC GAMMLN: GAMMLN_V
GMODPROC GAMMLN: GAMMLN_S
GMODPROC SPRSDIAG: SPRSDIAG_DP
GMODPROC SPRSDIAG: SPRSDIAG_SP
GMODPROC RAN: RAN
GMODPROC RAN0: RAN0_V
GMODPROC RAN0: RAN0_S
GMODPROC RAN1: RAN1_V
GMODPROC RAN1: RAN1_S
GMODPROC INDEXX: INDEXX_I4B
GMODPROC INDEXX: INDEXX_SP
GMODPROC RAN2: RAN2_V
GMODPROC RAN2: RAN2_S
GMODPROC GASDEV: GASDEV_V
GMODPROC GASDEV: GASDEV_S
GMODPROC SVBKSB: SVBKSB_DP
GMODPROC SVBKSB: SVBKSB_SP
GMODPROC SVDCMP: SVDCMP_DP
GMODPROC SVDCMP: SVDCMP_SP
GMODPROC PYTHAG: PYTHAG_DP
GMODPROC PYTHAG: PYTHAG_SP
GMODPROC SPRSIN: SPRSIN_DP
GMODPROC SPRSIN: SPRSIN_SP
PROC GAUSSJ,2,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:SWAP=>SWAP,OUTERPROD=>OUTERPROD,OUTERAND=>OUTERAND,NRERROR=>NRERROR,ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR A,3,,: 2,1,4,0,1,183,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,3
VAR B,3,,: 2,1,4,0,1,183,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,3
ENDPROC
PROC BCUINT,13,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:NRERROR=>NRERROR
USE NRTYPE 2
VAR Y,3,,: 2,1,4,0,1,3,0 (1,1,4: 1,1,1,4),0,1000,1
VAR Y1,3,,: 2,1,4,0,1,3,0 (1,1,4: 1,1,1,4),0,1000,1
VAR Y2,3,,: 2,1,4,0,1,3,0 (1,1,4: 1,1,1,4),0,1000,1
VAR Y12,3,,: 2,1,4,0,1,3,0 (1,1,4: 1,1,1,4),0,1000,1
VAR X1L,3,,: 2,1,4,0,0,103,0,0,0,1
VAR X1U,3,,: 2,1,4,0,0,103,0,0,0,1
VAR X2L,3,,: 2,1,4,0,0,103,0,0,0,1
VAR X2U,3,,: 2,1,4,0,0,103,0,0,0,1
VAR X1,3,,: 2,1,4,0,0,103,0,0,0,1
VAR X2,3,,: 2,1,4,0,0,103,0,0,0,1
VAR ANSY,3,,: 2,1,4,0,0,183,0,0,0,2
VAR ANSY1,3,,: 2,1,4,0,0,183,0,0,0,2
VAR ANSY2,3,,: 2,1,4,0,0,183,0,0,0,2
ENDPROC
PROC SORT_BYRESHAPE,1,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:SWAP=>SWAP
USE NRTYPE 2
VAR ARR,3,,: 2,1,4,0,1,83,0 (1,5,0: 5,3,1,0),0,0,3
ENDPROC
PROC SPRSIN_DP,3,8,0,17,0: 8,0,0,0,0,40200,1,0,40000,0
USE NRUTIL 2,ONLY:ASSERT_EQ=>ASSERT_EQ,ARTH=>ARTH
USE NRTYPE 2
VAR A,3,,: 2,2,5,0,1,3,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,1
VAR THRESH,3,,: 2,2,5,0,0,103,0,0,0,1
VAR SA,3,,: 7,SPRS2_DP,b,0,0,183,0,0,0,2
ENDPROC
PROC INDEX_BYPACK,3,8,0,17,0: 8,0,0,0,0,60a00,1,0,0,0
USE NRUTIL 2,ONLY:ASSERT_EQ=>ASSERT_EQ,ARTH=>ARTH,ARRAY_COPY=>ARRAY_COPY
USE NRTYPE 2
VAR ARR,12,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,1000,1
VAR INDEX,3,,: 1,3,3,0,1,183,0 (1,5,0: 5,3,1,0),0,0,3
VAR PARTIAL,3,,: 1,3,3,0,0,13,0,0,0,1
ENDPROC
PROC SELECT,2,8,0,17,0: 2,1,4,0,0,60281,1,0,0,0
USE NRUTIL 2,ONLY:SWAP=>SWAP,ASSERT=>ASSERT
USE NRTYPE 2
VAR K,3,,: 1,3,3,0,0,103,0,0,0,1
VAR ARR,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,0,3
ENDPROC
PROC SPLIN2,6,8,0,17,0: 2,1,4,0,0,40281,1,0,0,0
USE NRUTIL 2,ONLY:ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR X1A,3,,: 2,1,4,0,1,3,0 (1,5,0: 5,3,1,0),0,1000,1
VAR X2A,3,,: 2,1,4,0,1,3,0 (1,5,0: 5,3,1,0),0,1000,1
VAR YA,3,,: 2,1,4,0,1,103,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,1
VAR Y2A,3,,: 2,1,4,0,1,103,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,1
VAR X1,3,,: 2,1,4,0,0,3,0,0,1000,1
VAR X2,3,,: 2,1,4,0,0,3,0,0,1000,1
ENDPROC
PROC EULSUM,3,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:REALLOCATE=>REALLOCATE,POLY_TERM=>POLY_TERM
USE NRTYPE 2
VAR SUM,3,,: 2,1,4,0,0,183,0,0,0,3
VAR TERM,3,,: 2,1,4,0,0,103,0,0,0,1
VAR JTERM,3,,: 1,3,3,0,0,103,0,0,0,1
ENDPROC
PROC SPRSAX_SP,3,8,0,17,0: 8,0,0,0,0,40200,1,0,40000,0
USE NRUTIL 2,ONLY:SCATTER_ADD=>SCATTER_ADD,ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR SA,3,,: 7,SPRS2_SP,b,0,0,103,0,0,0,1
VAR X,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR B,3,,: 2,1,4,0,1,83,0 (1,5,0: 5,3,1,0),0,1000,2
ENDPROC
PROC SPRSTX_DP,3,8,0,17,0: 8,0,0,0,0,40200,1,0,40000,0
USE NRUTIL 2,ONLY:SCATTER_ADD=>SCATTER_ADD,ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR SA,3,,: 7,SPRS2_DP,b,0,0,103,0,0,0,1
VAR X,3,,: 2,2,5,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR B,3,,: 2,2,5,0,1,83,0 (1,5,0: 5,3,1,0),0,1000,2
ENDPROC
PROC SELECT_INPLACE,2,8,0,17,0: 2,1,4,0,0,40281,1,0,0,0
USE NRTYPE 2
VAR K,3,,: 1,3,3,0,0,3,0,0,1000,1
VAR ARR,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC SPRSAX_DP,3,8,0,17,0: 8,0,0,0,0,40200,1,0,40000,0
USE NRUTIL 2,ONLY:SCATTER_ADD=>SCATTER_ADD,ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR SA,3,,: 7,SPRS2_DP,b,0,0,103,0,0,0,1
VAR X,3,,: 2,2,5,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR B,3,,: 2,2,5,0,1,83,0 (1,5,0: 5,3,1,0),0,1000,2
ENDPROC
PROC INDEXX_I4B,2,8,0,17,0: 8,0,0,0,0,40200,1,0,40000,0
USE NRUTIL 2,ONLY:SWAP=>SWAP,NRERROR=>NRERROR,ASSERT_EQ=>ASSERT_EQ,ARTH=>ARTH
USE NRTYPE 2
VAR IARR,12,,: 1,3,3,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR INDEX,3,,: 1,3,3,0,1,183,0 (1,5,0: 5,3,1,0),0,0,2
ENDPROC
PROC RATINT,5,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:NRERROR=>NRERROR,IMINLOC=>IMINLOC,ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR XA,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR YA,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR X,3,,: 2,1,4,0,0,103,0,0,0,1
VAR Y,3,,: 2,1,4,0,0,183,0,0,0,2
VAR DY,3,,: 2,1,4,0,0,183,0,0,0,2
ENDPROC
PROC CYCLIC,7,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:ASSERT_EQ=>ASSERT_EQ,ASSERT=>ASSERT
USE NRTYPE 2
VAR A,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR B,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR C,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR ALPHA,3,,: 2,1,4,0,0,103,0,0,0,1
VAR BETA,3,,: 2,1,4,0,0,103,0,0,0,1
VAR R,3,,: 2,1,4,0,1,3,0 (1,5,0: 5,3,1,0),0,1000,1
VAR X,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,1000,2
ENDPROC
PROC SPRSTX_SP,3,8,0,17,0: 8,0,0,0,0,40200,1,0,40000,0
USE NRUTIL 2,ONLY:SCATTER_ADD=>SCATTER_ADD,ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR SA,3,,: 7,SPRS2_SP,b,0,0,103,0,0,0,1
VAR X,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR B,3,,: 2,1,4,0,1,83,0 (1,5,0: 5,3,1,0),0,1000,2
ENDPROC
PROC SELECT_BYPACK,2,8,0,17,0: 2,1,4,0,0,40281,1,0,0,0
USE NRUTIL 2,ONLY:SWAP=>SWAP,ASSERT=>ASSERT,ARRAY_COPY=>ARRAY_COPY
USE NRTYPE 2
VAR K,3,,: 1,3,3,0,0,103,0,0,0,1
VAR ARR,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,0,3
ENDPROC
PROC TRIDAG_PAR,5,8,0,17,0: 8,0,0,0,0,60a00,1,0,0,0
USE NRUTIL 2,ONLY:NRERROR=>NRERROR,ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR A,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,1000,1
VAR B,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,1000,1
VAR C,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,1000,1
VAR R,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,1000,1
VAR U,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,1000,2
ENDPROC
PROC MPROVE,5,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR A,3,,: 2,1,4,0,1,3,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,1
VAR ALUD,3,,: 2,1,4,0,1,3,0 (2,5,0: 5,3,1,0,5,3,1,0),0,1000,1
VAR INDX,3,,: 1,3,3,0,1,3,0 (1,5,0: 5,3,1,0),0,1000,1
VAR B,3,,: 2,1,4,0,1,3,0 (1,5,0: 5,3,1,0),0,0,1
VAR X,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,0,3
ENDPROC
PROC GAMMLN_S,1,8,0,17,0: 2,1,4,0,0,40281,1,0,40000,0
USE NRUTIL 2,ONLY:ASSERT=>ASSERT,ARTH=>ARTH
USE NRTYPE 2
VAR XX,3,,: 2,1,4,0,0,103,0,0,0,1
ENDPROC
PROC GAMMLN_V,1,8,0,17,0: 2,1,4,0,1,40281,1 (1,2,0: 1,2,1,2),0,40000,0
USE NRUTIL 2,ONLY:ASSERT=>ASSERT
USE NRTYPE 2
VAR XX,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC POLCOF,2,8,0,17,0: 2,1,4,0,1,40301,1 (1,2,0: 1,2,1,2),0,0,0
USE NRUTIL 2,ONLY:IMINLOC=>IMINLOC,ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR XA,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR YA,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC ATIMES,3,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE XLINBCG_DATA 2
USE NRUTIL 2,ONLY:ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR X,3,,: 2,2,5,0,1,3,0 (1,5,0: 5,3,1,0),0,1000,1
VAR R,3,,: 2,2,5,0,1,3,0 (1,5,0: 5,3,1,0),0,1000,2
VAR ITRNSP,3,,: 1,3,3,0,0,103,0,0,0,1
ENDPROC
PROC SPRSDIAG_SP,2,8,0,17,0: 8,0,0,0,0,40200,1,0,40000,0
USE NRUTIL 2,ONLY:ASSERT_EQ=>ASSERT_EQ,ARRAY_COPY=>ARRAY_COPY
USE NRTYPE 2
VAR SA,3,,: 7,SPRS2_SP,b,0,0,103,0,0,0,1
VAR B,3,,: 2,1,4,0,1,83,0 (1,5,0: 5,3,1,0),0,0,2
ENDPROC
PROC SORT_SHELL,1,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRTYPE 2
VAR ARR,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,0,3
ENDPROC
PROC SNRM,2,8,0,17,0: 2,2,5,0,0,40281,1,0,0,0
USE NRTYPE 2
VAR SX,3,,: 2,2,5,0,1,3,0 (1,5,0: 5,3,1,0),0,0,1
VAR ITOL,3,,: 1,3,3,0,0,103,0,0,0,1
ENDPROC
PROC RAN0_S,1,8,0,17,0: 8,0,0,0,0,40200,1,0,40000,0
USE RAN_STATE 2,ONLY:RANS=>RANS,NRAN0=>NRAN0,KRAN0=>KRAN0,JRAN0=>JRAN0,IRAN0=>IRAN0,RAN_INIT=>RAN_INIT,LENRAN=>LENRAN,AMM=>AMM
USE NRTYPE 2
VAR HARVEST,3,,: 2,1,4,0,0,83,0,0,0,2
ENDPROC
PROC SPLINE,5,8,0,17,0: 8,0,0,0,0,60200,1,0,0,0
USE NRUTIL 2,ONLY:ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR X,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR Y,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR YP1,3,,: 2,1,4,0,0,103,0,0,0,1
VAR YPN,3,,: 2,1,4,0,0,103,0,0,0,1
VAR Y2,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,2
ENDPROC
PROC LUDCMP,3,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:SWAP=>SWAP,OUTERPROD=>OUTERPROD,NRERROR=>NRERROR,IMAXLOC=>IMAXLOC,ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR A,3,,: 2,1,4,0,1,183,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,3
VAR INDX,3,,: 1,3,3,0,1,83,0 (1,5,0: 5,3,1,0),0,0,2
VAR D,3,,: 2,1,4,0,0,183,0,0,0,2
ENDPROC
PROC SPRSDIAG_DP,2,8,0,17,0: 8,0,0,0,0,40200,1,0,40000,0
USE NRUTIL 2,ONLY:ASSERT_EQ=>ASSERT_EQ,ARRAY_COPY=>ARRAY_COPY
USE NRTYPE 2
VAR SA,3,,: 7,SPRS2_DP,b,0,0,103,0,0,0,1
VAR B,3,,: 2,2,5,0,1,83,0 (1,5,0: 5,3,1,0),0,0,2
ENDPROC
PROC RAN0_V,1,8,0,17,0: 8,0,0,0,0,40200,1,0,40000,0
USE RAN_STATE 2,ONLY:RANV=>RANV,NRAN=>NRAN,KRAN=>KRAN,JRAN=>JRAN,IRAN=>IRAN,RAN_INIT=>RAN_INIT,LENRAN=>LENRAN,AMM=>AMM,K4B=>K4B
USE NRTYPE 2
VAR HARVEST,3,,: 2,1,4,0,1,83,0 (1,5,0: 5,3,1,0),0,0,2
ENDPROC
PROC BANBKS,6,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:SWAP=>SWAP,ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR A,3,,: 2,1,4,0,1,103,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,1
VAR M1,3,,: 1,3,3,0,0,103,0,0,1000,1
VAR M2,3,,: 1,3,3,0,0,103,0,0,0,1
VAR AL,3,,: 2,1,4,0,1,103,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,1
VAR INDX,3,,: 1,3,3,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR B,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,0,3
ENDPROC
PROC RAN1_S,1,8,0,17,0: 8,0,0,0,0,40200,1,0,40000,0
USE RAN_STATE 2,ONLY:RANS=>RANS,MRAN0=>MRAN0,NRAN0=>NRAN0,KRAN0=>KRAN0,JRAN0=>JRAN0,IRAN0=>IRAN0,RAN_INIT=>RAN_INIT,LENRAN=>LENRAN,AMM=>AMM,K4B=>K4B
USE NRTYPE 2
VAR HARVEST,3,,: 2,1,4,0,0,83,0,0,0,2
ENDPROC
PROC INDEXX_SP,2,8,0,17,0: 8,0,0,0,0,40200,1,0,40000,0
USE NRUTIL 2,ONLY:SWAP=>SWAP,NRERROR=>NRERROR,ASSERT_EQ=>ASSERT_EQ,ARTH=>ARTH
USE NRTYPE 2
VAR ARR,12,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR INDEX,3,,: 1,3,3,0,1,183,0 (1,5,0: 5,3,1,0),0,0,2
ENDPROC
PROC CHOLDC,2,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:NRERROR=>NRERROR,ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR A,3,,: 2,1,4,0,1,183,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,3
VAR P,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,0,2
ENDPROC
PROC RAN1_V,1,8,0,17,0: 8,0,0,0,0,40200,1,0,40000,0
USE RAN_STATE 2,ONLY:RANV=>RANV,MRAN=>MRAN,NRAN=>NRAN,KRAN=>KRAN,JRAN=>JRAN,IRAN=>IRAN,LENRAN=>LENRAN,AMM=>AMM,K4B=>K4B
USE NRTYPE 2
VAR HARVEST,3,,: 2,1,4,0,1,83,0 (1,5,0: 5,3,1,0),0,0,2
ENDPROC
PROC SORT,1,8,0,17,0: 8,0,0,0,0,60200,1,0,0,0
USE NRUTIL 2,ONLY:NRERROR=>NRERROR,SWAP=>SWAP
USE NRTYPE 2
VAR ARR,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,0,3
ENDPROC
PROC RAN2_S,1,8,0,17,0: 8,0,0,0,0,40200,1,0,40000,0
USE RAN_STATE 2,ONLY:RANS=>RANS,MRAN0=>MRAN0,NRAN0=>NRAN0,KRAN0=>KRAN0,JRAN0=>JRAN0,IRAN0=>IRAN0,RAN_INIT=>RAN_INIT,LENRAN=>LENRAN,AMM=>AMM,K4B=>K4B
USE NRTYPE 2
VAR HARVEST,3,,: 2,1,4,0,0,83,0,0,0,2
ENDPROC
PROC SPLINT,4,8,0,17,0: 2,1,4,0,0,60281,1,0,0,0
USE NRUTIL 2,ONLY:NRERROR=>NRERROR,ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR XA,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,1000,1
VAR YA,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR Y2A,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR X,3,,: 2,1,4,0,0,103,0,0,1000,1
ENDPROC
PROC RAN2_V,1,8,0,17,0: 8,0,0,0,0,40200,1,0,40000,0
USE RAN_STATE 2,ONLY:RANV=>RANV,MRAN=>MRAN,NRAN=>NRAN,KRAN=>KRAN,JRAN=>JRAN,IRAN=>IRAN,RAN_INIT=>RAN_INIT,LENRAN=>LENRAN,AMM=>AMM,K4B=>K4B
USE NRTYPE 2
VAR HARVEST,3,,: 2,1,4,0,1,83,0 (1,5,0: 5,3,1,0),0,0,2
ENDPROC
PROC CHOLSL,4,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR A,3,,: 2,1,4,0,1,103,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,1
VAR P,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR B,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR X,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,0,3
ENDPROC
PROC GASDEV_S,1,8,0,17,0: 8,0,0,0,0,40200,1,0,40000,0
USE NRTYPE 2
VAR HARVEST,3,,: 2,1,4,0,0,83,0,0,0,2
ENDPROC
PROC BANMUL,5,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:ARTH=>ARTH,ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR A,3,,: 2,1,4,0,1,103,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,1
VAR M1,3,,: 1,3,3,0,0,103,0,0,0,1
VAR M2,3,,: 1,3,3,0,0,103,0,0,0,1
VAR X,3,,: 2,1,4,0,1,3,0 (1,5,0: 5,3,1,0),0,0,1
VAR B,3,,: 2,1,4,0,1,83,0 (1,5,0: 5,3,1,0),0,0,2
ENDPROC
PROC QRDCMP,4,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:VABS=>VABS,OUTERPROD=>OUTERPROD,ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR A,3,,: 2,1,4,0,1,183,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,3
VAR C,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,0,2
VAR D,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,0,2
VAR SING,3,,: 4,3,7,0,0,83,0,0,0,2
ENDPROC
PROC GASDEV_V,1,8,0,17,0: 8,0,0,0,0,40200,1,0,40000,0
USE NRUTIL 2,ONLY:ARRAY_COPY=>ARRAY_COPY
USE NRTYPE 2
VAR HARVEST,3,,: 2,1,4,0,1,83,0 (1,5,0: 5,3,1,0),0,0,2
ENDPROC
PROC ECLAZZ,2,8,0,17,0: 1,3,3,0,1,40381,1 (1,2,0: 1,2,1,2),0,0,0
USE NRUTIL 2,ONLY:ARTH=>ARTH
USE NRTYPE 2
PROC EQUIV,2,10,0,17,0: 4,3,7,0,0,60003,0,0,0,0
USE NRTYPE 2
VAR I,3,,: 1,3,3,0,0,3,0,0,0,1
VAR J,3,,: 1,3,3,0,0,3,0,0,0,1
ENDPROC
VAR N,3,,: 1,3,3,0,0,103,0,0,1000,1
ENDPROC
PROC LOCATE,2,8,0,17,0: 1,3,3,0,0,40281,1,0,0,0
USE NRTYPE 2
VAR XX,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR X,3,,: 2,1,4,0,0,103,0,0,0,1
ENDPROC
PROC BCUCOF,7,8,0,17,0: 8,0,0,0,0,60200,1,0,0,0
USE NRTYPE 2
VAR Y,3,,: 2,1,4,0,1,103,0 (1,1,4: 1,1,1,4),0,0,1
VAR Y1,3,,: 2,1,4,0,1,103,0 (1,1,4: 1,1,1,4),0,0,1
VAR Y2,3,,: 2,1,4,0,1,103,0 (1,1,4: 1,1,1,4),0,0,1
VAR Y12,3,,: 2,1,4,0,1,103,0 (1,1,4: 1,1,1,4),0,0,1
VAR D1,3,,: 2,1,4,0,0,103,0,0,0,1
VAR D2,3,,: 2,1,4,0,0,103,0,0,0,1
VAR C,3,,: 2,1,4,0,1,83,0 (2,1,16: 1,1,1,4,1,1,1,4),0,0,2
ENDPROC
PROC SPRSTP,1,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRTYPE 2
VAR SA,3,,: 7,SPRS2_SP,b,0,0,83,0,0,0,3
ENDPROC
PROC SVBKSB_SP,5,8,0,17,0: 8,0,0,0,0,40200,1,0,40000,0
USE NRUTIL 2,ONLY:ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR U,3,,: 2,1,4,0,1,3,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,1
VAR W,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR V,3,,: 2,1,4,0,1,3,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,1
VAR B,3,,: 2,1,4,0,1,3,0 (1,5,0: 5,3,1,0),0,0,1
VAR X,3,,: 2,1,4,0,1,83,0 (1,5,0: 5,3,1,0),0,0,2
ENDPROC
PROC LUBKSB,3,8,0,17,0: 8,0,0,0,0,60200,1,0,0,0
USE NRUTIL 2,ONLY:ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR A,3,,: 2,1,4,0,1,103,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,1
VAR INDX,3,,: 1,3,3,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR B,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,0,3
ENDPROC
PROC DDPOLY,3,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:POLY_TERM=>POLY_TERM,CUMPROD=>CUMPROD,ARTH=>ARTH
USE NRTYPE 2
VAR C,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR X,3,,: 2,1,4,0,0,3,0,0,1000,1
VAR PD,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,0,2
ENDPROC
PROC SORT2,2,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR ARR,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,1000,3
VAR SLAVE,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,0,3
ENDPROC
PROC QRSOLV,4,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR A,3,,: 2,1,4,0,1,103,0 (2,5,0: 5,3,1,0,5,3,1,0),0,1000,1
VAR C,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR D,3,,: 2,1,4,0,1,3,0 (1,5,0: 5,3,1,0),0,1000,1
VAR B,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,1000,3
ENDPROC
PROC SVBKSB_DP,5,8,0,17,0: 8,0,0,0,0,40200,1,0,40000,0
USE NRUTIL 2,ONLY:ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR U,3,,: 2,2,5,0,1,3,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,1
VAR W,3,,: 2,2,5,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR V,3,,: 2,2,5,0,1,3,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,1
VAR B,3,,: 2,2,5,0,1,3,0 (1,5,0: 5,3,1,0),0,0,1
VAR X,3,,: 2,2,5,0,1,83,0 (1,5,0: 5,3,1,0),0,0,2
ENDPROC
PROC SPLIE2,4,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR X1A,3,,: 2,1,4,0,1,3,0 (1,5,0: 5,3,1,0),0,0,1
VAR X2A,3,,: 2,1,4,0,1,3,0 (1,5,0: 5,3,1,0),0,1000,1
VAR YA,3,,: 2,1,4,0,1,103,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,1
VAR Y2A,3,,: 2,1,4,0,1,103,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,2
ENDPROC
PROC RSOLV,3,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR A,3,,: 2,1,4,0,1,103,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,1
VAR D,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR B,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,0,3
ENDPROC
PROC SVDCMP_SP,3,8,0,17,0: 8,0,0,0,0,40200,1,0,40000,0
USE NRUTIL 2,ONLY:OUTERPROD=>OUTERPROD,NRERROR=>NRERROR,ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR A,3,,: 2,1,4,0,1,183,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,3
VAR W,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,0,2
VAR V,3,,: 2,1,4,0,1,183,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,2
ENDPROC
PROC TRAPZD,5,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:ARTH=>ARTH
USE NRTYPE 2
PROC FUNC,1,10,0,17,0: 2,1,4,0,1,60003,0 (1,2,0: 1,2,1,2),0,0,0
USE NRTYPE 2
VAR X,3,,: 2,1,4,0,1,3,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
VAR A,3,,: 2,1,4,0,0,103,0,0,0,1
VAR B,3,,: 2,1,4,0,0,103,0,0,0,1
VAR S,3,,: 2,1,4,0,0,183,0,0,0,3
VAR N,3,,: 1,3,3,0,0,103,0,0,0,1
ENDPROC
PROC VANDER,2,8,0,17,0: 2,2,5,0,1,40381,1 (1,2,0: 1,2,1,2),0,0,0
USE NRUTIL 2,ONLY:OUTERDIFF=>OUTERDIFF,ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR X,3,,: 2,2,5,0,1,103,0 (1,5,0: 5,3,1,0),0,1000,1
VAR Q,3,,: 2,2,5,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC SORT3,3,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR ARR,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,1000,3
VAR SLAVE1,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,0,3
VAR SLAVE2,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,0,3
ENDPROC
PROC HUNT,3,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRTYPE 2
VAR XX,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR X,3,,: 2,1,4,0,0,103,0,0,0,1
VAR JLO,3,,: 1,3,3,0,0,183,0,0,0,3
ENDPROC
PROC RAN,1,8,0,17,0: 2,1,4,0,0,40281,1,0,40000,0
VAR IDUM,3,,: 1,3,3,0,0,183,0,0,0,3
ENDPROC
PROC SVDCMP_DP,3,8,0,17,0: 8,0,0,0,0,40200,1,0,40000,0
USE NRUTIL 2,ONLY:OUTERPROD=>OUTERPROD,NRERROR=>NRERROR,ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR A,3,,: 2,2,5,0,1,183,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,3
VAR W,3,,: 2,2,5,0,1,183,0 (1,5,0: 5,3,1,0),0,0,2
VAR V,3,,: 2,2,5,0,1,183,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,2
ENDPROC
PROC TRIDAG_SER,5,8,0,17,0: 8,0,0,0,0,60200,1,0,0,0
USE NRUTIL 2,ONLY:NRERROR=>NRERROR,ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR A,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR B,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR C,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR R,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR U,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,0,2
ENDPROC
PROC BANDEC,6,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:ARTH=>ARTH,SWAP=>SWAP,IMAXLOC=>IMAXLOC,ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR A,3,,: 2,1,4,0,1,183,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,3
VAR M1,3,,: 1,3,3,0,0,103,0,0,1000,1
VAR M2,3,,: 1,3,3,0,0,103,0,0,0,1
VAR AL,3,,: 2,1,4,0,1,83,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,2
VAR INDX,3,,: 1,3,3,0,1,83,0 (1,5,0: 5,3,1,0),0,0,2
VAR D,3,,: 2,1,4,0,0,183,0,0,0,2
ENDPROC
PROC SORT_PICK,1,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRTYPE 2
VAR ARR,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,0,3
ENDPROC
PROC SORT_BYPACK,1,8,0,17,0: 8,0,0,0,0,60a00,1,0,0,0
USE NRUTIL 2,ONLY:SWAP=>SWAP,ARRAY_COPY=>ARRAY_COPY
USE NRTYPE 2
VAR ARR,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,0,3
ENDPROC
PROC RANK,1,8,0,17,0: 1,3,3,0,1,40281,1 (1,2,0: 1,2,1,2),0,0,0
USE NRUTIL 2,ONLY:ARTH=>ARTH
USE NRTYPE 2
VAR INDEX,3,,: 1,3,3,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC POLCOE,2,8,0,17,0: 2,1,4,0,1,40381,1 (1,2,0: 1,2,1,2),0,0,0
USE NRUTIL 2,ONLY:OUTERDIFF=>OUTERDIFF,ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR X,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,1000,1
VAR Y,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC LINBCG,7,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:NRERROR=>NRERROR,ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR B,3,,: 2,2,5,0,1,103,0 (1,5,0: 5,3,1,0),0,1000,1
VAR X,3,,: 2,2,5,0,1,183,0 (1,5,0: 5,3,1,0),0,1000,3
VAR ITOL,3,,: 1,3,3,0,0,103,0,0,1000,1
VAR TOL,3,,: 2,2,5,0,0,103,0,0,0,1
VAR ITMAX,3,,: 1,3,3,0,0,103,0,0,0,1
VAR ITER,3,,: 1,3,3,0,0,183,0,0,0,2
VAR ERR,3,,: 2,2,5,0,0,183,0,0,0,2
ENDPROC
PROC SORT_HEAP,1,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:SWAP=>SWAP
USE NRTYPE 2
VAR ARR,12,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,0,3
ENDPROC
PROC TOEPLZ,2,8,0,17,0: 2,1,4,0,1,40381,1 (1,2,0: 1,2,1,2),0,0,0
USE NRUTIL 2,ONLY:NRERROR=>NRERROR,ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR R,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR Y,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC PYTHAG_SP,2,8,0,17,0: 2,1,4,0,0,40281,1,0,40000,0
USE NRTYPE 2
VAR A,3,,: 2,1,4,0,0,3,0,0,0,1
VAR B,3,,: 2,1,4,0,0,3,0,0,0,1
ENDPROC
PROC ROTATE,5,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR R,3,,: 2,1,4,0,1,1183,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,3
VAR QT,3,,: 2,1,4,0,1,1183,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,3
VAR I,3,,: 1,3,3,0,0,103,0,0,0,1
VAR A,3,,: 2,1,4,0,0,103,0,0,0,1
VAR B,3,,: 2,1,4,0,0,103,0,0,0,1
ENDPROC
PROC QRUPDT,4,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:IFIRSTLOC=>IFIRSTLOC,ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR R,3,,: 2,1,4,0,1,183,0 (2,5,0: 5,3,1,0,5,3,1,0),0,1000,3
VAR QT,3,,: 2,1,4,0,1,3,0 (2,5,0: 5,3,1,0,5,3,1,0),0,1000,3
VAR U,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,0,3
VAR V,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
ENDPROC
PROC ECLASS,3,8,0,17,0: 1,3,3,0,1,40381,1 (1,2,0: 1,2,1,2),0,0,0
USE NRUTIL 2,ONLY:ASSERT_EQ=>ASSERT_EQ,ARTH=>ARTH
USE NRTYPE 2
VAR LISTA,3,,: 1,3,3,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR LISTB,3,,: 1,3,3,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR N,3,,: 1,3,3,0,0,103,0,0,1000,1
ENDPROC
PROC PYTHAG_DP,2,8,0,17,0: 2,2,5,0,0,40281,1,0,40000,0
USE NRTYPE 2
VAR A,3,,: 2,2,5,0,0,3,0,0,0,1
VAR B,3,,: 2,2,5,0,0,3,0,0,0,1
ENDPROC
PROC SELECT_HEAP,2,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:SWAP=>SWAP,NRERROR=>NRERROR
USE NRTYPE 2
VAR ARR,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR HEAP,3,,: 2,1,4,0,1,183,0 (1,5,0: 5,3,1,0),0,1000,2
ENDPROC
PROC POLINT,5,8,0,17,0: 8,0,0,0,0,60200,1,0,0,0
USE NRUTIL 2,ONLY:NRERROR=>NRERROR,IMINLOC=>IMINLOC,ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR XA,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR YA,3,,: 2,1,4,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR X,3,,: 2,1,4,0,0,103,0,0,0,1
VAR Y,3,,: 2,1,4,0,0,183,0,0,0,2
VAR DY,3,,: 2,1,4,0,0,183,0,0,0,2
ENDPROC
PROC ASOLVE,3,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE XLINBCG_DATA 2
USE NRUTIL 2,ONLY:NRERROR=>NRERROR,ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR B,3,,: 2,2,5,0,1,103,0 (1,5,0: 5,3,1,0),0,0,1
VAR X,3,,: 2,2,5,0,1,183,0 (1,5,0: 5,3,1,0),0,1000,2
VAR ITRNSP,3,,: 1,3,3,0,0,3,0,0,0,1
ENDPROC
PROC SORT_RADIX,1,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:NRERROR=>NRERROR,ARRAY_COPY=>ARRAY_COPY
USE NRTYPE 2
VAR ARR,3,,: 2,1,4,0,1,83,0 (1,5,0: 5,3,1,0),0,0,3
ENDPROC
PROC POLIN2,7,8,0,17,0: 8,0,0,0,0,40200,1,0,0,0
USE NRUTIL 2,ONLY:ASSERT_EQ=>ASSERT_EQ
USE NRTYPE 2
VAR X1A,3,,: 2,1,4,0,1,3,0 (1,5,0: 5,3,1,0),0,1000,1
VAR X2A,3,,: 2,1,4,0,1,3,0 (1,5,0: 5,3,1,0),0,1000,1
VAR YA,3,,: 2,1,4,0,1,103,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,1
VAR X1,3,,: 2,1,4,0,0,3,0,0,1000,1
VAR X2,3,,: 2,1,4,0,0,3,0,0,1000,1
VAR Y,3,,: 2,1,4,0,0,3,0,0,1000,2
VAR DY,3,,: 2,1,4,0,0,3,0,0,1000,2
ENDPROC
PROC SPRSIN_SP,3,8,0,17,0: 8,0,0,0,0,40200,1,0,40000,0
USE NRUTIL 2,ONLY:ASSERT_EQ=>ASSERT_EQ,ARTH=>ARTH
USE NRTYPE 2
VAR A,3,,: 2,1,4,0,1,3,0 (2,5,0: 5,3,1,0,5,3,1,0),0,0,1
VAR THRESH,3,,: 2,1,4,0,0,103,0,0,0,1
VAR SA,3,,: 7,SPRS2_SP,b,0,0,183,0,0,0,2
ENDPROC
END