#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "init.h"

#define NX 600
#define NY 520
#define DX 50.0
#define DY 50.0
#define DT 0.05
#define TOTAL_STEPS 144000
#define SAVE_EACH 300

#define G 9.81
#define MANNING_N 0.035
#define H_UPSTREAM 130.5
#define H_DOWNSTREAM 68.0
#define DAM_X_POS 0
#define DAM_HEIGHT 73.0
#define SLOPE 0.00015

#define BREACH_WIDTH 200.0
#define BRACH_START_TIME 5.0
#define BREACH_FORMATION_TIME 60.0

#define USE_REAL_TOPOGRAPHY 0