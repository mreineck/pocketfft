#ifndef POCKETFFT_H
#define POCKETFFT_H

#include <stdlib.h>

struct cfft_plan_i;
typedef struct cfft_plan_i * cfft_plan;
cfft_plan make_cfft_plan (size_t length);
void destroy_cfft_plan (cfft_plan plan);
size_t cfft_worksize (cfft_plan plan);
void cfft_backward(cfft_plan plan, double c[]);
void cfft_forward(cfft_plan plan, double c[]);
//void cfft_backward_with_work(cfft_plan plan, double c[], double w[]);
//void cfft_forward_with_work(cfft_plan plan, double c[], double w[]);

struct rfft_plan_i;
typedef struct rfft_plan_i * rfft_plan;
rfft_plan make_rfft_plan (size_t length);
void destroy_rfft_plan (rfft_plan plan);
size_t rfft_worksize (cfft_plan plan);
void rfft_backward(rfft_plan plan, double c[]);
void rfft_forward(rfft_plan plan, double c[]);
//void rfft_backward_with_work(rfft_plan plan, double c[], double w[]);
//void rfft_forward_with_work(rfft_plan plan, double c[], double w[]);

void cfft_backward_noplan(double c[], size_t length);
void cfft_forward_noplan(double c[], size_t length);
void rfft_backward_noplan(double c[], size_t length);
void rfft_forward_noplan(double c[], size_t length);

#endif
