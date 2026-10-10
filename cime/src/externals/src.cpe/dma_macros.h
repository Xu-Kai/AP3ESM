// Author: Xiaohui Duan
// A set of Macros with DMA and athread_get/put implementation.
// #define DMA_FAST to use DMA implementation, by default, athread
// implementation is used.
// PLEASE DO NOT USE DMA_FAST at begining.
// DMA_FAST can be tried when athread_get/put implementation can
// produce excepted results.

#ifndef DMA_MACROS_H
#define DMA_MACROS_H

#include <crts.h>
// Pretending to be initializing, actually, it is defining reply, count.
// No descriptor required, setting stride, bsize and mask are done by
// setting variables.
// REMINDER: please clear masks, strides and bsizes after using strided or
// rank mode, and, currently it is not so compatible with DMA_FAST since
// all modes share a same set of masks.
#define dma_init()                              \
  volatile long reply_shadow = 0;                \
  long count_shadow = 0;                         \
  long pe_get_stride_shadow = 0;                \
  long pe_put_stride_shadow = 0;                \
  long row_get_stride_shadow = 0;               \
  long row_put_stride_shadow = 0;               \
  long rank_get_stride_shadow = 0;              \
  long rank_put_stride_shadow = 0;              \
  long brow_get_stride_shadow = 0;              \
  long bcast_get_stride_shadow = 0;             \
                                                \
  long pe_get_bsize_shadow = 0;                 \
  long pe_put_bsize_shadow = 0;                 \
  long row_get_bsize_shadow = 0;                \
  long row_put_bsize_shadow = 0;                \
  long rank_get_bsize_shadow = 0;               \
  long rank_put_bsize_shadow = 0;               \
  long brow_get_bsize_shadow = 0;               \
  long bcast_get_bsize_shadow = 0;              \
                                                \
  int bcast_get_mask_shadow = 0xff;             \
  int brow_get_mask_shadow = 0xff;              \

// The comments for the following macros can be found in DMA_FAST mode.
#define dma_set_stride(mode, stride, bsize) {   \
    mode ## _stride_shadow = (long)(stride);    \
    mode ## _bsize_shadow = (long)(bsize);      \
  }

// #define pe_get(mem, ldm, size) {                                        \
//     athread_get(PE_MODE, (mem), (ldm), (size), (void*)&reply_shadow     \
//         , 0, pe_get_stride_shadow, pe_get_bsize_shadow);                \
//     count_shadow ++;                                                    \
//   }

// #define pe_put(mem, ldm, size) {                                        \
//     athread_put(PE_MODE, (ldm), (mem), (size), (void*)&reply_shadow     \
//         , pe_put_stride_shadow, pe_put_bsize_shadow);                   \
//     count_shadow ++;                                                    \
//   }

#define pe_get(mem, ldm, size) {                                        \
    CRTS_dma_iget((ldm), (mem), (size), (void*)&reply_shadow);          \
    count_shadow ++;                                                    \
  }

#define pe_put(mem, ldm, size) {                                        \
    CRTS_dma_iput((mem), (ldm), (size), (void*)&reply_shadow);          \
    count_shadow ++;                                                    \
  }

#define dma_reset_stride(mode) dma_set_stride(mode, 0, 0)

// Sync the reply and count. "memb" is always a good choice to avoid
// god's punishment.
// #define dma_syn() {                                     \
//     while (reply_shadow != count_shadow) {              \
//     };                                                  \
//     asm volatile("memb\n\t");                           \
//   }

#define dma_syn() {                                     \
    CRTS_dma_wait_value(&reply_shadow, count_shadow);   \
  }                                                     \

#endif