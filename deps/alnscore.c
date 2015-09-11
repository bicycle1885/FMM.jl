#include "stdio.h"
#include "stdlib.h"
#include "stdint.h"
#include "emmintrin.h"

// 16-byte-aligned memory allocation (copied from Julia src/gc.c)
static inline void *malloc_a16(size_t sz)
{
    void *ptr;
    if (posix_memalign(&ptr, 16, sz))
        return NULL;
    return ptr;
}

// alignment score type
typedef int16_t score_t;
#define SCORE_MIN INT16_MIN

// DNA encoding:
//   A: 0x00
//   C: 0x01
//   G: 0x02
//   T: 0x03
//   N: 0x04

const int n_nucs = 5;

typedef struct dnaseq_s
{
    uint8_t* seq;
    size_t len;
} dnaseq_t;

size_t seqlen(dnaseq_t seq)
{
    return seq.len;
}

typedef struct seqprof_s
{
    dnaseq_t seq;
    score_t* profile;
} seqprof_t;

typedef struct params_s
{
    score_t matching_score;
    score_t mismatching_score;
    score_t gap_open_penalty;
    score_t gap_ext_penalty;
} params_t;

const params_t default_params = {
    .matching_score     =  0,
    .mismatching_score  = -6,
    .gap_open_penalty   =  5,
    .gap_ext_penalty    =  3
};

// align 8 sequences in parallel (= 128-bit register / 16-bit score type)
const int N_PAR = 8;

int alignment_score(const params_t* params, const dnaseq_t query, const dnaseq_t* seqs, score_t* scores)
{
    // read params
    if (params == NULL)
        params = &default_params;
    score_t matching_score = params->matching_score;
    score_t mismatching_score = params->mismatching_score;
    score_t gap_open_penalty = params->gap_open_penalty;
    score_t gap_ext_penalty = params->gap_ext_penalty;

    // check length
    size_t len = seqlen(seqs[0]);
    for (size_t n = 1; n < N_PAR; n++)
        if (seqlen(seqs[n]) != len)
            return 1;

    // allocate working space
    size_t qlen = seqlen(query);
    size_t column_bytes = sizeof(__m128i) * (qlen + 1);
    size_t profile_bytes = sizeof(__m128i) * n_nucs * len;
    __m128i* columnH = malloc_a16(column_bytes * 2 + profile_bytes);
    __m128i* columnE = columnH + (qlen + 1);
    __m128i* profile = columnE + (qlen + 1);
    if (columnH == NULL)
        return 1;

    // initialize profile
    for (size_t j = 0; j < len; j++) {
        for (uint8_t nt = 0; nt < n_nucs; nt++) {
            score_t P[N_PAR];
            for (size_t n = 0; n < N_PAR; n++)
                P[n] = nt == seqs[n].seq[j] ? matching_score : mismatching_score;
            // set profile
            _mm_store_si128(
                profile + n_nucs * j + nt,
                _mm_set_epi16(P[7], P[6], P[5], P[4], P[3], P[2], P[1], P[0])
            );
        }
    }

    // initialize columnH and columnE, and set other vectors
    columnH[0] = _mm_set1_epi16(0);
    columnE[0] = _mm_set1_epi16(SCORE_MIN);
    for (size_t i = 1; i <= qlen; i++) {
        columnH[i] = _mm_set1_epi16(-gap_open_penalty - gap_ext_penalty * i);
        columnE[i] = _mm_set1_epi16(SCORE_MIN);
    }
    __m128i gap_open = _mm_set1_epi16(gap_open_penalty + gap_ext_penalty);
    __m128i gap_ext  = _mm_set1_epi16(gap_ext_penalty);
    __m128i maxH = _mm_set1_epi16(SCORE_MIN);

    // run dynamic programming
    uint8_t* qseq = query.seq;
    for (size_t j = 1; j <= len; j++) {
        __m128i H, E, F, diagH;
        F = _mm_set1_epi16(SCORE_MIN);
        diagH = columnH[0];
        columnH[0] = _mm_set1_epi16(-gap_open_penalty - gap_ext_penalty * j);
        __m128i* prof = profile + (j - 1) * n_nucs;
        for (size_t i = 1; i <= qlen; i++) {
            E = _mm_max_epi16(
                _mm_subs_epi16(columnE[i], gap_ext),
                _mm_subs_epi16(columnH[i], gap_open)
            );
            F = _mm_max_epi16(
                _mm_subs_epi16(F, gap_ext),
                _mm_subs_epi16(columnH[i-1], gap_open)
            );
            __m128i P = *(prof + qseq[i-1]);
            H = _mm_max_epi16(E, F);
            H = _mm_max_epi16(H, _mm_adds_epi16(diagH, P));
            // update
            diagH = columnH[i];
            columnE[i] = E;
            columnH[i] = H;
        }
        // update maximum score
        maxH = _mm_max_epi16(maxH, columnH[qlen]);
    }

    // save best scores
    _mm_store_si128((__m128i*)scores, maxH);

    // free working space
    free(columnH);

    return 0;
}
