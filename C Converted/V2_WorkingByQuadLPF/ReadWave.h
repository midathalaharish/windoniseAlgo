//
//  ReadWave.h
//  ShivaKumarChikkala
//
//  Created by Midathala Harish on 05.10.21.
//

#ifndef ReadWave_h
#define ReadWave_h

#include <stdio.h>
#include "ByquadFilter.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <err.h>
#include <unistd.h>
#define FTAP_LEN 3
#if 1
//const double lpf_hpf_ba[][FTAP_LEN] = {{0.2534,0.2534},
//                                       {1.0000,-0.4931},
//                                       {0.7466,-0.7466},
//                                       {1.0000,-0.4931}};
const double lpf_hpf_ba[][FTAP_LEN] = {{0.072231,0.144462,0.072231},
{1.0000,-1.1092,0.3982},
{0.6268,-1.2537,0.6268},
{1.0000,-1.1092,0.3982}};
#else
//const double lpf_hpf_ba[][FTAP_LEN] = {{5.5418e-03,2.2167e-02,3.3251e-02,2.2167e-02,5.5418e-03},
//                                {1.0000,-2.3024,2.2091,-0.9922,0.1742},
//                                {0.4174,-1.6695,2.5042,-1.6695,0.4174},
//                                {1.0000,-2.3024,2.2091,-0.9922,0.1742}};
#endif

typedef struct {
    char     chunk_id[4];
    uint32_t chunk_size;
    char     format[4];

    char     fmtchunk_id[4];
    uint32_t fmtchunk_size;
    uint16_t audio_format;
    uint16_t num_channels;
    uint32_t sample_rate;
    uint32_t byte_rate;
    uint16_t block_align;
    uint16_t bps;

    char     datachunk_id[4];
    uint32_t datachunk_size;
}WavHeader;
#endif /* ReadWave_h */
