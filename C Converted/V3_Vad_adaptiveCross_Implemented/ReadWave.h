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
void wavread(char *file_name, int16_t **samples, uint32_t *samplelength);
void wavwrite(char *file_name, int16_t *samples);
#endif /* ReadWave_h */
