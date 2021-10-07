//
//  ReadWave.c
//  ShivaKumarChikkala
//
//  Created by Midathala Harish on 05.10.21.
//

#include "ReadWave.h"
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

WavHeader *header;

void wavread(char *file_name, int16_t **samples)
{
    int fd;

    if (!file_name)
        errx(1, "Filename not specified");
    
    if ((fd = open(file_name, O_RDONLY)) < 0)
        errx(1, "Error opening file %s", file_name);

    if (!header)
        header = (WavHeader*)malloc(sizeof(WavHeader));

    if (read(fd, header, sizeof(WavHeader)) < sizeof(WavHeader))
        errx(1, "File broken: header");

    if (strncmp(header->chunk_id, "RIFF", 4) ||
        strncmp(header->format, "WAVE", 4))
        errx(1, "Not a wav file");

    if (header->audio_format != 1)
        errx(1, "Only PCM encoding supported");

//    if (*samples) free(*samples);
    *samples = (int16_t*)malloc(header->datachunk_size);
    printf("The sample length == %d \n",header->datachunk_size);
    if (!*samples)
        errx(1, "Error allocating memory");

    if (read(fd, *samples, header->datachunk_size) < header->datachunk_size)
        errx(1, "File broken: samples");

    close(fd);
}

void wavwrite(char *file_name, int16_t *samples)
{
    int fd;

    if (!file_name)
        errx(1, "Filename not specified");

    if (!samples)
        errx(1, "Samples buffer not specified");

    if ((fd = creat(file_name, 0666)) < 1)
        errx(1, "Error creating file");

    if (write(fd, header, sizeof(WavHeader)) < sizeof(WavHeader))
        errx(1, "Error writing header");

    if (write(fd, samples, header->datachunk_size) < header->datachunk_size)
        errx(1, "Error writing samples");

    close(fd);
}


int main(int argc, const char * argv[]) {

    BQF BQF_local;

    char *infileName = "/Users/harry/Desktop/PyLearn/Dev/WindNoiseSupress/ShivaKumarChikkala/ShivaKumarChikkala/OTOG2_Speech_silent_wind3_wind5_FanB_mic7.wav";
    int16_t *samples;
    BQF_local.type = bq_type_lowpass;
    BQF_local.a0 = 1.0;
    BQF_local.a1 = BQF_local.a2 = BQF_local.b1 = BQF_local.b2 = 0.0;
    BQF_local.Fc = 0.10;
    BQF_local.Q = 0.707;
    BQF_local.peakGain = 0.0;
    BQF_local.z1 = BQF_local.z2 = 0.0;
    calcBiquad(&BQF_local);
    GetBEQ(&BQF_local);
    wavread(infileName, &samples);
    free(samples);
    return 0;
}

