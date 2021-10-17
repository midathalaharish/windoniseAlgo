//
//  ReadWave.c
//  ShivaKumarChikkala
//
//  Created by Midathala Harish on 05.10.21.
//

#include "ReadWave.h"

WavHeader *header;

void wavread(char *file_name, int16_t **samples, uint32_t *samplelength)
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

    if (header->sample_rate != 48000)
        printf("ERROR: incorrect sampling rate\n");
    
    *samples = (int16_t *)malloc(header->datachunk_size);
    printf("The wave information \n");
    printf("num_channels = %d \n",header->num_channels);
    printf("Sampling rate= %d \n",header->sample_rate);
    printf("The sample length == %d \n",header->datachunk_size);
    printf("The byte_rate == %d \n",header->byte_rate);
    printf("The block_align == %d \n",header->block_align);
    printf("The bits per sample == %d \n",header->bps);

    *samplelength = header->datachunk_size;
    if (!*samples)
        errx(1, "Error allocating memory");

    if (read(fd, *samples, header->datachunk_size) < header->datachunk_size)
        errx(1, "File broken: samples");

    free(header);
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
