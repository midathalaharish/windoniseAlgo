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

//    if (*samples) free(*samples);
    *samples = (int16_t *)malloc(header->datachunk_size);
    printf("The sample length == %d \n",header->datachunk_size);
    *samplelength = header->datachunk_size;
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


int main(int argc, const char * argv[])
{
    /*Var*/
    BQF BQF_local;
    uint32_t i=0;
    uint32_t samplelength=0;
    char *infileName = "/Users/harry/Desktop/PyLearn/Dev/WindNoiseSupress/ShivaKumarChikkala/ShivaKumarChikkala/chirp_20_3sec.wav";
    char *outfileName = "/Users/harry/Desktop/PyLearn/Dev/WindNoiseSupress/ShivaKumarChikkala/ShivaKumarChikkala/out.wav";
    int16_t *samples,*out;
    
    /*LoadFilter PArams*/
    
    BQF_local.type = 0; //Check Enum
    
    for(i=0;i<FTAP_LEN;i++)
    {
        BQF_local.lpf_b[i]= lpf_hpf_ba[0][i];
        BQF_local.lpf_a[i]= lpf_hpf_ba[1][i];
        BQF_local.hpf_b[i]= lpf_hpf_ba[2][i];
        BQF_local.hpf_a[i]= lpf_hpf_ba[3][i];
    }
    BQF_local.Fc = 0.4;
    BQF_local.Q = 0.707;
    BQF_local.peakGain = 30;
    BQF_local.lz1 = BQF_local.lz2 = 0.0;
    BQF_local.hz1 = BQF_local.hz2 = 0.0;
    
    
    wavread(infileName, &samples, &samplelength);
    printf("length = %d \n",samplelength);
    out = (int16_t *)malloc(samplelength);
    
    //Process Byquad samples
    for(i=0;i<samplelength;i++)
    {
        *(out+i) = processlpf(&BQF_local, *(samples+i));
    }
    
    //Write output if needeed
    wavwrite(outfileName, out);
    free(samples);
    free(out);
    return 0;
}

