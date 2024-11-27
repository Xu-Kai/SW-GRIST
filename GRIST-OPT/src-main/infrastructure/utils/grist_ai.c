#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <athread.h>
#include <omp.h>
#include <math.h>

#define NCO 128
#define CO_PE 8
#define CO_PP (NCO/CO_PE)
#define NCI 128
#define LEV 32
#define NBATCH 16
#define BA_PE (64/CO_PE)
#define VWID 8
#define CIBLK 8
#define INPUT_NCI 5
#define INPUT_L 30
#define KERNEL_SIZE 3
#define OUTPUT_NCO 2

#define GB 512
#define GM 512
#define GN 512
#define B 16
#define M 8
#define N 512

struct conv1d3_128x128x32_param {
  float (*conv_out)[NCO][LEV];
  float (*conv_in)[NCI][LEV];
  float (*kern)[NCI][3][VWID];
  float *bias;
  int nbatch;
};

struct Conv1D{
    __attribute__((aligned(256))) float weight[NCO/VWID][NCI][3][VWID]; 
    __attribute__((aligned(256))) float bias[NCO];
    int in_channel;
    int out_channel;
    int kernel_size;
};


#define IN_DIM 512
#define OUT_DIM 512
#define INPUT_DIM 92
#define OUTPUT_DIM 3


struct Linear
{
    /* data */
    __attribute__((aligned(256))) float weight[OUT_DIM/M][IN_DIM][M]; 
    __attribute__((aligned(256))) float bias[OUT_DIM];
    int in_dim;
    int out_dim;
};

float *ra_weight_buffer;
float *cu_weight_buffer;
struct Conv1D cu_weights[24];
struct Linear ra_weights[24];
int ra_batch_size;
int cu_batch_size;

__attribute__((aligned(256))) float (*ra_out1)[OUT_DIM][B]; 
__attribute__((aligned(256))) float (*ra_out2)[OUT_DIM][B]; 
__attribute__((aligned(256))) float (*ra_residual)[OUT_DIM][B]; 


__attribute__((aligned(256))) float (*cu_out1)[NCO][LEV]; 
__attribute__((aligned(256))) float (*cu_out2)[NCO][LEV]; 
__attribute__((aligned(256))) float (*cu_residual)[NCO][LEV];





void read_weight(const char* file_name, float **weight_buffer)
{
    FILE *fp = fopen(file_name, "r");
    long file_size = -1;
    long cur_offset = ftell(fp);
    fseek(fp, 0, SEEK_END);
    file_size = ftell(fp);
    fseek(fp, cur_offset, SEEK_SET);
    *weight_buffer = (float*)malloc(file_size);
    fread(*weight_buffer, 4, file_size/4, fp);
    // printf("file size is %lld\n", file_size);
    fclose(fp);

};

void init_ra_weight_(int *batch_size)
{   

    char * file_name = NULL;
    if(file_name == NULL)
    {
        // printf("batch size is %d\n", *batch_size);
        // printf("file name is not set. set default path: /home/export/online1/mdt00/shisuan/swgbcm/yc_data/weight/ra_weight_np.bin\n");
        read_weight("/home/export/online1/mdt00/shisuan/swgbcm/yc_data/weight/ra_weight_np.bin", &ra_weight_buffer);
    }
    else
    {
        read_weight(file_name, &ra_weight_buffer);
    }


    long offset = 0;

    for(int i = 0; i < OUT_DIM;  i++)
    {
        for(int j = 0; j < IN_DIM; j ++)
        {
            ra_weights[0].weight[i/M][j][i%M] = j < INPUT_DIM?-ra_weight_buffer[offset + i * INPUT_DIM + j]:0;
        }
        ra_weights[0].bias[i] = ra_weight_buffer[offset + INPUT_DIM * OUT_DIM + i];
    }
    offset += INPUT_DIM * OUT_DIM + OUT_DIM;


    for(int k = 1; k < 14; k ++)
    {
        for(int i = 0; i < OUT_DIM;  i++)
        {
            for(int j = 0; j < IN_DIM; j ++)
            {

                ra_weights[k].weight[i/M][j][i%M] = -ra_weight_buffer[offset + i * IN_DIM + j];

            }
            ra_weights[k].bias[i] = ra_weight_buffer[offset + IN_DIM * OUT_DIM + i];
        }
        offset += IN_DIM * OUT_DIM + OUT_DIM;
    }
    for(int i = 0; i < OUT_DIM;  i++)
    {
        for(int j = 0; j < IN_DIM; j ++)
        {
            ra_weights[14].weight[i/M][j][i%M] = i < OUTPUT_DIM?-ra_weight_buffer[offset + i * IN_DIM + j]:0;
        }
        ra_weights[14].bias[i] = i < OUTPUT_DIM ? ra_weight_buffer[offset + OUTPUT_DIM * IN_DIM + i]:0;
    }

    ra_batch_size = *batch_size;
    int aligned_size = ((ra_batch_size +  B - 1)/B)*B;
    // printf("ra batch size is %d %d\n", ra_batch_size, aligned_size);
    ra_out1 = malloc(aligned_size * IN_DIM * sizeof(float));
    ra_out2 = malloc(aligned_size * IN_DIM * sizeof(float));
    ra_residual = malloc(aligned_size * IN_DIM * sizeof(float));
    // printf("aligned size is %d\n", aligned_size);

};
extern void convert_weight(float (*out)[NCI][KERNEL_SIZE][VWID], float (*in)[NCI][KERNEL_SIZE]);



void init_cu_weight_(int *batch_size)
{   
    char *file_name = NULL;
    if(file_name == NULL)
    {
        // printf("batch size is %d\n", *batch_size);
        // printf("file name is not set. set default path: /home/export/online1/mdt00/shisuan/swgbcm/yc_data/weight/cu_weight_np.bin\n");
        read_weight("/home/export/online1/mdt00/shisuan/swgbcm/yc_data/weight/cu_weight_np.bin", &cu_weight_buffer);
    }
    else
    {
        read_weight(file_name, &cu_weight_buffer);
    }


    float conv1d_in[NCO][NCI][KERNEL_SIZE];
    long offset = 0;

    for(int i = 0; i < NCO; i ++)
    {
        for(int j = 0; j < NCI; j ++)
        {
            if(j < INPUT_NCI)
            {
                conv1d_in[i][j][0] = cu_weight_buffer[offset + i*INPUT_NCI*KERNEL_SIZE + j*KERNEL_SIZE + 0];
                conv1d_in[i][j][1] = cu_weight_buffer[offset + i*INPUT_NCI*KERNEL_SIZE + j*KERNEL_SIZE + 1];;
                conv1d_in[i][j][2] = cu_weight_buffer[offset + i*INPUT_NCI*KERNEL_SIZE + j*KERNEL_SIZE + 2];;
            }
            else
            {
                conv1d_in[i][j][0] = 0;
                conv1d_in[i][j][1] = 0;
                conv1d_in[i][j][2] = 0;
            }
        }
        cu_weights[0].bias[i] = cu_weight_buffer[offset + NCO * INPUT_NCI * KERNEL_SIZE + i];
    }


    offset += NCO * INPUT_NCI * KERNEL_SIZE + NCO;
    convert_weight(cu_weights[0].weight, conv1d_in);

    for(int k = 1; k < 11; k ++)
    {
        float *tmp = &cu_weight_buffer[offset];
        convert_weight(cu_weights[k].weight, tmp);
        for(int i = 0; i < NCO; i ++)
        {
            cu_weights[k].bias[i] = cu_weight_buffer[offset + NCO * NCI * KERNEL_SIZE + i];
        }
        offset += NCO * NCI * KERNEL_SIZE + NCO;
    }

    for(int i = 0; i < NCO; i ++)
    {
        for(int j = 0; j < NCI; j ++)
        {
            if(i<OUTPUT_NCO)
            {
                conv1d_in[i][j][0] = cu_weight_buffer[offset + i*NCI*KERNEL_SIZE + j*KERNEL_SIZE + 0];
                conv1d_in[i][j][1] = cu_weight_buffer[offset + i*NCI*KERNEL_SIZE + j*KERNEL_SIZE + 1];;
                conv1d_in[i][j][2] = cu_weight_buffer[offset + i*NCI*KERNEL_SIZE + j*KERNEL_SIZE + 2];;
            }
            else
            {
                conv1d_in[i][j][0] = 0;
                conv1d_in[i][j][1] = 0;
                conv1d_in[i][j][2] = 0;
            }
        }
         cu_weights[11].bias[i] = i < OUTPUT_NCO? cu_weight_buffer[offset + OUTPUT_NCO * NCI * KERNEL_SIZE + i]:0;

    }
    convert_weight(cu_weights[11].weight, conv1d_in);


    cu_batch_size = *batch_size;
    cu_out1 = malloc(cu_batch_size*NCO*LEV*sizeof(float));
    cu_out2 = malloc(cu_batch_size*NCO*LEV*sizeof(float));
    cu_residual = malloc(cu_batch_size*NCO*LEV*sizeof(float));

};


extern void slave_conv1d3_128x128(struct conv1d3_128x128x32_param *pm);



void cu_resnet_forward_(float *input, float *output, int *batch_size)

{
    if(*batch_size > cu_batch_size)
    {
        printf("cu resnet wrong batch %d %d\n", *batch_size, cu_batch_size);
    }
    // printf("weight is %f input %f\n", ra_weights[0].weight[0][2], input[5]);

    int batch = *batch_size;
    #pragma omp target parallel for collapse(2)
    for(int i = 0; i < batch; i ++)
    {
        for(int j = 0; j < NCI; j ++)
        {
            if(j < INPUT_NCI)
            {
                for(int k = 0; k < LEV; k ++)
                {
                    cu_out1[i][j][k] = k < INPUT_L?input[i*INPUT_NCI*INPUT_L + j * INPUT_L + k]:0;
                }
            }
            else
            {
                for(int k = 0; k < LEV; k ++)
                {     
                    cu_out1[i][j][k] = 0;
                }
            }
        }
    }

    // printf("num is %f\n", cu_out1[0][0][3]);

    struct conv1d3_128x128x32_param pm = {cu_residual, cu_out1, cu_weights[0].weight, cu_weights[0].bias, batch};

    __real_athread_spawn(slave_conv1d3_128x128, &pm, 1);
    athread_join();
     for(int i = 0; i < 5; i ++)
     {
        #pragma omp target parallel for collapse(2)
        for(int b = 0; b < batch; b ++)
        {
            for(int j = 0; j < NCO; j ++)
            {
                for(int k = 0; k < LEV; k ++)
                {
                    cu_out1[b][j][k] = cu_residual[b][j][k] < 0?0:cu_residual[b][j][k];

                }
            }

        }

        struct conv1d3_128x128x32_param pm1 = {cu_out2, cu_out1, cu_weights[2*i + 1].weight, cu_weights[2*i + 1].bias, batch};
        __real_athread_spawn(slave_conv1d3_128x128, &pm1, 1);
        athread_join();
        #pragma omp target parallel for collapse(2)
        for(int b = 0; b < batch; b ++)
        {
            for(int j = 0; j < NCO; j ++)
            {
                for(int k = 0; k < LEV; k ++)
                {
                    if(cu_out2[b][j][k] < 0)
                    {
                        cu_out2[b][j][k] = 0;
                    }

                }
            }
        }

        struct conv1d3_128x128x32_param pm2 = {cu_out1, cu_out2, cu_weights[2*i + 2].weight, cu_weights[2*i + 2].bias, batch};
        __real_athread_spawn(slave_conv1d3_128x128, &pm2, 1);
        athread_join();
        
        #pragma omp target parallel for collapse(2)
        for(int b = 0; b < batch; b ++)
        {
            for(int j = 0; j < NCO; j ++)
            {
                for(int k = 0; k < LEV; k ++)
                {
                     cu_residual[b][j][k] += cu_out1[b][j][k];
                }
            }
        }


     }

    #pragma omp target parallel for  collapse(2)
    for(int b = 0; b < batch; b ++)
    {
        for(int j = 0; j < NCO; j ++)
        {
            for(int k = 0; k < LEV; k ++)
            {

                cu_out1[b][j][k] = cu_residual[b][j][k] < 0?0:cu_residual[b][j][k];
            }
        }
    }

    struct conv1d3_128x128x32_param pm1 = {cu_out2, cu_out1, cu_weights[11].weight, cu_weights[11].bias, batch};
    __real_athread_spawn(slave_conv1d3_128x128, &pm1, 1);
    athread_join();
    
    #pragma omp target parallel for collapse(2)
    for(int b = 0; b < batch; b ++)
    {
        for(int j = 0; j < OUTPUT_NCO; j ++)
        {
            for(int k = 0; k < INPUT_L; k ++)
            {
                float t1 = exp(cu_out2[b][j][k]);
                float t2 = exp(-cu_out2[b][j][k]);
                output[b*OUTPUT_NCO*INPUT_L + j*INPUT_L + k]= (t1 - t2)/(t1 + t2);

            }
        }
    }
}



// void linear(float *a, float *b, float *c, float *bias, int m, int n, int k)
// {
//     //#pragma omp target parallel for collapse(2)
//     for(int i = 0; i < m; i ++){
//         for(int j = 0; j < n; j ++)
//         {
//             float sum = 0;
//             for (int t = 0; t < k; t++){
//                 sum += a[i * k + t] * b[j * k + t];
//             }
//             c[i*n + j] = sum + bias[j];
//         }
//     }
// }

struct mm_par_t {
  float (*gout)[GM][B], (*gin)[GN][B], (*gw)[GN][M], *gbias;
  int nB;
};
extern void slave_linear(void *);
void ra_forward_(float *input, float *output,  int *batch_size)
{
    if(ra_batch_size < *batch_size)
    {
        printf("wrong in ra_forward %d %d\n", ra_batch_size, *batch_size);
    }

    int m = ((*batch_size + B - 1) / B) * B;

    #pragma omp target parallel for
    for(int i = 0; i < m; i ++)
    {
        for(int j = 0; j < IN_DIM; j ++)
        {
            if( j < INPUT_DIM && i < ra_batch_size)
            {
                ra_residual[i/B][j][i%B] = input[i*INPUT_DIM + j];
            }
            else{
                ra_residual[i/B][j][i%B]  = 0;
            }
        }
    }

        for(int i = 0; i < 7; i ++)
        { 
            struct mm_par_t param1 = {ra_out1, ra_residual,  ra_weights[2*i].weight,  ra_weights[2*i].bias, m/B};
            __real_athread_spawn(slave_linear, &param1, 1);
            athread_join();

            // linear(&(ra_residual[0][0]), &(ra_weights[2*i].weight[0][0]), &(ra_out1[0][0]), ra_weights[2*i].bias, m, OUT_DIM, IN_DIM);
            #pragma omp target parallel for collapse(2)
            for(int j = 0; j < m/B; j ++)
            {
                for(int k = 0; k < OUT_DIM; k ++)
                {
                    for(int b = 0; b < B; b ++)
                    {
                        if (ra_out1[j][k][b] < 0)
                            ra_out1[j][k][b] = 0; 
                    }
                }
            }

            // linear(&(ra_out1[0][0]), &(ra_weights[2*i + 1].weight[0][0]), &(ra_out2[0][0]), ra_weights[2*i + 1].bias, m, OUT_DIM, IN_DIM);
            struct mm_par_t param2 = {ra_out2, ra_out1,  ra_weights[2*i + 1].weight,  ra_weights[2*i + 1].bias, m/B};
            __real_athread_spawn(slave_linear, &param2, 1);
            athread_join();     
    
            #pragma omp target parallel for collapse(2)
            for(int j = 0; j < m/B; j ++)
            {
                for(int k = 0; k < OUT_DIM; k ++)
                {
                    for(int b = 0; b < B; b ++)
                    {
                        if(ra_out2[j][k][b] < 0)
                        {
                            ra_out2[j][k][b] = 0;
                        }
                        ra_residual[j][k][b] += ra_out2[j][k][b];   
                    }
                }
            }
        }

    struct mm_par_t param2 = {ra_out1, ra_residual,  ra_weights[14].weight,  ra_weights[14].bias, m/B};
    __real_athread_spawn(slave_linear, &param2, 1);
    athread_join();
    #pragma omp target parallel for
    for(int i = 0; i < m/B; i ++)
    {
        for(int j = 0; j <  OUTPUT_DIM; j ++)
        {
            for(int b = 0; b < B; b ++)
            {
                if(i * B + b < ra_batch_size)
                    output[(i*B + b)*OUTPUT_DIM + j] = ra_out1[i][j][b];
            }
        }
    }
}



#ifdef AI_TEST

int main()
{

    athread_init();

     int cu_batch = 16;
     float *cu_input = NULL;
     float *cu_output = NULL;
     init_cu_weight_(&cu_batch);

     read_weight("./cu_input.bin",  &cu_input);
     read_weight("./cu_output.bin", &cu_output);
     float *cu_result = (float*)malloc(cu_batch*NCO*30*sizeof(float));
     cu_resnet_forward_(cu_input, cu_result, &cu_batch);
     float maxrel =  1e-5;
     for(int i = 0; i < 16; i ++)
     {
         for(int j = 0; j < 2; j ++)
        {
             for(int k = 0; k < 30; k++)
             {
                float ref = cu_output[i* OUTPUT_NCO * 30  + j *30 + k];
                float out = cu_result[i*OUTPUT_NCO*INPUT_L + j * INPUT_L + k];
                 float diff = fabs(ref - out);
                 float max = fabs(ref) > fabs(out) ? fabs(ref) : fabs(out);
                 if (diff / max > maxrel) {
                    // maxrel = diff/max;
                    printf("%d %d %d %f %f %f\n", i, j, k, maxrel, ref, out);
                 }
             }
         }
     }
   
    int ra_batch = 16;
    init_ra_weight_(&ra_batch);
    float *ra_input = NULL;
    float *ra_output = NULL;
    read_weight("./ra_input.bin",  &ra_input);
    read_weight("./ra_output.bin", &ra_output);

    float *ra_result = malloc(sizeof(float) * ra_batch * OUT_DIM);
    maxrel = 1e-5;
    ra_forward_(ra_input, ra_result, &ra_batch);
    for(int i = 0; i < ra_batch; i ++)
    {
        for(int j = 0; j < OUTPUT_DIM; j ++)
        {
            float ref = ra_output[i*OUTPUT_DIM  + j];
            float out = ra_result[i*OUTPUT_DIM + j];
            float diff = fabs(ref - out);
            float max = fabs(ref) > fabs(out) ? fabs(ref) : fabs(out);
            if (diff / max > maxrel) {
                // maxrel = diff/max;
                printf("%d %d %f %f %f\n", i, j, maxrel, ref, out);
            }
        }
    }

    return 0;
}



#endif
