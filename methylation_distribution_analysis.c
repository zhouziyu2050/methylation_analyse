#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <zlib.h>
#include <glob.h>
#include <time.h>

#define MAX_LINE_LENGTH 256
#define MAX_PERCENTAGE 101 // 0-100%
#define MAX_CONTEXTS 10    // 预设最大上下文数量
#define MAX_FILES 200      // 预设最大文件数量

typedef struct
{
    char context[10];
    long long percentage_counts[MAX_PERCENTAGE];
    long long readsM_sums[MAX_PERCENTAGE];
    long long readsN_sums[MAX_PERCENTAGE];
} ContextStats;

void process_file(const char *filename, ContextStats *context_stats, int *context_count, int max_contexts)
{
    FILE *file;
    gzFile gz_file = NULL;
    int is_gz = (strstr(filename, ".gz") != NULL);

    if (is_gz)
    {
        gz_file = gzopen(filename, "rb");
        if (!gz_file)
        {
            perror("Error opening gzipped file");
            return;
        }
    }
    else
    {
        file = fopen(filename, "r");
        if (!file)
        {
            perror("Error opening file");
            return;
        }
    }

    char line[MAX_LINE_LENGTH];
    while (is_gz ? gzgets(gz_file, line, sizeof(line)) : fgets(line, sizeof(line), file))
    {
        char chr[50];
        int pos;
        char strand;
        int readsM;
        int readsN;
        int non_methylation;
        char context[10];
        char sequence[10];

        // 解析每一行
        if (sscanf(line, "%s %d %c %d %d %s %s", chr, &pos, &strand, &readsM, &non_methylation, context, sequence) == 7)
        {
            readsN = readsM + non_methylation;
            if (readsN > 0)
            {
                // 计算甲基化率并转换为百分比
                double methylation_rate = (double)readsM /readsN * 100;
                int percentage = (int)round(methylation_rate);

                // 查找上下文是否已存在
                int context_index = -1;
                for (int i = 0; i < *context_count; i++)
                {
                    if (strcmp(context_stats[i].context, context) == 0)
                    {
                        context_index = i;
                        break;
                    }
                }

                // 如果上下文不存在，则添加新上下文
                if (context_index == -1)
                {
                    if (*context_count >= max_contexts)
                    {
                        fprintf(stderr, "Error: Exceeded maximum context limit.\n");
                        if (is_gz)
                            gzclose(gz_file);
                        else
                            fclose(file);
                        return;
                    }
                    strcpy(context_stats[*context_count].context, context);
                    memset(context_stats[*context_count].percentage_counts, 0, sizeof(context_stats[*context_count].percentage_counts));
                    memset(context_stats[*context_count].readsM_sums, 0, sizeof(context_stats[*context_count].readsM_sums));
                    memset(context_stats[*context_count].readsN_sums, 0, sizeof(context_stats[*context_count].readsN_sums));
                    context_index = (*context_count)++;
                }

                // 更新上下文的百分比计数和 readsM, readsN 总和
                if (percentage >= 0 && percentage < MAX_PERCENTAGE)
                {
                    context_stats[context_index].percentage_counts[percentage]++;
                    context_stats[context_index].readsM_sums[percentage] += readsM;
                    context_stats[context_index].readsN_sums[percentage] += readsN;

                    // printf("readsN: %d\n",context_stats[context_index].readsN_sums[percentage]);
                }
            }
        }
    }
    if (is_gz)
        gzclose(gz_file);
    else
        fclose(file);
}

void calculate_methylation_rates(const char *input_pattern, const char *output_filename)
{
    glob_t glob_result;
    int context_count = 0;
    ContextStats context_stats[MAX_CONTEXTS];
    int max_files = MAX_FILES;

    if (glob(input_pattern, 0, NULL, &glob_result) != 0)
    {
        perror("Error matching file pattern");
        return;
    }

    int total_files = glob_result.gl_pathc;
    time_t start_time = time(NULL);

    for (size_t i = 0; i < total_files; i++)
    {
        time_t file_start_time = time(NULL);
        process_file(glob_result.gl_pathv[i], context_stats, &context_count, MAX_CONTEXTS);
        time_t file_end_time = time(NULL);

        double elapsed = difftime(file_end_time, file_start_time);
        double total_elapsed = difftime(file_end_time, start_time);

        printf("Processed file %zu/%d: %s (%.2f seconds)\n", i + 1, total_files, glob_result.gl_pathv[i], elapsed);
        printf("Total elapsed time: %.2f seconds\n", total_elapsed);
    }

    globfree(&glob_result);

    // 输出统计结果到文件
    FILE *output_file = fopen(output_filename, "w");
    if (!output_file)
    {
        perror("Error opening output file");
        return;
    }

    fprintf(output_file, "context\tmethylation_level\tcount\treadsM\treadsN\n");
    for (int i = 0; i < context_count; i++)
    {
        for (int j = 0; j < MAX_PERCENTAGE; j++)
        {
            if (context_stats[i].percentage_counts[j] > 0)
            {
                fprintf(output_file, "%s\t%d\t%d\t%d\t%d\n", context_stats[i].context, j, context_stats[i].percentage_counts[j],
                        context_stats[i].readsM_sums[j], context_stats[i].readsN_sums[j]);
            }
        }
    }

    fclose(output_file);
}

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        fprintf(stderr, "Usage: %s <input_pattern> <output_file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    calculate_methylation_rates(argv[1], argv[2]);

    return EXIT_SUCCESS;
}
