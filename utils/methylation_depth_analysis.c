#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <zlib.h>
#include <glob.h>
#include <time.h>

#define MAX_LINE_LENGTH 256
#define MAX_DEPTH 200     // 设定最大覆盖深度为 200，超过200的按200统计
#define MAX_CONTEXTS 10   // 预设最大上下文数量
#define MAX_FILES 200     // 预设最大文件数量

typedef struct {
    char context[10];
    long long depth_counts[MAX_DEPTH + 1];  // 加1以便储存超过200的统计
} ContextStats;

void process_file(const char *filename, ContextStats *context_stats, int *context_count, int max_contexts) {
    FILE *file;
    gzFile gz_file = NULL;
    int is_gz = (strstr(filename, ".gz") != NULL);

    if (is_gz) {
        gz_file = gzopen(filename, "rb");
        if (!gz_file) {
            perror("Error opening gzipped file");
            return;
        }
    } else {
        file = fopen(filename, "r");
        if (!file) {
            perror("Error opening file");
            return;
        }
    }

    char line[MAX_LINE_LENGTH];
    while (is_gz ? gzgets(gz_file, line, sizeof(line)) : fgets(line, sizeof(line), file)) {
        char chr[50];
        int pos;
        char strand;
        int readsM;
        int readsN;
        char context[10];
        char sequence[10];

        // 解析每一行
        if (sscanf(line, "%s %d %c %d %d %s %s", chr, &pos, &strand, &readsM, &readsN, context, sequence) == 7) {
            int depth = readsM + readsN;

            // 确保 depth 在有效范围内
            if (depth > 0) {
                if (depth > MAX_DEPTH) {
                    depth = MAX_DEPTH; // 超过200的覆盖深度都归入200这一类
                }

                // 查找上下文是否已存在
                int context_index = -1;
                for (int i = 0; i < *context_count; i++) {
                    if (strcmp(context_stats[i].context, context) == 0) {
                        context_index = i;
                        break;
                    }
                }

                // 如果上下文不存在，则添加新上下文
                if (context_index == -1) {
                    if (*context_count >= max_contexts) {
                        fprintf(stderr, "Error: Exceeded maximum context limit.\n");
                        if (is_gz) gzclose(gz_file);
                        else fclose(file);
                        return;
                    }
                    strcpy(context_stats[*context_count].context, context);
                    memset(context_stats[*context_count].depth_counts, 0, sizeof(context_stats[*context_count].depth_counts));
                    context_index = (*context_count)++;
                }

                // 更新上下文的覆盖深度计数
                context_stats[context_index].depth_counts[depth]++;
            }
        }
    }

    if (is_gz) gzclose(gz_file);
    else fclose(file);
}

void calculate_methylation_rates(const char *input_pattern, const char *output_filename) {
    glob_t glob_result;
    int context_count = 0;
    ContextStats context_stats[MAX_CONTEXTS];
    int max_files = MAX_FILES;

    if (glob(input_pattern, 0, NULL, &glob_result) != 0) {
        perror("Error matching file pattern");
        return;
    }

    int total_files = glob_result.gl_pathc;
    time_t start_time = time(NULL);

    for (size_t i = 0; i < total_files; i++) {
        time_t file_start_time = time(NULL);
        process_file(glob_result.gl_pathv[i], context_stats, &context_count, MAX_CONTEXTS);
        time_t file_end_time = time(NULL);

        double elapsed = difftime(file_end_time, file_start_time);
        double total_elapsed = difftime(file_end_time, start_time);

        printf("Processed file %zu/%d: %s (%.2f seconds)\n", i + 1, total_files, glob_result.gl_pathv[i], elapsed);
        printf("Total elapsed time: %.2f seconds\n", total_elapsed);
    }

    globfree(&glob_result);

    // 输出二维表格到文件
    FILE *output_file = fopen(output_filename, "w");
    if (!output_file) {
        perror("Error opening output file");
        return;
    }

    // 输出表头
    fprintf(output_file, "Depth");
    for (int i = 0; i < context_count; i++) {
        fprintf(output_file, "\t%s", context_stats[i].context);
    }
    fprintf(output_file, "\n");

    // 输出每个深度的计数
    for (int depth = 1; depth <= MAX_DEPTH; depth++) {
        fprintf(output_file, "%d", depth);
        for (int i = 0; i < context_count; i++) {
            fprintf(output_file, "\t%lld", context_stats[i].depth_counts[depth]);
        }
        fprintf(output_file, "\n");
    }

    fclose(output_file);
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input_pattern> <output_file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    calculate_methylation_rates(argv[1], argv[2]);

    return EXIT_SUCCESS;
}
