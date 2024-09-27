#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <glob.h>

#define MAX_LINE_LENGTH 256
#define HASH_TABLE_SIZE 1000
#define MAX_CONTEXTS 3

typedef struct
{
    char context[4]; // 甲基化上下文 (CG, CHG, CHH)
    int count;       // 出现次数
    int covered;     // 覆盖次数
    int totalReadsM; // readsM 的总和
    int totalReadsN; // readsN 的总和
} ContextStats;

typedef struct
{
    char chromosome[50];
    ContextStats contexts[MAX_CONTEXTS];
} ChromosomeStats;

typedef struct HashTableEntry
{
    ChromosomeStats data;
    struct HashTableEntry *next;
} HashTableEntry;

// 函数声明
unsigned int hashFunction(const char *str);
ChromosomeStats *findOrInsertChromosome(HashTableEntry **hashTable, const char *chromosome);
void processLine(const char *line, HashTableEntry **hashTable);
void printResultsAsTable(HashTableEntry **hashTable, FILE *outputFile);
void processFile(const char *filename, HashTableEntry **hashTable);
void processFilesWithWildcard(const char *pattern, HashTableEntry **hashTable);

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        fprintf(stderr, "Usage: %s <input_file_pattern> <output_file>\n", argv[0]);
        return 1;
    }

    HashTableEntry *hashTable[HASH_TABLE_SIZE] = {NULL};

    // 处理文件（支持通配符）
    processFilesWithWildcard(argv[1], hashTable);

    // 打开输出文件
    FILE *outputFile = fopen(argv[2], "w");
    if (outputFile == NULL)
    {
        perror("Error opening output file");
        return 1;
    }

    // 打印结果为表格格式并输出到文件
    printResultsAsTable(hashTable, outputFile);

    // 关闭输出文件
    fclose(outputFile);

    return 0;
}

unsigned int hashFunction(const char *str)
{
    unsigned int hash = 0;
    while (*str)
    {
        hash = (hash << 5) + *str++;
    }
    return hash % HASH_TABLE_SIZE;
}

ChromosomeStats *findOrInsertChromosome(HashTableEntry **hashTable, const char *chromosome)
{
    unsigned int index = hashFunction(chromosome);
    HashTableEntry *entry = hashTable[index];

    while (entry != NULL)
    {
        if (strcmp(entry->data.chromosome, chromosome) == 0)
        {
            return &entry->data;
        }
        entry = entry->next;
    }

    // 如果没有找到染色体，则创建一个新的条目
    entry = malloc(sizeof(HashTableEntry));
    if (!entry)
    {
        fprintf(stderr, "Memory allocation error\n");
        exit(1);
    }
    strcpy(entry->data.chromosome, chromosome);
    strcpy(entry->data.contexts[0].context, "CG");
    strcpy(entry->data.contexts[1].context, "CHG");
    strcpy(entry->data.contexts[2].context, "CHH");
    for (int i = 0; i < MAX_CONTEXTS; i++)
    {
        entry->data.contexts[i].count = 0;
        entry->data.contexts[i].covered = 0;
        entry->data.contexts[i].totalReadsM = 0;
        entry->data.contexts[i].totalReadsN = 0;
    }
    entry->next = hashTable[index];
    hashTable[index] = entry;

    return &entry->data;
}

void processFile(const char *filename, HashTableEntry **hashTable)
{
    FILE *file = NULL;
    gzFile gzfile = NULL;
    char line[MAX_LINE_LENGTH];

    // 判断文件扩展名并打开文件
    if (strstr(filename, ".gz") != NULL)
    {
        gzfile = gzopen(filename, "r");
        if (gzfile == NULL)
        {
            perror("Error opening gz file");
            exit(1);
        }

        while (gzgets(gzfile, line, sizeof(line)))
        {
            processLine(line, hashTable);
        }
        gzclose(gzfile);
    }
    else if (strstr(filename, ".txt") != NULL)
    {
        file = fopen(filename, "r");
        if (file == NULL)
        {
            perror("Error opening txt file");
            exit(1);
        }

        while (fgets(line, sizeof(line), file))
        {
            processLine(line, hashTable);
        }
        fclose(file);
    }
    else
    {
        fprintf(stderr, "Unsupported file format. Please provide a .txt or .gz file.\n");
        exit(1);
    }
}

void processFilesWithWildcard(const char *pattern, HashTableEntry **hashTable)
{
    glob_t glob_result;
    memset(&glob_result, 0, sizeof(glob_result));

    // 使用glob函数查找所有匹配的文件
    int return_value = glob(pattern, 0, NULL, &glob_result);
    if (return_value != 0)
    {
        globfree(&glob_result);
        fprintf(stderr, "Error matching files with pattern: %s\n", pattern);
        exit(1);
    }

    // 处理所有匹配的文件
    for (size_t i = 0; i < glob_result.gl_pathc; i++)
    {
        processFile(glob_result.gl_pathv[i], hashTable);
    }

    // 释放glob分配的内存
    globfree(&glob_result);
}

void processLine(const char *line, HashTableEntry **hashTable)
{
    char chromosome[50];
    int position;
    char strand;
    int readsM, readsN, non_methylation;
    char context[4];
    char trinucleotide[4];

    // 解析每一行的数据
    sscanf(line, "%s\t%d\t%c\t%d\t%d\t%s\t%s", chromosome, &position, &strand, &readsM, &non_methylation, context, trinucleotide);

    readsN = readsM + non_methylation;

    // 判断上下文
    int contextIndex = -1;
    if (strcmp(context, "CG") == 0)
        contextIndex = 0;
    else if (strcmp(context, "CHG") == 0)
        contextIndex = 1;
    else if (strcmp(context, "CHH") == 0)
        contextIndex = 2;

    if (contextIndex == -1)
        return; // 不匹配的上下文

    // 查找或插入染色体
    ChromosomeStats *chromStats = findOrInsertChromosome(hashTable, chromosome);

    // 更新统计信息
    chromStats->contexts[contextIndex].count++;
    chromStats->contexts[contextIndex].totalReadsM += readsM;
    chromStats->contexts[contextIndex].totalReadsN += readsN;
    if (readsN >= 1)
    {
        chromStats->contexts[contextIndex].covered++;
    }
}

void printResultsAsTable(HashTableEntry **hashTable, FILE *outputFile)
{
    // 打印表头
    fprintf(outputFile, "Chromosome\tContext\tCount\tcovered\ttotalReadsM\ttotalReadsN\n");

    // 打印每个染色体的统计数据
    for (int i = 0; i < HASH_TABLE_SIZE; i++)
    {
        HashTableEntry *entry = hashTable[i];
        while (entry != NULL)
        {
            for (int j = 0; j < MAX_CONTEXTS; j++)
            {
                fprintf(outputFile, "%s\t%s\t%d\t%d\t%d\t%d\n",
                        entry->data.chromosome,
                        entry->data.contexts[j].context,
                        entry->data.contexts[j].count,
                        entry->data.contexts[j].covered,
                        entry->data.contexts[j].totalReadsM,
                        entry->data.contexts[j].totalReadsN);
            }
            entry = entry->next;
        }
    }
}
