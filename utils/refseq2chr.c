#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#define MAX_LINE_LEN 1024
#define HASH_SIZE 23

typedef struct HashNode {
    char refseq[20];
    char chr[10];
    struct HashNode *next;
} HashNode;

typedef struct {
    HashNode *table[HASH_SIZE];
} HashTable;

unsigned int hash(const char *str) {
    unsigned int hash = 0;
    while (*str) {
        hash = (hash << 5) + *str++;
    }
    return hash % HASH_SIZE;
}

HashTable *createHashTable() {
    HashTable *hashTable = (HashTable *)malloc(sizeof(HashTable));
    for (int i = 0; i < HASH_SIZE; i++) {
        hashTable->table[i] = NULL;
    }
    return hashTable;
}

void insert(HashTable *hashTable, const char *refseq, const char *chr) {
    unsigned int index = hash(refseq);
    HashNode *newNode = (HashNode *)malloc(sizeof(HashNode));
    strcpy(newNode->refseq, refseq);
    strcpy(newNode->chr, chr);
    newNode->next = hashTable->table[index];
    hashTable->table[index] = newNode;
}

const char *search(HashTable *hashTable, const char *refseq) {
    unsigned int index = hash(refseq);
    HashNode *node = hashTable->table[index];
    while (node != NULL) {
        if (strcmp(node->refseq, refseq) == 0) {
            return node->chr;
        }
        node = node->next;
    }
    return NULL;
}

void freeHashTable(HashTable *hashTable) {
    for (int i = 0; i < HASH_SIZE; i++) {
        HashNode *node = hashTable->table[i];
        while (node != NULL) {
            HashNode *temp = node;
            node = node->next;
            free(temp);
        }
    }
    free(hashTable);
}

int isGzipped(const char *filename) {
    return (strstr(filename, ".gz") != NULL);
}

FILE *openFile(const char *filename, const char *mode, int isGzipped) {
    if (isGzipped) {
        gzFile gz = gzopen(filename, mode);
        if (!gz) {
            perror("Error opening gzipped file");
            return NULL;
        }
        return (FILE *)gz;
    } else {
        FILE *file = fopen(filename, mode);
        if (!file) {
            perror("Error opening file");
            return NULL;
        }
        return file;
    }
}

int closeFile(FILE *file, int isGzipped) {
    if (isGzipped) {
        return gzclose((gzFile)file);
    } else {
        return fclose(file);
    }
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input_file> <output_file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    int inputGzipped = isGzipped(argv[1]);
    int outputGzipped = isGzipped(argv[2]);

    HashTable *hashTable = createHashTable();
    insert(hashTable, "NC_000067.7", "chr1");
    insert(hashTable, "NC_000068.8", "chr2");
    insert(hashTable, "NC_000069.7", "chr3");
    insert(hashTable, "NC_000070.7", "chr4");
    insert(hashTable, "NC_000071.7", "chr5");
    insert(hashTable, "NC_000072.7", "chr6");
    insert(hashTable, "NC_000073.7", "chr7");
    insert(hashTable, "NC_000074.7", "chr8");
    insert(hashTable, "NC_000075.7", "chr9");
    insert(hashTable, "NC_000076.7", "chr10");
    insert(hashTable, "NC_000077.7", "chr11");
    insert(hashTable, "NC_000078.7", "chr12");
    insert(hashTable, "NC_000079.7", "chr13");
    insert(hashTable, "NC_000080.7", "chr14");
    insert(hashTable, "NC_000081.7", "chr15");
    insert(hashTable, "NC_000082.7", "chr16");
    insert(hashTable, "NC_000083.7", "chr17");
    insert(hashTable, "NC_000084.7", "chr18");
    insert(hashTable, "NC_000085.7", "chr19");
    insert(hashTable, "NC_000086.8", "chrX");
    insert(hashTable, "NC_000087.8", "chrY");
    insert(hashTable, "NC_005089.1", "chrMT");

    FILE *inputFile = openFile(argv[1], "r", inputGzipped);
    if (!inputFile) {
        freeHashTable(hashTable);
        return EXIT_FAILURE;
    }

    FILE *outputFile = openFile(argv[2], "w", outputGzipped);
    if (!outputFile) {
        closeFile(inputFile, inputGzipped);
        freeHashTable(hashTable);
        return EXIT_FAILURE;
    }

    char line[MAX_LINE_LEN];
    while (fgets(line, sizeof(line), inputFile)) {
        char *refseq = strtok(line, "\t");
        const char *chr = search(hashTable, refseq);
        if (chr) {
            fprintf(outputFile, "%s", chr);
        } else {
            fprintf(outputFile, "%s", refseq); // 保留未匹配行的第一列
        }
        char *remaining = strtok(NULL, "\n");
        if (remaining) {
            fprintf(outputFile, "\t%s\n", remaining); // 保留其余列
        } else {
            fprintf(outputFile, "\n");
        }
    }

    closeFile(inputFile, inputGzipped);
    closeFile(outputFile, outputGzipped);
    freeHashTable(hashTable);

    return EXIT_SUCCESS;
}
