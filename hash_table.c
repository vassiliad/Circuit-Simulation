#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "hash_table.h"

#define HASHES_SIZE 1023

typedef struct HASH_T
{
  char *str;
  int id;
} hash_t;

typedef struct HASH_BUCKET_T
{
  hash_t *hashes;
  int size;
} hash_bucket_t;

static hash_bucket_t hashes[HASHES_SIZE];

int unique_hash=0;


//-------------------
void hash_initialize()
{
  memset(hashes, 0, sizeof(hash_bucket_t)*HASHES_SIZE);
  unique_hash = 1;
}

void hash_cleanup()
{
  int i = 0;
  int num = 0;
  int j = 0;
  for ( i=0; i<HASHES_SIZE; i++ ) {
    if ( hashes[i].size ) {
      for (j=0; j<hashes[i].size; j ++ ) {
        free(hashes[i].hashes[j].str);
        num++;
      }
      free(hashes[i].hashes);
    } 
  }
}

int hash(unsigned char *str)
{
  assert(str != 0 && "Passed NULL to hashing function");

  unsigned long h= 5381;
  int c;

  while ((c = *(str++)))
      h= ((h<< 5) + h) + c; /* hash * 33 + c */

  return h%HASHES_SIZE;
}

int hash_get(char *str)
{
  int h = hash((unsigned char*)str);
  int i;
  
  for (i=0; i< hashes[h].size; i++ ) {
    if ( strcasecmp(hashes[h].hashes[i].str, str) == 0 ) {
      free(str);
      return hashes[h].hashes[i].id;
    }
  }
  

  hashes[h].hashes = (hash_t*) realloc(hashes[h].hashes, sizeof(hash_t)*(hashes[h].size+1));
  hashes[h].hashes[ (hashes[h].size) ].str = str;
  hashes[h].hashes[ (hashes[h].size)++ ].id= unique_hash;
#ifdef VERBOSE_HASH
  printf("translating %s to %d\n", str, unique_hash);
#endif

  return unique_hash++;
}


