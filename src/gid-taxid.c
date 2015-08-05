#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <string.h>
#include "buffers.h"

#include "avl.h"

#define MAX_LINE_BUFFER 2048L


/*#define DEBUG_ME*/

/*---------------------------------------------------------------------------------------------------- */

struct  avl_table *	nameTagAVL  = NULL;
const char UsageString[] = "Usage: gid-taxid <list of genbank ids> <GenBank file mapping gids to taxids (gi_taxid files from ftp://ftp.ncbi.nih.gov/pub/taxonomy/)>\n";
int digits[256];

/*---------------------------------------------------------------------------------------------------- */

struct		bufferedString *globalNameBuffer;

/*---------------------------------------------------------------------------------------------------- */

struct		bufferedTags
{
  long *gID,
  *startTag,
  *lengthTag;
  
  long sLength,
  saLength;
}
*globalTagBuffer;

/*---------------------------------------------------------------------------------------------------- */

struct bufferedTags *allocateNewBufferedTag (void)
{
  struct bufferedTags *newS			= (struct bufferedTags*)malloc (sizeof (struct bufferedTags));
  check_pointer (allocateNewBufferedTag);
  
  check_pointer (newS->gID		= (long*)malloc (sizeof(long)*DEFAULT_STRING_ALLOC));
  check_pointer (newS->startTag	= (long*)malloc (sizeof(long)*DEFAULT_STRING_ALLOC));
  check_pointer (newS->lengthTag	= (long*)malloc (sizeof(long)*DEFAULT_STRING_ALLOC));
  
  newS->sLength  = 0;
  newS->saLength = DEFAULT_STRING_ALLOC;
  return newS;
}

/*---------------------------------------------------------------------------------------------------- */

void appendTagToBuffer (struct bufferedTags * d, const long id, const long start, const long length)
{
  long addThis;
		
  if (d->saLength == d->sLength)
  {
    addThis = d->saLength / 8;
    
    if (DEFAULT_STRING_ALLOC > addThis)
      addThis = DEFAULT_STRING_ALLOC;
    
    d->saLength += addThis;
    check_pointer (d->gID		= (long*)realloc (d->gID,sizeof(long)*d->saLength));
    check_pointer (d->startTag	= (long*)realloc (d->startTag,sizeof(long)*d->saLength));
    check_pointer (d->lengthTag	= (long*)realloc (d->lengthTag,sizeof(long)*d->saLength));
  }
  
  d->gID	   [d->sLength]   = id;
  d->startTag[d->sLength]   = start;
  d->lengthTag[d->sLength++] = length;
  
}


/*---------------------------------------------------------------------------------------------------- */
/*---------------------------------------------------------------------------------------------------- */

struct		storedNameTag
{
  long	taxID,
  startIndex,
  mappedTag,
  length;
};

/*---------------------------------------------------------------------------------------------------- */

struct storedNameTag * allocateNameTag (void)
{
  struct storedNameTag *newS = (struct storedNameTag*)malloc (sizeof (struct storedNameTag));
  check_pointer (newS);
  newS->taxID			= 0;
  newS->startIndex	= 0;
  newS->length		= 0;
  newS->mappedTag		= -1;
  
  return newS;
}

/*---------------------------------------------------------------------------------------------------- */

int compare_tags (const void *avl_a, const void *avl_b, void * xtra)
{
  long n1 = ((struct storedNameTag*)avl_a)->taxID;
  long n2 = ((struct storedNameTag*)avl_b)->taxID;
  
  return n1 - n2;
  /*if (n1 > n2) return 1;
   if (n1 < n2) return -1;
   return 0;*/
}

/*---------------------------------------------------------------------------------------------------- */

struct bufferedString * nameByID (long taxID)
{
  static struct			storedNameTag queryTag;
  struct storedNameTag *	res;
  struct bufferedString * resStr = NULL;
		
  queryTag.taxID = taxID;
  res = (struct storedNameTag *)avl_find(nameTagAVL, &queryTag);
  if (res)
  {
    resStr = allocateNewString();
    appendRangeToString(resStr,globalNameBuffer,res->startIndex,res->startIndex+res->length-1);
  }
  return resStr;
}

/*---------------------------------------------------------------------------------------------------- */

struct storedNameTag * tagByID (long taxID)
{
  static struct	storedNameTag queryTag;
  queryTag.taxID = taxID;
  return  (struct storedNameTag *)avl_find(nameTagAVL, &queryTag);
}




/*---------------------------------------------------------------------------------------------------- */

int main (int argc, const char * argv[])
{
  FILE							*inFile,
  *mapFile;
  
  struct	bufferedString			**currentBuffers;
  
  struct  storedNameTag			*aTag,
  *aTag2;
  
  char	automatonState			= 0,
  currentField			= 0,
  currentChar				= 0,
  line_buffer       [MAX_LINE_BUFFER];
  
  long	currentLineID			= 1,
  expectedFields			= 2,
  mappedID				,
  indexer;
  
  size_t  current_buffer_index    = 0L,
  current_buffer_size     = 0L;
  
  
 
  globalNameBuffer = allocateNewString();
  globalTagBuffer	 = allocateNewBufferedTag();
  
  currentBuffers   = (struct	bufferedString**)malloc (expectedFields*sizeof (struct	bufferedString*));
  nameTagAVL		 = avl_create (compare_tags, NULL, NULL);
  
  for (indexer = 0; indexer < expectedFields; indexer++)
    currentBuffers[indexer] = allocateNewString();
  
		
  if (argc != 3)
  {
    fprintf (stderr,"%s\n", UsageString);
    return 1;
  }
		
  inFile		= fopen (argv[1], "rb");
  
  if (!inFile)
  {
	   fprintf (stderr,"Failed to open input file: %s\n", argv[1]);
	   return 1;
  }
  
  mapFile		= fopen (argv[2], "rb");
  
  if (!mapFile)
  {
	   fprintf (stderr,"Failed to open gid-taxid map file: %s\n", argv[2]);
	   return 1;
  }
  
  for (indexer = 0; indexer < 256; indexer ++) {
    digits [indexer] = 0;
  }
  
  for (indexer = (long)'0'; indexer <= (long)'9'; indexer ++) {
    digits [indexer] = 1;
  }
  
  for (;;) {
    if (current_buffer_index >= current_buffer_size) {
      current_buffer_size = fread (line_buffer, 1L, MAX_LINE_BUFFER, inFile);
      if (current_buffer_size == 0L) {
        break;
      }
      current_buffer_index = 1L;
      currentChar = line_buffer[0];
    } else {
      currentChar = line_buffer[current_buffer_index++];
    }
  
    switch (automatonState)
    {
      case 0: /* start of the line; expecting numbers */
        if (digits[currentChar])
        {
          automatonState = 1; /* reading sequence ID */
          appendCharacterToString(currentBuffers[currentField],currentChar);
        }
        else
          if (!(currentChar == '\n' || currentChar == '\r'))
            reportErrorLine ("Could not find a valid sequence ID to start the line",currentLineID);
        break;
        
      case 1: /* reading file ID */
        if (currentChar >= '0' && currentChar <='9')
          appendCharacterToString(currentBuffers[currentField],currentChar);
        else
          if (currentChar == '\t')
          {
            currentField ++;
            automatonState = 2;
            /* reading file ID */
          }
          else
            reportErrorLine ("Expected a tab following the gi ID",currentLineID);
        break;
        
      case 2: /* read a field */
        if (isalnum(currentChar) || currentChar == '/' || currentChar == '.' || currentChar == '_' || currentChar == '-')
          appendCharacterToString(currentBuffers[currentField],currentChar);
        else
          if (currentChar == '\n' || currentChar == '\r')
          {
            automatonState = 0;
            currentLineID ++;
            currentField = 0;
            
            aTag			 = allocateNameTag();
            aTag->taxID      = atoi(currentBuffers[0]->sData);
            aTag->startIndex = globalNameBuffer->sLength;
            aTag->length	 = appendRangeToString(globalNameBuffer,currentBuffers[1],0,currentBuffers[1]->sLength-1);
            
            appendTagToBuffer (globalTagBuffer, aTag->taxID, aTag->startIndex, aTag->length);
            
            if (aTag->length <= 0)
              reportErrorLine ("Empty name tag",currentLineID);
            if (*avl_probe(nameTagAVL,aTag) != aTag)
            {
              free			(aTag);
                //reportErrorLine ("Duplicate name tag",currentLineID);
            }
            for (indexer = 0; indexer < expectedFields; indexer++)
              clear_buffered_string(currentBuffers[indexer]);
          }
          else
            reportErrorLine ("Unexpected character",currentLineID);
        break;
        
        
    }
  }
  
  fclose (inFile);
		
  aTag	    = allocateNameTag();
  fprintf (stderr, "Reading gid-taxid mapping file...\n");
  currentLineID = 0;
  current_buffer_index = current_buffer_size = 0L;

  for (;;) {
    if (current_buffer_index >= current_buffer_size) {
      current_buffer_size = fread (line_buffer, 1L, MAX_LINE_BUFFER, mapFile);
      if (current_buffer_size == 0L) {
        break;
      }
      current_buffer_index = 1L;
      currentChar = line_buffer[0];
    } else {
      currentChar = line_buffer[current_buffer_index++];
    }
    if (currentLineID > 1 && currentLineID % 5000000 == 0 && automatonState == 0) {
      fprintf (stderr, "Read %ld lines\n", currentLineID);
    }
    switch (automatonState)
    {
      case 0: /* start of the line; expecting numbers */
        if (currentChar >= '0' && currentChar <='9')
        {
          automatonState = 1; /* reading sequence ID */
          appendCharacterToString(currentBuffers[currentField],currentChar);
        }
        else
          if (!(currentChar == '\n' || currentChar == '\r'))
            reportErrorLine ("Could not find a valid sequence ID to start the line",currentLineID);
        break;
        
      case 1: /* reading file ID */
        if (digits[currentChar])
          appendCharacterToString(currentBuffers[currentField],currentChar);
        else
          if (currentChar == '\t')
          {
            currentField ++;
            automatonState = 2;
            /* reading file ID */
          }
          else
            reportErrorLine ("Expected a tab following the tax ID",currentLineID);
        break;
        
      case 2: /* read a field */
        if (digits[currentChar])
          appendCharacterToString(currentBuffers[currentField],currentChar);
        else
          if (currentChar == '\n' || currentChar == '\r')
          {
            automatonState = 0;
            currentLineID ++;
            currentField = 0;
            
            aTag->taxID      = atoi(currentBuffers[0]->sData);
            mappedID		 = atoi(currentBuffers[1]->sData);
            
            aTag2 = avl_find(nameTagAVL,aTag);
            if (aTag2)
              aTag2->mappedTag = mappedID;
            for (indexer = 0; indexer < expectedFields; indexer++)
              clear_buffered_string(currentBuffers[indexer]);
          }
          else
            reportErrorLine ("Unexpected character",currentLineID);
        break;
        
        
    }
  }
  
  fclose (mapFile);
  
  for (indexer = 0L; indexer < expectedFields; indexer++)
    destroy_string(currentBuffers[indexer]);
  
  free(currentBuffers);
  
  for (indexer = 0L; indexer < globalTagBuffer->sLength; indexer++) {
    aTag->taxID      = globalTagBuffer->gID[indexer];
    aTag2 = avl_find(nameTagAVL,aTag);
    if (aTag2) {
      printf ("%ld\t%ld\t", aTag->taxID,aTag2->mappedTag);
    }
    else {
      printf ("%ld\t-1\t", aTag->taxID);
    }
    for (expectedFields = globalTagBuffer->startTag[indexer]; 
         expectedFields < globalTagBuffer->startTag[indexer] + globalTagBuffer->lengthTag[indexer];
         expectedFields ++)
    {
      printf ("%c",globalNameBuffer->sData[expectedFields]);
    }
    printf ("\n");
  }
  
    //printf ("-1,NULL");
  free (aTag);
  return 0;
}
