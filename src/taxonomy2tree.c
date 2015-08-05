
#include "avl.h"
#include "buffers.h"
#include "taxinfo.h"
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#define	NUMBER_OF_FIELDS		  24
#define AVL_THRESHOLD			  64


static char *const Usage = "taxonomy2tree \
\n\t<tax dump file(tab separated taxonomy dump)>\
\n\t<max tree level (<=0 to show all)>\
\n\t<write tree file here>\
\n\t<write tab separated summary here>\
\n\tinclude empty nodes; 0 or 1; default 0)]\n";

/*---------------------------------------------------------------------------------------------------- */

void			check_pointer	(void*);

char			validTaxonNameChar[256];

long			currentLineID			= 1;

/*---------------------------------------------------------------------------------------------------- */

struct  avl_table *	  idTagAVL     = NULL,
*   uniqueTags   = NULL,
**  byLevel	   = NULL;

/*---------------------------------------------------------------------------------------------------- */

struct		bufferedString *globalNameBuffer;

/*---------------------------------------------------------------------------------------------------- */


struct treeNode
{
  struct treeNode * parent;
  long   startIndex,
  hitCount,
  length
  /*,beenhere*/;
  
  struct vector   * children;
  struct avl_table* cachedChildren;
  
} *globalTreeRoot;

/*---------------------------------------------------------------------------------------------------- */

struct treeNode * allocateNewTreeNode (void)
{
  struct treeNode *newN = (struct treeNode*)malloc (sizeof (struct treeNode));
  check_pointer (newN);
  newN->parent = NULL;
  newN->startIndex     = 0;
  newN->length	     = 0;
  newN->hitCount       = 0;
  //newN->beenhere		 = 0;
  newN->cachedChildren = NULL;
  newN->children		 = allocateNewVector();
  check_pointer (newN->children);
  return newN;
}

/*---------------------------------------------------------------------------------------------------- */

int compare_tree_nodes (const void *avl_a, const void *avl_b, void * xtra)
{
  long  t1 = ((struct treeNode*)avl_a)->startIndex,
  t2 = ((struct treeNode*)avl_b)->startIndex;
  
  if (t1 < t2)
    return -1;
  if (t1 > t2)
    return 1;
  
  return 0;
}

/*---------------------------------------------------------------------------------------------------- */

void	destroyTreeNode (struct treeNode* n)
{
  if (n->cachedChildren)
    avl_destroy (n->cachedChildren, NULL);
  free (n->children);
  free (n);
}


/*---------------------------------------------------------------------------------------------------- */

void	traverseTree (FILE *summaryFile, struct treeNode* n, long maxDepth, long currentDepth)
{
  long i = 0,
  printedPars = 0;
  char c;
  
  if (currentDepth <= maxDepth)
  {
    //echoNode (n,currentDepth);
    if (n->children->vLength && maxDepth > currentDepth)
    {
      if (n->startIndex < 0)
      {
        for (i=0; i>n->startIndex && currentDepth - i < maxDepth; i--, printedPars++)
          fputc ('(', summaryFile);
        if (currentDepth - i < maxDepth)
          i = 0;
        else
          fprintf (summaryFile, "n:%ld", n->hitCount);
        
      }
      else
        fputc ('(',summaryFile);
      
      if (i==0)
        for (; i<n->children->vLength; i++)
        {
          traverseTree (summaryFile, ((struct treeNode**)n->children->vData)[i], maxDepth, currentDepth+1);
          if (i<n->children->vLength-1)
            fputc (',',summaryFile);
        }
      
      if (n->startIndex < 0)
        for (;printedPars>1;printedPars--)
          fprintf (summaryFile,")n:%ld", n->hitCount);
      fputc (')',summaryFile);
    }
    
    if (n->startIndex>= 0)
    {
      for (i=n->startIndex; i<n->startIndex+n->length;i++)
      {
        c = globalNameBuffer->sData[i];
        if (validTaxonNameChar [c])
          fputc (c,summaryFile);
        else
          fputc ('_',summaryFile);
      }
      fprintf (summaryFile,":%ld", n->hitCount);
    }
    else
      fprintf (summaryFile, "n:%ld", n->hitCount);
  }
}

/*---------------------------------------------------------------------------------------------------- */
struct treeNode * addAChild (struct treeNode* p, struct treeNode * c, char killIfD)
{
  struct treeNode * c2			= NULL;
  long			  i				= 0;
  char			  addNode		= 0;
  
  if (p->cachedChildren)
  {
    if ((c2 = ((struct treeNode*)avl_find (p->cachedChildren, c))) == NULL)
      addNode = 1;
    /*else
     {
     if (c2->startIndex > globalNameBuffer->sLength || c2->startIndex < 0)
     printf ("Reject node add %x %d %d\n", c2, c2->startIndex, c->startIndex);
     }*/
  }
  else
  {
    for (; i<p->children->vLength; i++)
      if (((struct treeNode**)p->children->vData)[i]->startIndex == c->startIndex)
        break;
    addNode = (p->children->vLength == i);
  }
  
  if (addNode)
  {
    if (c->parent && c->parent != p)
    {
      c2 = allocateNewTreeNode();
      c2->startIndex = c->startIndex;
      c2->length	   = c->length;
      c2->hitCount   = c->hitCount;
      //if (p->cachedChildren)
      //	fprintf (stdout, "%x %x\n", p, c->parent);
      c = c2;
      
    }
    c->parent  = p;
    appendValueToVector (p->children, (long)c);
    if (p->children->vLength > AVL_THRESHOLD) {
      if (p->cachedChildren == NULL) {
        //fprintf (stdout, "Switch %d\n", currentLineID);
        // switch over the the avl representation of nodes
        p->cachedChildren = avl_create (compare_tree_nodes, NULL, NULL);
        for (i=0; i<p->children->vLength; i++)
          avl_probe(p->cachedChildren, ((struct node**)p->children->vData)[i]);
      }
      else {
        avl_probe (p->cachedChildren, c);
      }
    }
  }
  else
  {
    if (killIfD)
      destroyTreeNode (c);
    if (c2)
      return c2;
    else
      c = ((struct treeNode**)p->children->vData)[i];
    
  }
  return c;
}


/*---------------------------------------------------------------------------------------------------- */
/*---------------------------------------------------------------------------------------------------- */

struct		storedIDTag
{
  long		taxID,
  hit_count;
  struct		treeNode * tNode;
};

/*---------------------------------------------------------------------------------------------------- */

struct storedIDTag *allocateIDTag (void)
{
  struct storedIDTag *newS = (struct storedIDTag*)malloc (sizeof (struct storedIDTag));
  check_pointer (newS);
  newS->taxID		  = -1;
  newS->hit_count   = 0;
  newS->tNode		  = NULL;
  return newS;
}

/*---------------------------------------------------------------------------------------------------- */

int compare_id_tags (const void *avl_a, const void *avl_b, void * xtra)
{
  long  t1 = ((struct storedIDTag*)avl_a)->taxID,
  t2 = ((struct storedIDTag*)avl_b)->taxID;
  
  if (t1 < t2)
    return -1;
  if (t1 > t2)
    return 1;
  
  return 0;
}


/*---------------------------------------------------------------------------------------------------- */
/*---------------------------------------------------------------------------------------------------- */

struct		storedSequenceTag
{
  long		startIndex,
  length;
  
  struct		treeNode * refNode;
};

/*---------------------------------------------------------------------------------------------------- */

struct storedSequenceTag *allocateStringTag (void)
{
  struct storedSequenceTag *newS = (struct storedSequenceTag*)malloc (sizeof (struct storedSequenceTag));
  check_pointer (newS);
  newS->startIndex  = -1;
  newS->length	  = 0;
  newS->refNode	  = NULL;
  return newS;
}


/*---------------------------------------------------------------------------------------------------- */
/*
 struct storedSequenceTag * storedSequenceTag (void)
 {
	struct storedSequenceTag *newS = (struct storedSequenceTag*)malloc (sizeof (struct storedSequenceTag));
	check_pointer		(newS);
	newS->startIndex	= 0;
	newS->length		= 0;
	return newS;
 }*/

/*---------------------------------------------------------------------------------------------------- */

int compare_tags (const void *avl_a, const void *avl_b, void * xtra)
{
  char* n1 = globalNameBuffer->sData+((struct storedSequenceTag*)avl_a)->startIndex,
		* n2 = globalNameBuffer->sData+((struct storedSequenceTag*)avl_b)->startIndex;
		
  long  l1 = ((struct storedSequenceTag*)avl_a)->length,
  l2 = ((struct storedSequenceTag*)avl_b)->length,
  i;
  
  signed char c;
  
  for (i = 0; i < l1 && i < l2; i++)
  {
    c = n1[i] - n2[i];
    if (c < 0)
      return -1;
    if (c>0)
      return 1;
  }
  
  if (l1 < l2)
    return -1;
  if (l1 > l2)
    return 1;
  
  return 0;
}


/*---------------------------------------------------------------------------------------------------- */

int main (int argc, const char * argv[])
{
		
  struct	bufferedString			**currentBuffers;
  
  struct  storedIDTag				*aTag = allocateIDTag(),
  *aTag2;
  
  struct  storedSequenceTag		*sTag = allocateStringTag(),
  *sTag2,
  *sTag3;
  
  struct  avl_traverser			avlt;
  
  struct  treeNode				*currentParent = NULL,
  *currentNode;
  
  char	automatonState			= 0,
  currentField			= 0,
  currentChar				= 0,
  showEmptyNodes			= 0;
  
  long	expectedFields			= NUMBER_OF_FIELDS,
  indexer,
  indexer2,
  indexer3,
  maxTreeLevel			= 0,
  setEOF					= 0,
  nRunCounter				= 0,
  lineHitCount            = 0;
  
  FILE  	* treeFile				= NULL,
  * summaryFile			= NULL,
  * inFile				= NULL;
  
  globalNameBuffer  = allocateNewString();
  globalTreeRoot	  = allocateNewTreeNode();
  currentBuffers    = (struct	bufferedString**)malloc (expectedFields*sizeof (struct	bufferedString*));
  byLevel			  = (struct avl_table**)malloc (sizeof (struct avl_table*) * NUMBER_OF_TAX_FIELDS);
  check_pointer	  (uniqueTags		  = avl_create (compare_tags, NULL, NULL));
  check_pointer	  (idTagAVL		      = avl_create (compare_id_tags, NULL, NULL));
  
  for (indexer = 0; indexer < NUMBER_OF_TAX_FIELDS; indexer++)
    check_pointer (byLevel[indexer] = avl_create (compare_tags, NULL, NULL));
		
  for (indexer = 0; indexer < expectedFields; indexer++)
    currentBuffers[indexer] = allocateNewString();
  
  if (argc != 6 && argc != 5)
  {
    fprintf (stderr,"%s",Usage);
    return 1;
  }
  
  maxTreeLevel = atoi (argv[2]);
		
  inFile		   = fopen (argv[1], "rb");
  if (!inFile)
  {
    fprintf (stderr, "Failed to open the input file %s\n", argv[1]);
    return 1;
  }
  
  treeFile		   = fopen (argv[3], "wb");
  if (!treeFile)
  {
    fprintf (stderr, "Failed to open the tree output file %s\n", argv[3]);
    return 1;
  }
  
  summaryFile		   = fopen (argv[4], "wb");
  if (!summaryFile)
  {
    fprintf (stderr, "Failed to open the summary output file %s\n", argv[4]);
    return 1;
  }
  
  if (argc == 6)
    showEmptyNodes = atoi (argv[5]);
  
  currentChar  = fgetc(inFile);
  currentField = 0;
  while (setEOF < 2)
  {
    switch (automatonState)
    {
      case 0: /* start of the line; expecting numbers or charcters */
        if (isalnum(currentChar))
        {
          automatonState = 1; /* reading sequence ID */
          appendCharacterToString(currentBuffers[currentField],toupper(currentChar));
        }
        else
          if (!(currentChar == '\n' || currentChar == '\r'))
          {
            reportSkippedLine ("Could not find a valid gid number to start the line",currentLineID);
            automatonState = 6;
            continue;
          }
        break;
        
      case 1: /* reading sequence ID */
        if (currentChar == '\t')
        {
          automatonState = 2;
          continue;
        }
        else
        {
          if (currentChar == '\n' || currentChar == '\r')
          {
            reportSkippedLine ("Expected a tab following the gid",currentLineID);
            automatonState = 6;
            continue;
          }
          else
            appendCharacterToString(currentBuffers[currentField],toupper(currentChar));
          
        }
        break;
        
      case 2: /* looking for a \t or a \n|\r*/
        if (currentChar == '\t')
        {
          currentField   ++;
          if (currentField < expectedFields)
            automatonState = 3;
          
          //reportSkippedLine ("Too many fields",currentLineID);
          //automatonState = 6;
          //continue;
        }
        else
          if (currentChar == '\n' || currentChar == '\r')
            // finish the read
          {
            if (currentField < expectedFields-1)
              reportSkippedLine ("Too few fields",currentLineID);
            else {
              aTag->taxID     = atoi (currentBuffers[1]->sData);
              lineHitCount    = atoi (currentBuffers[0]->sData);
              aTag2 = *(struct storedIDTag**)avl_probe(idTagAVL, aTag);
              //fprintf (stdout, "%d %d %x %x\n", currentLineID,aTag->taxID,aTag,aTag2);
              if (aTag == aTag2) // new taxID
              {
                // process fields
                currentParent = globalTreeRoot;
                nRunCounter   = 0;
                for (indexer = NUMBER_OF_FIELDS-NUMBER_OF_TAX_FIELDS; indexer < NUMBER_OF_FIELDS; indexer++)
                {
                  indexer2 = strlen(currentBuffers[indexer]->sData);
                  //fprintf (stdout, "%d %d\n", indexer, nRunCounter,indexer2);
                  if ((currentBuffers[indexer])->sData[0] == 'n' && indexer2 == 1)
                    nRunCounter++;
                  else
                  {
                    indexer3 = globalNameBuffer->sLength;
                    appendRangeToString (globalNameBuffer,currentBuffers[indexer],0,indexer2-1);
                    sTag->startIndex = indexer3;
                    sTag->length	 = indexer2;
                    sTag2			 = *(struct storedSequenceTag**)avl_probe (uniqueTags,sTag);
                    if (sTag == sTag2) // added
                      sTag = allocateStringTag();
                    else
                      globalNameBuffer->sLength = indexer3;
                    
                    sTag->startIndex = sTag2->startIndex;
                    sTag->length	 = sTag2->length;
                    
                    sTag3 = *(struct storedSequenceTag**)avl_probe (byLevel[indexer-(NUMBER_OF_FIELDS-NUMBER_OF_TAX_FIELDS)],sTag);
                    if (sTag == sTag3) // new node
                    {
                      //fprintf (stderr, "Add node level %d %d\n",indexer-(NUMBER_OF_FIELDS-NUMBER_OF_TAX_FIELDS), sTag2->startIndex);
                      sTag3 = sTag;
                      sTag3->refNode					= allocateNewTreeNode();
                      sTag3->refNode->length			= sTag2->length;
                      sTag3->refNode->startIndex		= sTag2->startIndex;
                      sTag3->refNode->parent			= NULL;
                      sTag = allocateStringTag();
                    }
                    //else
                    //fprintf (stderr, "Have node %x %x %c%c %d\n",sTag3->refNode, (currentBuffers[indexer])->sData[0],(currentBuffers[indexer])->sData[1], currentParent->children->vLength);
                    
                    sTag3->refNode->hitCount += lineHitCount;
                    
                    if (showEmptyNodes && nRunCounter>0)
                    {
                      currentNode             = allocateNewTreeNode();
                      currentNode->startIndex = -nRunCounter;
                      currentNode->hitCount	= lineHitCount;
                      if (currentParent)
                        currentParent = addAChild (currentParent,currentNode,1);
                      else
                        reportError ("Attempting to attach an empty node a null parent.");
                      
                    }
                    
                    if (currentParent)
                      currentParent = addAChild (currentParent,sTag3->refNode,0);
                    else
                      currentParent = sTag3->refNode;
                    
                    nRunCounter = 0;
                  }
                }
                aTag->tNode = currentParent;
                aTag = allocateIDTag();
              }
              
            else
            {
              aTag2->hit_count+=lineHitCount;
              currentParent = aTag2->tNode;
              while (currentParent)
              {
                currentParent->hitCount+=lineHitCount;
                currentParent = currentParent->parent;
              }
              //printf ("Duplicate tag %d %d\n", aTag2->hit_count, aTag2->taxID);
            }
            //traverseCheck (globalTreeRoot);
            }
            automatonState = 5;
            continue;
          }
          else
          {
            if (currentField < expectedFields)
            {
              reportSkippedLine ("Expected a tab following a field",currentLineID);
              automatonState = 6;
              continue;
            }
          }
        break;
        
      case 3: /* read a field */
        if (currentChar == '\t')
        {
          automatonState = 2;
          continue;
        }
        else
          if (currentChar == '\n' || currentChar == '\r')
          {
            currentField++;
            automatonState = 2;
            continue;
          }
          else
          {
            if (currentChar == '\'')
              automatonState = 4;
            else
              appendCharacterToString(currentBuffers[currentField],currentChar);
          }
        break;
        
      case 4: /* inside ' ' */
        if (currentChar == '\'')
          automatonState = 3;
        else
          if (currentChar == '\n' || currentChar == '\r')
          {
            reportSkippedLine ("Premature line end while reading a literal",currentLineID);
            ungetc (currentChar,inFile);
            automatonState = 5;
          }
          else
            appendCharacterToString(currentBuffers[currentField],currentChar);
        break;
        
      case 5:
        automatonState = 0;
        currentLineID ++;
        currentField = 0;
        for (indexer = 0; indexer < expectedFields; indexer++)
          clear_buffered_string(currentBuffers[indexer]);
        break;
        
      case 6:
        if (currentChar == '\n' || currentChar == '\r')
          automatonState = 5;
        break;
    }
    currentChar = fgetc(inFile);
    if (feof (inFile))
    {
      setEOF		++;
      currentChar = '\n';
    }
  }
  
  fprintf (stderr, "Read %ld unique taxIDs\n", (long)avl_count (idTagAVL));
  
  fclose  (inFile);
  
  indexer3 = ((maxTreeLevel>0)?(maxTreeLevel+1):0xfffffL);
  if (indexer3 > NUMBER_OF_TAX_FIELDS)
    indexer3 = NUMBER_OF_TAX_FIELDS;
  for (indexer = 0; indexer < indexer3; indexer++)
  {
    avl_t_init (&avlt, byLevel[indexer]);
    sTag = (struct storedSequenceTag*)avl_t_first (&avlt, byLevel[indexer]);
    while (sTag)
    {
      fprintf (summaryFile, "%s\t", rankLabels[indexer]);
      for (indexer2 = 0; indexer2 < sTag->length; indexer2++)
        fputc(*(globalNameBuffer->sData+sTag->startIndex+indexer2),summaryFile);
      fprintf (summaryFile, "\t%ld\n", sTag->refNode->hitCount);
      sTag = (struct storedSequenceTag*)avl_t_next (&avlt);
    }
    //fprintf (stderr, "Level %d unique IDs = %d\n", indexer, avl_count (byLevel[indexer]));
  }
  
  for (indexer = 0; indexer < 256; indexer ++)
    validTaxonNameChar[indexer] = 1;
		
  validTaxonNameChar[','] = 0;
  validTaxonNameChar[')'] = 0;
  validTaxonNameChar['('] = 0;
  validTaxonNameChar[':'] = 0;
  
  
  traverseTree (treeFile,((struct treeNode**)globalTreeRoot->children->vData)[0],
                (maxTreeLevel>0)?(maxTreeLevel):0xfffffL,
                0);
  fclose (treeFile);
  return 0;
}

