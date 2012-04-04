#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include "buffers.h"

/*---------------------------------------------------------------------------------------------------- */

struct bufferedString *allocateNewString (void)
{
	struct bufferedString *newS = (struct bufferedString*)malloc (sizeof (struct bufferedString));
	check_pointer (newS);
	check_pointer (newS->sData = (char*)malloc (DEFAULT_STRING_ALLOC+1));
	newS->sLength  = 0;
	newS->saLength = DEFAULT_STRING_ALLOC;
	newS->sData[0] = 0;
	return newS;
}

/*---------------------------------------------------------------------------------------------------- */

void clear_buffered_string (struct bufferedString* theString)
{
	theString->sLength = 0;
}

/*---------------------------------------------------------------------------------------------------- */

void appendCharacterToString (struct bufferedString * s, const char c)
{
	long addThis; 
	if (s->saLength == s->sLength)
	{
		addThis = s->saLength / 8;
		if (DEFAULT_STRING_ALLOC > addThis)
			addThis = DEFAULT_STRING_ALLOC;
		s->saLength += addThis;
		check_pointer (s->sData = realloc (s->sData,s->saLength+1));
	}
	s->sData[s->sLength] = c;
	s->sData[++s->sLength] = 0;
}

/*---------------------------------------------------------------------------------------------------- */

long appendRangeToString (struct bufferedString * d, struct bufferedString *s, long from, long to)
{
	long addThis,
		 pl = to-from+1;
	
	if (pl<=0)
		return -1;
		
	if (d->saLength < d->sLength + pl)
	{
		addThis = d->saLength / 8;
		
		if (DEFAULT_STRING_ALLOC > addThis)
			addThis = DEFAULT_STRING_ALLOC;
		if (addThis < pl)
			addThis = pl;
			
		d->saLength += addThis;
		check_pointer (d->sData = realloc (d->sData,d->saLength+1));
	}
	for (addThis = from; addThis <=to; addThis++)
		d->sData[d->sLength++] = s->sData[addThis];

	d->sData[d->sLength] = 0;
	
	return pl;
}

/*---------------------------------------------------------------------------------------------------- */

void appendCharRangeToString (struct bufferedString * d, char * buffer)
{
	long addThis,
		 pl = strlen(buffer);
	
	if (pl<=0)
		return;
		
	if (d->saLength < d->sLength + pl)
	{
		addThis = d->saLength / 8;
		
		if (DEFAULT_STRING_ALLOC > addThis)
			addThis = DEFAULT_STRING_ALLOC;
		if (addThis < pl)
			addThis = pl;
			
		d->saLength += addThis;
		check_pointer (d->sData = realloc (d->sData,d->saLength+1));
	}
	for (addThis = 0; addThis <pl; addThis++)
		d->sData[d->sLength++] = buffer[addThis];

	d->sData[d->sLength] = 0;
	
}

/*---------------------------------------------------------------------------------------------------- */

long appendCharBufferToString (struct bufferedString * d, const char * b)
{
	long addThis,
		 pl = strlen (b);
	
	if (pl<=0)
		return -1;
		
	if (d->saLength < d->sLength + pl)
	{
		addThis = d->saLength / 8;
		
		if (DEFAULT_STRING_ALLOC > addThis)
			addThis = DEFAULT_STRING_ALLOC;
		if (addThis < pl)
			addThis = pl;
			
		d->saLength += addThis;
		check_pointer (d->sData = realloc (d->sData,d->saLength+1));
	}
	for (addThis = 0; addThis <pl; addThis++)
		d->sData[d->sLength++] = b[addThis];

	d->sData[d->sLength] = 0;
	return pl;
}

/*---------------------------------------------------------------------------------------------------- */

void	check_pointer (void * p)
{
    if (p == NULL)
    {
        fprintf (stderr,"Memory allocation error\n");
        exit (1);
    }
}

/*---------------------------------------------------------------------------------------------------- */

int		compare_strings (const struct bufferedString * s1, const struct bufferedString * s2)
{
	long upTo,
		 i;
		 
	if  (s1->sLength>s2->sLength)
		upTo = s2->sLength;
	else
		upTo = s1->sLength;

	for (i=0; i<upTo; i++)
	{
		int res = (s1->sData[i]-s2->sData[i]);
	 	if (res < 0)
	 		return -1;
	 	else
	 		if (res>0)
	 			return 1;
	}
	
	if (s1->sLength == s2->sLength)
		return 0;

	return 1-2*(s1->sLength<s2->sLength);
}

/*---------------------------------------------------------------------------------------------------- */

int		compare_string_and_char (const struct bufferedString * s1, const char * s2)
{
	long upTo,
		 i,
         s2l = strlen (s2);
		 
	if  (s1->sLength>s2l)
		upTo = s2l;
	else
		upTo = s1->sLength;

	for (i=0; i<upTo; i++)
	{
		int res = (s1->sData[i]-s2[i]);
	 	if (res < 0)
	 		return -1;
	 	else
	 		if (res>0)
	 			return 1;
	}
	
	if (s1->sLength == s2l)
		return 0;

	return 1-2*(s1->sLength<s2l);
}

/*---------------------------------------------------------------------------------------------------- */

void destroy_string (struct bufferedString* aStr)
{
	free (aStr->sData);
	free (aStr);
}


/*---------------------------------------------------------------------------------------------------- */

void reportError (char * theMessage)
{
	fprintf (stderr, "\nERROR: %s\n", theMessage);
	exit (1);
}

/*---------------------------------------------------------------------------------------------------- */

void reportErrorLine (char * theMessage, long lineID)
{
	fprintf (stderr, "\nERROR in line %ld: %s\n", lineID, theMessage);
	exit (1);
}

/*---------------------------------------------------------------------------------------------------- */

void reportSkippedLine (char * theMessage, long lineID)
{
	fprintf (stderr, "SKIPPED line %ld: %s\n", lineID, theMessage);
}

/*---------------------------------------------------------------------------------------------------- */

struct vector *allocateNewVector (void)
{
	struct vector *newS = (struct vector*)malloc (sizeof (struct vector));
	check_pointer (newS);
	check_pointer (newS->vData = (long*)malloc (DEFAULT_STRING_ALLOC*sizeof(long)));
	newS->vLength  = 0;
	newS->vaLength = DEFAULT_STRING_ALLOC;
	return newS;
}

/*---------------------------------------------------------------------------------------------------- */

void clear_vector (struct vector* v)
{
	v->vLength = 0;
}

/*---------------------------------------------------------------------------------------------------- */

void appendValueToVector (struct vector * v, long c)
{
	long addThis; 
	if (v->vaLength == v->vLength)
	{
		addThis = v->vaLength / 8;
		if (DEFAULT_STRING_ALLOC > addThis)
			addThis = DEFAULT_STRING_ALLOC;
		v->vaLength += addThis;
		check_pointer (v->vData = realloc (v->vData,sizeof(long)*v->vaLength));
	}
	v->vData[v->vLength++] = c;
}

/*---------------------------------------------------------------------------------------------------- */

long popValueFromVector (struct vector * v)
{
	if (v->vLength)
		return v->vData[--v->vLength];

	return 0;
}
