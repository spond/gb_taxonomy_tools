

#define         DEFAULT_STRING_ALLOC	  16L


struct		bufferedString
{
	char *sData;
	long sLength,
		 saLength;
};

struct		vector
{
	long *vData;

	long vLength,
		 vaLength;
};



struct bufferedString *allocateNewString (void);
void   clear_buffered_string (struct bufferedString* theString);
void   appendCharacterToString (struct bufferedString * s, const char c);
long   appendRangeToString (struct bufferedString * d, struct bufferedString *s, long from, long to);
void   appendCharRangeToString (struct bufferedString * d, char * buffer);
long   appendCharBufferToString (struct bufferedString * d, const char * b);
int    compare_strings (const struct bufferedString * s1, const struct bufferedString * s2);
int    compare_string_and_char (const struct bufferedString * s1, const char * s2);
void   destroy_string (struct bufferedString* aStr);
void   appendValueToVector (struct vector * v, long c);


void   reportError (char * theMessage);
void   reportErrorLine (char * theMessage, long lineID);
void   check_pointer	(void*);
void   reportSkippedLine (char * theMessage, long lineID);

struct vector *allocateNewVector (void);
void   clear_vector (struct vector* v);
long   popValueFromVector (struct vector * v);

