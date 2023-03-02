# ifndef __LAMA_RUNTIME__
# define __LAMA_RUNTIME__

# ifndef __cplusplus

# include <stdio.h>
# include <stdio.h>
# include <string.h>
# include <stdarg.h>
# include <stdlib.h>
# include <sys/mman.h>
# include <assert.h>
# include <errno.h>
# include <regex.h>
# include <time.h>
# include <limits.h>
# include <ctype.h>

# endif // __cplusplus

# ifdef __cplusplus
extern "C"
{
# endif // __cplusplus

# define UNBOXED(x)  (((int) (x)) &  0x0001)
# define UNBOX(x)    (((int) (x)) >> 1)
# define BOX(x)      ((((int) (x)) << 1) | 0x0001)

# define STRING_TAG  0x00000001
# define ARRAY_TAG   0x00000003
# define SEXP_TAG    0x00000005
# define CLOSURE_TAG 0x00000007
# define UNBOXED_TAG 0x00000009 // Not actually a tag; used to return from LkindOf

# define WORD_SIZE (CHAR_BIT * sizeof(int))

void failure (char *s, ...);
void printValue (void *p);

int Lread ();
int Lwrite (int n);

void* LmakeArray (int length);

int Barray_patt (void *d, int n);

void* Belem (void *p, int i);

void* Bsta (void *v, int i, void *x);

int Llength (void *p);

void* BsexpTag (int bn, int tag);

int Btag (void *d, int t, int n);

int LtagHash (char*);

# ifdef __cplusplus
}
# endif // __cplusplus

# endif
