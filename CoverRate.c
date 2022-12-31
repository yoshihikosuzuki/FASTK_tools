#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>

#include "libfastk.h"

static char *Usage = "[-l<int(1)>] [-r<int(32767)>] <source_root>[.prof]";

int main(int argc, char *argv[])
{ Profile_Index *P;
  int MULT_LOW, MULT_HGH;

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    (void) flags;

    ARG_INIT("CoverRate");

    MULT_LOW = 1;
    MULT_HGH = 32767;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("")
            break;
          case 'l':
            ARG_POSITIVE(MULT_LOW,"Minimum k-mer multiplicity")
            break;
          case 'r':
            ARG_POSITIVE(MULT_HGH,"Maximum k-mer multiplicity")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -l: Minimum multiplicity of k-mers to be counted\n");
        fprintf(stderr,"      -r: Maximum multiplicity of k-mers to be counted\n");
        exit (1);
      }
  }

  P = Open_Profiles(argv[1]);
  if (P == NULL)
    { fprintf(stderr,"%s: Cannot open %s\n",Prog_Name,argv[1]);
      exit (1);
    }

  { int64   id;
    char   *eptr;
    uint16 *profile;
    int     plen, tlen;
    int     count;

    plen    = 20000;
    profile = Malloc(plen*sizeof(uint16),"Profile array");

    printf("RID\tN_TOT_KMER\tN_FOCAL_KMER\tCOV_RATE\n");
    for (id = 1; id <= P->nreads; id++)
      { tlen = Fetch_Profile(P,id-1,plen,profile);
        if (tlen > plen)
          { plen    = 1.2*tlen + 1000;
            profile = Realloc(profile,plen*sizeof(uint16),"Profile array");
            Fetch_Profile(P,id-1,plen,profile);
          }
        count = 0;
        for (int i = 0; i < tlen; i++)
          if (MULT_LOW <= profile[i] && profile[i] <= MULT_HGH) count++;
        printf("%lld\t%d\t%d\t%.3lf\n", id, tlen, count, (double)count/tlen);
      }
    free(profile);
  }

  Free_Profiles(P);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}
