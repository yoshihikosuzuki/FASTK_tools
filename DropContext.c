#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include "libfastk.h"
#include "DB.h"

static char *Usage = "-e<int> -g<int>:<int> <source_root>";

int main(int argc, char *argv[])
{ Profile_Index *P;
  DAZZ_DB        _db, *db = &_db;
  char          *seq, *pseq;
  int            nreads;
  int            ERROR, GOOD_LOW, GOOD_HGH;

  { int    i, j, k;
    int    flags[128];
    char  *eptr, *fptr;

    (void) flags;

    ARG_INIT("DropContext");

    ERROR    = -1;
    GOOD_LOW = -1;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("")
            break;
          case 'e':
            ERROR = strtol(argv[i]+2,&eptr,10);
            if (eptr > argv[i]+2 && *eptr == '\0')
              { if (ERROR < 1 || ERROR > 0x7fff)
                  { fprintf(stderr,"%s: Error threshold %d is out of range\n",
                                   Prog_Name,ERROR);
                    exit (1);
                  }
                break;
              }
            fprintf(stderr,"%s: Syntax of -e option invalid -e<int>\n",Prog_Name);
            exit (1);
          case 'g':
            GOOD_LOW = strtol(argv[i]+2,&eptr,10);
            if (eptr > argv[i]+2)
              { if (GOOD_LOW < 1 || GOOD_LOW > 0x7fff)
                  { fprintf(stderr,"%s: Minimum valid count %d is out of range\n",
                                   Prog_Name,GOOD_LOW);
                    exit (1);
                  }
                if (*eptr == ':')
                  { GOOD_HGH = strtol(eptr+1,&fptr,10);
                    if (fptr > eptr+1 && *fptr == '\0')
                      { if (GOOD_HGH < 1 || GOOD_HGH > 0x7fff)
                          { fprintf(stderr,"%s: Maximum valid count %d is out of range\n",
                                           Prog_Name,GOOD_HGH);
                            exit (1);
                          }
                        if (GOOD_LOW > GOOD_HGH)
                          { fprintf(stderr,"%s: Good count range is invalid\n",Prog_Name);
                            exit (1);
                          }
                        break;
                      }
                  }
              }
            fprintf(stderr,"%s: Syntax of -g option invalid -g<int>:<int>\n",Prog_Name);
            exit (1);
        }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc != 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -e: Counts <= this value are considered errors.\n");
        fprintf(stderr,"      -g: Counts in this range are considered correct.\n");
        exit (1);
      }

    if (ERROR < 0)
      { fprintf(stderr,"%s: Must give error count threshold -e\n",Prog_Name);
        exit (1);
      }
    if (GOOD_LOW < 0)
      { fprintf(stderr,"%s: Must give good count range -g\n",Prog_Name);
        exit (1);
      }
  }

  //  Open DB or DAM, and if a DAM open also .hdr file
  { int status;

    status = Open_DB(argv[1],db);
    if (status < 0)
      { fprintf(stderr,"%s: Cannot open %s.db\n",Prog_Name,argv[1]);
        exit (1);
      }
    if (db->part > 0)
      { fprintf(stderr,"%s: Cannot be called on a block\n",Prog_Name);
        exit (1);
      }
    Trim_DB(db);
    seq = New_Read_Buffer(db);
    nreads = db->nreads;
  }

  // Load Profile for each read
  P = Open_Profiles(argv[1]);
  if (P == NULL)
    { fprintf(stderr,"%s: Cannot open %s.prof\n",Prog_Name,argv[1]);
      exit (1);
    }
  if (nreads != P->nreads)
    { fprintf(stderr,"# of reads must be same\n");
      exit(1);
    }

  { int     id;
    uint16 *profile;
    int     plen, tlen;
    int     km1;
    int     homo, di, tri;

    km1 = P->kmer - 1;
    plen    = 20000;
    profile = Malloc(plen*sizeof(uint16),"Profile array");
    for (id = 1; id <= nreads; id++)
      { tlen = Fetch_Profile(P,(int64) id-1,plen,profile);
        if (tlen > plen)
          { plen    = 1.2*tlen + 1000;
            profile = Realloc(profile,plen*sizeof(uint16),"Profile array");
            Fetch_Profile(P,(int64) id-1,plen,profile);
          }

        // Ignore reads with repeats
        int flag_rep = 0;
        for (int i = 0; i < tlen; i++) {
          if (profile[i] > GOOD_HGH) {
            flag_rep = 1;
            break;
          }
        }
        if (flag_rep == 1) continue;

        Load_Read(db,id-1,seq,1);
        pseq = seq + km1;

        // NOTE: check if a drop happens at position i+1
        for (int i = 0; i < tlen - 1; i++)
          { if (GOOD_LOW <= profile[i] && profile[i] <= GOOD_HGH
                && GOOD_LOW <= profile[i+1] && profile[i+1] <= GOOD_HGH
                && 3 <= profile[i]-profile[i+1] && profile[i]-profile[i+1] <= 10)
              { homo = di = tri = 0;
                if (pseq[i-1] == pseq[i])
                  homo = 1;
                if (pseq[i-1] != pseq[i]
                    && pseq[i-3] == pseq[i-1] && pseq[i-2] == pseq[i])
                  di = 1;
                if ((pseq[i-2] != pseq[i-1] || pseq[i-1] != pseq[i])
                    && pseq[i-5] == pseq[i-2] && pseq[i-4] == pseq[i-1] && pseq[i-3] == pseq[i])
                  tri = 1;
                printf("Read %8d, @ %8d, %3d -> %3d (-%2d) ...%.*s|%.*s... %c %c %c\n",
                       id,km1+i+1,profile[i],profile[i+1],profile[i]-profile[i+1],
                       20,pseq+i-20+1,
                       20,pseq+i+1,
                       homo?'1':' ',
                       di?'2':' ',
                       tri?'3':' ');
              }
            
          }
      }
    free(profile);
  }

  Free_Profiles(P);

  Close_DB(db);

  Catenate(NULL,NULL,NULL,NULL);
  Numbered_Suffix(NULL,0,NULL);
  free(Prog_Name);

  exit (0);
}