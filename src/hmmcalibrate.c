/************************************************************
 * HMMER - Biological sequence analysis with profile-HMMs
 * Copyright (C) 1992-1998 Washington University School of Medicine
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* hmmcalibrate.c
 * SRE, Fri Oct 31 09:25:21 1997 [St. Louis]
 * 
 * Score an HMM against random sequence data sets;
 * set histogram fitting parameters.
 * 
 * RCS $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>

#ifdef HMMER_THREADS
#include <pthread.h>
#endif

#include "squid.h"
#include "config.h"
#include "structs.h"
#include "funcs.h"
#include "version.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

#include "globals.h"

static char banner[] = "hmmcalibrate -- calibrate HMM search statistics";

static char usage[] = "\
Usage: hmmcalibrate [-options] <hmmfile>\n\
Available options are:\n\
  -h             : print short usage and version info, then exit\n\
";

static char experts[] = "\
  --benchmark    : run hmmcalibrate in benchmarking mode\n\
  --cpu <n>      : run <n> threads in parallel (if threaded)\n\
  --fixed <n>    : fix random sequence length at <n>\n\
  --histfile <f> : save histogram(s) to file <f>\n\
  --mean <x>     : set random seq length mean at <x> [350]\n\
  --num <n>      : set number of sampled seqs to <n> [5000]\n\
  --sd <x>       : set random seq length std. dev to <x> [350]\n\
  --seed <n>     : set random seed to <n> [time()]\n\
";

static struct opt_s OPTIONS[] = {
   { "-h", TRUE, sqdARG_NONE  },
   { "--benchmark",FALSE, sqdARG_NONE }, 
   { "--cpu",      FALSE, sqdARG_INT },
   { "--fixed",    FALSE, sqdARG_INT   },
   { "--histfile", FALSE, sqdARG_STRING },
   { "--mean",     FALSE, sqdARG_FLOAT },
   { "--num",      FALSE, sqdARG_INT   },
   { "--sd",       FALSE, sqdARG_FLOAT },   
   { "--seed",     FALSE, sqdARG_INT}, 
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv)
{
  char    *hmmfile;             /* HMM file to open                */
  char    *tmpfile;             /* temporary calibrated HMM file   */
  HMMFILE *hmmfp;               /* opened hmm file pointer         */
  FILE    *outfp;               /* for writing HMM(s) into tmpfile */
  char    *mode;                /* write mode, "w" or "wb"         */
  struct plan7_s     *hmm;      /* the hidden Markov model         */
  struct histogram_s *hist;	/* score histogram                 */
  int     idx;			/* counter over sequences          */
  char   *seq;			/* a random sequence               */
  char   *dsq;			/* seq, digitized for alignment    */
  float   randomseq[MAXABET];	/* random sequence model           */
  float   p1;			/* random sequence model p1        */
  float   score;		/* score of an alignment           */
  float   max;			/* maximum score                   */
  int     sqlen;		/* length of sampled sequences     */
  sigset_t blocksigs;		/* list of signals to protect from */
  int     fitok;		/* TRUE if the EVD fit went ok     */

  int     nsample;		/* number of random seqs to sample */
  int     seed;			/* random number seed              */
  int     fixedlen;		/* fixed length, or 0 if unused    */
  float   lenmean;		/* mean of length distribution     */
  float   lensd;		/* std dev of length distribution  */
  int     do_benchmark;		/* TRUE to go into benchmark mode  */
  char   *histfile;             /* histogram save file             */
  FILE   *hfp;                  /* open file pointer for histfile  */
  time_t  sttime, endtime;	/* start, end time for run         */
  clock_t stcpu, endcpu;	/* start, end CPU time for run     */

  char *optname;		/* name of option found by Getopt() */
  char *optarg;			/* argument found by Getopt()       */
  int   optind;		        /* index in argv[]                  */

#ifdef HMMER_THREADS
  struct vpool_s *vpool;        /* pool of worker threads  */
  char           *odsq;         /* digitized seq returned by thread to be free'd */
#endif
  int   num_threads;            /* number of worker threads */   

#ifdef MEMDEBUG
  unsigned long histid1, histid2, orig_size, current_size;
  orig_size = malloc_inuse(&histid1);
  fprintf(stderr, "[... memory debugging is ON ...]\n");
#endif

  /***********************************************
   * Parse the command line
   ***********************************************/

  nsample      = 5000;
  fixedlen     = 0;
  lenmean      = 325.;
  lensd        = 200.;
  seed         = (int) time ((time_t *) NULL);
  do_benchmark = FALSE;
  histfile     = NULL;
#ifdef HMMER_THREADS
  num_threads  = ThreadNumber(); /* only matters if we're threaded */
#endif

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
		&optind, &optname, &optarg))
    {
      if      (strcmp(optname, "--benchmark")== 0) do_benchmark = TRUE;
      else if (strcmp(optname, "--cpu")      == 0) num_threads  = atoi(optarg);
      else if (strcmp(optname, "--fixed")    == 0) fixedlen = atoi(optarg);
      else if (strcmp(optname, "--histfile") == 0) histfile = optarg;
      else if (strcmp(optname, "--mean")     == 0) lenmean  = atof(optarg); 
      else if (strcmp(optname, "--num")      == 0) nsample  = atoi(optarg); 
      else if (strcmp(optname, "--sd")       == 0) lensd    = atof(optarg); 
      else if (strcmp(optname, "--seed")     == 0) seed     = atoi(optarg);
      else if (strcmp(optname, "-h") == 0)
	{
	  Banner(stdout, banner);
	  puts(usage);
	  puts(experts);
	  exit(0);
	}
    }

  if (argc - optind != 1) Die("Incorrect number of arguments.\n%s\n", usage);
  hmmfile = argv[optind++];

  sre_srandom(seed);

  /* Benchmark mode: timing of HMMER on different platforms/compilers.
   * overrides several defaults and options. 
   * Sets --seed 0, --fixed 300, --num 5000.
   * Will not show banner; shows an alternative output that
   * is the benchmark result. Does not save results in hmm file.
   * Expects to be called with the Demos/rrm.hmm file.
   */
  if (do_benchmark) 
    {
      fixedlen = 300;
      nsample  = 5000;
      seed     = 0;
      sttime   = time(NULL);
      stcpu    = clock();
      printf("HMMER Viterbi benchmark \n");
      printf("  (will take about 1-4 minutes...)\n"); 
    }

  /***********************************************
   * Open our i/o file pointers, make sure all is well
   ***********************************************/

  if ((hmmfp = HMMFileOpen(hmmfile, "HMMERDB")) == NULL)
    Die("failed to open HMM file %s for reading.", hmmfile);

  /* Generate calibrated HMM(s) in a tmp file in the current
   * directory. When we're finished, we delete the original
   * HMM file and rename() this one. That way, the worst
   * effect of a catastrophic failure should be that we
   * leave a tmp file lying around, but the original HMM
   * file remains uncorrupted. tmpnam() doesn't work here,
   * because it'll put the file in /tmp and we won't
   * necessarily be able to rename() it from there.
   */
  tmpfile = MallocOrDie(strlen(hmmfile) + 5);
  strcpy(tmpfile, hmmfile);
  strcat(tmpfile, ".xxx");	/* could be more inventive here... */
  if (FileExists(tmpfile))
    Die("temporary file %s already exists; please delete it first", tmpfile);
  if (hmmfp->is_binary) mode = "wb";
  else                  mode = "w"; 
  if ((outfp = fopen(tmpfile, mode)) == NULL)
    Die("temporary file %s couldn't be opened for writing", tmpfile); 

				/* histogram file */
  if (histfile != NULL)
    {
      if ((hfp = fopen(histfile, "w")) == NULL)
	Die("Failed to open histogram save file %s for writing\n", histfile);
    }

  /*********************************************** 
   * Show the banner
   ***********************************************/

  if (! do_benchmark) 
    {
      Banner(stdout, banner);
      printf("HMM file:                 %s\n", hmmfile);
      if (fixedlen) 
	printf("Length fixed to:          %d\n", fixedlen);
      else {
	printf("Length distribution mean: %.0f\n", lenmean);
	printf("Length distribution s.d.: %.0f\n", lensd);
      }
      printf("Number of samples:        %d\n", nsample);
      printf("random seed:              %d\n", seed);
      printf("histogram(s) saved to:    %s\n",
	     histfile != NULL ? histfile : "[not saved]");
      printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
    }

#ifdef HMMER_THREADS
  vpool = VpoolInit(FALSE, FALSE, num_threads, num_threads*2);
#endif

  /***********************************************
   * Calibrate each model in turn
   ***********************************************/

  while (HMMFileRead(hmmfp, &hmm)) 
    {	
      if (hmm == NULL) 
	Die("HMM file %s may be corrupt or in incorrect format; parse failed", hmmfile);
      P7Logoddsify(hmm, TRUE);
				/* we could use the null model in the HMM? */
      P7DefaultNullModel(randomseq, &p1);
      
      hist = AllocHistogram(-200, 200, 100);
      max = -FLT_MAX;
				/* deliberately generating one too many,
				   because of how we handle the threads */
      idx = 0;
      for (;;)
	{
	  if (idx < nsample) {
				/* choose length of random sequence */
	    if (fixedlen) sqlen = fixedlen;
	    else do sqlen = (int) Gaussrandom(lenmean, lensd); while (sqlen < 1);
				/* generate it */
	    seq = RandomSequence(Alphabet, randomseq, Alphabet_size, sqlen);
	    dsq = DigitizeSequence(seq, sqlen);
	  }

#ifdef HMMER_THREADS
	  /* If we're done preloading num_threads worth of input, 
	   * then block while waiting for output
	   */
	  if (vpool->shutdown || idx >= num_threads)
	    {
	      SQD_DPRINTF2(("Boss has %d in hand; looking for output\n", idx));
	      if (!VpoolGetResults(vpool, NULL, &odsq, NULL, NULL, &score, NULL))
		break;		/* no more output */
	    }
				/* add work, or shut down when finished */
	  if (idx < nsample) 
	    {
	      SQD_DPRINTF2(("Boss adding sequence %d as work\n", idx));
	      VpoolAddWork(vpool, hmm, dsq, NULL, sqlen);
	      free(seq);
	    }
	  else if (idx == nsample) VpoolShutdown(vpool);

				/* process  results */
	  if (vpool->shutdown || idx >= num_threads)
	    {
	      AddToHistogram(hist, score);
	      if (score > max) max = score;
	      free(odsq); 
	    }
#else /* unthreaded version */

	  if (idx == nsample) break; /* terminate serial version */

	  if (P7ViterbiSize(sqlen, hmm->M) <= RAMLIMIT)
	    score = P7Viterbi(dsq, sqlen, hmm, NULL);
	  else
	    score = P7SmallViterbi(dsq, sqlen, hmm, NULL);

	  AddToHistogram(hist, score);
	  if (score > max) max = score;
	  free(dsq); 
	  free(seq);
#endif

	  idx++;
	}


      /* Fit an EVD to the observed histogram.
       * The TRUE left-censors and fits only the right slope of the histogram.
       * The 9999. is an arbitrary high number that means we won't trim outliers
       * on the right.
       */
      fitok = ExtremeValueFitHistogram(hist, TRUE, 9999.);

      /* Set HMM EVD parameters 
       */
      if (fitok)
	{
	  hmm->mu      = hist->param[EVD_MU];
	  hmm->lambda  = hist->param[EVD_LAMBDA];
	  hmm->flags  |= PLAN7_STATS;
	}
      else
	printf(" -- fit failed; -n may be set too small?\n");

      /* Record command line in comlog
       */
      Plan7ComlogAppend(hmm, argc, argv);

      /* Save HMM to tmpfile
       */
      if (hmmfp->is_binary) WriteBinHMM(outfp, hmm);
      else                  WriteAscHMM(outfp, hmm); 

      /* Output results
       */
      if (! do_benchmark)
	{
	  printf("HMM    : %s\n", hmm->name);
	  if (fitok)
	    {
	      printf("mu     : %12f\n", hmm->mu);
	      printf("lambda : %12f\n", hmm->lambda);
	    }
	  else
	    {
	      printf("mu     : [undetermined]\n");
	      printf("lambda : [undetermined]\n");
	    }
	  printf("max    : %12f\n", max);
	  printf("//\n");
	}

      if (histfile != NULL) 
	{
	  fprintf(hfp, "HMM: %s\n", hmm->name);
	  PrintASCIIHistogram(hfp, hist);
	  fprintf(hfp, "//\n");
	}

      FreePlan7(hmm);
      FreeHistogram(hist);
    }

  /* Now, carefully remove original file and replace it
   * with the tmpfile. Note the protection from signals;
   * we wouldn't want a user to ctrl-C just as we've deleted
   * their HMM but before the new one is moved.
   */
  HMMFileClose(hmmfp);
  if (fclose(outfp)   != 0) PANIC;

  if (! do_benchmark) 
    {
      if (sigemptyset(&blocksigs) != 0) PANIC;
      if (sigaddset(&blocksigs, SIGINT) != 0) PANIC;
      if (sigprocmask(SIG_BLOCK, &blocksigs, NULL) != 0)   PANIC;
      if (remove(hmmfile) != 0)                            PANIC;
      if (rename(tmpfile, hmmfile) != 0)                   PANIC;
      if (sigprocmask(SIG_UNBLOCK, &blocksigs, NULL) != 0) PANIC;
    }
  else
    {
      remove(tmpfile);
      endtime = time(NULL);
      endcpu  = clock();
    }

  /***********************************************
   * Benchmark output.
   ***********************************************/

  /* Benchmark output assumes that rrm.hmm (72 positions)
   * was used: 72 x 300 x 5000 = 108 million cells
   * 264.922 sec is a baseline for woozle, an SGI Indy R4600PC/133.
   */
  if (do_benchmark)
    {
      double clocksec;
      double cpusec;

      if (sttime == -1 || endtime == -1 || stcpu == -1 || endcpu == -1)
	Die("Time was not available for some reason.");

      clocksec = difftime(endtime,sttime);
      cpusec   = (double) (endcpu - stcpu) / (double) CLOCKS_PER_SEC;

      printf("   Wall clock: %.1f sec\n", clocksec);
      printf("   CPU time:   %.1f sec\n", cpusec);
      printf("   Mcells/sec: %.3f\n", 108./cpusec);
      printf("   Relative:   %.1f\n", 264.922/cpusec);
    }

  /***********************************************
   * Exit
   ***********************************************/
  free(tmpfile);
  SqdClean();
#ifdef MEMDEBUG
  current_size = malloc_size(&histid2);
  if (current_size != orig_size)
    malloc_list(2, histid1, histid2);
  else
    fprintf(stderr, "[No memory leaks]\n");
#endif

  return 0;
}


