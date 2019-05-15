/* probify: reverse engineer scorematrix to get implicit probabilistic model.
 * This may fail! Essentially copied from esl_scorematrix's Example program.
 * Added option of providing bg probabilities
 */

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_fileparser.h"
#include "esl_getopts.h"
#include "esl_dmatrix.h"
#include "esl_vectorops.h"
#include "esl_scorematrix.h"
#include "esl_scorematrix.h"
#include "esl_composition.h"
#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name             type          default  env  range    toggles          reqs incomp  help                                       docgroup*/
  { "-h",          eslARG_NONE,       FALSE,  NULL, NULL,  NULL,             NULL, NULL, "show brief help on version and usage",        0 },
  { "--dna",       eslARG_NONE,       FALSE,  NULL, NULL,  "--dna,--amino",  NULL, NULL, "use DNA alphabet",                            0 },
  { "--amino",     eslARG_NONE,      "TRUE",  NULL, NULL,  "--dna,--amino",  NULL, NULL, "use protein alphabet",                        0 },
  { "--bg_file",   eslARG_INFILE,      NULL,  NULL, NULL,  NULL,             NULL, NULL, "name the (optional) background file <s>",                0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <mxfile>";
static char banner[] = "estimate parameters underlying score matrix";


int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go        = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char            *scorefile = esl_opt_GetArg(go, 1);
  ESL_ALPHABET    *abc       = NULL;
  ESL_FILEPARSER  *efp       = NULL;
  ESL_SCOREMATRIX *S         = NULL;
  ESL_DMATRIX     *P1        = NULL; /* implicit probability basis, bg unknown */
  ESL_DMATRIX     *P2        = NULL; /* implicit probability basis, bg known   */
  double          *fi        = NULL;
  double          *fj        = NULL;
  double           lambda, D, E;
  int              vstatus, status;

  if      (esl_opt_GetBoolean(go, "--dna"))   abc = esl_alphabet_Create(eslDNA);
  else if (esl_opt_GetBoolean(go, "--amino")) abc = esl_alphabet_Create(eslAMINO);

  /* Input a score matrix from a file. */
  if ( esl_fileparser_Open(scorefile, NULL, &efp) != eslOK) esl_fatal("failed to open score file %s",         scorefile);
  if ( esl_scorematrix_Read(efp, abc, &S)         != eslOK) esl_fatal("failed to read matrix from %s:\n  %s", scorefile, efp->errbuf);
  esl_fileparser_Close(efp);

  /* Try to reverse engineer it to get implicit probabilistic model. This may fail! */
  if (esl_opt_IsOn(go, "--bg_file")) {
      P7_BG  *bg    = NULL;
      bg            = p7_bg_Create(abc);
      vstatus = p7_bg_Read(esl_opt_GetString(go, "--bg_file"), bg, NULL);
      if (vstatus != eslOK) p7_Fail("Trouble reading bgfile\n");


      ESL_ALLOC(fi, sizeof(double)*4);
      esl_vec_F2D(bg->f,4,fi);

      vstatus = esl_scorematrix_ProbifyGivenBG(S, fi, fi, &lambda, &P1  );

      if (vstatus == eslOK)
      { /* Print some info, and the joint probabilities. */
          esl_scorematrix_RelEntropy   (S, fi, fi, lambda,   &D);
          esl_scorematrix_ExpectedScore(S, fi, fi,           &E);

          printf("By solving for lambda from given background frequencies:\n");
          printf("Lambda           = %.4f\n",      lambda);
          printf("Relative entropy = %.4f bits\n", D);
          printf("Expected score   = %.4f bits\n", E * lambda * eslCONST_LOG2R);

          printf("p_ij's are:\n");   esl_dmatrix_Dump(stdout, P1, abc->sym, abc->sym);
          printf("fi's are:\n");     esl_vec_DDump(stdout, fi, S->K, abc->sym);
          printf("============================================================\n\n");

      }
      p7_bg_Destroy(bg);
  }
  else
  {
      vstatus = esl_scorematrix_Probify(S, &P1, &fi, &fj, &lambda);

      if (vstatus == eslOK)
      { /* Print some info, and the joint probabilities. */

          esl_scorematrix_RelEntropy   (S, fi, fj, lambda, &D);
          esl_scorematrix_ExpectedScore(S, fi, fj,         &E);

          printf("By Yu/Altschul (2003,2005) procedure:\n");
          printf("Lambda           = %.4f\n",      lambda);
          printf("Relative entropy = %.4f bits\n", D);
          printf("Expected score   = %.4f bits\n", E * lambda * eslCONST_LOG2R);

          printf("p_ij's are:\n");  esl_dmatrix_Dump(stdout, P1, abc->sym, abc->sym);
          printf("fi's are:\n");    esl_vec_DDump(stdout, fi, S->K, abc->sym);
          printf("fj's are:\n");    esl_vec_DDump(stdout, fj, S->K, abc->sym);
          printf("============================================================\n\n");
      }
      else
      {
          printf("Yu/Altschul procedure FAILS to find a valid implicit probability basis!\n");
          printf("Lambda  = %.4f\n",      lambda);
          printf("p_ij's are:\n");  esl_dmatrix_Dump(stdout, P1, abc->sym, abc->sym);
          printf("fi's are:\n");    esl_vec_DDump(stdout, fi, S->K, abc->sym);
          printf("fj's are:\n");    esl_vec_DDump(stdout, fj, S->K, abc->sym);
          printf("============================================================\n\n");

          esl_composition_BL62(fi); esl_composition_BL62(fj);
      }

      /* Now reverse engineer it again, this time using "known" background probs */
      esl_scorematrix_ProbifyGivenBG(S, fi, fj, &lambda, &P2);
      esl_scorematrix_RelEntropy   (S, fi, fj, lambda,   &D);
      esl_scorematrix_ExpectedScore(S, fi, fj,           &E);

      printf("By solving for lambda from given background frequencies:\n");
      printf("Lambda           = %.4f\n",      lambda);
      printf("Relative entropy = %.4f bits\n", D);
      printf("Expected score   = %.4f bits\n", E * lambda * eslCONST_LOG2R);

      printf("p_ij's are:\n");   esl_dmatrix_Dump(stdout, P2, abc->sym, abc->sym);
      printf("fi's are:\n");     esl_vec_DDump(stdout, fi, S->K, abc->sym);
      printf("fj's are:\n");     esl_vec_DDump(stdout, fj, S->K, abc->sym);
      printf("============================================================\n\n");

      /* Now recalculate a score matrix from the probabilistic basis */
      printf("Before:\n");
      esl_scorematrix_Write(stdout, S);
      printf("After:\n");
      esl_scorematrix_SetFromProbs(S, lambda, P2, fi, fj);
      esl_scorematrix_Write(stdout, S);

      free(fi);
      free(fj);
  }


  if (P1 != NULL) {esl_dmatrix_Destroy(P1);}
  if (P2 != NULL) {esl_dmatrix_Destroy(P2);}
  esl_scorematrix_Destroy(S);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);

  return 0;

ERROR:
  exit(1);
}
