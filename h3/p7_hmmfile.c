/* Input/output of HMMs.
 * 
 * Contents:
 *     1. The P7_HMMFILE object.
 *     2. Writing HMMs to save files.
 *     3. Reading HMMs from save files.
 *     4. Some private functions.
 *     5. Unit tests.
 *     6. Test driver.
 *     7. Copyright and license.
 * 
 * SRE, Wed Jan  3 13:48:12 2007 [Janelia]
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdio.h>
#include "easel.h"
#include "esl_ssi.h" 		/* this gives us esl_byteswap */
#include "p7_hmmfile.h"

/* Magic numbers identifying binary formats.
 * Do not change the old magics! Necessary for backwards compatibility.
 */
static uint32_t  v10magic = 0xe8ededb1; /* v1.0 binary: "hmm1" + 0x80808080 */
static uint32_t  v10swap  = 0xb1edede8; /* byteswapped v1.0                 */
static uint32_t  v11magic = 0xe8ededb2; /* v1.1 binary: "hmm2" + 0x80808080 */
static uint32_t  v11swap  = 0xb2edede8; /* byteswapped v1.1                 */
static uint32_t  v17magic = 0xe8ededb3; /* v1.7 binary: "hmm3" + 0x80808080 */
static uint32_t  v17swap  = 0xb3edede8; /* byteswapped v1.7                 */
static uint32_t  v19magic = 0xe8ededb4; /* V1.9 binary: "hmm4" + 0x80808080 */
static uint32_t  v19swap  = 0xb4edede8; /* V1.9 binary, byteswapped         */ 
static uint32_t  v20magic = 0xe8ededb5; /* V2.0 binary: "hmm5" + 0x80808080 */
static uint32_t  v20swap  = 0xb5edede8; /* V2.0 binary, byteswapped         */
static uint32_t  v30magic = 0xe8ededb6; /* V3.0 binary: "hmm6" + 0x80808080 */
static uint32_t  v30swap  = 0xb6edede8; /* V3.0 binary, byteswapped         */

static int write_bin_string(FILE *fp, char *s);
static int read_bin_string (FILE *fp, int doswap, char **ret_s);


/*****************************************************************
 * 1. The P7_HMMFILE object.
 *****************************************************************/

/* Function:  esl_hmmfile_Open()
 * Incept:    SRE, Wed Jan  3 18:38:10 2007 [Casa de Gatos]
 *
 * Purpose:   Open an HMM file or HMM database in the file <filename>
 *            and prepare to read the first HMM.
 *            
 *            We look for <filename> relative to the current working
 *            directory. Additionally, if we don't find it in the cwd
 *            and <env> is non-NULL, we will look for <filename>
 *            relative to one or more directories in a colon-delimited
 *            list obtained from the environment variable <env>. For
 *            example, if we had <setenv HMMERDB
 *            /misc/db/Pfam:/misc/db/Rfam> in the environment, a
 *            profile HMM application might pass "HMMERDB" as <env>.
 *            
 *            As a special case, if <filename> is "-", then HMMs will
 *            be read from <stdin>. In this case, <env> has no effect.
 *            [NOT YET IMPLEMENTED]
 *            
 *            As another special case, if <filename> ends in a <.gz>
 *            suffix, the file is assumed to be compressed by GNU
 *            <gzip>, and it is opened for reading from a pipe with
 *            <gunzip -dc>. This feature is only available on
 *            POSIX-compliant systems that have a <popen()> call, and
 *            <HAVE_POPEN> is defined by the configure script at
 *            compile time. 
 *            [NOT YET IMPLEMENTED]
 *            
 * Args:      
 *
 * Returns:   <eslOK> on success, and the open <ESL_HMMFILE> is returned
 *            in <*ret_hfp>.
 *            
 *            <eslENOTFOUND> if <filename> can't be opened for
 *            reading, even after the list of directories in <env> (if
 *            any) is checked.

 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_hmmfile_Open(char *filename, char *env, P7_HMMFILE **ret_hfp)
{
  P7_HMMFILE *hfp = NULL;
  int         status;

  ESL_ALLOC(hfp, sizeof(P7_HMMFILE));
  hfp->f      = NULL;
  hfp->parser = NULL;

  if ((hfp->f = fopen(filename, "r")) != NULL) 
    ;
  else if (esl_FileEnvOpen(filename, env, &(hfp->f), &envfile) == eslOK)
    ;
  else
    { status = eslENOTFOUND; goto ERROR; }
  
  hfp->parser = read_bin30hmm;

  *ret_hfp = hfp;
  return eslOK;

 ERROR:
  if (hfp != NULL) p7_hmmfile_Close(hfp);
  *ret_hfp = NULL;
  return status;
}

/* Function:  p7_hmmfile_Close()
 * Incept:    SRE, Wed Jan  3 18:48:44 2007 [Casa de Gatos]
 *
 * Purpose:   Closes an open HMM file <hfp>.
 *
 * Returns:   (void)
 */
void
p7_hmmfile_Close(P7_HMMFILE *hfp)
{
  if (hfp != NULL)
    {
      if (hfp->f != NULL) fclose(afp->f);
    }
  return;
}



/* Function:  p7_hmmfile_Write()
 * Incept:    SRE, Wed Jan  3 13:50:26 2007 [Janelia]			
 * 
 * Purpose:   Writes an HMM to a file in binary format.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_hmmfile_Write(FILE *fp, P7_HMM *hmm)
{
  int k;

  /* ye olde magic number */
  fwrite((char *) &(v30magic), sizeof(uint32_t), 1, fp);

  /* info necessary for sizes of things
   */
  fwrite((char *) &(hmm->flags),      sizeof(int),  1,   fp);
  fwrite((char *) &(hmm->M),          sizeof(int),  1,   fp);
  fwrite((char *) &(hmm->abc->type),  sizeof(int),  1,   fp);

  /* The core model probabilities
   */
  for (k = 1; k <= hmm->M; k++)
    fwrite((char *) hmm->mat[k], sizeof(float), hmm->abc->K, fp);
  for (k = 1; k < hmm->M; k++)
    fwrite((char *) hmm->ins[k], sizeof(float), hmm->abc->K, fp);
  for (k = 1; k < hmm->M; k++)
    fwrite((char *) hmm->t[k], sizeof(float), 7, fp);

  /* annotation section
   */
  write_bin_string(fp, hmm->name);
  if (hmm->flags & p7_ACC)  write_bin_string(fp, hmm->acc);
  if (hmm->flags & p7_DESC) write_bin_string(fp, hmm->desc);
  if (hmm->flags & p7_RF)   fwrite((char *) hmm->rf,  sizeof(char), hmm->M+1, fp);
  if (hmm->flags & p7_CS)   fwrite((char *) hmm->cs,  sizeof(char), hmm->M+1, fp);
  if (hmm->flags & p7_CA)   fwrite((char *) hmm->ca,  sizeof(char), hmm->M+1, fp);
  write_bin_string(fp, hmm->comlog);
  fwrite((char *) &(hmm->nseq), sizeof(int),  1,   fp);
  write_bin_string(fp, hmm->ctime);
  if (hmm->flags & p7_MAP)  fwrite((char *) hmm->map, sizeof(int), hmm->M+1, fp);
  fwrite((char *) &(hmm->checksum), sizeof(int),  1,   fp);

  /* Pfam cutoffs section
   */
  if (hmm->flags & p7_GA) {
    fwrite((char *) &(hmm->ga1), sizeof(float), 1, fp);
    fwrite((char *) &(hmm->ga2), sizeof(float), 1, fp);
  }
  if (hmm->flags & p7_TC) {
    fwrite((char *) &(hmm->tc1), sizeof(float), 1, fp);
    fwrite((char *) &(hmm->tc2), sizeof(float), 1, fp);
  }
  if (hmm->flags & p7_NC) {
    fwrite((char *) &(hmm->nc1), sizeof(float), 1, fp);
    fwrite((char *) &(hmm->nc2), sizeof(float), 1, fp);
  }
  return eslOK;
}


int
p7_hmmfile_Read(P7_HMMFILE *hfp, ESL_ALPHABET **ret_abc,  P7_HMM **ret_hmm)
{
  /* A call to SSI to remember file position will go here.
   */
  return (*hfp->parser)(hfp, ret_abc, ret_hmm);
}



/*****************************************************************
 * 2. Private functions for reading/parsing various save file formats
 *****************************************************************/

/* Binary save files from HMMER 3.x
 */
static int
read_bin30hmm(P7_HMMFILE *hmmfp, ESL_ALPHABET **ret_abc, P7_HMM **ret_hmm)
{
  int     status;
  P7_HMM *hmm = NULL;
  uint32_t magic;
  int     flags;
  int     M;
  int     alphabet_type;
  int     k,x;

  /* Check magic. */
  if (feof(hmmfp->f)) { status = eslEOD; goto ERROR; }
  if (! fread((char *) &magic, sizeof(uint32_t), 1, hmmfp->f))    { status = eslEOD;       goto ERROR; }
  if (magic != v30magic)                                          { status = eslEFORMAT;   goto ERROR; }

  /* Get sizes of things */
  if (! fread((char *) &flags,         sizeof(int), 1, hmmfp->f)) { status = eslEOD;       goto ERROR; }
  if (! fread((char *) &M,             sizeof(int), 1, hmmfp->f)) { status = eslEOD;       goto ERROR; }
  if (! fread((char *) &alphabet_type, sizeof(int), 1, hmmfp->f)) { status = eslEOD;       goto ERROR; }
  
  /* Set or verify alphabet. */
  if (*ret_abc == NULL)	{	/* still unknown: set it, pass control of it back to caller */
    if ((*ret_abc = esl_alphabet_Create(alphabet_type)) == NULL)  { status = eslEMEM;      goto ERROR; }
  } else {			/* already known: check it */
    if ((*ret_abc)->type != alphabet_type)                        { status = eslEINCOMPAT; goto ERROR; }
  }

  /* Allocate the new HMM. */
  if ((hmm = p7_hmm_Create(M, *ret_abc)) == NULL)                 { status = eslEMEM;      goto ERROR; }  
  
  /* Core model probabilities. */
  for (k = 1; k <= hmm->M; k++)
    if (! fread((char *) hmm->mat[k], sizeof(float), hmm->abc->K, hfp->f)) { status = eslEOD; goto ERROR; }
  for (k = 1; k < hmm->M; k++)
    if (! fread((char *) hmm->ins[k], sizeof(float), hmm->abc->K, hfp->f)) { status = eslEOD; goto ERROR; }
  for (k = 1; k < hmm->M; k++)
    if (! fread((char *) hmm->t[k],   sizeof(float), 7,           hfp->f)) { status = eslEOD; goto ERROR; }
  
  /* Annotations. */
  if ((status = read_bin_string(hmmfp->f, hmmfp->byteswap, &(hmm->name))) != eslOK) goto ERROR;
      


}


/*****************************************************************
 * 3. Private functions involved in i/o
 *****************************************************************/

/* Function: write_bin_string()
 * Date:     SRE, Wed Oct 29 13:49:27 1997 [TWA 721 over Canada]
 * 
 * Purpose:  Write a string in binary save format: an integer
 *           for the string length (including \0), followed by
 *           the string.
 */
static int
write_bin_string(FILE *fp, char *s)
{
  int len;
  if (s != NULL) 
    {
      len = strlen(s) + 1;
      fwrite((char *) &len, sizeof(int),  1,   fp);
      fwrite((char *) s,    sizeof(char), len, fp);
    }
  else
    {
      len = 0;
      fwrite((char *) &len, sizeof(int), 1, fp);
    }
  return eslOK;
}

/* Function: read_bin_string()
 * Date:     SRE, Wed Oct 29 14:03:23 1997 [TWA 721]
 * 
 * Purpose:  Read in a string from a binary file, where
 *           the first integer is the length (including '\0').
 *           
 * Args:     fp       - FILE to read from
 *           ret_s    - string to read into
 *                             
 * Return:   <eslOK> on success. ret_s is malloc'ed here.
 */                            
static int
read_bin_string(FILE *fp, char **ret_s)
{
  int   status;
  char *s = NULL;
  int   len;

  if (! fread((char *) &len, sizeof(int), 1, fp)) { status = eslEOD; goto FAILURE; }
  ESL_ALLOC(s,  (sizeof(char) * len));
  if (! fread((char *) s, sizeof(char), len, fp)) { status = eslEOD; goto FAILURE; }
  *ret_s = s;
  return eslOK;

 ERROR:
  if (s != NULL) free(s);
  *ret_s = NULL;
  return status;
}


/*****************************************************************
 * @LICENSE@
 *****************************************************************/