
/**
 * @file fna.c
 *
 * @brief FASTA / FASTQ reader implementation.
 *
 * @author Hajime Suzuki
 * @date 2015/05/23
 * @license Apache v2.
 *
 * @detail
 * Supported formats:
 *   FASTA (raw and gzipped)
 *   FASTQ (raw and gzipped)
 *   FAST5 (unsupported for now!!!)
 *
 * List of APIs:
 *   Basic readers:
 *     fna_t *fna_init(char const *path, fna_params_t *params);
 *     fna_seq_t *fna_read(fna_t const *fna, fna_seq_t *seq);
 *     void fna_seq_free(fna_seq_t *seq);
 *     void fna_close(fna_t *fna);
 *
 *   Sequence duplicators:
 *     fna_seq_t *fna_duplicate(fna_seq_t const *seq);
 *     fna_seq_t *fna_revcomp(fna_seq_t const *seq);
 *
 *   Sequence modifiers:
 *     void fna_append(fna_seq_t *dst, fna_seq_t const *src);
 *     void fna_append_revcomp(fna_seq_t *dst, fna_seq_t const *src);
 *
 * Types and members:
 *   fna_t (alias to struct _fna): sequence reader instance container.
 *     path: path to the file.
 *   fna_seq_t (alias to struct _fna_seq): sequence container.
 *     name: sequence name container (kvec_t(char) instance)
 *       name.a: pointer to the sequence name (null-terminated ASCII).
 *     seq: sequence container (kvec_t(uint8_t) or kpvec_t(uint8_t) instance)
 *       seq.a: pointer to the sequence (null-terminated when fna->encode == FNA_RAW)
 */
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdint.h>
#include "fna.h"
#include "kopen.h"
#include "kvec.h"
#include "sassert.h"

#define UNITTEST_UNIQUE_ID			33

#ifdef TEST
/* use auto-generated main function to run tests */
#define UNITTEST 					1
#define UNITTEST_ALIAS_MAIN			1
#endif

#include "unittest.h"

/**
 * @struct fna_context_s
 *
 * @brief a struct for fna context
 */
struct fna_context_s {
	char *path;
	int32_t format;
	int32_t encode;
	int32_t status;
	int fd;					/** file descriptor */
	void *ko;				/** kopen file pointer */
	gzFile fp;				/** zlib file pointer */
	int32_t state;
	uint8_t *buf;
	uint64_t blen;
	uint16_t head_margin;	/** margin at the head of fna_seq_t */
	uint16_t tail_margin;
	struct fna_seq_intl_s *(*read)(struct fna_context_s *fna);
};
_static_assert_offset(struct fna_s, path, struct fna_context_s, path, 0);
_static_assert_offset(struct fna_s, file_format, struct fna_context_s, format, 0);
_static_assert_offset(struct fna_s, seq_encode, struct fna_context_s, encode, 0);
_static_assert_offset(struct fna_s, status, struct fna_context_s, status, 0);

/**
 * @struct fna_seq_intl_s
 *
 * @brief a struct which contains individual sequence.
 */
struct fna_seq_intl_s {
	kvec_t(char) name;
	kvec_t(uint8_t) seq;
	// char *name;
	// uint8_t *seq;
	int64_t len;
	int32_t encode;			/** one of _fna_flag_encode */
	uint16_t head_margin;	/** margin at the head of fna_seq_t */
	uint16_t tail_margin;
};
_static_assert_offset(struct fna_seq_s, name, struct fna_seq_intl_s, name.a, 0);
_static_assert_offset(struct fna_seq_s, seq, struct fna_seq_intl_s, seq.a, 0);
_static_assert_offset(struct fna_seq_s, len, struct fna_seq_intl_s, len, 0);
_static_assert_offset(struct fna_seq_s, seq_encode, struct fna_seq_intl_s, encode, 0);

/**
 * @enum _fna_state
 *
 * @brief state flag constant
 */
enum _fna_state {
	FNA_SEQ_NAME	= 1,
	FNA_SEQ 		= 2,
	FNA_QUAL_NAME	= 3,
	FNA_QUAL 		= 4
};

/* function delcarations */
static struct fna_seq_intl_s *fna_read_fasta(struct fna_context_s *fna);
static struct fna_seq_intl_s *fna_read_fastq(struct fna_context_s *fna);
static struct fna_seq_intl_s *fna_read_fast5(struct fna_context_s *fna);

/**
 * @fn fna_init
 *
 * @brief create a sequence reader context
 *
 * @param[in] path : a path to a file to open.
 *
 * @return a pointer to the context
 */
fna_t *fna_init(
	char const *path,
	fna_params_t const *params)
{
	int64_t i;
	struct fna_context_s *fna = NULL;
	size_t len;
	size_t const blen = 16 * 1024;

	/* extension determination */
	struct _ext {
		char const *ext;
		int32_t format;
	};
	struct _ext const ext[] = {
		{".fasta", FNA_FASTA},
		{".fas",   FNA_FASTA},
		{".seq",   FNA_FASTA},
		{".fna",   FNA_FASTA},
		{".ffn",   FNA_FASTA},
		{".fa",    FNA_FASTA},
		{".fastq", FNA_FASTQ},
		{".fq",    FNA_FASTQ},
		{".fast5", FNA_FAST5},
		{".f5",    FNA_FAST5},
		{NULL,     0}
	};
	struct _ext const *ep;

	/* read functions */
	struct fna_seq_intl_s *(*read[4])(
		struct fna_context_s *fna) = {
		[FNA_FASTA] = fna_read_fasta,
		[FNA_FASTQ] = fna_read_fastq,
		[FNA_FAST5] = fna_read_fast5
	};

	/* default params */
	struct fna_params_s default_params = {
		.seq_encode = FNA_ASCII,
		.file_format = FNA_UNKNOWN,
		.head_margin = 0
	};

	if(path == NULL) { return NULL; }
	if(params == NULL) { params = &default_params; }

	if((fna = (struct fna_context_s *)malloc(sizeof(struct fna_context_s))) == NULL) {
		goto _fna_init_error_handler;
	}
	fna->path = NULL;
	fna->encode = params->seq_encode;	/** encode sequence to 2-bit if encode == FNA_2BITPACKED */
	fna->fd = -1;
	fna->fp = NULL;
	fna->ko = NULL;
	fna->format = params->file_format;	/** format (see enum FNA_FORMAT) */
	fna->state = 0;						/** initial state (see enum FNA_STATE in fna.h) */
	fna->head_margin = params->head_margin;
	fna->tail_margin = params->tail_margin;
	if((fna->buf = (uint8_t *)malloc(sizeof(uint8_t) * blen)) == NULL) {
		fna->status = FNA_ERROR_OUT_OF_MEM;
		goto _fna_init_error_handler;
	}
	fna->blen = blen;

	/**
	 * restore defaults
	 */
	if(fna->encode == 0) {
		fna->encode = FNA_ASCII;
	}

	/**
	 * check if the file is gzipped
	 */
	len = strlen(path);
	if(strncmp(path + len - 3, ".gz", 3) == 0) { len -= 3; }
//	if(strncmp(path + len - 3, ".xz", 3) == 0) { len -= 3; }	/** currently unsupported */

	/**
	 * if fna->format is not specified...
	 * 1. determine file format from the path extension
	 */
	if(fna->format == 0) {
		for(ep = ext; ep->ext != NULL; ep++) {
			if(strncmp(path + len - strlen(ep->ext), ep->ext, strlen(ep->ext)) == 0) {
				fna->format = ep->format; break;
			}
		}
	}

	/**
	 * 2. determine file format from the content of the file
	 */
	if(fna->format == 0) {
		fna->ko = kopen(path, &(fna->fd));
		if((fna->fp = gzdopen(fna->fd, "rb")) == NULL) {
			fna->status = FNA_ERROR_FILE_OPEN;
			goto _fna_init_error_handler;
		}
		for(i = 0; i < 4; i++) {
			switch(gzgetc(fna->fp)) {
				case '>': fna->format = FNA_FASTA; break;
				case '@': fna->format = FNA_FASTQ; break;
			}
		}
		gzclose(fna->fp); fna->fp = NULL;
		kclose(fna->ko); fna->ko = NULL;
	}
	if(fna->format == 0) {
		// log_error("Couldn't determine file format `%s'.\n", path);
		fna->status = FNA_ERROR_UNKNOWN_FORMAT;
		goto _fna_init_error_handler;
	}
	#ifndef HAVE_HDF5
		if(fna->format == FNA_FAST5) {
			// log_error("Fast5 file format is not supported in this build.\n");
			fna->status = FNA_ERROR_UNKNOWN_FORMAT;
			goto _fna_init_error_handler;
		}
	#endif
	fna->read = read[fna->format];
//	printf("%d\n", fna->format);
	fna->path = strdup(path);

	/**
	 * open file and set buffer size
	 */
	fna->ko = kopen(path, &(fna->fd));
	if((fna->fp = gzdopen(fna->fd, "rb")) == NULL) {
		fna->status = FNA_ERROR_FILE_OPEN;
		goto _fna_init_error_handler;
	}
	#if ZLIB_VERNUM >= 0x1240
		/** resize read buffer if supported */
		gzbuffer(fna->fp, 512 * 1024);		/** 512 kilobytes */
	#endif

	return((struct fna_s *)fna);

_fna_init_error_handler:
	if(fna != NULL) {
		if(fna->fp != NULL) {
			gzclose(fna->fp); fna->fp = NULL;
		}
		if(fna->ko != NULL) {
			kclose(fna->ko); fna->ko = NULL;
		}
		if(fna->path != NULL) {
			free(fna->path); fna->path = NULL;
		}
		if(fna->buf != NULL) {
			free(fna->buf); fna->buf = NULL;
		}
		// free(fna); fna = NULL;
		return((struct fna_s *)fna);
	}
	return NULL;
}

/**
 * @fn fna_close
 *
 * @brief clean up sequence reader context
 */
void fna_close(fna_t *ctx)
{
	struct fna_context_s *fna = (struct fna_context_s *)ctx;

	if(fna != NULL) {
		if(fna->fp != NULL) {
			gzclose(fna->fp); fna->fp = NULL;
		}
		if(fna->ko != NULL) {
			kclose(fna->ko); fna->ko = NULL;
		}
		if(fna->path != NULL) {
			free(fna->path); fna->path = NULL;
		}
		if(fna->buf != NULL) {
			free(fna->buf); fna->buf = NULL;
		}
		free(fna); fna = NULL;
	}
	return;
}

/**
 * miscellaneous tables and functions
 */
/**
 * @fn fna_encode_base
 * @brief (internal) encode ascii base to 2-bit base.
 */
static inline
uint8_t fna_encode_base(char c)
{
	uint8_t const table[256] = {
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
	};
	return(table[(unsigned char)c]);
}

/**
 * @fn fna_incr
 * @brief (internal) encode ascii base to 2-bit base.
 */
static inline
uint8_t fna_type(char c)
{
	uint8_t const table[256] = {
		0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
		0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
		0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 2,  0, 0, 0, 0,
		0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 2, 0, 
		2, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, 
		1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 0,  0, 0, 0, 0, 
		0, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, 
		1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 0,  0, 0, 0, 0, 
		0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
		0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
		0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
		0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
		0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
		0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
		0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 
		0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0
	};
	return(table[(unsigned char)c]);
}

/**
 * @fn fna_read_fasta
 *
 * @brief (internal) fasta format parser
 */
static
struct fna_seq_intl_s *fna_read_fasta(
	struct fna_context_s *fna)
{
	char c;

	/** check state consistency */
	if(fna->state != 0) {
		// log_error("Inconsistent fasta file format.");
		fna->status = FNA_ERROR_BROKEN_FORMAT;
		return NULL;
	}

	/** malloc mem */
	void *ptr = malloc(sizeof(struct fna_seq_intl_s)
		+ fna->head_margin + fna->tail_margin);
	struct fna_seq_intl_s *seq = (struct fna_seq_intl_s *)(ptr + fna->head_margin);
	seq->head_margin = fna->head_margin;		/** set head margin size */
	seq->tail_margin = fna->tail_margin;

	/** init vector */
	kv_init(seq->name);
	kv_init(seq->seq);

	/**
	 * read sequence name
	 */
	while((c = gzgetc(fna->fp)) != '>') {		/** search seqname token */
		if(gzeof(fna->fp)) {
			fna->status = FNA_EOF;
			goto _fna_read_fasta_error_handler;	/** no sequences found */
		}
	}
	fna->state = FNA_SEQ_NAME;
	kv_clear(seq->name);
	int head_spaces = 1;
	while((c = gzgetc(fna->fp)) != '\n') {
		/** header without sequence */
		if(gzeof(fna->fp)) {
			fna->status = FNA_EOF;
			goto _fna_read_fasta_error_handler;
		}
		/** skip contiguous spaces at the head */
		if(head_spaces && c == ' ') { continue; } else { head_spaces = 0; }
		kv_push(seq->name, c);
	}
	kv_push(seq->name, '\0');					/** terminator */

	/**
	 * read sequence
	 */
	fna->state = FNA_SEQ;
	kv_clear(seq->seq);
	while((c = gzgetc(fna->fp)) != '>') {
		if(gzeof(fna->fp)) {
			fna->status = FNA_EOF;				/** the last sequence */
			break;
		}
		if(fna_type(c) == 1) {
			switch(fna->encode) {
				case FNA_RAW: kv_push(seq->seq, c); break;
				case FNA_2BIT: kv_push(seq->seq, fna_encode_base(c)); break;
				case FNA_2BITPACKED: kpv_push(seq->seq, fna_encode_base(c)); break;
			}
		}
	}
	if(fna->encode == FNA_2BITPACKED) {
		seq->len = kpv_size(seq->seq);
	} else {
		seq->len = kv_size(seq->seq);
		if(fna->encode == FNA_RAW) {
			kv_push(seq->seq, '\0'); 			/** terminator */
		}
	}
	seq->encode = fna->encode;					/** set packing format */

	gzungetc(c, fna->fp);						/** push back '>' */
	fna->state = 0;								/** initial state */

	return(seq);

_fna_read_fasta_error_handler:
	if(seq != NULL) {
		if(kv_ptr(seq->name) != NULL) {
			kv_destroy(seq->name);
		}
		if(kv_ptr(seq->seq) != NULL) {
			kv_destroy(seq->seq);
		}
		free(seq); seq = NULL;
	}
	return NULL;
}

/**
 * @fn fna_read_fastq
 *
 * @brief (internal) fastq format parser
 */
static
struct fna_seq_intl_s *fna_read_fastq(
	struct fna_context_s *fna)
{
	char c;

	/** check state consistency */
	if(fna->state != 0) {
		// log_error("Inconsistent fastq file format.");
		fna->status = FNA_ERROR_BROKEN_FORMAT;
		return NULL;
	}

	void *ptr = malloc(sizeof(struct fna_seq_intl_s)
		+ fna->head_margin + fna->tail_margin);
	struct fna_seq_intl_s *seq = (struct fna_seq_intl_s *)(ptr + fna->head_margin);
	seq->head_margin = fna->head_margin;		/** set head margin size */
	seq->tail_margin = fna->tail_margin;

	/** init vector */
	kv_init(seq->name);
	kv_init(seq->seq);

	/**
	 * read sequence name
	 */
	while((c = gzgetc(fna->fp)) != '@') {		/** search seqname token */
		if(gzeof(fna->fp)) {
			fna->status = FNA_EOF;
			goto _fna_read_fastq_error_handler;
		}
	}
	fna->state = FNA_SEQ_NAME;
	kv_clear(seq->name);
	int head_spaces = 1;
	while((c = gzgetc(fna->fp)) != '\n') {
		if(gzeof(fna->fp)) {
			fna->status = FNA_EOF;
			goto _fna_read_fastq_error_handler;
		}
		if(head_spaces && c == ' ') { continue; } else { head_spaces = 0; }
		kv_push(seq->name, c);
	}
	kv_push(seq->name, '\0');						/** terminator */

	/**
	 * read sequence
	 */
	fna->state = FNA_SEQ;
	kv_clear(seq->seq);
	while((c = gzgetc(fna->fp)) != '+') {
		if(gzeof(fna->fp)) {
			fna->status = FNA_EOF;
			break;
		}
		if(fna_type(c) == 1) {
			switch(fna->encode) {
				case FNA_RAW: kv_push(seq->seq, c); break;
				case FNA_2BIT: kv_push(seq->seq, fna_encode_base(c)); break;
				case FNA_2BITPACKED: kpv_push(seq->seq, fna_encode_base(c)); break;
			}
		}
	}
	if(fna->encode == FNA_2BITPACKED) {
		seq->len = kpv_size(seq->seq);
	} else {
		seq->len = kv_size(seq->seq);
		if(fna->encode == FNA_RAW) {
			kv_push(seq->seq, '\0'); 				/** terminator */
		}
	}

	seq->encode = fna->encode;					/** set packing format */
	while((c = gzgetc(fna->fp)) != '@') {		/** skip quality info */
		if(gzeof(fna->fp)) {
			fna->status = FNA_EOF;
			break;
		}
	}
	gzungetc(c, fna->fp);						/** push back '@' */
	fna->state = 0;								/** initial state */
	return(seq);

_fna_read_fastq_error_handler:
	if(seq != NULL) {
		if(kv_ptr(seq->name) != NULL) {
			kv_destroy(seq->name);
		}
		if(kv_ptr(seq->seq) != NULL) {
			kv_destroy(seq->seq);
		}
		free(seq); seq = NULL;
	}
	return NULL;
}

/**
 * @fn fna_read_fast5
 *
 * @brief (internal) fast5 format parser
 */
static
struct fna_seq_intl_s *fna_read_fast5(
	struct fna_context_s *fna)
{
#ifdef HAVE_HDF5
	/** not implemented yet */
#endif /* HAVE_HDF5 */
	return(NULL);
}

/**
 * @fn fna_read
 *
 * @brief read a sequence
 *
 * @param[in] fna : a pointer to the context
 *
 * @return a pointer to a sequence object
 */
fna_seq_t *fna_read(fna_t *ctx)
{
	struct fna_context_s *fna = (struct fna_context_s *)ctx;
	if(fna == NULL) { return NULL; }
	return((fna_seq_t *)fna->read(fna));
}

/**
 * @fn fna_seq_free
 *
 * @brief clean up sequence object
 */
void fna_seq_free(fna_seq_t *seq)
{
	struct fna_seq_intl_s *s = (struct fna_seq_intl_s *)seq;
	if(s != NULL) {
		if(kv_ptr(s->name) != NULL) { kv_destroy(s->name); }
		if(kv_ptr(s->seq) != NULL) { kv_destroy(s->seq); }
		free((void *)s - s->head_margin);
	}
	return;
}

/**
 * @fn fna_base_comp
 * @brief (internal) make complemented base
 */
static inline
uint8_t fna_base_comp(uint8_t c)
{
	switch(c) {
		case 'a': return 't';
		case 'c': return 'g';
		case 'g': return 'c';
		case 't': return 'a';
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
		case 'T': return 'A';
		default: return 'N';
	}
}

/**
 * @fn fna_append
 *
 * @brief concatenate src sesquence after dst sequence
 */
void fna_append(
	fna_seq_t *_dst,
	fna_seq_t const *_src)
{
	struct fna_seq_intl_s *dst = (struct fna_seq_intl_s *)_dst;
	struct fna_seq_intl_s *src = (struct fna_seq_intl_s *)_src;

	int64_t i;
	int64_t size = src->len;

	if(src->encode != dst->encode || src->encode != src->encode) {
		return;
	}
	if(dst->encode == FNA_RAW) {
		(void)kv_pop(dst->seq);
	}
	dst->len += size;
	if(src->encode == FNA_RAW) {
		kv_reserve(dst->seq, dst->len + 1);
		for(i = 0; i < size; i++) {
			kv_push(dst->seq, kv_at(src->seq, i));
		}
		kv_push(dst->seq, '\0');
	} else if(src->encode == FNA_2BIT) {
		kv_reserve(dst->seq, dst->len);
		for(i = 0; i < size; i++) {
			kv_push(dst->seq, kv_at(src->seq, i));
		}
	} else if(src->encode == FNA_2BITPACKED) {
		kpv_reserve(dst->seq, dst->len);
		for(i = 0; i < size; i++) {
			kpv_push(dst->seq, kpv_at(src->seq, i));
		}
	}
	return;
}

/**
 * @fn fna_duplicate
 *
 * @brief duplicate sequence
 */
fna_seq_t *fna_duplicate(
	fna_seq_t const *_seq)
{
	struct fna_seq_intl_s *seq = (struct fna_seq_intl_s *)_seq;

	void *ptr = malloc(sizeof(struct fna_seq_intl_s)
		+ seq->head_margin + seq->tail_margin);
	struct fna_seq_intl_s *dup = (struct fna_seq_intl_s *)(ptr + seq->head_margin);
	dup->head_margin = seq->head_margin;		/** set head margin size */
	dup->tail_margin = seq->tail_margin;

	kv_init(dup->name);
	kv_init(dup->seq);
	dup->len = 0;
	if(seq->encode == FNA_RAW) {
		kv_push(dup->seq, '\0');
	}

	fna_append((fna_seq_t *)dup, (fna_seq_t *)seq);

	dup->name = seq->name;
	dup->name.a = strdup(seq->name.a);
	dup->encode = seq->encode;
	return((fna_seq_t *)dup);
}

/**
 * @fn fna_append_revcomp
 *
 * @brief append reverse complement of src after dst sequence
 */
void fna_append_revcomp(
	fna_seq_t *_dst,
	fna_seq_t const *_src)
{
	struct fna_seq_intl_s *dst = (struct fna_seq_intl_s *)_dst;
	struct fna_seq_intl_s *src = (struct fna_seq_intl_s *)_src;

	int64_t i;
	int64_t size = src->len;

	if(src->encode != dst->encode || src->encode != src->encode) {
		return;
	}

	if(dst->encode == FNA_RAW) {
		(void)kv_pop(dst->seq);
	}
	dst->len += size;
	if(src->encode == FNA_RAW) {
		kv_reserve(dst->seq, dst->len + 1);
		for(i = 0; i < size; i++) {
			kv_push(dst->seq, fna_base_comp(kv_at(src->seq, size - i - 1)));
		}
		kv_push(dst->seq, '\0');
	} else if(src->encode == FNA_2BIT) {
		kv_reserve(dst->seq, dst->len);
		for(i = 0; i < size; i++) {
			kv_push(dst->seq, 0x03 - kv_at(src->seq, size - i - 1));
		}
	} else if(src->encode == FNA_2BITPACKED) {
		kv_reserve(dst->seq, dst->len);
		for(i = 0; i < size; i++) {
			kpv_push(dst->seq, 0x03 - kpv_at(src->seq, size - i - 1));
		}
	}
	return;
}

/**
 * @fn fna_revcomp
 * @brief make reverse complemented sequence
 */
fna_seq_t *fna_revcomp(
	fna_seq_t const *_seq)
{
	struct fna_seq_intl_s *seq = (struct fna_seq_intl_s *)_seq;

	void *ptr = malloc(sizeof(struct fna_seq_intl_s)
		+ seq->head_margin + seq->tail_margin);
	struct fna_seq_intl_s *rev = (struct fna_seq_intl_s *)(ptr + seq->head_margin);
	rev->head_margin = seq->head_margin;		/** set head margin size */
	rev->tail_margin = seq->tail_margin;

	kv_init(rev->name);
	kv_init(rev->seq);
	rev->len = 0;
	if(seq->encode == FNA_RAW) {
		kv_push(rev->seq, '\0');
	}

	fna_append_revcomp((fna_seq_t *)rev, (fna_seq_t *)seq);

	rev->name = seq->name;
	rev->name.a = strdup(seq->name.a);
	rev->encode = seq->encode;
	return((fna_seq_t *)rev);
}

#ifdef TEST
/**
 * unittests
 */
#include <sys/stat.h>

unittest_config(
	.name = "fna"
);

/**
 * @fn fdump
 * @brief dump string to file, returns 1 if succeeded
 */
static inline
int fdump(
	char const *filename,
	char const *content)
{
	FILE *fp = fopen(filename, "w");
	int64_t l = fprintf(fp, "%s", content);
	fclose(fp);
	return(l == strlen(content));
}

/**
 * @fn fcmp
 * @brief compare file, returns zero if the content is equivalent to arr
 */
static inline
int fcmp(char const *filename, int64_t size, uint8_t const *arr)
{
	int64_t res;
	struct stat st;
	FILE *fp = NULL;
	uint8_t *buf = NULL;

	if((fp = fopen(filename, "rb")) == NULL) { return(1); }
	fstat(fileno(fp), &st);
	buf = malloc(sizeof(uint8_t) * st.st_size);
	fread(buf, st.st_size, sizeof(uint8_t), fp);
	res = memcmp(buf, arr, size);
	free(buf);
	return(res == 0);
}

/**
 * basic FASTA parsing
 */
unittest()
{
	/**
	 * create file
	 * test0: valid FASTA format
	 * test1: a space in header
	 * test2: two spaces in header and two \n's between header and content
	 */
	char const *fasta_filename = "test_fna_0.fa";
	char const *fasta_content =
		">test0\nAAAA\n"
		"> test1\nATAT\nCGCG\n"
		">  test2\n\nAAAA\n";
	assert(fdump(fasta_filename, fasta_content));
	assert(fcmp(fasta_filename, strlen(fasta_content), (uint8_t const *)fasta_content));

	fna_t *fna = fna_init(fasta_filename, NULL);
	assert(fna != NULL, "fna(%p)", fna);

	/* format detection (internal) */
	assert(fna->file_format == FNA_FASTA, "fna->file_format(%d)", fna->file_format);

	/* test0 */
	fna_seq_t *seq = fna_read(fna);
	assert(strcmp(seq->name, "test0") == 0, "name(%s)", seq->name);
	assert(strcmp((char const *)seq->seq, "AAAA") == 0, "seq(%s)", (char const *)seq->seq);
	assert(seq->len == 4, "len(%lld)", seq->len);
	fna_seq_free(seq);

	/* test1 */
	seq = fna_read(fna);
	assert(strcmp(seq->name, "test1") == 0, "name(%s)", seq->name);
	assert(strcmp((char const *)seq->seq, "ATATCGCG") == 0, "seq(%s)", (char const *)seq->seq);
	assert(seq->len == 8, "len(%lld)", seq->len);
	fna_seq_free(seq);

	/* test2 */
	seq = fna_read(fna);
	assert(strcmp(seq->name, "test2") == 0, "name(%s)", seq->name);
	assert(strcmp((char const *)seq->seq, "AAAA") == 0, "seq(%s)", (char const *)seq->seq);
	assert(seq->len == 4, "len(%lld)", seq->len);
	fna_seq_free(seq);

	fna_close(fna);

	/** cleanup files */
	remove(fasta_filename);
}

/**
 * basic FASTQ parsing
 */
unittest()
{
	/**
	 * create file
	 * test0: valid FASTQ format
	 * test1: a space in header
	 * test2: two spaces in header and two \n's between header and content
	 */
	char const *fastq_filename = "test_fna_10.fq";
	char const *fastq_content =
		"@test0\nAAAA\n+test0\nNNNN\n"
		"@ test1\nATAT\nCGCG\n+ test1\nNNNN\nNNNN\n"
		"@  test2\n\nAAAA\n+  test2\n\nNNNN\n";
	assert(fdump(fastq_filename, fastq_content));
	assert(fcmp(fastq_filename, strlen(fastq_content), (uint8_t const *)fastq_content));

	fna_t *fna = fna_init(fastq_filename, NULL);
	assert(fna != NULL, "fna(%p)", fna);

	/* format detection (internal) */
	assert(fna->file_format == FNA_FASTQ, "fna->file_format(%d)", fna->file_format);

	/* test0 */
	fna_seq_t *seq = fna_read(fna);
	assert(strcmp(seq->name, "test0") == 0, "name(%s)", seq->name);
	assert(strcmp((char const *)seq->seq, "AAAA") == 0, "seq(%s)", (char const *)seq->seq);
	assert(seq->len == 4, "len(%lld)", seq->len);
	fna_seq_free(seq);

	/* test1 */
	seq = fna_read(fna);
	assert(strcmp(seq->name, "test1") == 0, "name(%s)", seq->name);
	assert(strcmp((char const *)seq->seq, "ATATCGCG") == 0, "seq(%s)", (char const *)seq->seq);
	assert(seq->len == 8, "len(%lld)", seq->len);
	fna_seq_free(seq);

	/* test2 */
	seq = fna_read(fna);
	assert(strcmp(seq->name, "test2") == 0, "name(%s)", seq->name);
	assert(strcmp((char const *)seq->seq, "AAAA") == 0, "seq(%s)", (char const *)seq->seq);
	assert(seq->len == 4, "len(%lld)", seq->len);
	fna_seq_free(seq);

	fna_close(fna);

	/** cleanup files */
	remove(fastq_filename);
	return;
}

/**
 * format detection from the content of the file (FASTA)
 */
unittest()
{
	char const *fasta_filename = "test_fna_20.txt";
	char const *fasta_content =
		">test0\nAAAA\n"
		"> test1\nATAT\nCGCG\n"
		">  test2\n\nAAAA\n";
	assert(fdump(fasta_filename, fasta_content));
	assert(fcmp(fasta_filename, strlen(fasta_content), (uint8_t const *)fasta_content));
	fna_t *fna = fna_init(fasta_filename, NULL);
	assert(fna != NULL, "fna(%p)", fna);

	/* format detection (internal) */
	assert(fna->file_format == FNA_FASTA, "fna->file_format(%d)", fna->file_format);

	fna_close(fna);

	/** cleanup file */
	remove(fasta_filename);
	return;
}

/**
 * format detection (FASTQ)
 */
unittest()
{
	char const *fastq_filename = "test_fna_30.txt";
	char const *fastq_content =
		"@test0\nAAAA\n+test0\nNNNN\n"
		"@ test1\nATAT\nCGCG\n+ test1\nNNNN\nNNNN\n"
		"@  test2\n\nAAAA\n+  test2\n\nNNNN\n";
	assert(fdump(fastq_filename, fastq_content));
	assert(fcmp(fastq_filename, strlen(fastq_content), (uint8_t const *)fastq_content));

	fna_t *fna = fna_init(fastq_filename, NULL);
	assert(fna != NULL, "fna(%p)", fna);

	/* format detection (internal) */
	assert(fna->file_format == FNA_FASTQ, "fna->file_format(%d)", fna->file_format);

	fna_close(fna);

	/** cleanup file */
	remove(fastq_filename);
	return;
}

/**
 * sequence handling
 */
unittest()
{
	char const *fasta_filename = "test_fna_40.fa";
	char const *fasta_content = ">test0\nAACA\n";
	int32_t const margin = 32;
	char const *magic[3] = {
		"The quick brown fox jumps over the lazy dog.",
		"Lorem ipsum dolor sit amet, consectetur adipisicing elit,",
		"ETAOIN SHRDLU CMFWYP VBGKQJ XZ  "
	};
	assert(fdump(fasta_filename, fasta_content));
	assert(fcmp(fasta_filename, strlen(fasta_content), (uint8_t const *)fasta_content));

	fna_t *fna = fna_init(fasta_filename,
		FNA_PARAMS(
			.head_margin = margin,
			.tail_margin = margin
		));
	assert(fna != NULL, "fna(%p)", fna);

	fna_seq_t *seq = fna_read(fna);
	/* fill margin with magic */
	memcpy((void *)seq - margin, magic[0], margin);
	memcpy((void *)(seq + 1), magic[0], margin);

	/* duplicate */
	fna_seq_t *dup = fna_duplicate(seq);
	memcpy((void *)dup - margin, magic[1], margin);
	memcpy((void *)(dup + 1), magic[1], margin);
	assert(strcmp(dup->name, "test0") == 0, "name(%s)", dup->name);
	assert(strcmp((char const *)dup->seq, "AACA") == 0, "dup(%s)", (char const *)dup->seq);
	assert(dup->len == 4, "len(%lld)", dup->len);

	/* generate reverse complement */
	fna_seq_t *rev = fna_revcomp(dup);
	memcpy((void *)rev - margin, magic[2], margin);
	memcpy((void *)(rev + 1), magic[2], margin);
	assert(strcmp(rev->name, "test0") == 0, "name(%s)", rev->name);
	assert(strcmp((char const *)rev->seq, "TGTT") == 0, "rev(%s)", (char const *)rev->seq);
	assert(rev->len == 4, "len(%lld)", rev->len);

	/* append */
	fna_append(dup, dup);
	assert(strcmp(dup->name, "test0") == 0, "name(%s)", dup->name);
	assert(strcmp((char const *)dup->seq, "AACAAACA") == 0, "dup(%s)", (char const *)dup->seq);
	assert(dup->len == 8, "len(%lld)", dup->len);

	/* append reverse complement */
	fna_append_revcomp(rev, rev);
	assert(strcmp(rev->name, "test0") == 0, "name(%s)", rev->name);
	assert(strcmp((char const *)rev->seq, "TGTTAACA") == 0, "rev(%s)", (char const *)rev->seq);
	assert(rev->len == 8, "len(%lld)", rev->len);

	/* check margin */
	assert(((struct fna_seq_intl_s *)seq)->head_margin == margin,
		"margin(%d)", ((struct fna_seq_intl_s *)seq)->head_margin);
	assert(((struct fna_seq_intl_s *)dup)->head_margin == margin,
		"margin(%d)", ((struct fna_seq_intl_s *)seq)->head_margin);
	assert(((struct fna_seq_intl_s *)rev)->head_margin == margin,
		"margin(%d)", ((struct fna_seq_intl_s *)seq)->head_margin);

	assert(strncmp((char const *)seq - margin, magic[0], margin) == 0, "");
	assert(strncmp((char const *)dup - margin, magic[1], margin) == 0, "");
	assert(strncmp((char const *)rev - margin, magic[2], margin) == 0, "");

	assert(((struct fna_seq_intl_s *)seq)->tail_margin == margin,
		"margin(%d)", ((struct fna_seq_intl_s *)seq)->tail_margin);
	assert(((struct fna_seq_intl_s *)dup)->tail_margin == margin,
		"margin(%d)", ((struct fna_seq_intl_s *)seq)->tail_margin);
	assert(((struct fna_seq_intl_s *)rev)->tail_margin == margin,
		"margin(%d)", ((struct fna_seq_intl_s *)seq)->tail_margin);

	assert(strncmp((char const *)(seq + 1), magic[0], margin) == 0, "");
	assert(strncmp((char const *)(dup + 1), magic[1], margin) == 0, "");
	assert(strncmp((char const *)(rev + 1), magic[2], margin) == 0, "");

	/* cleanup */
	fna_seq_free(seq);
	fna_seq_free(dup);
	fna_seq_free(rev);

	fna_close(fna);
	remove(fasta_filename);
}
#endif /* #ifdef TEST */

/**
 * end of fna.c
 */
