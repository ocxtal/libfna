
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
#include "zf/zf.h"
#include "fna.h"
#include "kvec.h"
#include "sassert.h"
#include "log.h"

#define UNITTEST_UNIQUE_ID			33

#ifdef TEST
/* use auto-generated main function to run tests */
#define UNITTEST 					1
#define UNITTEST_ALIAS_MAIN			1
#endif

#include "unittest.h"
unittest_config(
	.name = "fna"
);

/* inline directive */
#define _force_inline				inline

/* roundup */
#define roundup(x, base)			( (((x) + (base) - 1) / (base)) * (base) )

/* type aliasing for returning values */
typedef kvec_t(uint8_t) kvec_uint8_t;

/**
 * @struct fna_context_s
 *
 * @brief a struct for fna context
 */
struct fna_context_s {
	char *path;
	uint8_t format;
	uint8_t encode;
	uint16_t options;
	int32_t status;
	zf_t *fp;				/** zf context (file pointer) */
	uint16_t head_margin;	/** margin at the head of fna_seq_t */
	uint16_t tail_margin;
	uint16_t seq_head_margin;	/** margin at the head of seq buffer */
	uint16_t seq_tail_margin;	/** margin at the tail of seq buffer */

	/* file format specific parser */
	struct fna_seq_intl_s *(*read)(struct fna_context_s *fna);

	/* output sequence format specific parser */
	int64_t (*read_seq)(struct fna_context_s *fna, kvec_uint8_t *v, char delim);
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
	char *name;
	int64_t name_len;
	uint8_t *seq;
	int64_t seq_len;
	uint8_t *qual;
	int64_t qual_len;
	uint8_t encode;			/** one of _fna_flag_encode */
	uint8_t reserved;
	uint16_t options;
	uint16_t head_margin;	/** margin at the head of fna_seq_t */
	uint16_t tail_margin;
	uint16_t seq_head_margin;	/** margin at the head of seq buffer */
	uint16_t seq_tail_margin;	/** margin at the tail of seq buffer */
};
_static_assert_offset(struct fna_seq_s, name, struct fna_seq_intl_s, name, 0);
_static_assert_offset(struct fna_seq_s, seq, struct fna_seq_intl_s, seq, 0);
_static_assert_offset(struct fna_seq_s, seq_len, struct fna_seq_intl_s, seq_len, 0);
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
static void fna_read_head_fasta(struct fna_context_s *fna);
static void fna_read_head_fastq(struct fna_context_s *fna);
static void fna_read_head_fast5(struct fna_context_s *fna);

static struct fna_seq_intl_s *fna_read_fasta(struct fna_context_s *fna);
static struct fna_seq_intl_s *fna_read_fastq(struct fna_context_s *fna);
static struct fna_seq_intl_s *fna_read_fast5(struct fna_context_s *fna);

static int64_t fna_read_seq_ascii(struct fna_context_s *fna, kvec_uint8_t *v, char delim);
static int64_t fna_read_seq_2bit(struct fna_context_s *fna, kvec_uint8_t *v, char delim);
static int64_t fna_read_seq_2bitpacked(struct fna_context_s *fna, kvec_uint8_t *v, char delim);
static int64_t fna_read_seq_4bit(struct fna_context_s *fna, kvec_uint8_t *v, char delim);
static int64_t fna_read_seq_4bitpacked(struct fna_context_s *fna, kvec_uint8_t *v, char delim);

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
	struct fna_context_s *fna = NULL;

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
	void (*read_head[])(
		struct fna_context_s *fna) = {
		[FNA_FASTA] = fna_read_head_fasta,
		[FNA_FASTQ] = fna_read_head_fastq,
		[FNA_FAST5] = fna_read_head_fast5
	};
	struct fna_seq_intl_s *(*read[])(
		struct fna_context_s *fna) = {
		[FNA_FASTA] = fna_read_fasta,
		[FNA_FASTQ] = fna_read_fastq,
		[FNA_FAST5] = fna_read_fast5
	};

	/* pack functions */
	int64_t (*read_seq[])(
		struct fna_context_s *fna,
		kvec_uint8_t *v,
		char delim) = {
		[FNA_ASCII] = fna_read_seq_ascii,
		[FNA_2BIT] = fna_read_seq_2bit,
		[FNA_2BITPACKED] = fna_read_seq_2bitpacked,
		[FNA_4BIT] = fna_read_seq_4bit,
		[FNA_4BITPACKED] = fna_read_seq_4bitpacked,
	};

	/* default params */
	struct fna_params_s default_params = {
		.seq_encode = FNA_ASCII,
		.file_format = FNA_UNKNOWN,
		.options = 0,
		.head_margin = 0,
		.tail_margin = 0,
		.seq_head_margin = 0,
		.seq_tail_margin = 0
	};

	if(path == NULL) { return NULL; }
	if(params == NULL) { params = &default_params; }

	if((fna = (struct fna_context_s *)malloc(sizeof(struct fna_context_s))) == NULL) {
		goto _fna_init_error_handler;
	}
	fna->path = NULL;
	fna->fp = NULL;

	/* copy params */
	fna->encode = params->seq_encode;	/** encode sequence to 2-bit if encode == FNA_2BITPACKED */
	fna->format = params->file_format;	/** format (see enum FNA_FORMAT) */
	fna->options = params->options;
	fna->head_margin = roundup(params->head_margin, 16);
	fna->tail_margin = roundup(params->tail_margin, 16);
	fna->seq_head_margin = roundup(params->seq_head_margin, 16);
	fna->seq_tail_margin = roundup(params->seq_tail_margin, 16);

	/* restore defaults */
	if(fna->encode == 0) { fna->encode = FNA_ASCII; }

	/* open file */
	fna->fp = zfopen(path, "r");
	if(fna->fp == NULL) { goto _fna_init_error_handler; }

	/**
	 * if fna->format is not specified...
	 * 1. determine file format from the path extension
	 */
	if(fna->format == 0) {
		char const *path_tail = fna->fp->path + strlen(fna->fp->path);
		for(ep = ext; ep->ext != NULL; ep++) {
			if(strncmp(path_tail - strlen(ep->ext), ep->ext, strlen(ep->ext)) == 0) {
				fna->format = ep->format; break;
			}
		}
	}

	/**
	 * 2. determine file format from the content of the file
	 */
	if(fna->format == 0) {
		/* peek the head of the file */
		char buf[32];
		uint64_t len = zfpeek(fna->fp, buf, 32);
		for(int64_t i = 0; i < len; i++) {
			switch(buf[i]) {
				case '>': fna->format = FNA_FASTA; break;
				case '@': fna->format = FNA_FASTQ; break;
			}
		}
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
	fna->read_seq = read_seq[fna->encode];
	fna->path = strdup(path);

	/* parse header */
	read_head[fna->format](fna);
	return((struct fna_s *)fna);

_fna_init_error_handler:
	if(fna != NULL) {
		if(fna->fp != NULL) {
			zfclose(fna->fp); fna->fp = NULL;
		}
		if(fna->path != NULL) {
			free(fna->path); fna->path = NULL;
		}
		return((struct fna_s *)fna);
	}
	return(NULL);
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
			zfclose(fna->fp); fna->fp = NULL;
		}
		if(fna->path != NULL) {
			free(fna->path); fna->path = NULL;
		}
		free(fna); fna = NULL;
	}
	return;
}

/**
 * miscellaneous tables and functions
 */

/**
 * @macro _fna_pack_seq
 */
#define _fna_pack_seq(_fna, _type, _name, _id, _seq, _len) ({ \
	void *ptr = malloc(sizeof(struct fna_seq_intl_s) \
		+ (_fna)->head_margin + (_fna)->tail_margin); \
	struct fna_seq_intl_s *seq = (struct fna_seq_intl_s *)(ptr + (_fna)->head_margin); \
	seq->head_margin = (_fna)->head_margin; \
	seq->tail_margin = (_fna)->tail_margin; \
	/* copy content */ \
	seq->type = (_type); \
	seq->name = (_name); \
	seq->id = (_id); \
	seq->seq = (_seq); \
	seq->seq_len = (_len); \
	seq->encode = (_fna)->encode; \
	seq; \
})

/**
 * @fn fna_seq_make_margin
 */
static _force_inline
void fna_seq_make_margin(
	kvec_uint8_t *v,
	int64_t len)
{
	for(int64_t i = 0; i < len; i++) {
		kv_push(*v, 0);
	}
	return;
}

/**
 * @fn fna_parse_version_string
 * @brief parse version string in ("%d.%d.%d", major, minor, patch) format,
 * return 0x10000 * major + 0x100 * minor + patch
 */
static _force_inline
uint64_t fna_parse_version_string(
	char const *str)
{
	char buf[256];
	uint64_t v[3] = { 0, 0, 0 };

	for(int64_t i = 0; i < 3; i++) {
		int64_t j = 0;
		while(*str != '\0') {
			if(*str == '.') { break; }
			buf[j++] = *str++;
		}
		buf[j] = '\0';
		v[i] = (uint64_t)strtoll(buf, NULL, 10);
		str++;
	}

	return(0x10000 * v[0] + 0x100 * v[1] + v[2]);
}

/**
 * @fn fna_read_ascii
 * @brief read ascii until delim
 */
static _force_inline
int64_t fna_read_ascii(
	struct fna_context_s *fna,
	kvec_uint8_t *v,
	char delim)
{
	int64_t len = 0;
	int c;

	/* strip spaces at the head */
	while((c = zfgetc(fna->fp)) != EOF && isspace(c)) {}
	if(c == EOF) {
		debug("reached eof");
		kv_push(*v, '\0');
		return(len);
	}

	/* read line until delim */
	kv_push(*v, c); len++;
	while((c = zfgetc(fna->fp)) != EOF && (char)c != delim) {
		debug("%c, %d", c, c);
		kv_push(*v, c); len++;
	}

	/* strip spaces at the tail */
	while(len-- > 0 && isspace(c = kv_pop(*v))) {
		debug("%c, %d", c, c);
	}
	kv_push(*v, c); len++;	/* push back last non-space char */

	kv_push(*v, '\0');		/* push null terminator */
	debug("finished, len(%lld)", len);
	return(len);
}

/**
 * @fn fna_read_skip
 */
static _force_inline
int64_t fna_read_skip(
	struct fna_context_s *fna,
	char delim)
{
	int c;
	while((c = zfgetc(fna->fp)) != EOF && (char)c != delim) {}
	return(0);
}


/**
 * @fn fna_type
 * @brief (internal) encode ascii base to 2-bit base.
 */
static _force_inline
uint8_t fna_type(int c)
{
	uint8_t const table[256] = {
		0,		/* EOF */
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
		0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0
	};
	return(table[(uint8_t)(c + 1)]);
}

/**
 * @fn fna_read_seq_ascii
 * @brief read seq until delim, with conv table
 */
static
int64_t fna_read_seq_ascii(
	struct fna_context_s *fna,
	kvec_uint8_t *v,
	char delim)
{
	int64_t len = 0;
	debug("%c, %d", delim, delim);
	while(1) {
		int c = zfgetc(fna->fp);
		if((char)c == delim || c == EOF) { break; }
		if(fna_type(c) != 1) { continue; }
		debug("%c, %d", c, c);
		kv_push(*v, (uint8_t)c); len++;
	}
	kv_push(*v, '\0');
	debug("finished, len(%lld)", len);

	fna->status = zfeof(fna->fp) ? FNA_EOF : FNA_SUCCESS;
	return(len);
}

/**
 * @fn fna_encode_2bit
 * @brief mapping IUPAC amb. to 4bit encoding
 */
static _force_inline
uint8_t fna_encode_2bit(
	int c)
{
	/* convert to upper case and subtract offset by 0x40 */
	#define _b(x)	( (x) & 0x1f )

	/* conversion tables */
	enum bases {
		A = 0x00, C = 0x01, G = 0x02, T = 0x03
	};
	uint8_t const table[] = {
		[_b('A')] = A,
		[_b('C')] = C,
		[_b('G')] = G,
		[_b('T')] = T,
		[_b('U')] = T,
		[_b('N')] = A,		/* treat 'N' as 'A' */
		[_b('_')] = 0		/* sentinel */
	};
	return(table[_b((uint8_t)c)]);

	#undef _b
}

/**
 * @fn fna_read_seq_2bit
 * @brief read seq until delim, with conv table
 */
static
int64_t fna_read_seq_2bit(
	struct fna_context_s *fna,
	kvec_uint8_t *v,
	char delim)
{
	int64_t len = 0;
	while(1) {
		int c = zfgetc(fna->fp);
		if((char)c == delim || c == EOF) { break; }
		if(fna_type(c) != 1) { continue; }
		kv_push(*v, fna_encode_2bit(c)); len++;
	}

	fna->status = zfeof(fna->fp) ? FNA_EOF : FNA_SUCCESS;
	return(len);
}

/**
 * @fn fna_read_seq_2bitpacked
 * @brief read seq until delim, with conv table
 */
static
int64_t fna_read_seq_2bitpacked(
	struct fna_context_s *fna,
	kvec_uint8_t *v,
	char delim)
{
	int64_t len = 0;
	uint64_t rem = 8;
	uint8_t arr = 0;
	/* 4x unrolled loop */
	while(1) {
		#define _fetch(_fna) ({ \
			char _c; \
			while(fna_type(_c = zfgetc(_fna->fp)) != 1) { \
				if((char)_c == delim || _c == EOF) { \
					goto _fna_read_seq_2bitpacked_finish; \
				} \
			} \
			_c; \
		}) \

		rem = 8;
		arr = (arr>>2) | (fna_encode_2bit(_fetch(fna))<<6); rem -= 2;
		arr = (arr>>2) | (fna_encode_2bit(_fetch(fna))<<6); rem -= 2;
		arr = (arr>>2) | (fna_encode_2bit(_fetch(fna))<<6); rem -= 2;
		arr = (arr>>2) | (fna_encode_2bit(_fetch(fna))<<6);
		kv_push(*v, arr); len += 4;

		#undef _fetch
	}
_fna_read_seq_2bitpacked_finish:;
	kv_push(*v, arr>>rem); len += (8 - rem) / 2;

	fna->status = zfeof(fna->fp) ? FNA_EOF : FNA_SUCCESS;
	return(len);
}

/**
 * @fn fna_encode_4bit
 * @brief mapping IUPAC amb. to 4bit encoding
 */
static _force_inline
uint8_t fna_encode_4bit(
	int c)
{
	/* convert to upper case and subtract offset by 0x40 */
	#define _b(x)	( (x) & 0x1f )

	/* conversion tables */
	enum bases {
		A = 0x01, C = 0x02, G = 0x04, T = 0x08
	};
	uint8_t const table[] = {
		[_b('A')] = A,
		[_b('C')] = C,
		[_b('G')] = G,
		[_b('T')] = T,
		[_b('U')] = T,
		[_b('R')] = A | G,
		[_b('Y')] = C | T,
		[_b('S')] = G | C,
		[_b('W')] = A | T,
		[_b('K')] = G | T,
		[_b('M')] = A | C,
		[_b('B')] = C | G | T,
		[_b('D')] = A | G | T,
		[_b('H')] = A | C | T,
		[_b('V')] = A | C | G,
		[_b('N')] = 0,		/* treat 'N' as a gap */
		[_b('_')] = 0		/* sentinel */
	};
	return(table[_b((uint8_t)c)]);

	#undef _b
}

/**
 * @fn fna_read_seq_4bit
 * @brief read seq until delim, with conv table
 */
static
int64_t fna_read_seq_4bit(
	struct fna_context_s *fna,
	kvec_uint8_t *v,
	char delim)
{
	int64_t len = 0;
	while(1) {
		int c = zfgetc(fna->fp);
		if((char)c == delim || c == EOF) { break; }
		if(fna_type(c) != 1) { continue; }
		kv_push(*v, fna_encode_4bit(c)); len++;
	}

	fna->status = zfeof(fna->fp) ? FNA_EOF : FNA_SUCCESS;
	return(len);
}

/**
 * @fn fna_read_seq_4bitpacked
 * @brief read seq until delim, with conv table
 */
static
int64_t fna_read_seq_4bitpacked(
	struct fna_context_s *fna,
	kvec_uint8_t *v,
	char delim)
{
	int64_t len = 0;
	uint64_t rem = 8;
	uint8_t arr = 0;
	/* 2x unrolled loop */
	while(1) {
		#define _fetch(_fna) ({ \
			char _c; \
			while(fna_type(_c = zfgetc(_fna->fp)) != 1) { \
				if((char)_c == delim || _c == EOF) { \
					goto _fna_read_seq_4bitpacked_finish; \
				} \
			} \
			_c; \
		}) \

		rem = 8;
		arr = (arr>>2) | (fna_encode_4bit(_fetch(fna))<<4); rem -= 4;
		arr = (arr>>2) | (fna_encode_4bit(_fetch(fna))<<4);
		kv_push(*v, arr); len += 2;

		#undef _fetch
	}
_fna_read_seq_4bitpacked_finish:;
	kv_push(*v, arr>>rem); len += (8 - rem) / 4;

	fna->status = zfeof(fna->fp) ? FNA_EOF : FNA_SUCCESS;
	return(len);
}

/**
 * @fn fna_read_head_fasta
 */
static
void fna_read_head_fasta(
	struct fna_context_s *fna)
{
	/* eat '>' at the head */
	int c;
	while((c = zfgetc(fna->fp)) != EOF && c != '>') {}
	return;
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
	kvec_uint8_t v;
	kv_init(v);

	/* make margin at the head of seq */
	fna_seq_make_margin(&v, fna->head_margin);

	kv_pusha(struct fna_seq_intl_s, v, ((struct fna_seq_intl_s){
		.head_margin = fna->head_margin,
		.tail_margin = fna->tail_margin,
		.seq_head_margin = fna->seq_head_margin,
		.seq_tail_margin = fna->seq_tail_margin
	}));
	dump(kv_ptr(v), 64);

	/* parse name */
	int64_t name_len = fna_read_ascii(fna, &v, '\n');	/* fasta header line must ends with '\n' */
	dump(kv_ptr(v), 64);

	/* parse seq */
	fna_seq_make_margin(&v, fna->seq_head_margin);
	int64_t seq_len = fna->read_seq(fna, &v, '>');
	dump(kv_ptr(v), 64);

	/* check termination */
	if(name_len == 0 && seq_len == 0) {
		kv_destroy(v);
		return(NULL);
	}

	/* make margin at the tail */
	fna_seq_make_margin(&v, fna->seq_tail_margin);
	kv_push(v, '\0');
	fna_seq_make_margin(&v, fna->tail_margin);

	/* finished, build links */
	struct fna_seq_intl_s *r = (struct fna_seq_intl_s *)(
		kv_ptr(v) + fna->head_margin);
	r->name = (char *)(r + 1);
	r->name_len = name_len;
	r->seq = (uint8_t *)(r->name + (name_len + 1) + r->seq_head_margin);
	r->seq_len = seq_len;
	r->qual = (uint8_t *)(r->seq + seq_len + 1);
	r->qual_len = 0;
	debug("%s, %s, %s", r->name, r->seq, r->qual);
	dump(kv_ptr(v), 64);

	return(r);
}

/**
 * @fn fna_read_head_fastq
 */
static
void fna_read_head_fastq(
	struct fna_context_s *fna)
{
	/* eat '@' at the head */
	char c;
	while((c = zfgetc(fna->fp)) != '@') {}
	return;
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
	kvec_uint8_t v;
	kv_init(v);

	/* make margin at the head of seq */
	fna_seq_make_margin(&v, fna->head_margin);

	kv_pusha(struct fna_seq_intl_s, v, ((struct fna_seq_intl_s){
		.head_margin = fna->head_margin,
		.tail_margin = fna->tail_margin,
		.seq_head_margin = fna->seq_head_margin,
		.seq_tail_margin = fna->seq_tail_margin
	}));

	/* parse name */
	int64_t name_len = fna_read_ascii(fna, &v, '\n');	/* fasta header line must ends with '\n' */

	/* parse seq */
	fna_seq_make_margin(&v, fna->seq_head_margin);
	int64_t seq_len = fna->read_seq(fna, &v, '+');
	fna_seq_make_margin(&v, fna->seq_tail_margin);

	/* skip name */
	fna_read_skip(fna, '\n');

	/* parse qual */
	int64_t qual_len = ((fna->options & FNA_SKIP_QUAL) == 0)
		? fna->read_seq(fna, &v, '@')
		: fna_read_skip(fna, '@');

	/* check termination */
	if(name_len == 0 && seq_len == 0) {
		kv_destroy(v);
		return(NULL);
	}

	/* make margin at the tail */
	fna_seq_make_margin(&v, fna->tail_margin);


	/* finished, build links */
	struct fna_seq_intl_s *r = (struct fna_seq_intl_s *)(
		kv_ptr(v) + fna->head_margin);
	r->name = (char *)(r + 1);
	r->name_len = name_len;
	r->seq = (uint8_t *)(r->name + name_len + 1 + r->seq_head_margin);
	r->seq_len = seq_len;
	r->qual = (uint8_t *)(r->seq + seq_len + 1 + r->seq_tail_margin);
	r->qual_len = qual_len;

	return(r);
}

/**
 * @fn fna_read_head_fast5
 */
static
void fna_read_head_fast5(
	struct fna_context_s *fna)
{
#ifdef HAVE_HDF5
	/** not implemented yet */
#endif /* HAVE_HDF5 */
	return;
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
		/* free if external mem is used */
		char const *name_base = (char *)(s + 1);
		if(s->name != name_base) { free(s->name); }
		
		uint8_t const *seq_base = (uint8_t const *)(
			name_base + s->seq_head_margin + s->name_len + 1);
		if(s->seq != seq_base) { free(s->seq); }
		
		uint8_t const *qual_base = (uint8_t const *)(
			seq_base + s->seq_tail_margin + s->seq_len + 1);
		if(s->qual != qual_base) { free(s->qual); }

		/* free context */
		free((void *)((uint8_t *)s - s->head_margin));
	}
	return;
}

/**
 * @fn fna_base_comp
 * @brief (internal) make complemented base
 */
static _force_inline
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

#if 0
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
#endif

#ifdef TEST
/**
 * unittests
 */
#include <sys/stat.h>

/**
 * @fn fdump
 * @brief dump string to file, returns 1 if succeeded
 */
static _force_inline
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
static _force_inline
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

/* unittest for parse_version_string */
unittest()
{
	#define assert_parse_version_string(_str, _num) \
		assert(fna_parse_version_string(_str) == (_num), \
			"%u", \
			fna_parse_version_string(_str));

	assert_parse_version_string("0.0.0", 0x000000);
	assert_parse_version_string("0.0.1", 0x000001);
	assert_parse_version_string("0.1.0", 0x000100);
	assert_parse_version_string("1.2.3", 0x010203);
	assert_parse_version_string("100.200.50", 0x64c832);
	assert_parse_version_string("0.0.01", 0x000001);
	assert_parse_version_string("0.0.10", 0x00000a);
	assert_parse_version_string("0.0.15", 0x00000f);

	#undef assert_parse_version_string
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
	char const *fasta_filename = "test_fna.fa";
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
	assert(seq->seq_len == 4, "len(%lld)", seq->seq_len);
	fna_seq_free(seq);

	/* test1 */
	seq = fna_read(fna);
	assert(strcmp(seq->name, "test1") == 0, "name(%s)", seq->name);
	assert(strcmp((char const *)seq->seq, "ATATCGCG") == 0, "seq(%s)", (char const *)seq->seq);
	assert(seq->seq_len == 8, "len(%lld)", seq->seq_len);
	fna_seq_free(seq);

	/* test2 */
	seq = fna_read(fna);
	assert(strcmp(seq->name, "test2") == 0, "name(%s)", seq->name);
	assert(strcmp((char const *)seq->seq, "AAAA") == 0, "seq(%s)", (char const *)seq->seq);
	assert(seq->seq_len == 4, "len(%lld)", seq->seq_len);
	fna_seq_free(seq);

	/* test eof */
	seq = fna_read(fna);
	assert(seq == NULL, "seq(%p)", seq);
	assert(fna->status == FNA_EOF, "status(%d)", fna->status);
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
	char const *fastq_filename = "test_fna.fq";
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
	assert(seq->seq_len == 4, "len(%lld)", seq->seq_len);
	assert(strcmp((char const *)seq->qual, "NNNN") == 0, "seq(%s)", (char const *)seq->qual);
	assert(seq->qual_len == 4, "len(%lld)", seq->qual_len);
	fna_seq_free(seq);

	/* test1 */
	seq = fna_read(fna);
	assert(strcmp(seq->name, "test1") == 0, "name(%s)", seq->name);
	assert(strcmp((char const *)seq->seq, "ATATCGCG") == 0, "seq(%s)", (char const *)seq->seq);
	assert(seq->seq_len == 8, "len(%lld)", seq->seq_len);
	assert(strcmp((char const *)seq->qual, "NNNNNNNN") == 0, "seq(%s)", (char const *)seq->qual);
	assert(seq->qual_len == 8, "len(%lld)", seq->qual_len);
	fna_seq_free(seq);

	/* test2 */
	seq = fna_read(fna);
	assert(strcmp(seq->name, "test2") == 0, "name(%s)", seq->name);
	assert(strcmp((char const *)seq->seq, "AAAA") == 0, "seq(%s)", (char const *)seq->seq);
	assert(seq->seq_len == 4, "len(%lld)", seq->seq_len);
	assert(strcmp((char const *)seq->qual, "NNNN") == 0, "seq(%s)", (char const *)seq->qual);
	assert(seq->qual_len == 4, "len(%lld)", seq->qual_len);
	fna_seq_free(seq);

	/* test eof */
	seq = fna_read(fna);
	assert(seq == NULL, "seq(%p)", seq);
	assert(fna->status == FNA_EOF, "status(%d)", fna->status);
	fna_seq_free(seq);

	fna_close(fna);

	/** cleanup files */
	remove(fastq_filename);
	return;
}

/* fastq with qual skipping */
unittest()
{
	/**
	 * create file
	 * test0: valid FASTQ format
	 * test1: a space in header
	 * test2: two spaces in header and two \n's between header and content
	 */
	char const *fastq_filename = "test_fna.fq";
	char const *fastq_content =
		"@test0\nAAAA\n+test0\nNNNN\n"
		"@ test1\nATAT\nCGCG\n+ test1\nNNNN\nNNNN\n"
		"@  test2\n\nAAAA\n+  test2\n\nNNNN\n";
	assert(fdump(fastq_filename, fastq_content));
	assert(fcmp(fastq_filename, strlen(fastq_content), (uint8_t const *)fastq_content));

	fna_t *fna = fna_init(fastq_filename, FNA_PARAMS(.options = FNA_SKIP_QUAL));
	assert(fna != NULL, "fna(%p)", fna);

	/* format detection (internal) */
	assert(fna->file_format == FNA_FASTQ, "fna->file_format(%d)", fna->file_format);

	/* test0 */
	fna_seq_t *seq = fna_read(fna);
	assert(strcmp(seq->name, "test0") == 0, "name(%s)", seq->name);
	assert(strcmp((char const *)seq->seq, "AAAA") == 0, "seq(%s)", (char const *)seq->seq);
	assert(seq->seq_len == 4, "len(%lld)", seq->seq_len);
	assert(strcmp((char const *)seq->qual, "") == 0, "seq(%s)", (char const *)seq->qual);
	assert(seq->qual_len == 0, "len(%lld)", seq->qual_len);
	fna_seq_free(seq);

	/* test1 */
	seq = fna_read(fna);
	assert(strcmp(seq->name, "test1") == 0, "name(%s)", seq->name);
	assert(strcmp((char const *)seq->seq, "ATATCGCG") == 0, "seq(%s)", (char const *)seq->seq);
	assert(seq->seq_len == 8, "len(%lld)", seq->seq_len);
	assert(strcmp((char const *)seq->qual, "") == 0, "seq(%s)", (char const *)seq->qual);
	assert(seq->qual_len == 0, "len(%lld)", seq->qual_len);
	fna_seq_free(seq);

	/* test2 */
	seq = fna_read(fna);
	assert(strcmp(seq->name, "test2") == 0, "name(%s)", seq->name);
	assert(strcmp((char const *)seq->seq, "AAAA") == 0, "seq(%s)", (char const *)seq->seq);
	assert(seq->seq_len == 4, "len(%lld)", seq->seq_len);
	assert(strcmp((char const *)seq->qual, "") == 0, "seq(%s)", (char const *)seq->qual);
	assert(seq->qual_len == 0, "len(%lld)", seq->qual_len);
	fna_seq_free(seq);

	/* test eof */
	seq = fna_read(fna);
	assert(seq == NULL, "seq(%p)", seq);
	assert(fna->status == FNA_EOF, "status(%d)", fna->status);
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
	char const *fasta_filename = "test_fna.txt";
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
	char const *fastq_filename = "test_fna.txt";
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

#if 0
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
#endif
#endif /* #ifdef TEST */

/**
 * end of fna.c
 */
