
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
 *       seq.a: pointer to the sequence (null-terminated when fna->seq_encode == FNA_RAW)
 */
#include <stdint.h>
#include <string.h>
#include "zf/zf.h"
#include "fna.h"
#include "kvec.h"
#include "sassert.h"
#include "log.h"

#define UNITTEST_UNIQUE_ID			33
#define UNITTEST 					1

#include "unittest.h"


/* inline directive */
#define _force_inline				inline

/* roundup */
#define _roundup(x, base)			( (((x) + (base) - 1) / (base)) * (base) )

/* type aliasing for returning values */
typedef kvec_t(uint8_t) kvec_uint8_t;

/**
 * @struct fna_read_ret_s
 */
struct fna_read_ret_s {
	int64_t len;
	char c;
};

/**
 * @struct fna_context_s
 *
 * @brief a struct for fna context
 */
struct fna_context_s {
	char *path;
	uint8_t file_format;
	uint8_t seq_encode;
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
	struct fna_read_ret_s (*read_seq)(struct fna_context_s *fna, kvec_uint8_t *v, uint8_t const *delim_table);
};
_static_assert_offset(struct fna_s, path, struct fna_context_s, path, 0);
_static_assert_offset(struct fna_s, file_format, struct fna_context_s, file_format, 0);
_static_assert_offset(struct fna_s, seq_encode, struct fna_context_s, seq_encode, 0);
_static_assert_offset(struct fna_s, status, struct fna_context_s, status, 0);

/**
 * @struct fna_seq_intl_s
 *
 * @brief a struct which contains individual sequence.
 */
struct fna_seq_intl_s {
	uint8_t type;			/** type, seq or link */
	uint8_t seq_encode;		/** one of _fna_flag_encode */
	uint16_t options;
	union {
		struct {
			char *name;					/** sequence name */
			int32_t name_len;
			int32_t reserved;
			uint8_t *seq;				/** sequence */
			int64_t seq_len;			/** sequence length */
			uint8_t *qual;
			int64_t qual_len;
		} segment;
		struct {
			char *from;
			int32_t from_len;
			int32_t from_ori;
			char *to;
			int32_t to_len;
			int32_t to_ori;
			char *cigar;				/** contains link cigar string */
			int64_t cigar_len;			/** contains link cigar string length (== strlen(seq)) */
		} link;
	};
	uint16_t head_margin;	/** margin at the head of fna_seq_t */
	uint16_t tail_margin;
	uint16_t seq_head_margin;	/** margin at the head of seq buffer */
	uint16_t seq_tail_margin;	/** margin at the tail of seq buffer */
};
_static_assert_offset(struct fna_seq_s, type, struct fna_seq_intl_s, type, 0);
_static_assert_offset(struct fna_seq_s, seq_encode, struct fna_seq_intl_s, seq_encode, 0);
_static_assert_offset(struct fna_seq_s, options, struct fna_seq_intl_s, options, 0);

/* segment members */
_static_assert_offset(struct fna_seq_s, segment.name, struct fna_seq_intl_s, segment.name, 0);
_static_assert_offset(struct fna_seq_s, segment.name_len, struct fna_seq_intl_s, segment.name_len, 0);
_static_assert_offset(struct fna_seq_s, segment.seq, struct fna_seq_intl_s, segment.seq, 0);
_static_assert_offset(struct fna_seq_s, segment.seq_len, struct fna_seq_intl_s, segment.seq_len, 0);
_static_assert_offset(struct fna_seq_s, segment.qual, struct fna_seq_intl_s, segment.qual, 0);
_static_assert_offset(struct fna_seq_s, segment.qual_len, struct fna_seq_intl_s, segment.qual_len, 0);

/* link members */
_static_assert_offset(struct fna_seq_s, link.from, struct fna_seq_intl_s, link.from, 0);
_static_assert_offset(struct fna_seq_s, link.from_len, struct fna_seq_intl_s, link.from_len, 0);
_static_assert_offset(struct fna_seq_s, link.to, struct fna_seq_intl_s, link.to, 0);
_static_assert_offset(struct fna_seq_s, link.to_len, struct fna_seq_intl_s, link.to_len, 0);
_static_assert_offset(struct fna_seq_s, link.cigar, struct fna_seq_intl_s, link.cigar, 0);
_static_assert_offset(struct fna_seq_s, link.cigar_len, struct fna_seq_intl_s, link.cigar_len, 0);


/* function delcarations */
static int fna_read_head_fasta(struct fna_context_s *fna);
static int fna_read_head_fastq(struct fna_context_s *fna);
static int fna_read_head_fast5(struct fna_context_s *fna);
static int fna_read_head_gfa(struct fna_context_s *fna);

static struct fna_seq_intl_s *fna_read_fasta(struct fna_context_s *fna);
static struct fna_seq_intl_s *fna_read_fastq(struct fna_context_s *fna);
static struct fna_seq_intl_s *fna_read_fast5(struct fna_context_s *fna);
static struct fna_seq_intl_s *fna_read_gfa(struct fna_context_s *fna);

static struct fna_read_ret_s fna_read_seq_ascii(struct fna_context_s *fna, kvec_uint8_t *v, uint8_t const *delim_table);
static struct fna_read_ret_s fna_read_seq_2bit(struct fna_context_s *fna, kvec_uint8_t *v, uint8_t const *delim_table);
static struct fna_read_ret_s fna_read_seq_2bitpacked(struct fna_context_s *fna, kvec_uint8_t *v, uint8_t const *delim_table);
static struct fna_read_ret_s fna_read_seq_4bit(struct fna_context_s *fna, kvec_uint8_t *v, uint8_t const *delim_table);
static struct fna_read_ret_s fna_read_seq_4bitpacked(struct fna_context_s *fna, kvec_uint8_t *v, uint8_t const *delim_table);

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
		int32_t file_format;
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
		{".gfa",   FNA_GFA},
		{NULL,     0}
	};
	struct _ext const *ep;

	/* read functions */
	int (*read_head[])(
		struct fna_context_s *fna) = {
		[FNA_FASTA] = fna_read_head_fasta,
		[FNA_FASTQ] = fna_read_head_fastq,
		[FNA_FAST5] = fna_read_head_fast5,
		[FNA_GFA]	= fna_read_head_gfa
	};
	struct fna_seq_intl_s *(*read[])(
		struct fna_context_s *fna) = {
		[FNA_FASTA] = fna_read_fasta,
		[FNA_FASTQ] = fna_read_fastq,
		[FNA_FAST5] = fna_read_fast5,
		[FNA_GFA]	= fna_read_gfa
	};

	/* pack functions */
	struct fna_read_ret_s (*read_seq[])(
		struct fna_context_s *fna,
		kvec_uint8_t *v,
		uint8_t const *delim_table) = {
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
	fna->seq_encode = params->seq_encode;	/** encode sequence to 2-bit if encode == FNA_2BITPACKED */
	fna->file_format = params->file_format;	/** format (see enum FNA_FORMAT) */
	fna->options = params->options;
	fna->head_margin = _roundup(params->head_margin, 16);
	fna->tail_margin = _roundup(params->tail_margin, 16);
	fna->seq_head_margin = _roundup(params->seq_head_margin, 16);
	fna->seq_tail_margin = _roundup(params->seq_tail_margin, 16);

	/* restore defaults */
	if(fna->seq_encode == 0) { fna->seq_encode = FNA_ASCII; }

	/* open file */
	fna->fp = zfopen(path, "r");
	if(fna->fp == NULL) { goto _fna_init_error_handler; }

	/**
	 * if fna->file_format is not specified...
	 * 1. determine file format from the path extension
	 */
	if(fna->file_format == 0) {
		char const *path_tail = fna->fp->path + strlen(fna->fp->path);
		for(ep = ext; ep->ext != NULL; ep++) {
			if(strncmp(path_tail - strlen(ep->ext), ep->ext, strlen(ep->ext)) == 0) {
				fna->file_format = ep->file_format; break;
			}
		}
	}

	/**
	 * 2. determine file format from the content of the file
	 */
	if(fna->file_format == 0) {
		/* peek the head of the file */
		char buf[33] = { 0 };
		uint64_t len = zfpeek(fna->fp, buf, 32);
		for(int64_t i = 0; i < len; i++) {
			switch(buf[i]) {
				case '>': fna->file_format = FNA_FASTA; break;
				case '@': fna->file_format = FNA_FASTQ; break;
				case 'H':
				if(buf[i + 1] == '\t') {
					fna->file_format = FNA_GFA; break;
				}
			}
		}
	}
	if(fna->file_format == 0) {
		// log_error("Couldn't determine file format `%s'.\n", path);
		fna->status = FNA_ERROR_UNKNOWN_FORMAT;
		goto _fna_init_error_handler;
	}
	#ifndef HAVE_HDF5
		if(fna->file_format == FNA_FAST5) {
			// log_error("Fast5 file format is not supported in this build.\n");
			fna->status = FNA_ERROR_UNKNOWN_FORMAT;
			goto _fna_init_error_handler;
		}
	#endif
	fna->read = read[fna->file_format];
	fna->read_seq = read_seq[fna->seq_encode];
	fna->path = strdup(path);

	/* parse header */
	if(read_head[fna->file_format](fna) != FNA_SUCCESS) {
		/* something is wrong */
		goto _fna_init_error_handler;
	}
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
 * @macro _fna_pack_segment
 */
#define _fna_pack_segment(_seg, _name, _name_len, _seq, _seq_len, _qual, _qual_len) ({ \
	(_seg)->segment.name = (_name); \
	(_seg)->segment.name_len = (_name_len); \
	(_seg)->segment.seq = (_seq); \
	(_seg)->segment.seq_len = (_seq_len); \
	(_seg)->segment.qual = (_qual); \
	(_seg)->segment.qual_len = (_qual_len); \
	(_seg); \
})

/**
 * @macro _fna_pack_link
 */
#define _fna_pack_link(_link, _from, _from_len, _from_ori, _to, _to_len, _to_ori, _cigar, _cigar_len) ({ \
	(_link)->link.from = (_from); \
	(_link)->link.from_len = (_from_len); \
	(_link)->link.from_ori = (_from_ori); \
	(_link)->link.to = (_to); \
	(_link)->link.to_len = (_to_len); \
	(_link)->link.to_ori = (_to_ori); \
	(_link)->link.cigar = (_cigar); \
	(_link)->link.cigar_len = (_cigar_len); \
	(_link); \
})

/**
 * @val ascii_table, alpha_table, space_table
 */
/* non-printable, first 32 elements of the ASCII table */
#define _non_printable(n) \
	(n),(n),(n),(n), (n),(n),(n),(n), (n),(n),(n),(n), (n),(n),(n),(n), \
	(n),(n),(n),(n), (n),(n),(n),(n), (n),(n),(n),(n), (n),(n),(n),(n)

static
uint8_t const space_table[256] = {
	['\0'] = 1,
	[' '] = 1,
	['\t'] = 1,
	['\v'] = 1,
	[0xff] = 0xff
};
static
uint8_t const delim_line[256] = {
	['\r'] = 1,
	['\n'] = 1,
	[0xff] = 0xff
};
static
uint8_t const delim_fasta_seq[256] = {
	_non_printable(2),
	['>'] = 1,
	[0xff] = 0xff
};
static
uint8_t const delim_fastq_seq[256] = {
	_non_printable(2),
	['@'] = 1,
	[0xff] = 0xff
};
static
uint8_t const delim_fastq_qual[256] = {
	_non_printable(2),
	['+'] = 1,
	[0xff] = 0xff
};
static
uint8_t const delim_gfa_field[256] = {
	['\t'] = 1,
	['\r'] = 1,
	['\n'] = 1,
	[0xff] = 0xff
};

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
struct fna_read_ret_s fna_read_ascii(
	struct fna_context_s *fna,
	kvec_uint8_t *v,
	uint8_t const *delim_table)
{
	int64_t len = 0;
	int c;

	/* strip spaces at the head */
	while(space_table[(uint8_t)(c = zfgetc(fna->fp))] == 1) {
		debug("%c, %d", c, c);
	}
	if(c == EOF) {
		debug("reached eof");
		kv_push(*v, '\0');
		return((struct fna_read_ret_s){
			.len = len,
			.c = (char)EOF
		});
	}

	/* read line until delim */
	debug("%c, %d", c, c);
	kv_push(*v, c); len++;
	while(delim_table[(uint8_t)(c = zfgetc(fna->fp))] == 0) {
		debug("%c, %d", c, c);
		kv_push(*v, c); len++;
	}
	char cret = (char)c;

	/* strip spaces at the tail */
	while(len-- > 0 && space_table[(uint8_t)(c = kv_pop(*v))] == 1) {
		debug("%c, %d", c, c);
	}
	kv_push(*v, c); len++;	/* push back last non-space char */

	kv_push(*v, '\0');		/* push null terminator */
	debug("finished, len(%lld)", len);
	return((struct fna_read_ret_s){
		.len = len,
		.c = cret
	});
}

/**
 * @fn fna_read_skip
 */
static _force_inline
struct fna_read_ret_s fna_read_skip(
	struct fna_context_s *fna,
	uint8_t const *delim_table)
{
	int c;
	while((delim_table[(uint8_t)(c = zfgetc(fna->fp))] & 0x01) == 0) {
		debug("%c, %d, %u", c, c, delim_table[(uint8_t)c]);
	}
	return((struct fna_read_ret_s){
		.len = 0,
		.c = (char)c
	});
}

/**
 * @fn fna_read_seq_ascii
 * @brief read seq until delim, with conv table
 */
static
struct fna_read_ret_s fna_read_seq_ascii(
	struct fna_context_s *fna,
	kvec_uint8_t *v,
	uint8_t const *delim_table)
{
	int c;
	int64_t len = 0;
	while(1) {
		c = zfgetc(fna->fp);
		uint8_t type = delim_table[(uint8_t)c];
		if(type & 0x01) { break; }
		if(type != 0) { continue; }
		debug("%c, %d", c, c);
		kv_push(*v, (uint8_t)c); len++;
	}
	kv_push(*v, '\0');
	debug("finished, len(%lld)", len);

	fna->status = zfeof(fna->fp) ? FNA_EOF : FNA_SUCCESS;
	return((struct fna_read_ret_s){
		.len = len,
		.c = (char)c
	});
}

/**
 * @fn fna_encode_2bit
 * @brief mapping IUPAC amb. to 2bit encoding
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
struct fna_read_ret_s fna_read_seq_2bit(
	struct fna_context_s *fna,
	kvec_uint8_t *v,
	uint8_t const *delim_table)
{
	int c;
	int64_t len = 0;
	while(1) {
		c = zfgetc(fna->fp);
		uint8_t type = delim_table[(uint8_t)c];
		if(type & 0x01) { break; }
		if(type != 0) { continue; }
		kv_push(*v, fna_encode_2bit(c)); len++;
	}

	fna->status = zfeof(fna->fp) ? FNA_EOF : FNA_SUCCESS;
	return((struct fna_read_ret_s){
		.len = len,
		.c = (char)c
	});
}

/**
 * @fn fna_read_seq_2bitpacked
 * @brief read seq until delim, with conv table
 */
static
struct fna_read_ret_s fna_read_seq_2bitpacked(
	struct fna_context_s *fna,
	kvec_uint8_t *v,
	uint8_t const *delim_table)
{
	int c;
	int64_t len = 0;
	uint64_t rem = 8;
	uint8_t arr = 0;
	/* 4x unrolled loop */
	while(1) {
		#define _fetch(_fna) ({ \
			char _c; \
			uint8_t _type; \
			while((_type = delim_table[(uint8_t)(_c = zfgetc(_fna->fp))]) != 0) { \
				if(_type & 0x01) { goto _fna_read_seq_2bitpacked_finish; } \
			} \
			_c; \
		}) \

		rem = 8;
		arr = (arr>>2) | (fna_encode_2bit(c = _fetch(fna))<<6); rem -= 2;
		arr = (arr>>2) | (fna_encode_2bit(c = _fetch(fna))<<6); rem -= 2;
		arr = (arr>>2) | (fna_encode_2bit(c = _fetch(fna))<<6); rem -= 2;
		arr = (arr>>2) | (fna_encode_2bit(c = _fetch(fna))<<6);
		kv_push(*v, arr); len += 4;

		#undef _fetch
	}
_fna_read_seq_2bitpacked_finish:;
	kv_push(*v, arr>>rem); len += (8 - rem) / 2;

	fna->status = zfeof(fna->fp) ? FNA_EOF : FNA_SUCCESS;
	return((struct fna_read_ret_s){
		.len = len,
		.c = (char)c
	});
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
struct fna_read_ret_s fna_read_seq_4bit(
	struct fna_context_s *fna,
	kvec_uint8_t *v,
	uint8_t const *delim_table)
{
	int c;
	int64_t len = 0;
	while(1) {
		c = zfgetc(fna->fp);
		uint8_t type = delim_table[(uint8_t)c];
		if(type & 0x01) { break; }
		if(type != 0) { continue; }
		kv_push(*v, fna_encode_4bit(c)); len++;
	}

	fna->status = zfeof(fna->fp) ? FNA_EOF : FNA_SUCCESS;
	return((struct fna_read_ret_s){
		.len = len,
		.c = (char)c
	});
}

/**
 * @fn fna_read_seq_4bitpacked
 * @brief read seq until delim, with conv table
 */
static
struct fna_read_ret_s fna_read_seq_4bitpacked(
	struct fna_context_s *fna,
	kvec_uint8_t *v,
	uint8_t const *delim_table)
{
	int c;
	int64_t len = 0;
	uint64_t rem = 8;
	uint8_t arr = 0;
	/* 2x unrolled loop */
	while(1) {
		#define _fetch(_fna) ({ \
			char _c; \
			uint8_t _type; \
			while((_type = delim_table[(uint8_t)(_c = zfgetc(_fna->fp))]) != 0) { \
				if(_type & 0x01) { goto _fna_read_seq_4bitpacked_finish; } \
			} \
			_c; \
		}) \

		rem = 8;
		arr = (arr>>2) | (fna_encode_4bit(c = _fetch(fna))<<4); rem -= 4;
		arr = (arr>>2) | (fna_encode_4bit(c = _fetch(fna))<<4);
		kv_push(*v, arr); len += 2;

		#undef _fetch
	}
_fna_read_seq_4bitpacked_finish:;
	kv_push(*v, arr>>rem); len += (8 - rem) / 4;

	fna->status = zfeof(fna->fp) ? FNA_EOF : FNA_SUCCESS;
	return((struct fna_read_ret_s){
		.len = len,
		.c = (char)c
	});
}

/**
 * @fn fna_read_head_fasta
 */
static
int fna_read_head_fasta(
	struct fna_context_s *fna)
{
	/* eat '>' at the head */
	fna_read_skip(fna, delim_fasta_seq);
	return(zfeof(fna->fp) ? FNA_EOF :  FNA_SUCCESS);
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
		.type = FNA_SEGMENT,
		.seq_encode = fna->seq_encode,
		.options = fna->options,
		.head_margin = fna->head_margin,
		.tail_margin = fna->tail_margin,
		.seq_head_margin = fna->seq_head_margin,
		.seq_tail_margin = fna->seq_tail_margin
	}));
	dump(kv_ptr(v), 64);

	/* parse name */
	int64_t name_len = fna_read_ascii(fna, &v, delim_line).len;	/* fasta header line must ends with '\n' */
	dump(kv_ptr(v), 64);

	/* parse seq */
	fna_seq_make_margin(&v, fna->seq_head_margin);
	int64_t seq_len = (fna->read_seq(fna, &v, delim_fasta_seq)).len;
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
	char *name = (char *)(r + 1);
	uint8_t *seq = (uint8_t *)(name + (name_len + 1) + r->seq_head_margin);
	uint8_t *qual = (uint8_t *)(seq + seq_len + 1);
	int64_t qual_len = 0;

	return(_fna_pack_segment(r,
		name, name_len,
		seq, seq_len,
		qual, qual_len));
}

/**
 * @fn fna_read_head_fastq
 */
static
int fna_read_head_fastq(
	struct fna_context_s *fna)
{
	/* eat '@' at the head */
	fna_read_skip(fna, delim_fastq_seq);
	return(zfeof(fna->fp) ? FNA_EOF :  FNA_SUCCESS);
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
		.type = FNA_SEGMENT,
		.seq_encode = fna->seq_encode,
		.options = fna->options,
		.head_margin = fna->head_margin,
		.tail_margin = fna->tail_margin,
		.seq_head_margin = fna->seq_head_margin,
		.seq_tail_margin = fna->seq_tail_margin
	}));

	/* parse name */
	int64_t name_len = fna_read_ascii(fna, &v, delim_line).len;

	/* parse seq */
	fna_seq_make_margin(&v, fna->seq_head_margin);
	int64_t seq_len = (fna->read_seq(fna, &v, delim_fastq_qual)).len;
	fna_seq_make_margin(&v, fna->seq_tail_margin);

	/* skip name */
	fna_read_skip(fna, delim_line);

	/* parse qual */
	int64_t qual_len = (((fna->options & FNA_SKIP_QUAL) == 0)
		? fna->read_seq(fna, &v, delim_fastq_seq)
		: fna_read_skip(fna, delim_fastq_seq)).len;

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
	char *name = (char *)(r + 1);
	uint8_t *seq = (uint8_t *)(name + (name_len + 1) + r->seq_head_margin);
	uint8_t *qual = (uint8_t *)(seq + seq_len + 1);

	return(_fna_pack_segment(r,
		name, name_len,
		seq, seq_len,
		qual, qual_len));
}

/**
 * @fn fna_read_head_fast5
 */
static
int fna_read_head_fast5(
	struct fna_context_s *fna)
{
#ifdef HAVE_HDF5
	/** not implemented yet */
#endif /* HAVE_HDF5 */
	return(FNA_EOF);
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
 * @fn fna_read_head_gfa
 */
static
int fna_read_head_gfa(
	struct fna_context_s *fna)
{
	kvec_uint8_t buf;
	kv_init(buf);

	debug("parse gfa header");

	/* read until '\n' */
	fna_read_ascii(fna, &buf, delim_line);

	/* check prefix */
	char const *prefix = "H\tVN:Z:";
	if(strncmp((char const *)kv_ptr(buf), prefix, strlen(prefix)) != 0) {
		debug("broken");
		return(fna->status = FNA_ERROR_BROKEN_FORMAT);
	}

	/* check version */
	uint64_t version = fna_parse_version_string(&((char const *)kv_ptr(buf))[strlen(prefix)]);

	debug("%llx", version);
	return((version >= 0x10000) ? FNA_SUCCESS : FNA_ERROR_UNSUPPORTED_VERSION);
}

/**
 * @fn fna_read_gfa_seq
 */
static _force_inline
struct fna_seq_intl_s *fna_read_gfa_seq(
	struct fna_context_s *fna)
{
	kvec_uint8_t v;
	kv_init(v);

	/* make margin at the head of seq */
	fna_seq_make_margin(&v, fna->head_margin);

	kv_pusha(struct fna_seq_intl_s, v, ((struct fna_seq_intl_s){
		.type = FNA_SEGMENT,
		.seq_encode = fna->seq_encode,
		.options = fna->options,
		.head_margin = fna->head_margin,
		.tail_margin = fna->tail_margin,
		.seq_head_margin = fna->seq_head_margin,
		.seq_tail_margin = fna->seq_tail_margin
	}));
	dump(kv_ptr(v), 64);

	/* parse name */
	int64_t name_len = fna_read_ascii(fna, &v, delim_gfa_field).len;
	dump(kv_ptr(v), 64);

	/* parse seq */
	fna_seq_make_margin(&v, fna->seq_head_margin);
	struct fna_read_ret_s ret = fna->read_seq(fna, &v, delim_gfa_field);
	int64_t seq_len = ret.len;
	dump(kv_ptr(v), 64);

	/* check if optional field remains */
	if(ret.c == '\t') {
		/* skip optional fields */
		fna_read_skip(fna, delim_line);
	}

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
	char *name = (char *)(r + 1);
	uint8_t *seq = (uint8_t *)(name + (name_len + 1) + r->seq_head_margin);
	uint8_t *qual = (uint8_t *)(seq + seq_len + 1);
	int64_t qual_len = 0;

	return(_fna_pack_segment(r,
		name, name_len,
		seq, seq_len,
		qual, qual_len));
}

/**
 * @fn fna_read_gfa_link
 */
static _force_inline
struct fna_seq_intl_s *fna_read_gfa_link(
	struct fna_context_s *fna)
{
	kvec_uint8_t v;
	kv_init(v);

	/* make margin at the head of seq */
	fna_seq_make_margin(&v, fna->head_margin);

	kv_pusha(struct fna_seq_intl_s, v, ((struct fna_seq_intl_s){
		.type = FNA_LINK,
		.seq_encode = fna->seq_encode,
		.options = fna->options,
		.head_margin = fna->head_margin,
		.tail_margin = fna->tail_margin,
		.seq_head_margin = fna->seq_head_margin,
		.seq_tail_margin = fna->seq_tail_margin
	}));

	/* parse from field */
	struct fna_read_ret_s ret_from = fna_read_ascii(fna, &v, delim_gfa_field);
	if(ret_from.c != '\t') {
		fna->status = FNA_ERROR_BROKEN_FORMAT;
		return(NULL);
	}

	/* direction */
	int64_t from_ori = (zfgetc(fna->fp) == '+') ? 1 : -1;
	if(zfgetc(fna->fp) != '\t') {
		fna->status = FNA_ERROR_BROKEN_FORMAT;
		return(NULL);
	}

	/* parse to field */
	struct fna_read_ret_s ret_to = fna_read_ascii(fna, &v, delim_gfa_field);
	if(ret_to.c != '\t') {
		fna->status = FNA_ERROR_BROKEN_FORMAT;
		return(NULL);
	}

	/* direction */
	int64_t to_ori = (zfgetc(fna->fp) == '+') ? 1 : -1;
	if(zfgetc(fna->fp) != '\t') {
		fna->status = FNA_ERROR_BROKEN_FORMAT;
		return(NULL);
	}

	/* parse to field */
	struct fna_read_ret_s ret_cig = fna_read_ascii(fna, &v, delim_gfa_field);

	/* check if optional field remains */
	if(ret_cig.c == '\t') {
		/* skip optional fields */
		fna_read_skip(fna, delim_line);
	}

	/* make margin at the tail */
	fna_seq_make_margin(&v, fna->tail_margin);

	/* finished, build links */
	struct fna_seq_intl_s *r = (struct fna_seq_intl_s *)(
		kv_ptr(v) + fna->head_margin);
	char *from = (char *)(r + 1);
	char *to = from + (ret_from.len + 1);
	char *cig = to + (ret_to.len + 1);

	return(_fna_pack_link(r,
		from, ret_from.len, from_ori,
		to, ret_to.len, to_ori,
		cig, ret_cig.len));
}

/**
 * @fn fna_read_gfa_cont
 */
static _force_inline
struct fna_seq_intl_s *fna_read_gfa_cont(
	struct fna_context_s *fna)
{
	/* fixme: ignoring containment information */
	return(NULL);
}

/**
 * @fn fna_read_gfa_path
 */
static _force_inline
struct fna_seq_intl_s *fna_read_gfa_path(
	struct fna_context_s *fna)
{
	/* fixme: ignoring path line */
	return(NULL);
}

/**
 * @fn fna_read_gfa
 *
 * @brief gfa parser driver
 */
static
struct fna_seq_intl_s *fna_read_gfa(
	struct fna_context_s *fna)
{
	int c;
	while((c = zfgetc(fna->fp)) != EOF) {

		/* eat tab after type character */
		if(zfgetc(fna->fp) != '\t') {
			fna->status = FNA_ERROR_BROKEN_FORMAT;
			return(NULL);
		}

		/* examine type */
		switch(c) {
			case 'S': return(fna_read_gfa_seq(fna));
			case 'L': return(fna_read_gfa_link(fna));

			case 'C':	/* fall throught to 'P' */
			case 'P':
			fna_read_skip(fna, delim_line);
			break;

			/*
			case 'C': return(fna_read_gfa_cont(fna));
			case 'P': return(fna_read_gfa_path(fna));
			*/
			
			default:
			/* broken broken broken */
			fna->status = FNA_ERROR_BROKEN_FORMAT;
			return(NULL);
		}

	}
	fna->status = FNA_EOF;
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
		if(s->type == FNA_SEGMENT) {

			/* segment */
			char const *name_base = (char *)(s + 1);
			if(s->segment.name != name_base) { free(s->segment.name); }
			
			uint8_t const *seq_base = (uint8_t const *)(
				name_base + s->seq_head_margin + s->segment.name_len + 1);
			if(s->segment.seq != seq_base) { free(s->segment.seq); }
			
			uint8_t const *qual_base = (uint8_t const *)(
				seq_base + s->seq_tail_margin + s->segment.seq_len + 1);
			if(s->segment.qual != qual_base) { free(s->segment.qual); }

		} else if(s->type == FNA_LINK) {

			/* link */
			char const *from_base = (char *)(s + 1);
			if(s->link.from != from_base) { free(s->link.from); }

			char const *to_base = from_base + s->link.from_len + 1;
			if(s->link.to != to_base) { free(s->link.to); }

			char const *cigar_base = to_base + s->link.to_len + 1;
			if(s->link.cigar != cigar_base) { free(s->link.cigar); }

		}

		/* free context */
		free((void *)((uint8_t *)s - s->head_margin));
	}
	return;
}

#if 0
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

	if(src->seq_encode != dst->seq_encode || src->seq_encode != src->seq_encode) {
		return;
	}
	if(dst->seq_encode == FNA_RAW) {
		(void)kv_pop(dst->seq);
	}
	dst->len += size;
	if(src->seq_encode == FNA_RAW) {
		kv_reserve(dst->seq, dst->len + 1);
		for(i = 0; i < size; i++) {
			kv_push(dst->seq, kv_at(src->seq, i));
		}
		kv_push(dst->seq, '\0');
	} else if(src->seq_encode == FNA_2BIT) {
		kv_reserve(dst->seq, dst->len);
		for(i = 0; i < size; i++) {
			kv_push(dst->seq, kv_at(src->seq, i));
		}
	} else if(src->seq_encode == FNA_2BITPACKED) {
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
	if(seq->segment.seq_encode == FNA_RAW) {
		kv_push(dup->seq, '\0');
	}

	fna_append((fna_seq_t *)dup, (fna_seq_t *)seq);

	dup->name = seq->segment.name;
	dup->name.a = strdup(seq->segment.name.a);
	dup->seq_encode = seq->segment.seq_encode;
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

	if(src->seq_encode != dst->seq_encode || src->seq_encode != src->seq_encode) {
		return;
	}

	if(dst->seq_encode == FNA_RAW) {
		(void)kv_pop(dst->seq);
	}
	dst->len += size;
	if(src->seq_encode == FNA_RAW) {
		kv_reserve(dst->seq, dst->len + 1);
		for(i = 0; i < size; i++) {
			kv_push(dst->seq, fna_base_comp(kv_at(src->seq, size - i - 1)));
		}
		kv_push(dst->seq, '\0');
	} else if(src->seq_encode == FNA_2BIT) {
		kv_reserve(dst->seq, dst->len);
		for(i = 0; i < size; i++) {
			kv_push(dst->seq, 0x03 - kv_at(src->seq, size - i - 1));
		}
	} else if(src->seq_encode == FNA_2BITPACKED) {
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
	if(seq->segment.seq_encode == FNA_RAW) {
		kv_push(rev->seq, '\0');
	}

	fna_append_revcomp((fna_seq_t *)rev, (fna_seq_t *)seq);

	rev->name = seq->segment.name;
	rev->name.a = strdup(seq->segment.name.a);
	rev->seq_encode = seq->segment.seq_encode;
	return((fna_seq_t *)rev);
}
#endif

/**
 * unittests
 */
#include <sys/stat.h>

unittest_config(
	.name = "fna",
	.depends_on = {"zf"}
);

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
		">  test2\n\nAAAA\n"
		">\ttest3\nACGT";
	assert(fdump(fasta_filename, fasta_content));
	assert(fcmp(fasta_filename, strlen(fasta_content), (uint8_t const *)fasta_content));

	fna_t *fna = fna_init(fasta_filename, NULL);
	assert(fna != NULL, "fna(%p)", fna);

	/* format detection (internal) */
	assert(fna->file_format == FNA_FASTA, "fna->file_format(%d)", fna->file_format);

	/* test0 */
	fna_seq_t *seq = fna_read(fna);
	assert(strcmp(seq->segment.name, "test0") == 0, "name(%s)", seq->segment.name);
	assert(strcmp((char const *)seq->segment.seq, "AAAA") == 0, "seq(%s)", (char const *)seq->segment.seq);
	assert(seq->segment.seq_len == 4, "len(%lld)", seq->segment.seq_len);
	fna_seq_free(seq);

	/* test1 */
	seq = fna_read(fna);
	assert(strcmp(seq->segment.name, "test1") == 0, "name(%s)", seq->segment.name);
	assert(strcmp((char const *)seq->segment.seq, "ATATCGCG") == 0, "seq(%s)", (char const *)seq->segment.seq);
	assert(seq->segment.seq_len == 8, "len(%lld)", seq->segment.seq_len);
	fna_seq_free(seq);

	/* test2 */
	seq = fna_read(fna);
	assert(strcmp(seq->segment.name, "test2") == 0, "name(%s)", seq->segment.name);
	assert(strcmp((char const *)seq->segment.seq, "AAAA") == 0, "seq(%s)", (char const *)seq->segment.seq);
	assert(seq->segment.seq_len == 4, "len(%lld)", seq->segment.seq_len);
	fna_seq_free(seq);

	/* test3 */
	seq = fna_read(fna);
	assert(strcmp(seq->segment.name, "test3") == 0, "name(%s)", seq->segment.name);
	assert(strcmp((char const *)seq->segment.seq, "ACGT") == 0, "seq(%s)", (char const *)seq->segment.seq);
	assert(seq->segment.seq_len == 4, "len(%lld)", seq->segment.seq_len);
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
		"@  test2\n\nAAAA\n+  test2\n\nNNNN\n"
		"@\ttest3\nACGT\n\n+\ttest3\nNNNN";
	assert(fdump(fastq_filename, fastq_content));
	assert(fcmp(fastq_filename, strlen(fastq_content), (uint8_t const *)fastq_content));

	fna_t *fna = fna_init(fastq_filename, NULL);
	assert(fna != NULL, "fna(%p)", fna);

	/* format detection (internal) */
	assert(fna->file_format == FNA_FASTQ, "fna->file_format(%d)", fna->file_format);

	/* test0 */
	fna_seq_t *seq = fna_read(fna);
	assert(strcmp(seq->segment.name, "test0") == 0, "name(%s)", seq->segment.name);
	assert(strcmp((char const *)seq->segment.seq, "AAAA") == 0, "seq(%s)", (char const *)seq->segment.seq);
	assert(seq->segment.seq_len == 4, "len(%lld)", seq->segment.seq_len);
	assert(strcmp((char const *)seq->segment.qual, "NNNN") == 0, "seq(%s)", (char const *)seq->segment.qual);
	assert(seq->segment.qual_len == 4, "len(%lld)", seq->segment.qual_len);
	fna_seq_free(seq);

	/* test1 */
	seq = fna_read(fna);
	assert(strcmp(seq->segment.name, "test1") == 0, "name(%s)", seq->segment.name);
	assert(strcmp((char const *)seq->segment.seq, "ATATCGCG") == 0, "seq(%s)", (char const *)seq->segment.seq);
	assert(seq->segment.seq_len == 8, "len(%lld)", seq->segment.seq_len);
	assert(strcmp((char const *)seq->segment.qual, "NNNNNNNN") == 0, "seq(%s)", (char const *)seq->segment.qual);
	assert(seq->segment.qual_len == 8, "len(%lld)", seq->segment.qual_len);
	fna_seq_free(seq);

	/* test2 */
	seq = fna_read(fna);
	assert(strcmp(seq->segment.name, "test2") == 0, "name(%s)", seq->segment.name);
	assert(strcmp((char const *)seq->segment.seq, "AAAA") == 0, "seq(%s)", (char const *)seq->segment.seq);
	assert(seq->segment.seq_len == 4, "len(%lld)", seq->segment.seq_len);
	assert(strcmp((char const *)seq->segment.qual, "NNNN") == 0, "seq(%s)", (char const *)seq->segment.qual);
	assert(seq->segment.qual_len == 4, "len(%lld)", seq->segment.qual_len);
	fna_seq_free(seq);

	/* test3 */
	seq = fna_read(fna);
	assert(strcmp(seq->segment.name, "test3") == 0, "name(%s)", seq->segment.name);
	assert(strcmp((char const *)seq->segment.seq, "ACGT") == 0, "seq(%s)", (char const *)seq->segment.seq);
	assert(seq->segment.seq_len == 4, "len(%lld)", seq->segment.seq_len);
	assert(strcmp((char const *)seq->segment.qual, "NNNN") == 0, "seq(%s)", (char const *)seq->segment.qual);
	assert(seq->segment.qual_len == 4, "len(%lld)", seq->segment.qual_len);
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
		"@  test2\n\nAAAA\n+  test2\n\nNNNN\n"
		"@\ttest3\nACGT\n\n+\ttest3\nNNNN";
	assert(fdump(fastq_filename, fastq_content));
	assert(fcmp(fastq_filename, strlen(fastq_content), (uint8_t const *)fastq_content));

	fna_t *fna = fna_init(fastq_filename, FNA_PARAMS(.options = FNA_SKIP_QUAL));
	assert(fna != NULL, "fna(%p)", fna);

	/* format detection (internal) */
	assert(fna->file_format == FNA_FASTQ, "fna->file_format(%d)", fna->file_format);

	/* test0 */
	fna_seq_t *seq = fna_read(fna);
	assert(strcmp(seq->segment.name, "test0") == 0, "name(%s)", seq->segment.name);
	assert(strcmp((char const *)seq->segment.seq, "AAAA") == 0, "seq(%s)", (char const *)seq->segment.seq);
	assert(seq->segment.seq_len == 4, "len(%lld)", seq->segment.seq_len);
	assert(strcmp((char const *)seq->segment.qual, "") == 0, "seq(%s)", (char const *)seq->segment.qual);
	assert(seq->segment.qual_len == 0, "len(%lld)", seq->segment.qual_len);
	fna_seq_free(seq);

	/* test1 */
	seq = fna_read(fna);
	assert(strcmp(seq->segment.name, "test1") == 0, "name(%s)", seq->segment.name);
	assert(strcmp((char const *)seq->segment.seq, "ATATCGCG") == 0, "seq(%s)", (char const *)seq->segment.seq);
	assert(seq->segment.seq_len == 8, "len(%lld)", seq->segment.seq_len);
	assert(strcmp((char const *)seq->segment.qual, "") == 0, "seq(%s)", (char const *)seq->segment.qual);
	assert(seq->segment.qual_len == 0, "len(%lld)", seq->segment.qual_len);
	fna_seq_free(seq);

	/* test2 */
	seq = fna_read(fna);
	assert(strcmp(seq->segment.name, "test2") == 0, "name(%s)", seq->segment.name);
	assert(strcmp((char const *)seq->segment.seq, "AAAA") == 0, "seq(%s)", (char const *)seq->segment.seq);
	assert(seq->segment.seq_len == 4, "len(%lld)", seq->segment.seq_len);
	assert(strcmp((char const *)seq->segment.qual, "") == 0, "seq(%s)", (char const *)seq->segment.qual);
	assert(seq->segment.qual_len == 0, "len(%lld)", seq->segment.qual_len);
	fna_seq_free(seq);

	/* test3 */
	seq = fna_read(fna);
	assert(strcmp(seq->segment.name, "test3") == 0, "name(%s)", seq->segment.name);
	assert(strcmp((char const *)seq->segment.seq, "ACGT") == 0, "seq(%s)", (char const *)seq->segment.seq);
	assert(seq->segment.seq_len == 4, "len(%lld)", seq->segment.seq_len);
	assert(strcmp((char const *)seq->segment.qual, "") == 0, "seq(%s)", (char const *)seq->segment.qual);
	assert(seq->segment.qual_len == 0, "len(%lld)", seq->segment.qual_len);
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
 * gfa parsing
 */
unittest()
{
	char const *gfa_filename = "test_gfa.gfa";
	char const *gfa_content =
		"H	VN:Z:1.0\n"
		"S	11	ACCTT\n"
		"S	12	TCAAGG\n"
		"S	13	CTTGATT\n"
		"L	11	+	12	-	4M\n"
		"L	12	-	13	+	5M\n"
		"L	11	+	13	+	3M\n"
		"P	14	11+,12-,13+	4M,5M\n"
		"S	15	CTTGATT\n";

	assert(fdump(gfa_filename, gfa_content));
	assert(fcmp(gfa_filename, strlen(gfa_content), (uint8_t const *)gfa_content));

	debug("start gfa");
	fna_t *fna = fna_init(gfa_filename, NULL);
	assert(fna != NULL, "fna(%p)", fna);

	/* format detection (internal) */
	assert(fna->file_format == FNA_GFA, "fna->file_format(%d)", fna->file_format);

	/* segment 11 */
	fna_seq_t *seq = fna_read(fna);
	assert(seq->type == FNA_SEGMENT, "type(%d)", seq->type);
	assert(strcmp(seq->segment.name, "11") == 0, "name(%s)", seq->segment.name);
	assert(strcmp((char const *)seq->segment.seq, "ACCTT") == 0, "seq(%s)", (char const *)seq->segment.seq);
	assert(seq->segment.seq_len == 5, "len(%lld)", seq->segment.seq_len);
	fna_seq_free(seq);

	/* segment 12 */
	seq = fna_read(fna);
	assert(seq->type == FNA_SEGMENT, "type(%d)", seq->type);
	assert(strcmp(seq->segment.name, "12") == 0, "name(%s)", seq->segment.name);
	assert(strcmp((char const *)seq->segment.seq, "TCAAGG") == 0, "seq(%s)", (char const *)seq->segment.seq);
	assert(seq->segment.seq_len == 6, "len(%lld)", seq->segment.seq_len);
	fna_seq_free(seq);

	/* segment 13 */
	seq = fna_read(fna);
	assert(seq->type == FNA_SEGMENT, "type(%d)", seq->type);
	assert(strcmp(seq->segment.name, "13") == 0, "name(%s)", seq->segment.name);
	assert(strcmp((char const *)seq->segment.seq, "CTTGATT") == 0, "seq(%s)", (char const *)seq->segment.seq);
	assert(seq->segment.seq_len == 7, "len(%lld)", seq->segment.seq_len);
	fna_seq_free(seq);

	/* link 11 to 12 */
	seq = fna_read(fna);
	assert(seq->type == FNA_LINK, "type(%d)", seq->type);
	assert(strcmp(seq->link.from, "11") == 0, "from(%s)", seq->link.from);
	assert(seq->link.from_ori == 1, "from_ori(%d)", seq->link.from_ori);
	assert(strcmp(seq->link.to, "12") == 0, "to(%s)", seq->link.to);
	assert(seq->link.to_ori == -1, "to_ori(%d)", seq->link.to_ori);
	assert(strcmp((char const *)seq->link.cigar, "4M") == 0, "cigar(%s)", (char const *)seq->link.cigar);
	fna_seq_free(seq);

	/* link 12 to 13 */
	seq = fna_read(fna);
	assert(seq->type == FNA_LINK, "type(%d)", seq->type);
	assert(strcmp(seq->link.from, "12") == 0, "from(%s)", seq->link.from);
	assert(seq->link.from_ori == -1, "from_ori(%d)", seq->link.from_ori);
	assert(strcmp(seq->link.to, "13") == 0, "to(%s)", seq->link.to);
	assert(seq->link.to_ori == 1, "to_ori(%d)", seq->link.to_ori);
	assert(strcmp((char const *)seq->link.cigar, "5M") == 0, "cigar(%s)", (char const *)seq->link.cigar);
	fna_seq_free(seq);

	/* link 11 to 13 */
	seq = fna_read(fna);
	assert(seq->type == FNA_LINK, "type(%d)", seq->type);
	assert(strcmp(seq->link.from, "11") == 0, "from(%s)", seq->link.from);
	assert(seq->link.from_ori == 1, "from_ori(%d)", seq->link.from_ori);
	assert(strcmp(seq->link.to, "13") == 0, "to(%s)", seq->link.to);
	assert(seq->link.to_ori == 1, "to_ori(%d)", seq->link.to_ori);
	assert(strcmp((char const *)seq->link.cigar, "3M") == 0, "cigar(%s)", (char const *)seq->link.cigar);
	fna_seq_free(seq);

	/* skip path line, not implemented yet */

	/* segment 15 */
	seq = fna_read(fna);
	assert(seq->type == FNA_SEGMENT, "type(%d)", seq->type);
	assert(strcmp(seq->segment.name, "15") == 0, "name(%s)", seq->segment.name);
	assert(strcmp((char const *)seq->segment.seq, "CTTGATT") == 0, "seq(%s)", (char const *)seq->segment.seq);
	assert(seq->segment.seq_len == 7, "len(%lld)", seq->segment.seq_len);
	fna_seq_free(seq);

	/* term */
	seq = fna_read(fna);
	assert(seq == NULL, "%p", seq);
	fna_seq_free(seq);

	fna_close(fna);

	/** cleanup files */
	remove(gfa_filename);
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

/* format detection */
unittest()
{
	char const *gfa_filename = "test_gfa.txt";
	char const *gfa_content =
		"H	VN:Z:1.0\n"
		"S	11	ACCTT\n"
		"S	12	TCAAGG\n"
		"S	13	CTTGATT\n"
		"L	11	+	12	-	4M\n"
		"L	12	-	13	+	5M\n"
		"L	11	+	13	+	3M\n"
		"P	14	11+,12-,13+	4M,5M\n"
		"S	15	CTTGATT\n";

	assert(fdump(gfa_filename, gfa_content));
	assert(fcmp(gfa_filename, strlen(gfa_content), (uint8_t const *)gfa_content));

	fna_t *fna = fna_init(gfa_filename, NULL);
	assert(fna != NULL, "fna(%p)", fna);

	/* format detection (internal) */
	assert(fna->file_format == FNA_GFA, "fna->file_format(%d)", fna->file_format);

	fna_close(fna);

	/** cleanup file */
	remove(gfa_filename);
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

/**
 * end of fna.c
 */
