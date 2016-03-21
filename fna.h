
/**
 * @file fna.h
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
 *     fna_t *fna_init(char const *path, int pack);
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
 *   fna_seq_t (alias to struct _seq): sequence container.
 *     name: sequence name container (kvec_t(char) instance)
 *       name.a: pointer to the sequence name (null-terminated ASCII).
 *     seq: sequence container (kvec_t(uint8_t) or kpvec_t(uint8_t) instance)
 *       seq.a: pointer to the sequence (null-terminated when fna->pack == FNA_RAW)
 */
#ifndef _FNA_H_INCLUDED
#define _FNA_H_INCLUDED

#include <stdint.h>

/**
 * @enum fna_flag_encode
 */
enum fna_flag_encode {
	FNA_RAW 		= 0,
	FNA_ASCII 		= 0,
	FNA_2BIT		= 1,
	FNA_2BITPACKED 	= 2
};

/**
 * @enum fna_format
 *
 * @brief format flag constant
 */
enum fna_format {
	FNA_UNKNOWN		= 0,
	FNA_FASTA		= 1,
	FNA_FASTQ		= 2,
	FNA_FAST5		= 3
};

/**
 * @enum fna_status
 */
enum fna_status {
	FNA_SUCCESS					= 0,
	FNA_ERROR_FILE_OPEN			= 1,
	FNA_ERROR_UNKNOWN_FORMAT	= 2,
	FNA_ERROR_BROKEN_FORMAT		= 3,
	FNA_ERROR_OUT_OF_MEM		= 4,
	FNA_EOF 					= 0
};

/**
 * @struct fna_params_s
 * @brief options
 */
struct fna_params_s {
	int32_t file_format;		/** see enum fna_format */
	int32_t seq_encode;			/** see enum fna_flag_encode */
	uint16_t head_margin;		/** margin at the head of fna_seq_t */
	uint16_t tail_margin;		/** margin at the tail of fna_seq_t	*/
};
typedef struct fna_params_s fna_params_t;

#define FNA_PARAMS(...)			( &((struct fna_params_s const) { __VA_ARGS__ }) )

/**
 * @struct fna_s
 *
 * @brief a struct for fna context
 */
struct fna_s {
	char *path;
	int32_t file_format;		/** see enum fna_format */
	int32_t seq_encode;			/** see enum fna_flag_encode */
	int32_t status;				/** see enum fna_status */
	uint32_t reserved[16];
};
typedef struct fna_s fna_t;

/**
 * @struct fna_seq_s
 *
 * @brief a struct to contain parsed sequence.
 */
struct fna_seq_s {
	uint64_t reserved1[2];
	char *name;					/** sequence name */
	uint64_t reserved2[2];
	uint8_t *seq;				/** sequence */
	int64_t len;				/** sequence length */
	int32_t seq_encode;			/** one of fna_flag_encode */
	uint32_t reserved3;
};
typedef struct fna_seq_s fna_seq_t;

/**
 * @fn fna_init
 *
 * @brief create a sequence reader context
 *
 * @param[in] path : a path to file to open.
 * @param[in] pack : see struct fna_params_s
 *
 * @return a pointer to the context, NULL if an error occurred (may be invalid path or invalid format)
 */
fna_t *fna_init(char const *path, fna_params_t const *params);

/**
 * @fn fna_close
 *
 * @brief clean up sequence reader context
 */
void fna_close(fna_t *fna);

/**
 * @fn fna_read
 *
 * @brief read a sequence
 *
 * @param[in] fna : a pointer to the context
 *
 * @return a pointer to a sequence object, NULL if the file pointer reached the end.
 */
fna_seq_t *fna_read(fna_t *fna);

/**
 * @fn fna_append
 *
 * @brief concatenate src sesquence after dst sequence
 */
void fna_append(fna_seq_t *dst, fna_seq_t const *src);

/**
 * @fn fna_duplicate
 *
 * @brief duplicate sequence
 */
fna_seq_t *fna_duplicate(fna_seq_t const *seq);

/**
 * @fn fna_append_revcomp
 *
 * @brief append reverse complemented sequence after the given sequence
 */
void fna_append_revcomp(fna_seq_t *seq, fna_seq_t const *src);

/**
 * @fn fna_revcomp
 *
 * @brief make reverse complemented sequence
 */
fna_seq_t *fna_revcomp(fna_seq_t const *seq);

/**
 * @fn fna_seq_free
 *
 * @brief clean up sequence object
 */
void fna_seq_free(fna_seq_t *seq);

#endif /* #ifndef _FNA_H_INCLUDED */
/**
 * end of fna.h
 */
