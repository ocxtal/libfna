# libfna

Libfna is a lightweight FASTA / FASTQ parser library written in pure C99, developed as a submodule of the comb aligner.

## Build

Python (2.7 or 3.x) is required to run the build script written in waf.

```
$ ./waf configure
$ ./waf build
```

## Usage

```
/* open file */
fna_t *fna = fna_init("sequence.fa", NULL);

/* read a sequence from the head */
fna_seq_t *seq = NULL;
while((seq = fna_read(fna)) != NULL) {
	printf("name = %s, seq = %s, len = %lld",
		seq->name, seq->seq, seq->len);
	
	/* cleanup the sequence */
	fna_seq_free(seq);
}

/* cleanup the context */
fna_close(fna);
```

## Lisence

MIT