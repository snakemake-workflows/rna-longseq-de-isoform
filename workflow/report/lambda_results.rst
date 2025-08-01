`Lambda <https://github.com/seqan/lambda>`_ was used to perform local alignment of the top differentially expressed genes, selected based on a log2FoldChange threshold of ``{{ snakemake.config["deseq2"]["lfc_null"] }}``.

Lambda outputs up to ``num_matches = {{ snakemake.config["protein_annotation"]["num_matches"] }}`` results per query. Each result includes twelve columns:

- Query Seq-id
- Subject Seq-id
- Percentage of identical matches
- Alignment length
- Number of mismatches
- Number of gap openings
- Start of alignment in query
- End of alignment in query
- Start of alignment in subject
- End of alignment in subject
- Expect value
- Bit score