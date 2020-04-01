# CtgRef-CNV
CNV calling method based on assembly and read depth of NGS data

1. De novo assembly of each accession based on the high depth (~50X) of NGS reads.

2. Mapping NGS reads to its assemblied genome (contig length >= 1000 bp) to calculate the depth along contigs.

3. Perfom the GC correction of the depth.

4. Alignment between the RefSeq and the assemblied genome.

5. Depth transformation depending on the results of Step 3 and 4: from contig-based depth to reference-based depth.

6. CNV calling and filtering.
