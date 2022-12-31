# FASTK_tools

Some programs using [FASTK](https://github.com/thegenemyers/FASTK) and [DAZZ_DB](https://github.com/thegenemyers/DAZZ_DB) libraries by Gene Myers.

```
CoverRate [-l<int(1)>] [-r<int(32767)>] <source_root>[.prof]
```

- For each read, `CoverRate` calculates the number of and the cover rate by k-mers whose multiplicities are within a specific range. `-l` and `-r` specify the minimum and the maximum value of the k-mer multiplicity, respectively.
- This can be used for:
  1. Estimation of per-read error rate using, for example, `-r5`. Note that the cover rate in this case represents the ratio of erroneous **k-mers**, not erroneous **bases**. One erroneous base results in ~k consecutive erroneous k-mers. If most erroneous k-mers do not overlap (like HiFi reads), the value of the erroneous k-mer rate divided by k is a good estimate of the error rate.
  2. Estimation of per-read repetitive content rate using, for example, `-l100`. Unlike sequencing errors, repeats should be defined on not single bases but k-mers, and thus the cover rate in this case directly represents the ratio of repetitive contents.
  3. Rough extraction of allele-specific reads using, for example, `-l15 -r25`, given the haploid sequencing depth is 20x. These values should be determined based on the histogram.

```
Classifier <source_root> <diplo_depth<int>>
```

- Both `<source_root>.db` and `<source_root>.prof` must exist, and these should be made from homopolymer compressed sequences. Current implementation requires that `<source_root>.db` is a database of homopolymer compressed reads, if one uses hoco reads.
- `diplo_depth` is the global avegarage depth per diploid. It will be automatically determined in the near future.
- Outputs `.<sequence_root>.class.anno` and `.<sequence_root>.class.data` as a DAZZ_DB track (0: Error, 1: Repeat, 2: Haplo, 3: Diplo).
