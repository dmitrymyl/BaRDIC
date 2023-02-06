# BaRDIC file formats

## dnah5

`dnah5` layout:
```
/
- chrom_sizes/
    - chrom | O
    - size | int64
- dna_parts/
    - rnaN/
        - chromN/
            - start | int64
            - end | int64
```

`dnah5` attributes:
|Field|Type|Description|
|-|:-:|-|
|are_binsizes_selected|bool|Whether bin sizes are selected for each RNA and corresponding data is recorded in the file|

`rnaN` attributes:

|Field|Type|Description|
|-|:-:|-|
|chrom|O|RNA gene chromosome|
|start|int64|RNA gene start|
|end|int64|RNA gene end|
|total_contacts|int64|Total number of RNA contacts|
|genic_contacts|int64|Number of RNA contacts inside its gene|
|cis_contacts|int64|Number of cis RNA contacts: on the RNA's origin chromosome but outside RNA gene|
|trans_contacts|int64|Number of trans RNA contacts: on all chromosomes except for RNA's origin one|
|eligible|bool|Whether this RNA has enough contacts for further processing|
|cis_factor|float64|Cis factor value for binning|
|cis_start|int64|Initial cis bin size|
|trans_bin_size|int64|Trans bin size|

## rdc

`rdc` layout:

```
/
- chrom_sizes/
    - chrom | O
    - size | int64
- background
    - chrN
        - start | int64
        - end | int64
        - count | int64
- pixels
    - rnaN
        - chrN
            - start | int64
            - end | int64
            - signal_count | float64
            - bg_count | float64
            - raw_bg_prob | float64
            - scaling_factor | float64
            - bg_prob | float64
            - impute | bool
            - signal_prob | float64
            - fc | float64
            - pvalue | float64
            - qvalue | float64
```

`rdc` attributes:
|Field|Type|Description|
|-|:-:|-|
|is_scaling_fitted|bool|Whether scaling is estimated and RNAs background levels of interactions are rescaled|
|are_peaks_estimated|bool|Whether p-values for peaks are estimated|


`rnaN` attributes:
|Field|Type|Description|
|-|:-:|-|
|chrom|O|RNA gene chromosome|
|start|int64|RNA gene start|
|end|int64|RNA gene end|
|total_contacts|int64|Total number of RNA contacts|
|genic_contacts|int64|Number of RNA contacts inside its gene|
|cis_contacts|int64|Number of cis RNA contacts: on the RNA's origin chromosome but outside RNA gene|
|trans_contacts|int64|Number of trans RNA contacts: on all chromosomes except for RNA's origin one|
|eligible|bool|Whether this RNA has enough contacts for further processing|
|cis_factor|float64|Cis factor value for binning|
|cis_start|int64|Initial cis bin size|
|trans_bin_size|int64|Trans bin size|
|scaling_spline_t|Float array|A vector of knots of a scaling B-spline|
|scaling_spline_c|Float array|A vector of scaling B-spline coefficients|
|scaling_spline_k|int64|A degree of a scaling B-spline|
