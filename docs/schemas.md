# Schemas


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
`rnaN` attributes:

|Field|Type|Description|
|-|:-:|-|
|chrom|O|abc|
|start|int64||
|end|int64||
|total_contacts|int64||
|genic_contacts|int64||
|cis_contacts|int64||
|trans_contacts|int64||
|eligible|bool||
|cis_factor|float||
|cis_start|int64||
|trans_bin_size|int64||

# rdc schema

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
            - signal_count | float
            - bg_count | float
            - raw_bg_prob | float
            - scaling_factor | float
            - bg_prob | float
            - impute | ?
            - signal_prob | float
            - fc | float
            - pvalue | float
            - qvalue | float
```
`rnaN` attributes
- coordinates
- contacts stats
- binning params
- scaling params

`rnaN/chrN` attributes
- bins: "linear" | "geom"
- type: "cis" | "trans"