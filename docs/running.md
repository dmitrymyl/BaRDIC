# How to run

The easiest way to run BaRDIC is to launch it from the command line. There will be an option for python API in the future.

## Input files

BaRDIC needs only 4 files as input:

|File|Explanation|Extension|
|-|-|-|
|dnaparts|DNA parts of RNA-DNA contacts as BED6 file with names of RNAs in the `name` field. `score` and `strand` fields are not relevant.|bed|
|annotation|Genome annotation that contains coordinates of RNAs genes in BED6 format.|bed|
|chromsizes|Either UCSC genome code name or two-column headerless tab-delimeted file with chromosome name and chromosome length in each row.|txt or genome code name|
|bg_rnas|List of RNAs names which DNA parts of contacts to use for background estimation, one RNA name per line.|txt|

## Command line interface

BaRDIC is a multi-command application. The list of commands is available with `bardic -h`:

```{bash}
usage: bardic [-h] SUBCOMMAND ...

Binomial RNA-DNA interaction caller.

optional arguments:
  -h, --help  show this help message and exit

Subcommands:
  SUBCOMMAND
    bed2h5    convert DNA parts to custom HDF5 format
    binsizes  Select bin size for each RNA and save it into dnah5 file
    background
              Create a bedGraph background track from DNA parts of selected RNAs
    makerdc   Create RDC file from dnah5 DNA parts and bedGraph background track.
    scaling   Estimate scaling by fitting splines and adjust background probabilities in RDC file.
    peaks     Estimate significance and fetch peaks at specified FDR level.
    run       Run pipeline with a single command.
```

The command `bardic run` will run the whole pipeline, while other commands launch single steps of the pipeline. 

### `bardic run` &mdash; run pipeline with single command

`bardic run` requires all input files as specified above. Additionally, algorithm can be tuned with a diverse set of parameters:

```{bash}
usage: bardic run [-h] [-f [{narrowPeak,bed}]] [-s [SCORE]] [-mcon [N_CONTACTS]] [-tmin [TRANS_MIN]] [-tmax [TRANS_MAX]] [-tstep [TRANS_STEP]]
                  [-cmin [CIS_MIN]] [-cmax [CIS_MAX]] [-cstep [CIS_STEP]] [-cstart [CIS_START]] [-tol [TOLERANCE]] [-w [W]] [-bs [BINSIZE]] [-i [IFACTOR]]
                  [-d [DEGREE]] [-mt [MAX_THRESHOLD]] [-nr] [-fv [FILL_VALUE]] [-q [QVAL_THRESHOLD]] [-c [N_CORES]]
                  dnaparts annotation chromsizes bg_rnas outdir

Run pipeline with a single command.

optional arguments:
  -h, --help            show this help message and exit

Input:
  dnaparts              BED6 file with coordinates of DNA parts. Names of corresponding RNAs are in the "name" column.
  annotation            RNA annotation in BED format.
  chromsizes            If filename, then it is a UCSC headerless chromsizes file; if genome abbreviation, then will fetch chromsizes from UCSC
  bg_rnas               A file with a list of RNAs with one RNA name per line.

Output:
  outdir                Output directory name.
  -f [{narrowPeak,bed}], --format [{narrowPeak,bed}]
                        Output peaks file format. (default: narrowPeak)
  -s [SCORE], --score [SCORE]
                        If --format=bed, which value to fill the score field with. If int, will fill every peak score with it; if str, will take
                        corresponding values from the column in RDC (choices: bg_count, raw_bg_prob, scaling_factor, bg_prob, signal_count, signal_prob,
                        impute, fc, pvalue, qvalue) (default: 0)

Binsize selection parameters:
  -mcon [N_CONTACTS], --min_contacts [N_CONTACTS]
                        Minimal number of contacts to consider an RNA. Any RNA with less contacts will be discarded from further processing. (default:
                        1000)
  -tmin [TRANS_MIN], --trans_min [TRANS_MIN]
                        Minimal trans bin size. (default: 10000)
  -tmax [TRANS_MAX], --trans_max [TRANS_MAX]
                        Maximal trans bin size. (default: 1000000)
  -tstep [TRANS_STEP], --trans_step [TRANS_STEP]
                        Step for increasing trans bin size. (default: 1000)
  -cmin [CIS_MIN], --cis_min [CIS_MIN]
                        Minimal cis factor. (default: 1.1)
  -cmax [CIS_MAX], --cis_max [CIS_MAX]
                        Maximal cis factor. (default: 2.0)
  -cstep [CIS_STEP], --cis_step [CIS_STEP]
                        Step for inreasing cis factor. (default: 0.01)
  -cstart [CIS_START], --cis_start [CIS_START]
                        Starting cis bin size. (default: 5000)
  -tol [TOLERANCE], --tolerance [TOLERANCE]
                        Maximal absolute difference between two consecutive cost function values to consider optimization converged. (default: 0.01)
  -w [W], --window [W]  Window size to average cost function values over. (default: 1)

Background parameters:
  -bs [BINSIZE], --binsize [BINSIZE]
                        Bin size of the background track. (default: 1000)

RDC creation parameters:
  -i [IFACTOR], --ifactor [IFACTOR]
                        Imputation factor: if background coverage of a bin is 0, this value is a multiplier of an average background coverage to impute
                        zero background coverage. (default: 0.01)

Scaling parameters:
  -d [DEGREE], --degree [DEGREE]
                        Spline degree. (default: 3)
  -mt [MAX_THRESHOLD], --max_threshold [MAX_THRESHOLD]
                        Maximal binomial test p-value to consider a point as an outlier in a spline refinement procedure. (default: 0.05)
  -nr, --no_refine      If included, do not apply a spline refinement procedure. (default: False)
  -fv [FILL_VALUE], --fill_value [FILL_VALUE]
                        Fold-change fill ratio in case of 0/0. (default: 1)

Peaks parameters:
  -q [QVAL_THRESHOLD], --qval_threshold [QVAL_THRESHOLD]
                        BH q-value threshold to consider bin a peak. (default: 0.05)

Processing:
  -c [N_CORES], --cores [N_CORES]
                        Maximal number of cores to use. (default: 1)
```

The output of `bardic run` will be put into a single directory. It consists of several files:

|File name|Description|File type|
|-|-|-|
|DnaDataset.dnah5|All DNA parts of RNA-DNA contacts compactly stored with binning parameters for each RNA. Refer to [schemas.md](./schemas.md) for the specification.|Custom HDF5 storage|
|background.bedGraph|A genomic track with background contacts.|bedGraph|
|contacts.rdc|Compact storage of binned contact profiles for each RNA along with the background track. Refer to [schemas.md](./schemas.md) for the specification.|Custom HDF5 storage|
|selection.tsv|A table with statistics on bin size selection for each RNA.|Tab delimited file with header|
|peaks.narrowPeak or peaks.bed|Resulting peaks.|Tab-delimited file (narrowPeak or bed)|