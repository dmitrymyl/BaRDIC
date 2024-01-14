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
usage: bardic run [-h] [-f [{narrowPeak,bed}]] [-s [score_field]] [-mcon [int]]
                  [-tmin [int]] [-tmax [int]] [-tstep [int]] [-cmin [float]]
                  [-cmax [float]] [-cstep [float]] [-cstart [int]]
                  [-tol [float]] [-w [float]] [-bs [int]] [-bt [{rnas,custom}]]
                  [-i [float]] [-d [int]] [-mt [float]] [-ns] [-nr]
                  [-fv [numeric]] [-q [float]] [-qt [{global,rna}]] [-c [int]]
                  dnaparts annotation chromsizes bgdata outdir

Run pipeline with a single command.

optional arguments:
  -h, --help            show this help message and exit

Input:
  dnaparts              BED6 file with coordinates of DNA parts. Names of
                        corresponding RNAs are in the "name" column.
  annotation            RNA annotation in BED format.
  chromsizes            If filename, then it is a UCSC headerless chromsizes
                        file; if genome abbreviation, then will fetch chromsizes
                        from UCSC
  bgdata                A file with data on background. If --bgtype="rnas", this
                        is a file with a list of RNAs with one RNA name per
                        line. If --bgtype="custom", this is a bedGraph file with
                        background signal in equally-sized bins.

Output:
  outdir                Output directory name.
  -f [{narrowPeak,bed}], --format [{narrowPeak,bed}]
                        Output peaks file format. (default: narrowPeak)
  -s [score_field], --score [score_field]
                        If --format=bed, which value to fill the score field
                        with. If int, will fill every peak score with it; if
                        str, will take corresponding values from the column in
                        RDC (choices: bg_count, raw_bg_prob, scaling_factor,
                        bg_prob, signal_count, signal_prob, impute, fc, pvalue,
                        qvalue) (default: 0)

Binsize selection parameters:
  -mcon [int], --min_contacts [int]
                        Minimal number of contacts to consider an RNA. Any RNA
                        with less contacts will be discarded from further
                        processing. (default: 1000)
  -tmin [int], --trans_min [int]
                        Minimal trans bin size. (default: 10000)
  -tmax [int], --trans_max [int]
                        Maximal trans bin size. (default: 1000000)
  -tstep [int], --trans_step [int]
                        Step for increasing trans bin size. (default: 1000)
  -cmin [float], --cis_min [float]
                        Minimal cis factor. (default: 1.1)
  -cmax [float], --cis_max [float]
                        Maximal cis factor. (default: 2.0)
  -cstep [float], --cis_step [float]
                        Step for inreasing cis factor. (default: 0.01)
  -cstart [int], --cis_start [int]
                        Starting cis bin size. (default: 5000)
  -tol [float], --tolerance [float]
                        Maximal absolute difference between two consecutive cost
                        function values to consider optimization converged.
                        (default: 0.01)
  -w [float], --window [float]
                        Window size to average cost function values over.
                        (default: 1)

Background parameters:
  -bs [int], --binsize [int]
                        Bin size of the background track. (default: 1000)
  -bt [{rnas,custom}], --bgtype [{rnas,custom}]
                        Type of backround. If "rnas", then will calculate
                        background from trans-contacts of RNAs supplied as
                        "bgdata". If "custom", will use bedgraph track provided
                        as "bgdata". (default: rnas)

RDC creation parameters:
  -i [float], --ifactor [float]
                        Imputation factor: if background coverage of a bin is 0,
                        this value is a multiplier of an average background
                        coverage to impute zero background coverage. (default:
                        0.01)

Scaling parameters:
  -d [int], --degree [int]
                        Spline degree. (default: 3)
  -mt [float], --max_threshold [float]
                        Maximal binomial test p-value to consider a point as an
                        outlier in a spline refinement procedure. (default:
                        0.05)
  -ns, --no_scaling     If included, do not estimate scaling. (default: False)
  -nr, --no_refine      If included, do not apply a spline refinement procedure.
                        (default: False)
  -fv [numeric], --fill_value [numeric]
                        Fold-change fill ratio in case of 0/0. (default: 1)

Peaks parameters:
  -q [float], --qval_threshold [float]
                        BH q-value threshold to consider bin a peak. (default:
                        0.05)
  -qt [{global,rna}], --qval_type [{global,rna}]
                        BH q-value type to use for peak calling. If "global"
                        (default), will use q-values calculated for all RNAs; if
                        "rna", will use q-values calculated for each RNA
                        separately. (default: global)

Processing:
  -c [int], --cores [int]
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