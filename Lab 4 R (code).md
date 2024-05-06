# Lab 4: plotting GO terms

This is a quick journal to show how to combine multiple GO-term files to create a single visualization that looks like the one found in Figure 6 of Van Weringh et al. (https://www.biorxiv.org/content/10.1101/2021.04.15.439991v1.abstract)

There are 8 files to import, all ending with the suffix: `_up_GO.txt`
- K60dvsw, K40dvsw, K20dvsw, Krwvsw, L60dvsw, L40dvsw, L20dvsw, Lrwvsw

Goals:
1. Determine the maximum and minimum value across all 8 sets. This will be used to set the range on the bubble/strip chart combination to come.
2. Create a bubble/strip chart. Y-axis GO terms, X axis each file set. Y-axis may be ordered by values in column 1 (data file "keep_formakinggofigures.txt")

Numbered steps below correspond to Lab Manual step numbers.

***
# 2. Loading packages


```R
library(tidyverse)
library(repr)
library(viridis)
```

    Warning message:
    “replacing previous import ‘ellipsis::check_dots_unnamed’ by ‘rlang::check_dots_unnamed’ when loading ‘tibble’”
    Warning message:
    “replacing previous import ‘ellipsis::check_dots_used’ by ‘rlang::check_dots_used’ when loading ‘tibble’”
    Warning message:
    “replacing previous import ‘ellipsis::check_dots_empty’ by ‘rlang::check_dots_empty’ when loading ‘tibble’”
    Registered S3 methods overwritten by 'tibble':
      method     from  
      format.tbl pillar
      print.tbl  pillar
    
    ── [1mAttaching packages[22m ─────────────────────────────────────── tidyverse 1.3.0 ──
    
    [32m✔[39m [34mggplot2[39m 3.3.0     [32m✔[39m [34mpurrr  [39m 1.0.1
    [32m✔[39m [34mtibble [39m 3.0.1     [32m✔[39m [34mdplyr  [39m 1.1.0
    [32m✔[39m [34mtidyr  [39m 1.3.0     [32m✔[39m [34mstringr[39m 1.5.0
    [32m✔[39m [34mreadr  [39m 1.3.1     [32m✔[39m [34mforcats[39m 0.5.0
    
    ── [1mConflicts[22m ────────────────────────────────────────── tidyverse_conflicts() ──
    [31m✖[39m [34mdplyr[39m::[32mfilter()[39m masks [34mstats[39m::filter()
    [31m✖[39m [34mdplyr[39m::[32mlag()[39m    masks [34mstats[39m::lag()
    
    Loading required package: viridisLite
    


***
# 3. Loading data

The first thing we'll do is load the data using the `lapply()` method to create a list of data frames. In this function, we'll use the file names and append the `_up_GO.txt` suffix.

Within our `lapply()` call, the function we'll apply to each data frame includes a `mutate()` step, which will add a new column containing the fileName prefix. We'll use this later to differentiate between our data sets.


```R
# Create a vector with your file names
fileNames <- c("K60dvsw", "K40dvsw", "K20dvsw", "Krwvsw", "L60dvsw", "L40dvsw", "L20dvsw", "Lrwvsw")

# While we import the files, we'll add a column to denote the dataset that it comes from. 
# This will help us put the data into a long format
data.list <- lapply(fileNames, FUN = function(x) mutate(read_tsv(paste0(x, "_up_GO.txt")), dataSet = x))
names(data.list) <- fileNames

# What is the structure of the data set
str(data.list, give.attr = FALSE)
# Note that L60dvsw_up_GO.txt has no data in it.
                    
# How many rows are there across all data sets?
sum(unlist(lapply(data.list, FUN = function(x) dim(x)[1])))
```

    Parsed with column specification:
    cols(
      GO_acc = [31mcol_character()[39m,
      term_type = [31mcol_character()[39m,
      Term = [31mcol_character()[39m,
      queryitem = [32mcol_double()[39m,
      querytotal = [32mcol_double()[39m,
      bgitem = [32mcol_double()[39m,
      bgtotal = [32mcol_double()[39m,
      pvalue = [32mcol_double()[39m,
      FDR = [32mcol_double()[39m,
      entries = [31mcol_character()[39m
    )
    
    Parsed with column specification:
    cols(
      GO_acc = [31mcol_character()[39m,
      term_type = [31mcol_character()[39m,
      Term = [31mcol_character()[39m,
      queryitem = [32mcol_double()[39m,
      querytotal = [32mcol_double()[39m,
      bgitem = [32mcol_double()[39m,
      bgtotal = [32mcol_double()[39m,
      pvalue = [32mcol_double()[39m,
      FDR = [32mcol_double()[39m,
      entries = [31mcol_character()[39m
    )
    
    Parsed with column specification:
    cols(
      GO_acc = [31mcol_character()[39m,
      term_type = [31mcol_character()[39m,
      Term = [31mcol_character()[39m,
      queryitem = [32mcol_double()[39m,
      querytotal = [32mcol_double()[39m,
      bgitem = [32mcol_double()[39m,
      bgtotal = [32mcol_double()[39m,
      pvalue = [32mcol_double()[39m,
      FDR = [32mcol_double()[39m,
      entries = [31mcol_character()[39m
    )
    
    Parsed with column specification:
    cols(
      GO_acc = [31mcol_character()[39m,
      term_type = [31mcol_character()[39m,
      Term = [31mcol_character()[39m,
      queryitem = [32mcol_double()[39m,
      querytotal = [32mcol_double()[39m,
      bgitem = [32mcol_double()[39m,
      bgtotal = [32mcol_double()[39m,
      pvalue = [32mcol_double()[39m,
      FDR = [32mcol_double()[39m,
      entries = [31mcol_character()[39m
    )
    
    Parsed with column specification:
    cols(
      GO_acc = [31mcol_character()[39m,
      term_type = [31mcol_character()[39m,
      Term = [31mcol_character()[39m,
      queryitem = [31mcol_character()[39m,
      querytotal = [31mcol_character()[39m,
      bgitem = [31mcol_character()[39m,
      bgtotal = [31mcol_character()[39m,
      pvalue = [31mcol_character()[39m,
      FDR = [31mcol_character()[39m,
      entries = [31mcol_character()[39m
    )
    
    Parsed with column specification:
    cols(
      GO_acc = [31mcol_character()[39m,
      term_type = [31mcol_character()[39m,
      Term = [31mcol_character()[39m,
      queryitem = [32mcol_double()[39m,
      querytotal = [32mcol_double()[39m,
      bgitem = [32mcol_double()[39m,
      bgtotal = [32mcol_double()[39m,
      pvalue = [32mcol_double()[39m,
      FDR = [32mcol_double()[39m,
      entries = [31mcol_character()[39m
    )
    
    Parsed with column specification:
    cols(
      GO_acc = [31mcol_character()[39m,
      term_type = [31mcol_character()[39m,
      Term = [31mcol_character()[39m,
      queryitem = [32mcol_double()[39m,
      querytotal = [32mcol_double()[39m,
      bgitem = [32mcol_double()[39m,
      bgtotal = [32mcol_double()[39m,
      pvalue = [32mcol_double()[39m,
      FDR = [32mcol_double()[39m,
      entries = [31mcol_character()[39m
    )
    
    Parsed with column specification:
    cols(
      GO_acc = [31mcol_character()[39m,
      term_type = [31mcol_character()[39m,
      Term = [31mcol_character()[39m,
      queryitem = [32mcol_double()[39m,
      querytotal = [32mcol_double()[39m,
      bgitem = [32mcol_double()[39m,
      bgtotal = [32mcol_double()[39m,
      pvalue = [32mcol_double()[39m,
      FDR = [32mcol_double()[39m,
      entries = [31mcol_character()[39m
    )
    


    List of 8
     $ K60dvsw: tbl_df [105 × 11] (S3: tbl_df/tbl/data.frame)
      ..$ GO_acc    : chr [1:105] "GO:0050896" "GO:0009611" "GO:0042221" "GO:0006950" ...
      ..$ term_type : chr [1:105] "P" "P" "P" "P" ...
      ..$ Term      : chr [1:105] "response to stimulus" "response to wounding" "response to chemical stimulus" "response to stress" ...
      ..$ queryitem : num [1:105] 30 7 18 18 6 6 9 7 8 9 ...
      ..$ querytotal: num [1:105] 94 94 94 94 94 94 94 94 94 94 ...
      ..$ bgitem    : num [1:105] 4057 197 2085 2320 229 ...
      ..$ bgtotal   : num [1:105] 37767 37767 37767 37767 37767 ...
      ..$ pvalue    : num [1:105] 2.7e-08 8.2e-07 3.7e-06 1.6e-05 2.8e-05 3.6e-05 3.4e-05 1.1e-04 1.4e-04 1.5e-04 ...
      ..$ FDR       : num [1:105] 8.8e-06 1.3e-04 4.0e-04 1.3e-03 1.7e-03 1.7e-03 1.7e-03 4.5e-03 4.9e-03 4.9e-03 ...
      ..$ entries   : chr [1:105] "// AT3G17790 // AT2G17840 // AT5G65140 // AT2G33380 // AT2G22330 // AT4G17615 // AT2G30020 // AT3G44260 // AT3G"| __truncated__ "// AT1G76650 // AT1G01470 // AT2G22330 // AT2G30020 // AT3G44260 // AT5G63770 // AT2G02990" "// AT5G01540 // AT2G17840 // AT5G65140 // AT3G61900 // AT1G01470 // AT2G35930 // AT3G17790 // AT2G33380 // AT1G"| __truncated__ "// AT1G61560 // AT3G17790 // AT2G17840 // AT3G25760 // AT1G76650 // AT1G12610 // AT2G35930 // AT5G57220 // AT2G"| __truncated__ ...
      ..$ dataSet   : chr [1:105] "K60dvsw" "K60dvsw" "K60dvsw" "K60dvsw" ...
     $ K40dvsw: tbl_df [249 × 11] (S3: tbl_df/tbl/data.frame)
      ..$ GO_acc    : chr [1:249] "GO:0050896" "GO:0006950" "GO:0010200" "GO:0042221" ...
      ..$ term_type : chr [1:249] "P" "P" "P" "P" ...
      ..$ Term      : chr [1:249] "response to stimulus" "response to stress" "response to chitin" "response to chemical stimulus" ...
      ..$ queryitem : num [1:249] 79 55 16 49 17 16 16 36 23 13 ...
      ..$ querytotal: num [1:249] 275 275 275 275 275 275 275 275 275 275 ...
      ..$ bgitem    : num [1:249] 4057 2320 151 2085 240 ...
      ..$ bgtotal   : num [1:249] 37767 37767 37767 37767 37767 ...
      ..$ pvalue    : num [1:249] 2.7e-16 1.3e-14 1.1e-13 6.5e-13 7.9e-12 ...
      ..$ FDR       : num [1:249] 2.4e-13 6.0e-12 3.4e-11 1.5e-10 1.4e-09 ...
      ..$ entries   : chr [1:249] "// AT5G47220 // AT2G06050 // AT2G17840 // AT2G30020 // AT1G25560 // AT1G32640 // AT2G34930 // AT2G47450 // AT5G"| __truncated__ "// AT5G47220 // AT2G06050 // AT2G17840 // AT5G64750 // AT1G32640 // AT2G34930 // AT5G54770 // AT1G29340 // AT1G"| __truncated__ "// AT5G47220 // AT1G66160 // AT3G52450 // AT4G17230 // AT1G32640 // AT5G67300 // AT2G35930 // AT4G31800 // AT4G"| __truncated__ "// AT5G47220 // AT2G17840 // AT5G64750 // AT1G32640 // AT4G34410 // AT2G14900 // AT2G23290 // AT1G61890 // AT1G"| __truncated__ ...
      ..$ dataSet   : chr [1:249] "K40dvsw" "K40dvsw" "K40dvsw" "K40dvsw" ...
     $ K20dvsw: tbl_df [541 × 11] (S3: tbl_df/tbl/data.frame)
      ..$ GO_acc    : chr [1:541] "GO:0050896" "GO:0006950" "GO:0042221" "GO:0009628" ...
      ..$ term_type : chr [1:541] "P" "P" "P" "P" ...
      ..$ Term      : chr [1:541] "response to stimulus" "response to stress" "response to chemical stimulus" "response to abiotic stimulus" ...
      ..$ queryitem : num [1:541] 232 144 131 97 43 109 111 111 109 109 ...
      ..$ querytotal: num [1:541] 1357 1357 1357 1357 1357 ...
      ..$ bgitem    : num [1:541] 4057 2320 2085 1471 485 ...
      ..$ bgtotal   : num [1:541] 37767 37767 37767 37767 37767 ...
      ..$ pvalue    : num [1:541] 3.4e-12 5.8e-10 1.8e-09 3.2e-08 2.7e-07 ...
      ..$ FDR       : num [1:541] 7.0e-09 5.9e-07 1.3e-06 1.6e-05 9.4e-05 9.4e-05 1.7e-04 1.7e-04 1.7e-04 1.7e-04 ...
      ..$ entries   : chr [1:541] "// AT1G80840 // AT1G33080 // AT4G22950 // AT2G17860 // AT1G78290 // AT1G73500 // AT4G35900 // AT5G65380 // AT5G"| __truncated__ "// AT1G80840 // AT4G22950 // AT1G78290 // AT1G73500 // AT1G73330 // AT5G42810 // AT5G63790 // AT4G37150 // AT1G"| __truncated__ "// AT1G80840 // AT1G33080 // AT3G21150 // AT1G73500 // AT1G73330 // AT5G65380 // AT5G58080 // AT5G63790 // AT1G"| __truncated__ "// AT4G22950 // AT1G78290 // AT1G73500 // AT4G35900 // AT5G27100 // AT4G25480 // AT4G28556 // AT1G20450 // AT2G"| __truncated__ ...
      ..$ dataSet   : chr [1:541] "K20dvsw" "K20dvsw" "K20dvsw" "K20dvsw" ...
     $ Krwvsw : tbl_df [795 × 11] (S3: tbl_df/tbl/data.frame)
      ..$ GO_acc    : chr [1:795] "GO:0006412" "GO:0042254" "GO:0009059" "GO:0034645" ...
      ..$ term_type : chr [1:795] "P" "P" "P" "P" ...
      ..$ Term      : chr [1:795] "translation" "ribosome biogenesis" "macromolecule biosynthetic process" "cellular macromolecule biosynthetic process" ...
      ..$ queryitem : num [1:795] 254 87 402 400 490 87 498 116 564 393 ...
      ..$ querytotal: num [1:795] 2076 2076 2076 2076 2076 ...
      ..$ bgitem    : num [1:795] 1445 241 3685 3661 4925 ...
      ..$ bgtotal   : num [1:795] 37767 37767 37767 37767 37767 ...
      ..$ pvalue    : num [1:795] 1.3e-53 2.5e-37 1.7e-37 1.9e-37 6.9e-37 ...
      ..$ FDR       : num [1:795] 4.2e-50 1.9e-34 1.9e-34 1.9e-34 4.4e-34 ...
      ..$ entries   : chr [1:795] "// AT1G61580 // AT4G36130 // AT3G43980 // AT1G04480 // AT3G09500 // AT5G60670 // AT1G09690 // AT2G37600 // AT3G"| __truncated__ "// AT3G09500 // AT5G60670 // AT3G55750 // AT3G44890 // AT3G21540 // AT5G02610 // AT4G22380 // AT4G27090 // AT3G"| __truncated__ "// AT1G61580 // AT3G43980 // AT1G60850 // AT1G09690 // AT1G17560 // AT3G55750 // AT3G44890 // AT1G56045 // AT5G"| __truncated__ "// AT1G61580 // AT3G43980 // AT1G60850 // AT1G09690 // AT1G17560 // AT3G55750 // AT3G44890 // AT1G56045 // AT5G"| __truncated__ ...
      ..$ dataSet   : chr [1:795] "Krwvsw" "Krwvsw" "Krwvsw" "Krwvsw" ...
     $ L60dvsw: tbl_df [0 × 11] (S3: tbl_df/tbl/data.frame)
      ..$ GO_acc    : chr(0) 
      ..$ term_type : chr(0) 
      ..$ Term      : chr(0) 
      ..$ queryitem : chr(0) 
      ..$ querytotal: chr(0) 
      ..$ bgitem    : chr(0) 
      ..$ bgtotal   : chr(0) 
      ..$ pvalue    : chr(0) 
      ..$ FDR       : chr(0) 
      ..$ entries   : chr(0) 
      ..$ dataSet   : chr(0) 
     $ L40dvsw: tbl_df [161 × 11] (S3: tbl_df/tbl/data.frame)
      ..$ GO_acc    : chr [1:161] "GO:0050896" "GO:0022621" "GO:0080090" "GO:0019222" ...
      ..$ term_type : chr [1:161] "P" "P" "P" "P" ...
      ..$ Term      : chr [1:161] "response to stimulus" "shoot system development" "regulation of primary metabolic process" "regulation of metabolic process" ...
      ..$ queryitem : num [1:161] 50 6 21 23 10 6 7 10 13 21 ...
      ..$ querytotal: num [1:161] 299 299 299 299 299 299 299 299 299 299 ...
      ..$ bgitem    : num [1:161] 4057 358 1952 2210 1228 ...
      ..$ bgtotal   : num [1:161] 37767 37767 37767 37767 37767 ...
      ..$ pvalue    : num [1:161] 0.0012 0.069 0.099 0.11 0.51 0.66 0.83 0.13 0.092 0.14 ...
      ..$ FDR       : num [1:161] 0.52 1 1 1 1 1 1 1 1 1 ...
      ..$ entries   : chr [1:161] "// AT4G12560 // AT4G04955 // AT5G49730 // AT1G01520 // AT3G61220 // AT5G47100 // AT1G76570 // AT3G22420 // AT5G"| __truncated__ "// AT4G01060 // AT1G64690 // AT3G61970 // AT2G23430 // AT3G25717 // AT3G50660" "// AT5G16600 // AT5G24120 // AT5G43700 // AT2G28200 // AT2G31230 // AT3G16220 // AT4G22950 // AT1G64860 // AT4G"| __truncated__ "// AT5G16600 // AT1G01520 // AT4G22950 // AT5G47140 // AT5G24120 // AT2G28200 // AT4G34060 // AT4G35270 // AT3G"| __truncated__ ...
      ..$ dataSet   : chr [1:161] "L40dvsw" "L40dvsw" "L40dvsw" "L40dvsw" ...
     $ L20dvsw: tbl_df [609 × 11] (S3: tbl_df/tbl/data.frame)
      ..$ GO_acc    : chr [1:609] "GO:0006915" "GO:0050896" "GO:0006950" "GO:0012501" ...
      ..$ term_type : chr [1:609] "P" "P" "P" "P" ...
      ..$ Term      : chr [1:609] "apoptosis" "response to stimulus" "response to stress" "programmed cell death" ...
      ..$ queryitem : num [1:609] 36 303 196 39 39 39 72 35 40 107 ...
      ..$ querytotal: num [1:609] 1904 1904 1904 1904 1904 ...
      ..$ bgitem    : num [1:609] 159 4057 2320 244 286 ...
      ..$ bgtotal   : num [1:609] 37767 37767 37767 37767 37767 ...
      ..$ pvalue    : num [1:609] 3.3e-12 1.5e-11 1.1e-11 3.6e-09 1.6e-07 ...
      ..$ FDR       : num [1:609] 7.7e-09 1.2e-08 1.2e-08 2.1e-06 6.2e-05 6.2e-05 5.2e-04 1.1e-03 7.9e-03 5.0e-02 ...
      ..$ entries   : chr [1:609] "// AT4G11280 // AT5G58120 // AT1G69550 // AT5G46490 // AT2G14080 // AT5G40060 // AT1G61180 // AT1G12060 // AT5G"| __truncated__ "// AT4G12560 // AT3G44480 // AT4G34135 // AT3G47950 // AT1G33080 // AT5G49730 // AT5G45380 // AT3G61220 // AT5G"| __truncated__ "// AT4G12560 // AT3G44480 // AT3G47950 // AT5G45380 // AT3G61220 // AT5G58120 // AT1G78290 // AT1G61180 // AT5G"| __truncated__ "// AT4G11280 // AT5G58120 // AT1G69550 // AT5G46490 // AT2G14080 // AT5G40060 // AT2G19860 // AT1G61180 // AT1G"| __truncated__ ...
      ..$ dataSet   : chr [1:609] "L20dvsw" "L20dvsw" "L20dvsw" "L20dvsw" ...
     $ Lrwvsw : tbl_df [456 × 11] (S3: tbl_df/tbl/data.frame)
      ..$ GO_acc    : chr [1:456] "GO:0050896" "GO:0042221" "GO:0008037" "GO:0048544" ...
      ..$ term_type : chr [1:456] "P" "P" "P" "P" ...
      ..$ Term      : chr [1:456] "response to stimulus" "response to chemical stimulus" "cell recognition" "recognition of pollen" ...
      ..$ queryitem : num [1:456] 163 93 7 7 41 43 7 14 41 18 ...
      ..$ querytotal: num [1:456] 1056 1056 1056 1056 1056 ...
      ..$ bgitem    : num [1:456] 4057 2085 32 32 766 ...
      ..$ bgtotal   : num [1:456] 37767 37767 37767 37767 37767 ...
      ..$ pvalue    : num [1:456] 2.6e-06 1.3e-05 7.8e-05 7.8e-05 1.2e-04 2.0e-04 2.3e-04 1.8e-04 1.5e-04 3.9e-04 ...
      ..$ FDR       : num [1:456] 0.0043 0.01 0.032 0.032 0.039 0.042 0.042 0.042 0.042 0.063 ...
      ..$ entries   : chr [1:456] "// AT1G64380 // AT5G26680 // AT4G34810 // AT5G45380 // AT1G66480 // AT1G12210 // AT2G16650 // AT2G30770 // AT4G"| __truncated__ "// AT1G64380 // AT1G52500 // AT2G16650 // AT4G23550 // AT3G14990 // AT1G08930 // AT2G37430 // AT1G07590 // AT4G"| __truncated__ "// AT1G61370 // AT1G11300 // AT1G61490 // AT4G27300 // AT1G61400 // AT1G61420 // AT1G61360" "// AT1G61370 // AT1G11300 // AT1G61490 // AT4G27300 // AT1G61400 // AT1G61420 // AT1G61360" ...
      ..$ dataSet   : chr [1:456] "Lrwvsw" "Lrwvsw" "Lrwvsw" "Lrwvsw" ...



2916


***
# 4. Combining the data frames

Now that all the data is in a list of data frames, we can do a quick `rbind()` call to combine it into a single data frame. After that, we'll create a new column "enrich" at the same time. We can combine the datasets at this point because we have created that additional `dataSet` column.

This will essentially create a "long-format" dataset.


```R
# Now that all the data is in a list, we can combine it into a single data frame
allData.df <- 
    do.call(rbind, data.list) %>% 
    mutate(enrich = queryitem/bgitem)

# What does the dataframe look like?
str(allData.df)
```

    tbl_df [2,916 × 12] (S3: tbl_df/tbl/data.frame)
     $ GO_acc    : chr [1:2916] "GO:0050896" "GO:0009611" "GO:0042221" "GO:0006950" ...
     $ term_type : chr [1:2916] "P" "P" "P" "P" ...
     $ Term      : chr [1:2916] "response to stimulus" "response to wounding" "response to chemical stimulus" "response to stress" ...
     $ queryitem : num [1:2916] 30 7 18 18 6 6 9 7 8 9 ...
     $ querytotal: num [1:2916] 94 94 94 94 94 94 94 94 94 94 ...
     $ bgitem    : num [1:2916] 4057 197 2085 2320 229 ...
     $ bgtotal   : num [1:2916] 37767 37767 37767 37767 37767 ...
     $ pvalue    : num [1:2916] 2.7e-08 8.2e-07 3.7e-06 1.6e-05 2.8e-05 3.6e-05 3.4e-05 1.1e-04 1.4e-04 1.5e-04 ...
     $ FDR       : num [1:2916] 8.8e-06 1.3e-04 4.0e-04 1.3e-03 1.7e-03 1.7e-03 1.7e-03 4.5e-03 4.9e-03 4.9e-03 ...
     $ entries   : chr [1:2916] "// AT3G17790 // AT2G17840 // AT5G65140 // AT2G33380 // AT2G22330 // AT4G17615 // AT2G30020 // AT3G44260 // AT3G"| __truncated__ "// AT1G76650 // AT1G01470 // AT2G22330 // AT2G30020 // AT3G44260 // AT5G63770 // AT2G02990" "// AT5G01540 // AT2G17840 // AT5G65140 // AT3G61900 // AT1G01470 // AT2G35930 // AT3G17790 // AT2G33380 // AT1G"| __truncated__ "// AT1G61560 // AT3G17790 // AT2G17840 // AT3G25760 // AT1G76650 // AT1G12610 // AT2G35930 // AT5G57220 // AT2G"| __truncated__ ...
     $ dataSet   : chr [1:2916] "K60dvsw" "K60dvsw" "K60dvsw" "K60dvsw" ...
     $ enrich    : num [1:2916] 0.00739 0.03553 0.00863 0.00776 0.0262 ...


***
# 5a. Pulling out the significant terms (optional)

Using the combined dataset, we can now filter and look for the significant terms. This yields about 118 unique terms but we'll actually use a premade set of terms in later steps. This is more just... to show we can do it.


```R
# Find the significant GO terms you're interested in
allSigTerms <-
allData.df %>% 
    # Filter the data
    filter(FDR <= 0.01, term_type == "P") %>% 
    # What are the terms left after filtering?
    pull(Term)

allSigTerms
# How many unique terms are there?
str(unique(allSigTerms))
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'response to stimulus'</li><li>'response to wounding'</li><li>'response to chemical stimulus'</li><li>'response to stress'</li><li>'response to water deprivation'</li><li>'response to water'</li><li>'response to biotic stimulus'</li><li>'response to external stimulus'</li><li>'response to other organism'</li><li>'multi-organism process'</li><li>'response to abiotic stimulus'</li><li>'response to stimulus'</li><li>'response to stress'</li><li>'response to chitin'</li><li>'response to chemical stimulus'</li><li>'response to carbohydrate stimulus'</li><li>'response to water deprivation'</li><li>'response to water'</li><li>'response to abiotic stimulus'</li><li>'response to biotic stimulus'</li><li>'response to wounding'</li><li>'response to organic substance'</li><li>'response to other organism'</li><li>'response to cold'</li><li>'multi-organism process'</li><li>'response to osmotic stress'</li><li>'response to fungus'</li><li>'response to temperature stimulus'</li><li>'cold acclimation'</li><li>'response to external stimulus'</li><li>'defense response'</li><li>'response to abscisic acid stimulus'</li><li>'di-, tri-valent inorganic cation transport'</li><li>'immune system process'</li><li>'transport'</li><li>'establishment of localization'</li><li>'immune response'</li><li>'response to endogenous stimulus'</li><li>'regulation of defense response'</li><li>'defense response to fungus'</li><li>'localization'</li><li>'response to jasmonic acid stimulus'</li><li>'multidrug transport'</li><li>'innate immune response'</li><li>'drug transport'</li><li>'response to drug'</li><li>'response to stimulus'</li><li>'response to stress'</li><li>'response to chemical stimulus'</li><li>'response to abiotic stimulus'</li><li>'response to temperature stimulus'</li><li>'regulation of transcription'</li><li>'regulation of biosynthetic process'</li><li>'regulation of cellular biosynthetic process'</li><li>'regulation of macromolecule biosynthetic process'</li><li>'regulation of nucleobase, nucleoside, nucleotide and nucleic acid metabolic process'</li><li>'starch metabolic process'</li><li>'regulation of cellular metabolic process'</li><li>'regulation of nitrogen compound metabolic process'</li><li>'response to water'</li><li>'regulation of primary metabolic process'</li><li>'carbohydrate biosynthetic process'</li><li>'transcription'</li><li>'response to water deprivation'</li><li>'cellular carbohydrate biosynthetic process'</li><li>'response to chitin'</li><li>'response to abscisic acid stimulus'</li><li>'response to hydrogen peroxide'</li><li>'regulation of gene expression'</li><li>'response to heat'</li><li>'regulation of transcription, DNA-dependent'</li><li>'flavonoid biosynthetic process'</li><li>'response to wounding'</li><li>'regulation of RNA metabolic process'</li><li>'cold acclimation'</li><li>'regulation of metabolic process'</li><li>'regulation of cellular process'</li><li>'cellular glucan metabolic process'</li><li>'RNA biosynthetic process'</li><li>'regulation of macromolecule metabolic process'</li><li>'response to organic substance'</li><li>'glucan metabolic process'</li><li>'transcription, DNA-dependent'</li><li>'cellular carbohydrate metabolic process'</li><li>'response to carbohydrate stimulus'</li><li>'response to cold'</li><li>'flavonoid metabolic process'</li><li>'response to oxidative stress'</li><li>'polysaccharide metabolic process'</li><li>'response to high light intensity'</li><li>'cellular polysaccharide metabolic process'</li><li>'response to reactive oxygen species'</li><li>'carbohydrate metabolic process'</li><li>'transmembrane transport'</li><li>'biological regulation'</li><li>'response to light intensity'</li><li>'translation'</li><li>'ribosome biogenesis'</li><li>'macromolecule biosynthetic process'</li><li>'cellular macromolecule biosynthetic process'</li><li>'cellular biosynthetic process'</li><li>'ribonucleoprotein complex biogenesis'</li><li>'biosynthetic process'</li><li>'cellular component biogenesis'</li><li>'cellular macromolecule metabolic process'</li><li>'gene expression'</li><li>'cellular protein metabolic process'</li><li>'cellular process'</li><li>'macromolecule metabolic process'</li><li>'cellular metabolic process'</li><li>'protein metabolic process'</li><li>'primary metabolic process'</li><li>'metabolic process'</li><li>'DNA replication'</li><li>'DNA metabolic process'</li><li>'response to stimulus'</li><li>'cell cycle'</li><li>'regulation of cell cycle'</li><li>'cellular nitrogen compound metabolic process'</li><li>'nitrogen compound metabolic process'</li><li>'response to auxin stimulus'</li><li>'DNA-dependent DNA replication'</li><li>'microtubule-based movement'</li><li>'response to stress'</li><li>'microtubule-based process'</li><li>'DNA replication initiation'</li><li>'response to hormone stimulus'</li><li>'post-embryonic development'</li><li>'cell cycle process'</li><li>'response to endogenous stimulus'</li><li>'nucleobase, nucleoside, nucleotide and nucleic acid metabolic process'</li><li>'response to chemical stimulus'</li><li>'cellular nitrogen compound biosynthetic process'</li><li>'lysine biosynthetic process via diaminopimelate'</li><li>'lysine biosynthetic process'</li><li>'DNA duplex unwinding'</li><li>'diaminopimelate metabolic process'</li><li>'DNA unwinding during replication'</li><li>'DNA geometric change'</li><li>'lysine metabolic process'</li><li>'response to other organism'</li><li>'cellular amino acid biosynthetic process'</li><li>'lipid localization'</li><li>'response to biotic stimulus'</li><li>'response to organic substance'</li><li>'apoptosis'</li><li>'response to stimulus'</li><li>'response to stress'</li><li>'programmed cell death'</li><li>'cell death'</li><li>'death'</li><li>'defense response'</li><li>'response to inorganic substance'</li><li>'response to abscisic acid stimulus'</li><li>'response to stimulus'</li><li>'response to chemical stimulus'</li></ol>



     chr [1:118] "response to stimulus" "response to wounding" ...


***
# 5b. Importing a list of GO terms

We'll use `read_file()` to import a text file of data terms. They are all separated by new line characters but it's not really a table. So we'll just import it as a single string and then separate based on the newline characters.

This will create `goTermKeep` which we will use to filter our data later.


```R
# Import the specific term list

goTermKeep <- read_file("keep_formakinggofigures.txt") %>% str_split(pattern = "\r\n") %>% unlist()

# What does it look like
goTermKeep

# What are the attributes
str(goTermKeep)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'multidrug transport'</li><li>'drug transport'</li><li>'transmembrane transport'</li><li>'ion transport'</li><li>'ion homeostasis'</li><li>'cell wall macromolecule metabolic process'</li><li>'cell wall modification'</li><li>'wax biosynthetic process'</li><li>'fatty acid biosynthetic process'</li><li>'lipid localization'</li><li>'lipid transport'</li><li>'glycine metabolic process'</li><li>'serine family amino acid catabolic process'</li><li>'lysine biosynthetic process'</li><li>'aspartate family amino acid biosynthetic process'</li><li>'glutamine family amino acid biosynthetic process'</li><li>'glucosinolate biosynthetic process'</li><li>'glycosinolate biosynthetic process'</li><li>'fat-soluble vitamin biosynthetic process'</li><li>'flavonoid biosynthetic process'</li><li>'flavonoid metabolic process'</li><li>'chlorophyll biosynthetic process'</li><li>'photorespiration'</li><li>'photosynthesis, dark reaction'</li><li>'photosynthetic electon transport chain'</li><li>'photosynthesis, light reaction'</li><li>'photosynthesis'</li><li>'starch metabolic process'</li><li>'polysaccharide metabolic process'</li><li>'glucan metabolic process'</li><li>'carbohydrate metabolic process'</li><li>'carbohydrate biosynthetic process'</li><li>'glycolysis'</li><li>'glucose metabolic process'</li><li>'monosaccharide metabolic process'</li><li>'apoptosis'</li><li>'programmed cell death'</li><li>'aging'</li><li>'cellular metabolic compound salvage'</li><li>'reproductive structure development'</li><li>'post-embryonic development'</li><li>'cell growth'</li><li>'cell division'</li><li>'regulation of cell cycle'</li><li>'translation'</li><li>'ribosome biogenesis'</li><li>'gene expression'</li><li>'nucleosome organization'</li><li>'regulation of gene expression'</li><li>'response to hydrogen peroxide'</li><li>'response to reactive oxygen species'</li><li>'systemic acquired resistance'</li><li>'immune response'</li><li>'defense response to fungus'</li><li>'defense response to bacterium'</li><li>'response to fungus'</li><li>'response to bacterium'</li><li>'response to other organism'</li><li>'response to ozone'</li><li>'response to UV'</li><li>'response to blue light'</li><li>'response to red light'</li><li>'response to high light intensity'</li><li>'response to salt stress'</li><li>'response to oxidative stress'</li><li>'response to osmotic stress'</li><li>'response to desiccation'</li><li>'response to water deprivation'</li><li>'response to water'</li><li>'cold acclimation'</li><li>'response to cold'</li><li>'response to heat'</li><li>'response to chitin'</li><li>'response to ethylene stiumlus'</li><li>'response to salicylic acid stimulus'</li><li>'response to auxin stimulus'</li><li>'response to jasmonic acid stimulus'</li><li>'response to abscisic acid stimulus'</li><li>'response to wounding'</li><li>'response to stress'</li></ol>



     chr [1:80] "multidrug transport" "drug transport" ...


***
# 6. Filtering data by GO terms

At this point we are ready to filter our data into a new set called `sigData.df`. However, we'll mutate a couple of variables first to make them into factors. In doing so, we'll artificially set the levels of those factors so that even if the expected values are missing (due to filtering), the memory of their existence as part of the set, will remain. This is mainly for the `dataSet` and `Term` variables.


```R
# What happens when we filter by these terms?
sigData.df <-
allData.df %>% 
    # This mutation step will return L60dvsw back as an x-axis value
    mutate(dataSet = factor(dataSet, levels = fileNames)) %>% 
    filter(Term %in% goTermKeep) %>% 
    mutate(Term = factor(Term, levels = goTermKeep)) %>% 
    # Filter the data again based on the lab notes
    filter(FDR <= 0.01, term_type == "P")

str(sigData.df)
```

    tbl_df [57 × 12] (S3: tbl_df/tbl/data.frame)
     $ GO_acc    : chr [1:57] "GO:0009611" "GO:0006950" "GO:0009414" "GO:0009415" ...
     $ term_type : chr [1:57] "P" "P" "P" "P" ...
     $ Term      : Factor w/ 80 levels "multidrug transport",..: 79 80 68 69 58 80 73 68 69 79 ...
     $ queryitem : num [1:57] 7 18 6 6 8 55 16 16 16 13 ...
     $ querytotal: num [1:57] 94 94 94 94 94 275 275 275 275 275 ...
     $ bgitem    : num [1:57] 197 2320 229 240 599 2320 151 229 240 197 ...
     $ bgtotal   : num [1:57] 37767 37767 37767 37767 37767 ...
     $ pvalue    : num [1:57] 8.2e-07 1.6e-05 2.8e-05 3.6e-05 1.4e-04 ...
     $ FDR       : num [1:57] 0.00013 0.0013 0.0017 0.0017 0.0049 ...
     $ entries   : chr [1:57] "// AT1G76650 // AT1G01470 // AT2G22330 // AT2G30020 // AT3G44260 // AT5G63770 // AT2G02990" "// AT1G61560 // AT3G17790 // AT2G17840 // AT3G25760 // AT1G76650 // AT1G12610 // AT2G35930 // AT5G57220 // AT2G"| __truncated__ "// AT2G17840 // AT1G01470 // AT2G35930 // AT2G33380 // AT4G17615 // AT3G25760" "// AT2G17840 // AT1G01470 // AT2G35930 // AT2G33380 // AT4G17615 // AT3G25760" ...
     $ dataSet   : Factor w/ 8 levels "K60dvsw","K40dvsw",..: 1 1 1 1 1 2 2 2 2 2 ...
     $ enrich    : num [1:57] 0.03553 0.00776 0.0262 0.025 0.01336 ...


***
# 7. Generating our figure

At the top of the code cell we'll create a few variables to help set limits on our scales:

- `max_dotsize`: sets how big our bubbles in the visualization will get
- `enrichMax`: sets the upper limit on our enrich values when setting bubble size scale.
- `min_FDR`: sets the lower limit of our FDR values when setting colour scale.

We'll use a Viridis colour scale but try out the grey colour scale by uncommenting the 5 lines before the `scale_colour_viridis()` line, and comment that line out. Something cool is that we can present __all levels__ of our data on the x and y axes by setting the `drop` parameter to `FALSE` in the `scale_x_discrete()` and `scale_y_discrete()` layers.

You can update the theme information as you see fit (bigger text, etc.)


```R
options(repr.plot.width = 16, repr.plot.height = 10)
# Set a maximum dot size
max_dotsize <-6
enrichMax <- max(sigData.df$enrich)
min_FDR = min(sigData.df$FDR)

# The hardest part of printing this is you want to maintain (as much of) the original goTermKeep values
# On the other hand you want to drop any data where FDR > 0.01 or the term_type != "P"

# If you filter the dataset this way, you'll lose y-axis values, even if converted to a factor with levels!
# To retain those missing values/levels from the factor, set scale_y_discrete(drop = FALSE)

ggplot(sigData.df) +
    # 2. Aesthetics and Theme
    aes(x = dataSet, y = Term, colour = FDR, size = enrich) + #Set the y and x axes

    labs(x = NULL, y = NULL) +
    theme(axis.title.x = element_text(size = rel(0.7))) +
    
    # 3. Scaling
    scale_size_continuous(range = c(1, max_dotsize), 
                          limits = c(0, enrichMax), 
                          breaks = seq(0, enrichMax, 0.1),
                          guide = guide_legend(order=2)
                         ) +
    scale_y_discrete(drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
#    scale_colour_gradient(limits = c(min_FDR, 0.01), 
#                           low="gray36", 
#                           high="gray87",
#                           guide = guide_colorbar(order=1)
#                          ) +
    scale_colour_viridis(guide = guide_colorbar(order=1)) +

    # 4. Geoms
    geom_point(aes(size = enrich, color = FDR))
    
```


![png](output_13_0.png)



```R
# In the above graph the spacing between the samples is a bit wide (depending on your screen)
# When we save this plot as a pdf, try playing around with the width and height parameters
# The ggsave function will automatically adjust column spacing, fonts etc.
# The following width and height settings result in the plot you see in the lab manual
ggsave("bubble_plot.pdf", width = 25,  height = 30,  units = "cm")
```


```R

```
