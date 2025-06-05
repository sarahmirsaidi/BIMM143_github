# Class 6: R functions
Sarah Mirsaidi Madjdabadi, A16890196

- [1. Function Basics](#1-function-basics)
- [2. Generate DNA Function](#2-generate-dna-function)
- [3. Generate Protein Function](#3-generate-protein-function)

## 1. Function Basics

Let’s start writing our first silly function to add some numbers:

Every r function has 3 things: - name (we get to pick this) - input
arguments (there can be loads of these separated by a comma) - the body
(the R code that does the work)

``` r
add <- function(x, y=10, z=0){
  x + y + z
}
```

I can just use this function like any other function as long as R knows
about it (i.e. run the code chunk)

``` r
add(1, 100)
```

    [1] 101

``` r
add( x=c(1,2,3,4), y=100)
```

    [1] 101 102 103 104

``` r
add(1)
```

    [1] 11

Functions can have “required” input arguments and “optional” input
arguments. The optional arguments are defined with an equals default
value (`y=10`) in the function definition.

``` r
add(x=1, y=100, z=10)
```

    [1] 111

> Q. Write a function to return a DNA sequence of a user specified
> length. Call it `generate_DNA()`

The `sample()` function can help here

``` r
#generate_dna <- function(size=5) {}

students <- c("jeff", "jeremy", "peter")

sample(students, size = 5, replace = TRUE)
```

    [1] "peter"  "peter"  "peter"  "peter"  "jeremy"

## 2. Generate DNA Function

Now work with `bases` rather than `students`

``` r
bases <- c("A", "C", "G", "T")

sample(bases, size = 10, replace = TRUE)
```

     [1] "A" "G" "A" "C" "A" "C" "A" "T" "C" "A"

Now I have a working `snippet` of code I can use this as the body of my
first function version here:

``` r
generate_dna <- function(size=5) {bases <- c("A", "C", "G", "T")

sample(bases, size = size, replace = TRUE)
  
}
```

``` r
generate_dna(100)
```

      [1] "C" "G" "A" "G" "G" "C" "C" "G" "C" "G" "G" "A" "A" "T" "T" "C" "T" "A"
     [19] "C" "A" "C" "C" "G" "C" "G" "A" "A" "C" "A" "A" "G" "C" "T" "T" "T" "G"
     [37] "A" "T" "A" "C" "A" "T" "A" "C" "T" "A" "C" "A" "A" "C" "A" "A" "T" "A"
     [55] "T" "C" "A" "G" "G" "G" "T" "T" "C" "T" "A" "G" "T" "G" "G" "G" "G" "A"
     [73] "C" "A" "A" "C" "C" "C" "G" "T" "T" "A" "C" "G" "A" "T" "C" "G" "C" "T"
     [91] "T" "G" "G" "T" "A" "G" "T" "G" "A" "T"

``` r
generate_dna()
```

    [1] "T" "A" "T" "T" "T"

I want the ability to return a sequence like “AGTACCTG” i.e. a one
element vector where the bases are all together.

``` r
generate_dna <- function(size=5, together=TRUE) {bases <- c("A", "C", "G", "T")

sequence <- sample(bases, size = size, replace = TRUE)
if(together) {
  sequence <- paste(sequence, collapse = "")
}
return(sequence)
}
```

``` r
generate_dna()
```

    [1] "AATCC"

``` r
generate_dna(together = F)
```

    [1] "A" "C" "G" "C" "A"

## 3. Generate Protein Function

> Q1. Write a protein sequence generating function that will return
> sequences of a user specified length?

> Q2.Generate random protein sequences of length 6 to 12 amino acids.

> Q3. Determine if these sequences can be found in nature or are they
> unique? Why or why not?

We can get the set of 20 natural amino acids from the **bio3d** package.

``` r
aa <- bio3d::aa.table$aa1[1:20]

aa
```

     [1] "A" "R" "N" "D" "C" "Q" "E" "G" "H" "I" "L" "K" "M" "F" "P" "S" "T" "W" "Y"
    [20] "V"

> Q1

``` r
generate_protein <- function(size=6, together=TRUE){aa
  sequence <- sample(aa, size=size, replace=TRUE)
  if(together){
   sequence <- paste(sequence, collapse = "") 
  }
  return(sequence)
}
```

``` r
generate_protein()
```

    [1] "KGMQYC"

> Q2

We can fix this inability to generate multiple sequences by either
editing and adding to the function body code (e.g. a for loop) or by
using the R **apply** family utility function.

``` r
sapply(6:12, generate_protein)
```

    [1] "NGWQLA"       "VDVHFEK"      "PYTAWAAR"     "PGCRGEGPH"    "GEVDSMNSST"  
    [6] "GMLFDAEEFRR"  "QLAAGQHTKCPK"

It would be cool and useful if I could get FASTA format output.

``` r
ans <- sapply(6:12, generate_protein)
ans
```

    [1] "DGAMNH"       "ASKSFLP"      "AWSIFHEY"     "NDFEYGVKH"    "FSQMVGDYKF"  
    [6] "HQCVYESMPPY"  "IAYPDQDNYECC"

``` r
cat(ans, sep="\n")
```

    DGAMNH
    ASKSFLP
    AWSIFHEY
    NDFEYGVKH
    FSQMVGDYKF
    HQCVYESMPPY
    IAYPDQDNYECC

I want this to look like FASTA format with an ID line

    >ID.6
    TPEKMV
    >ID.7
    KLCDMKR
    >ID.8
    RDYMASDR

The functions ‘paste()’ and ‘cat()’ can help us here…

``` r
cat(paste(">ID.", 7:12, "\n", ans, sep= ""), sep="\n")
```

    >ID.7
    DGAMNH
    >ID.8
    ASKSFLP
    >ID.9
    AWSIFHEY
    >ID.10
    NDFEYGVKH
    >ID.11
    FSQMVGDYKF
    >ID.12
    HQCVYESMPPY
    >ID.7
    IAYPDQDNYECC

``` r
id.line <- paste(">ID.", 6:12, sep="")
id.line
```

    [1] ">ID.6"  ">ID.7"  ">ID.8"  ">ID.9"  ">ID.10" ">ID.11" ">ID.12"

``` r
seq.line <- paste(id.line, ans, sep="\n")
cat(seq.line, sep="\n")
```

    >ID.6
    DGAMNH
    >ID.7
    ASKSFLP
    >ID.8
    AWSIFHEY
    >ID.9
    NDFEYGVKH
    >ID.10
    FSQMVGDYKF
    >ID.11
    HQCVYESMPPY
    >ID.12
    IAYPDQDNYECC

> Q3

I BLASTp searched my FASTA format sequences against NR and found that
length 6, 7, 8 are not unique and can be found in the databases with
100% coverage and 100% identity.

Random sequences of length 9 and above are unique and can’t be found in
the databases.

BLASTp has a window of 9.
