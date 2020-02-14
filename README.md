## Convenient functions for mapping orthologous genes using ENSEMBL BioMART

The package wraps around biomaRt package to provide a fast way of ortholog mapping.

[![DOI](https://zenodo.org/badge/240325005.svg)](https://zenodo.org/badge/latestdoi/240325005)

## Installation

```r
install.packages("BiocManager")
BiocManager::install(c("vitkl/orthologsBioMART"), dependencies=T)
```

## Usage

Mapping between any species:   
```r
findOrthologs = function(datasets_FROM_TO = loadBIOMARTdatasets(from = "hsapiens_gene_ensembl",
                                                                to = "mmusculus_gene_ensembl"),
                         from_filters = "hgnc_symbol", # type of input identifier
                         from_values = c("TP53", "TERT"), # gene identifiers to map
                         to_attributes = "external_gene_name", # type of output identifier
                         # column in the input database that maps to ENSEMBL id of the output database:
                         to_homolog_attribute = "mmusculus_homolog_ensembl_gene",  
                         # column with ENSEMBL gene id in the input database:  
                         from_gene_id_name = "human_ensembl_gene_id", 
                         # column with ENSEMBL gene id in the output database:  
                         to_gene_id_name = "mouse_ensembl_gene_id")
```

Mapping human to mouse and mouse to human:   
```r
findOrthologsHsMm(from_filters = "hgnc_symbol",
  from_values = c("TP53","TERT"), 
  to_attributes = "external_gene_name")

findOrthologsMmHs(from_filters = "ensembl_gene_id",
  from_values = c("ENSMUSG00000059552", "ENSMUSG00000021611"),
  to_attributes = "hgnc_symbol")
```

To find the correct `to_homolog_attribute` between species use:   
```r
attr = attributesFiltersFromTo(datasets_FROM_TO = loadBIOMARTdatasets())

# examine attr
attr

# find the desired column name
attr$from_attributes[grep("Mouse gene stable ID", attr, ignore.case = T),]
```

