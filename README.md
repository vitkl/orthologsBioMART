## Convenient functions for mapping orthologous genes using ENSEMBL BioMART

The package wraps around biomaRt package in R and pybiomart package in python to provide a concise way of ortholog mapping.

[![DOI](https://zenodo.org/badge/240325005.svg)](https://zenodo.org/badge/latestdoi/240325005)

## Installation R

```r
install.packages("BiocManager")
BiocManager::install(c("vitkl/orthologsBioMART"), dependencies=T)
```

## Usage R

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

## Installation python

```python
pip install git+https://github.com/vitkl/orthologsBioMART.git
```

## Usage python

Mapping between any species: 

```python
# find correct dataset name for your species
from pybiomart import Server
server = Server(host='http://www.ensembl.org')
server.marts['ENSEMBL_MART_ENSEMBL'].list_datasets()

from pyorthomap import FindOrthologs 
# then create the find orthogues object using correct datasets and attributes
# use help(FindOrthologs)
hs2mm = FindOrthologs(
          host = 'http://www.ensembl.org',
          mart = 'ENSEMBL_MART_ENSEMBL',
          from_dataset = 'hsapiens_gene_ensembl',
          to_dataset = 'mmusculus_gene_ensembl',
          from_filters = 'hgnc_symbol',
          from_values = ['TP53', 'TERT'],
          to_attributes = 'external_gene_name',
          to_homolog_attribute = 'mmusculus_homolog_ensembl_gene',
          from_gene_id_name = 'human_ensembl_gene_id',
          to_gene_id_name = 'mouse_ensembl_gene_id'
    )
    
hs2mm.map()
```

Mapping mouse to human and vice-versa:

```python
findOrthologsMmHs(from_filters = 'link_ensembl_gene_id',
                  from_values = ['ENSMUSG00000059552', 'ENSMUSG00000021611']).map()
findOrthologsMmHs(from_filters = 'external_gene_name',
                  from_values = ['Trp53', 'Tert']).map()

findOrthologsHsMm(from_filters = 'link_ensembl_gene_id',
                  from_values = ['ENSG00000141510', 'ENSG00000164362']).map()
findOrthologsHsMm(from_filters = 'hgnc_symbol',
                  from_values = ['TP53', 'TERT']).map()
```
