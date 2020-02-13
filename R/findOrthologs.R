##' map orthologous genes using ENSEMBL BioMART
##' @rdname findOrthologs
##' @name findOrthologs
##' @author Vitalii Kleshchevnikov
##' @description findOrthologs(): map orthologous genes using ENSEMBL BioMART (\link[biomaRt]{getBM}) from one species to another species. This function retrieves ids by 500 at a time to decrease the load on ENSEMBL servers.
##' @param datasets_FROM_TO environment containing "Mart" objects produced by \code{loadBIOMARTdatasets()}. It is convenient to specify datasets here rather than pre-loading. Full list of datasets in ENSEMBL BioMART: listDatasets("ensembl")
##' @param from_filters Filters (one or more) that should be used in the query in the "from" species database. A possible list of filters can be retrieved using \code{attributesFiltersFromTo(datasets_FROM_TO = loadBIOMARTdatasets())}.
##' @param from_values Values of the filter, e.g. vector of affy IDs in the "from" species database. If multiple filters are specified then the argument should be a list of vectors of which the position of each vector corresponds to the position of the filters in the filters argument. Please limit the number of values to less than 500.
##' @param to_attributes Attributes you want to retrieve from the "to" species database. A possible list of attributes can be retrieved using the function \code{attributesFiltersFromTo(datasets_FROM_TO = loadBIOMARTdatasets())}.
##' @param to_homolog_attribute Attribute that returns homolog gene id in the "from" database
##' @param from_gene_id_name name of the column containing the ENSEMBL gene ids of "from" species
##' @param to_gene_id_name name of the column containing the ENSEMBL gene ids of "from" species
##' @return findOrthologs(): data.frame containing the mapping of homologous genes from one species to another species including \code{to_attributes}
##' @import biomaRt
##' @export findOrthologs
##' @examples
##' findOrthologs()
findOrthologs = function(datasets_FROM_TO = loadBIOMARTdatasets(from = "hsapiens_gene_ensembl",
                                                                to = "mmusculus_gene_ensembl"),
                         from_filters = "hgnc_symbol",
                         from_values = c("TP53", "TERT"),
                         to_attributes = "external_gene_name",
                         to_homolog_attribute = "mmusculus_homolog_ensembl_gene",
                         from_gene_id_name = "human_ensembl_gene_id",
                         to_gene_id_name = "mouse_ensembl_gene_id"){

  ensembl_from = datasets_FROM_TO$ensembl_from
  ensembl_to = datasets_FROM_TO$ensembl_to

  # calculate where groups of 500 start and end
  n_from_values = length(from_values)
  starts = seq(1, n_from_values, 500)
  ends = starts + 499
  ends[ends == max(ends)] = n_from_values
  # download iteratively
  to_gene_id2whatever = lapply(1:length(ends), function(n){

    from_values_temp = from_values[starts[n]:ends[n]]
    from_whatever2gene_id = getBM(attributes = unique(c(from_filters, "ensembl_gene_id")),
                                  filters = from_filters,
                                  values = from_values_temp,
                                  mart = ensembl_from)
    from_whatever2gene_id[,from_gene_id_name] = from_whatever2gene_id$ensembl_gene_id
    from_whatever2gene_id$ensembl_gene_id = NULL

    from2to = getBM(attributes = c("ensembl_gene_id", to_homolog_attribute),
                    filters = "ensembl_gene_id",
                    values = from_whatever2gene_id[,from_gene_id_name],
                    mart = ensembl_from)
    from2to[, from_gene_id_name] = from2to$ensembl_gene_id
    from2to$ensembl_gene_id = NULL

    to_gene_id2whatever = getBM(attributes = unique(c(to_attributes, "ensembl_gene_id")),
                                filters = "ensembl_gene_id",
                                values = from2to[, to_homolog_attribute],
                                mart = ensembl_to)
    to_gene_id2whatever[, to_gene_id_name] = to_gene_id2whatever$ensembl_gene_id
    to_gene_id2whatever$ensembl_gene_id = NULL

    to_gene_id2whatever = merge(x = to_gene_id2whatever, y = from2to, by.x = to_gene_id_name, by.y = to_homolog_attribute)
    to_gene_id2whatever = merge(x = to_gene_id2whatever, y = from_whatever2gene_id, by.x = from_gene_id_name, by.y = from_gene_id_name)

    to_gene_id2whatever
  })
  # combine results
  to_gene_id2whatever = Reduce(rbind, to_gene_id2whatever)
  to_gene_id2whatever
}

##' @rdname findOrthologs
##' @name findOrthologsHsMm
##' @description findOrthologsHsMm(): map orthologous genes using ENSEMBL BioMART from human to mouse
##' @return findOrthologsHsMm(): data.frame containing the mapping of homologous genes from human to mouse including \code{to_attributes}
##' @import biomaRt
##' @export findOrthologsHsMm
##' @examples
##' findOrthologsHsMm()
findOrthologsHsMm = function(from_filters = "hgnc_symbol",
                             from_values = c("TP53", "TERT"),
                             to_attributes = "external_gene_name"){
  findOrthologs(datasets_FROM_TO = loadBIOMARTdatasets(from = "hsapiens_gene_ensembl",
                                                       to = "mmusculus_gene_ensembl"),
                from_filters = from_filters,
                from_values = from_values,
                to_attributes = to_attributes,
                to_homolog_attribute = "mmusculus_homolog_ensembl_gene",
                from_gene_id_name = "human_ensembl_gene_id",
                to_gene_id_name = "mouse_ensembl_gene_id")
}

##' @rdname findOrthologs
##' @name findOrthologsMmHs
##' @description findOrthologsMmHs(): map orthologous genes using ENSEMBL BioMART from mouse to human
##' @return findOrthologsMmHs(): data.frame containing the mapping of homologous genes from mouse to human including \code{to_attributes}
##' @import biomaRt
##' @export findOrthologsMmHs
##' @examples
##' findOrthologsMmHs()
findOrthologsMmHs = function(from_filters = "ensembl_gene_id",
                             from_values = c("ENSMUSG00000059552", "ENSMUSG00000021611"),
                             to_attributes = "hgnc_symbol"){
  findOrthologs(datasets_FROM_TO = loadBIOMARTdatasets(from = "mmusculus_gene_ensembl",
                                                       to = "hsapiens_gene_ensembl"),
                from_filters = from_filters,
                from_values = from_values,
                to_attributes = to_attributes,
                to_homolog_attribute = "hsapiens_homolog_ensembl_gene",
                from_gene_id_name = "mouse_ensembl_gene_id",
                to_gene_id_name = "human_ensembl_gene_id")
}

##' @rdname findOrthologs
##' @name loadBIOMARTdatasets
##' @description loadBIOMARTdatasets(): Load 2 ENSEMBL BioMART datasets for finding orthologs. Details \link[biomaRt]{useMart}
##' @param from ENSEMBL biomart dataset for a species whose identifiers you want to convert. Full list of possible options: listDatasets("ensembl")
##' @param to ENSEMBL biomart dataset for a species whose identifiers you want convert into.
##' @return loadBIOMARTdatasets(): environment containing 2 object of class 'Mart', one for each species dataset
##' @import biomaRt
##' @export loadBIOMARTdatasets
loadBIOMARTdatasets = function(from = "hsapiens_gene_ensembl", to = "mmusculus_gene_ensembl") {
  # load biomart ensembl datatsets
  # ensembl = useMart("ensembl")
  # listDatasets(ensembl)
  ensembl_from = useMart("ensembl", from)
  ensembl_to = useMart("ensembl", to)
  datasets_FROM_TO = new.env()
  datasets_FROM_TO$ensembl_from = ensembl_from
  datasets_FROM_TO$ensembl_to = ensembl_to
  datasets_FROM_TO
}

##' @rdname findOrthologs
##' @name attributesFiltersFromTo
##' @description attributesFiltersFromTo(): List filters and attributes available in each dataset. Details: \link[biomaRt]{listAttributes}
##' @param datasets_FROM_TO environment containing "Mart" objects produced by \code{loadBIOMARTdatasets()}
##' @return attributesFiltersFromTo(): list containing 4 vector of either filters or attributes for each dataset
##' @import biomaRt
##' @export attributesFiltersFromTo
##' @examples
##' attributes = attributesFiltersFromTo(datasets_FROM_TO = loadBIOMARTdatasets())
##' grep("Mouse gene stable ID", attributes$from_attributes$description, ignore.case = T, value = T)
##' attributes$from_attributes[which(attributes$from_attributes$description == "Mouse gene name"),]
##' grep("hgnc_symbol", attributes$from_attributes$name, value = T)
##' # find the desired column name
##' attributes$from_attributes[grep("Mouse gene stable ID", attributes$from_attributes$description, ignore.case = T),]
attributesFiltersFromTo = function(datasets_FROM_TO = loadBIOMARTdatasets()){
  ensembl_from = datasets_FROM_TO$ensembl_from
  ensembl_to = datasets_FROM_TO$ensembl_to
  # getting names of the things one can get from ENSEMBL
  list(from_attributes = listAttributes(ensembl_from),
       from_filters = listFilters(ensembl_from),
       to_attributes = listAttributes(ensembl_to),
       to_filters = listFilters(ensembl_to))
}

##' Map moleculer identifiers using ENSEMBL BioMart
##' @rdname mapIDs
##' @name mapIDs
##' @description mapIDs(): Map moleculer identifiers in a data.table using ENSEMBL BioMart - \link[biomaRt]{getBM}
##' @param DT data.table containing id to be converted
##' @param ids2convert_col column in \code{DT} containing identifiers to be converted
##' @param filters Filters (one, in ids2convert_col column) that should be used in the query. A possible list of filters can be retrieved using the function listFilters.
##' @param map_to Attributes you want to retrieve. A possible list of attributes can be retrieved using the function listAttributes.
##' @param map_to_name how to name the column in the output containing \code{map_to} identifiers
##' @param biomart_dataset which BioMart dataset to use for mapping. Full list of possible options: listDatasets("ensembl")
##' @return mapIDs(): data.table \code{DT} containing an additional column (\code{map_to_name})
##' @import biomaRt
##' @import data.table
##' @export mapIDs
##' @examples
##' # load TF regulon data
##' combined_file = "../regulatory_networks_by_cmap/data/validation_TF_regulons/all_simp_regulons_w_pos_effect.tsv"
##' regulons = fread(combined_file, header = T, stringsAsFactors = F)
##' regulons = mapIDs(regulons)
mapIDs = function(DT, ids2convert_col = "target", filters = "hgnc_symbol", map_to = "entrezgene", map_to_name =  "target_entrezgene", biomart_dataset = "hsapiens_gene_ensembl"){
  ## Map HGNC symbols to entrez gene id
  biomart = loadBIOMARTdatasets(from = biomart_dataset,
                                to = biomart_dataset)
  biomart = biomart$ensembl_from

  hgnc_symbol2entrezgene = getBM(attributes = c(map_to, filters),
                                 filters = filters,
                                 values = unique(DT[, ids2convert_col, with = F]),
                                 mart = biomart)

  hgnc_symbol2entrezgene = unique(as.data.table(hgnc_symbol2entrezgene))
  setnames(hgnc_symbol2entrezgene,map_to,map_to_name)
  unique(merge(x = DT, y = hgnc_symbol2entrezgene,
               by.x = ids2convert_col, by.y = filters,
               all.x = T, all.y = F, allow.cartesian = TRUE))

}
