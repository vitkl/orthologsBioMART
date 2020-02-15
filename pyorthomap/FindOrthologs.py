from pybiomart import Server
from tqdm.auto import tqdm
import pandas as pd

# Create a function called "chunks" with two arguments, l and n:
def chunks(l, n):
    r"""
    Chunk a long list into a list of smaller lists of length n
    :param l: list
    :param n: integer, the length of smaller lists
    """
    # For item i in a range that is a length of l,
    for i in range(0, len(l), n):
        # Create an index range for l of n items:
        yield l[i:i+n]

class FindOrthologs():
    r"""Map orthologous genes using ENSEMBL BioMart and pybiomart package.
    :param host: BioMart host server link
    :param mart: BioMart database version, e.g. 'ENSEMBL_MART_ENSEMBL'
    :param from_dataset: source BioMart dataset
    :param to_dataset: target BioMart dataset. For full list of datasets check:
            from pybiomart import Server
            server = Server(host='http://www.ensembl.org')
            server.marts['ENSEMBL_MART_ENSEMBL'].list_datasets()
    :param from_filters: Filters that should be used in the query in the source/"from" species database. To see the full list:
            from pybiomart import Server
            server = Server(host='http://www.ensembl.org')
            server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl'].list_filters()
            Some attributes can also be used as a filter (in which case from_attributes is ignored):
            server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl'].list_attributes()
    :param from_values: Values of the filter, e.g. list of gene names or ENSEMBL IDs
    :param from_attributes: Which additional attributes to collect from the source/"from" species database.
    :param to_attributes: Attributes you want to retrieve from the "to" species database. To see the full list:
            from pybiomart import Server
            server = Server(host='http://www.ensembl.org')
            server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl'].list_attributes()
    :param to_homolog_attribute: Attribute that corresponds to the homolog gene id in the "from" database
    :param from_gene_id_name: name of the column to assign to the ENSEMBL gene ids of "from" species
    :param to_gene_id_name: name of the column to assign the ENSEMBL gene ids of "to" species
    :param chunk_size: query the ENSEMBL BioMart webservice in chunks of this size
    """

    def __init__(
          self,
          host = 'http://www.ensembl.org',
          mart = 'ENSEMBL_MART_ENSEMBL',
          from_dataset = 'hsapiens_gene_ensembl',
          to_dataset = 'mmusculus_gene_ensembl',
          from_filters = 'hgnc_symbol',
          from_attributes = ['hgnc_symbol'],
          from_values = ['TP53', 'TERT'],
          to_attributes = ['external_gene_name'],
          to_homolog_attribute = 'mmusculus_homolog_ensembl_gene',
          from_gene_id_name = 'human_ensembl_gene_id',
          to_gene_id_name = 'mouse_ensembl_gene_id',
          chunk_size: int = 300
    ):

        # connect to server
        self.server = Server(host=host)
        self.ensembl_from = (self.server.marts[mart].datasets[from_dataset])
        self.ensembl_to = (self.server.marts[mart].datasets[to_dataset])

        # save parameters
        self.from_filters = from_filters
        self.from_values = from_values
        self.to_attributes = to_attributes
        self.to_homolog_attribute = to_homolog_attribute
        self.from_gene_id_name = from_gene_id_name
        self.to_gene_id_name = to_gene_id_name
        self.from_attributes = from_attributes
        self.chunk_size = chunk_size

    def map(self, from_values=None):
        r"""
        Map identifiers from from_values list. If not provided (None) uses the list provided when setting up the object.
        """

        if from_values is None:
            from_values = self.from_values

        self.to_gene_id2whatever = pd.DataFrame()

        for from_values_1 in tqdm(list(chunks(from_values, self.chunk_size))):

            # Map to ENSEMBL ID within source species
            if self.from_filters == 'link_ensembl_gene_id':
                from_whatever2gene_id = self.ensembl_from.query(\
                                            attributes=self.from_attributes + ['ensembl_gene_id'],
                                            filters={self.from_filters: from_values_1},
                                            use_attr_names=True)
            else:
                from_whatever2gene_id = self.ensembl_from.query(\
                                            attributes=[self.from_filters] + ['ensembl_gene_id'],
                                            use_attr_names=True)
                sel_ind = from_whatever2gene_id[self.from_filters].isin(from_values_1)
                from_whatever2gene_id = from_whatever2gene_id.loc[sel_ind,]

            from_whatever2gene_id = from_whatever2gene_id.rename(columns={'ensembl_gene_id':
                                                                          self.from_gene_id_name})

            # Map between species
            from2to = self.ensembl_from.query(attributes=['ensembl_gene_id', self.to_homolog_attribute],
                                              filters={'link_ensembl_gene_id': from_whatever2gene_id[self.from_gene_id_name].tolist()},
                                              use_attr_names=True)
            from2to = from2to.rename(columns={'ensembl_gene_id': self.from_gene_id_name})

            # Map to other ID within the target species
            to_gene_id2whatever = self.ensembl_to.query(attributes=self.to_attributes + ['ensembl_gene_id'],
                                                        filters={'link_ensembl_gene_id': from2to[self.to_homolog_attribute].tolist()},
                                                        use_attr_names=True)
            to_gene_id2whatever = to_gene_id2whatever.rename(columns={'ensembl_gene_id': self.to_gene_id_name})

            # Put the mapping together via merge
            to_gene_id2whatever = to_gene_id2whatever.merge(from2to, how='outer',
                                               left_on=self.to_gene_id_name, right_on=self.to_homolog_attribute)

            to_gene_id2whatever = to_gene_id2whatever.merge(from_whatever2gene_id, how='outer',
                                               left_on=self.from_gene_id_name, right_on=self.from_gene_id_name)
            to_gene_id2whatever = to_gene_id2whatever.drop(columns=self.to_homolog_attribute)

            self.to_gene_id2whatever = pd.concat([self.to_gene_id2whatever, to_gene_id2whatever], axis=0)

        return self.to_gene_id2whatever

class findOrthologsHsMm(FindOrthologs):

    def __init__(
          self,
          from_filters = 'hgnc_symbol',
          from_attributes: list = ['hgnc_symbol'],
          from_values: list = ['TP53', 'TERT', 'EPCAM'],
          to_attributes: list = ['external_gene_name']
      ):

        ############# Initialise parameters ################
          super().__init__(host = 'http://www.ensembl.org',
                mart = 'ENSEMBL_MART_ENSEMBL',
                from_dataset = 'hsapiens_gene_ensembl',
                to_dataset = 'mmusculus_gene_ensembl',
                from_filters = from_filters,
                from_attributes = from_attributes,
                from_values = from_values,
                to_attributes = to_attributes,
                to_homolog_attribute = 'mmusculus_homolog_ensembl_gene',
                from_gene_id_name = 'human_ensembl_gene_id',
                to_gene_id_name = 'mouse_ensembl_gene_id')


class findOrthologsMmHs(FindOrthologs):

    def __init__(
          self,
          from_filters = 'link_ensembl_gene_id',
          from_attributes: list = ['external_gene_name'],
          from_values: list = ['ENSMUSG00000059552', 'ENSMUSG00000021611'],
          to_attributes: list = ['hgnc_symbol']
      ):

        ############# Initialise parameters ################
          super().__init__(host = 'http://www.ensembl.org',
                mart = 'ENSEMBL_MART_ENSEMBL',
                from_dataset = 'mmusculus_gene_ensembl',
                to_dataset = 'hsapiens_gene_ensembl',
                from_filters = from_filters,
                from_attributes = from_attributes,
                from_values = from_values,
                to_attributes = to_attributes,
                to_homolog_attribute = 'hsapiens_homolog_ensembl_gene',
                from_gene_id_name = 'mouse_ensembl_gene_id',
                to_gene_id_name = 'human_ensembl_gene_id')
