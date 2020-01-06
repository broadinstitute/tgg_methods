import argparse
import logging


logging.basicConfig(format='%(asctime)s (%(name)s %(lineno)s): %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger('get_myoseq_gene_boundaries')
logger.setLevel(logging.INFO)


def get_genes(genes) -> set:
    '''
    Opens text file with  MyoSeq gene list and stores genes as list

    :param str genes: Path to text file with MyoSeq gene list
    :return: Set containing MyoSeq genes
    :rtype: set
    '''
    gene_list = []
    with open(genes) as g:
        for line in g:
            gene_list.append(line.strip())
    return set(gene_list)


def parse_gtf(gene_list, gtf) -> dict:
    '''
    Parses Gencode GTF for start and end positions of MyoSeq genes. 
    NOTE: Gencode format described at https://www.gencodegenes.org/pages/data_format.html

    :param set gene_list: Set containing MyoSeq genes
    :param str gtf: Path to Gencode GTF
    :return: Dictionary; key (str): gene, value (ints): (start, end) 
    :rtype: dict
    '''
    gene_boundaries = {}
    with open(gtf) as g:
        for line in g:
            line = line.strip().split('\t')

            # split last section in line; this is the part of the gtf that contains the gene name
            for item in line[-1].split(';'):

                if 'gene_name' in item:
                    gene = item.split(' ').replace('"', '')

                    if gene in gene_list:
                        start = int(line[3])
                        end = int(line[4])

                        # if gene is already in dictionary, take smallest start and largest end
                        if gene in gene_boundaries:
                            gene_boundaries[gene] = (min(start, gene_boundaries[gene][0]), max(end, gene_boundaries[gene][1]))
                        else:
                            gene_boundaries[gene] = (start, end)

    # sanity check that dictionary is same length as gene list
    if len(gene_boundaries) != len(gene_list):
        logger.warn('Length of gene boundaries dictionary is not the same as length of gene list. Make sure all genes are present in GTF.')

    return gene_boundaries 


def main(args):

    logger.info('Getting MyoSeq genes')
    gene_list = get_genes(args.genes)

    logger.info('Parsing GTF')
    gene_boundaries = parse_gtf(gene_list, args.gtf)

    logger.info('Writing out BEDs'

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Gets starts and ends of MyoSeq genes')
    parser.add_argument('-l', '--genes', help='Text file of MyoSeq genes, one gene per line', required=True)
    parser.add_argument('-g', '--gtf', help='Gencode GTF with gene boundaries', required=True)
    parser.add_argument('-o', '--out', help='Output directory', required=True)
    args = parser.parse_args()

    main(args)