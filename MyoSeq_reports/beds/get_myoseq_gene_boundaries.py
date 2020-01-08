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
    return sorted(set(gene_list))


def parse_gtf(gene_list, gtf) -> dict:
    '''
    Parses Gencode GTF for start and end positions of MyoSeq genes. 
    NOTE: Gencode format described at https://www.gencodegenes.org/pages/data_format.html

    :param set gene_list: Set containing MyoSeq genes
    :param str gtf: Path to Gencode GTF
    :return: Dictionary; key: gene, value: list of (chrom, start, end) 
    :rtype: dict
    '''
    gene_boundaries = {}
    with open(gtf) as g:
        for line in g:
            line = line.strip().split('\t')

            # split last section in line; this is the part of the gtf that contains the gene name
            for item in line[-1].split(';'):

                if 'gene_name' in item:
                    gene = item.replace('"', '').split(' ')[-1]

                    if gene in gene_list:
                        chrom = line[0]
                        start = int(line[3]) - 1
                        end = int(line[4])

                        if gene in gene_boundaries:
                            gene_boundaries[gene].append((chrom, start, end))
                            #gene_boundaries[gene] = (chrom, min(start, gene_boundaries[gene][1]), max(end, gene_boundaries[gene][2]))
                        else:
                            gene_boundaries[gene] = [(chrom, start, end)]

    # sanity check that dictionary is same length as gene list
    if len(gene_boundaries) != len(gene_list):
        logger.warn('Length of gene boundaries dictionary is not the same as length of gene list. Make sure all genes are present in GTF.')

    return gene_boundaries 


def write_beds(gene_boundaries, out) -> None:
    '''
    Writes BED files for MyoSeq genes

    :param dict gene_boundaries: Dictionary of MyoSeq genes (key; str) and their chromosome/starts/ends (value; str, int, int)
    :return: None
    :rtype: None
    '''
    for gene in gene_boundaries:
        out_bed = f'{out}/{gene}.bed'
        with open(out_bed, 'w') as o:
            lines = sorted(set(gene_boundaries[gene]))
            for line in lines:
                o.write(f'{line[0]}\t{line[1]}\t{line[2]}\n')
            #o.write(f'{gene_boundaries[gene][0]}\t{gene_boundaries[gene][1]}\t{gene_boundaries[gene][2]}\n')


def main(args):

    logger.info('Getting MyoSeq genes')
    gene_list = get_genes(args.genes)

    logger.info('Parsing GTF')
    gene_boundaries = parse_gtf(gene_list, args.gtf)

    logger.info('Writing out BEDs')
    write_beds(gene_boundaries, args.out)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Gets starts and ends of MyoSeq genes')
    parser.add_argument('-l', '--genes', help='Text file of MyoSeq genes, one gene per line', required=True)
    parser.add_argument('-g', '--gtf', help='Gencode GTF with gene boundaries', required=True)
    parser.add_argument('-o', '--out', help='Output directory', required=True)
    args = parser.parse_args()

    main(args)
