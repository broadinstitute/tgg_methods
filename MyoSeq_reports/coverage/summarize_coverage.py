import argparse
import gzip
import logging
import typing
from statistics import mean


logging.basicConfig(format='%(asctime)s (%(name)s %(lineno)s): %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger('summarize_gene_coverage')
logger.setLevel(logging.INFO)


def get_genes(genes) -> set:
    '''
    Opens text file with MyoSeq gene list and stores genes as list

    :param str genes: Path to text file with MyoSeq gene list
    :return: Set containing MyoSeq genes
    :rtype: set
    '''
    gene_list = []
    with open(genes) as g:
        for line in g:
            gene_list.append(line.strip())
    return sorted(set(gene_list))


def get_cram_order(cramslist) -> list:
    '''
    Opens text file with cram paths and stores samples, in order, in a list

    :param str cramslist: Path to text file with crams, one cram per line
    :return: List of samples (str)
    :rtype: list
    '''
    samples = []
    with open(cramslist) as c:
        for line in c:
            samples.append(line.strip().split('/')[-1].split('.')[0])
    return samples


def get_gene_size(bed: str) -> int:
    '''
    Opens bed file for gene and returns number of bases in gene

    :param str bed: Path to bedfile
    :return: Number of bases in gene
    :rtype: int
    '''
    total_bases = 0
    with open(bed) as b:
        chrom, start, end = b.readline().strip().split('\t')
        total_bases += int(end) - int(start)
    return total_bases


def summarize_coverage(gene_list: set, samples: list, covdir: str,  beddir: str) -> dict:
    '''
    Opens coverage files for processing. Summarizes coverage across gene across batch of samples.

    :param set gene_list: List of MyoSeq genes
    :param list samples: List of MyoSeq samples (order from cram list)
    :param str covdir: Directory containing all MyoSeq gene coverage files
    :param str beddir: Directory containing all MyoSeq genes
    :return: Dictionary of coverage metrics per gene: mean coverage, percent callable sites, and number uncallable sites, both across the batch and per sample
    :rtype: dict
    '''
    summary = {}
    for gene in gene_list:

        # set up variables for gene
        summary[gene] = {}
        tsv = f'{covdir}/{gene}.tsv.bgz'
        bed = f'{beddir}/{gene}.bed'
        total_bases = get_gene_size(bed)
        total_cov = 0
        callable_bases = 0
        logger.info(f'Working on {gene}...')
        logger.info(f'Number of bases: {total_bases}')

        # create dict to store coverage per sample
        sample_cov = {}
        sample_callable = {}
        for sample in samples:
            sample_cov[sample] = 0
            sample_callable[sample] = 0

        with gzip.open(tsv, 'rt') as t: 
            for line in t:
                line = line.strip().split('\t')
                chrom = line[0]
                pos = line[1]

                # check if position is considered "callable" (mean depth > 6, as defined by Monkol)
                mean_at_pos = mean(map(int, line[2:]))
                total_cov += mean_at_pos
                if mean_at_pos > 6:
                    callable_bases += 1

                # get sample coverage -- sample order in coverage file is same as cram order, with indices shifted by 2
                for i in range(2, len(line)):
                    cov = int(line[i])
                    if cov > 0:
                        sample_callable[samples[i-2]] += 1
                        sample_cov[samples[i-2]] += cov
                        sample_callable[samples[i-2]] += 1

        mean_cov = round(float(total_cov) / total_bases, 2)
        summary[gene]['mean'] = mean_cov
        summary[gene]['callable'] = round(float(callable_bases) / total_bases, 2)
        summary[gene]['uncallable'] = total_bases - callable_bases
        for sample in samples:
            logger.info(f'Working on {sample}...')
            summary[sample] = {}
            summary[sample][gene] = {}
            summary[sample][gene]['mean'] = round(float(sample_cov[sample]) / total_bases, 2)
            summary[sample][gene]['callable'] = round(float(sample_uncallable[sample]) / total_bases,  2)
            summary[sample][gene]['uncallable'] = total_bases - sample_uncallable[sample]
    return summary


def write_summary(samples: list, gene_list: set, summary: dict, out: str) -> None:
    '''
    Writes per sample summaries for all MyoSeq genes

    :param list summary: List of all samples
    :param set gene_list: List of MyoSeq genes
    :param dict summary: Summarized coverage information
    :param str out: Output directory
    :return: None
    :rtype: None
    '''
    for sample in samples:
        tsv = f'{out}/{sample}_coverage.tsv'

        with open(tsv, 'w') as t:

            # write header
            t.write('Gene\tCohort_Mean\tSample_Mean\tCohort_Callable\tSample_Callable\tCohort_Uncallable\tSample_Uncallable\n')

            # write summary per gene
            for gene in summary[sample]:
                cohort_mean = summary[gene]['mean']
                cohort_callable = summary[gene]['callable']
                cohort_uncallable = summary[gene]['uncallable']
                sample_mean = summary[sample][gene]['mean']
                sample_callable = summary[sample][gene]['callable']
                sample_uncallable = summary[sample][gene]['uncallable']
                t.write(f'{gene}\t{cohort_mean}\t{sample_mean}\t{cohort_callable}\t{sample_callable}\t{cohort_uncallable}\t{sample_uncallable}\n') 

def main(args):

    logger.info('Getting MyoSeq genes')
    gene_list = get_genes(args.genes)

    logger.info('Getting samples from list of crams')
    samples = get_cram_order(args.cramslist)

    logger.info('Summarizing cram coverage')
    summary = summarize_coverage(gene_list, samples, args.covdir, args.beddir)

    logger.info('Writing output files')
    write_summary(samples, gene_list, summary, args.out)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Summarizes coverage over MyoSeq genes')
    parser.add_argument('-g', '--genes', help='Text file of MyoSeq genes, one gene per line', required=True)
    parser.add_argument('-b', '--beddir', help='Directory containing BED files for all MyoSeq genes', required=True)
    parser.add_argument('-d', '--covdir', help='Directory containing coverage files', required=True)
    parser.add_argument('-c', '--cramslist', help='List of crams, one cram per line', required=True)
    parser.add_argument('-o', '--out', help='Output directory', required=True)
    args = parser.parse_args()

    main(args)
