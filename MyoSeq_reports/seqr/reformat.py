import argparse
import logging
import typing


logging.basicConfig(format='%(asctime)s (%(name)s %(lineno)s): %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger('reformat')
logger.setLevel(logging.INFO)


def get_sex(sex: str) -> dict:
    '''
    Extract inferred sex from sex tsv (created using sexcheck script)

    :param str sex: Input tsv with sex information
    :return: Dictionary; key: sample, value: inferred sex
    :rtype: dict
    '''
    samples = {}
    with open(sex) as s:
        for line in s:
            sample, reported_sex, inferred_sex, fstat, status = line.strip().split('\t')
            samples[sample] = inferred_sex
    return samples


def get_gnomad_info(tsv: str) -> dict:
    '''
    Extract gnomAD popmax AF, popmax pop, and  global AF from TSV (export using hail after joining variant loci/alleles with gnomAD tables)

    :param str tsv: Path to input tsv
    :return: Dictionary of variants; key: variants, value: global AF, popmax AF, popmax pop
    :rtype: dict
    '''
    variants = {}
    
    with open(tsv) as t:

        header = t.readline().strip().split('\t')

        for line in t:
            # extract necessary fields
            line = line.strip().split('\t')
            locus = line[header.index('locus')]
            alleles = line[header.index('alleles')].replace('[', '').replace(']', '').replace('"', '').split(',')
            global_af = line[header.index('gnomad_global_AF')]
            # correct global AF (hail reports 0.0000e+00 if AF was 0.0 from hail)
            if global_af == '0.0000e+00':
                global_af = '0.0'
            popmax_af = line[header.index('gnomad_exomes_popmax_AF')]
            # correct missing values for popmax AF
            if popmax_af == 'NA':
                popmax_af = '0.0'
            popmax_pop = line[header.index('gnomad_exomes_popmax_pop')]

            # add variant info to dict
            variant = f'{locus} {alleles[0]}>{alleles[1]}'
            if variant not in variants:
                variants[variant] = {}
            variants[variant]['global'] = global_af
            variants[variant]['popmax'] = popmax_af
            variants[variant]['pop'] = popmax_pop

    return variants


def get_output_fields(seqr: str, variants: dict, sex: dict) -> dict:
    """
    Reads in seqr export file and stores each sample's variants and necessary information for output file in dict

    :param str seqr: Path to seqr export tsv
    :param dict all_variants: Dictionary of variants and their gnomAD info
    :param dict sex: Dictionary of each sample and their inferred sex
    :return: Dictionary of samples (key) and their variants with frequencies (value)
    :rtype: dict
    """

    samples = {}
    with open(seqr) as s:

        # get header of seqr export file -- format should be stable but double check
        header = s.readline().strip().split('\t')
        logger.info(f'seqr format check: {header[0:6]}, {header[18:23]}')
        logger.info('expected: [chrom pos ref alt gene] [rsid hgvsc hgvsp clinvar_clinical_significance clinvar_gold_stars]')

        # cycle through seqr export file and pull relevant fields
        for line in s:
            line = line.strip().split('\t')
            chrom, pos, ref, alt, gene, worst_consequence = line[0:6]
            # add 'chr' to chrom
            if 'chr' not in chrom:
                chrom = f'chr{chrom}'
            rsid, hgvsc, hgvsp, clinvar_clinsig, clinvar_gold_stars = line[18:23]
            notes = line[26]

            # check if notes field has PMID; if not, don't include notes
            pmid = ''
            if 'PMID' in notes:
                notes = notes.split(' ')
                for i in range(len(notes)):
                    if 'PMID' in i:
                        pmid += f'{notes[i]}{notes[i+1]}'
            else:
                pmid = '.'

            # set clinvar values to '.' if they don't exist
            clinvar_clinsig = clinvar_clinsig.replace('_', ' ').lower()
            if len(clinvar_clinsig) == 0:
                clinvar_clinsig = '.'
            if len(clinvar_gold_stars) == 0:
                clinvar_gold_stars = '.'

            # get variant (format matches what will be output in reports)
            variant = '{}:{} {}>{}'.format(chrom, pos, ref, alt)

            # get gene name from column - check for seqr gene name bug
            if ',' in gene:
                gene = gene.replace("'","")

            # get popmax and global_af from files; if key doesn't exist, then use '0.0' (not seen in gnomAD)
            global_af = variants[variant]['global']
            popmax_af = variants[variant]['popmax']
            pop = variants[variant]['pop']

            # check all items in seqr download starting from index 27 (first index to potentially have sample information)
            for item in line[27:]:

                # sample_1:num_alt_alleles:gq:ab
                # skip if item doesn't have colon -- this means it doesn't have sample info
                if ':' not in item: 
                    continue

                # add sample to dict
                info = item.split(':')
                s = info[0]
                if s not in samples:
                    samples[s] = []

                # get sample genotype
                try:
                    geno = int(info[1])
                except:
                    logger.warning(f'{s} has non-int genotype. Skipping {s}')
                    continue

                # only keep heterozygous or homozygous variant genotypes
                if geno == 0:
                    continue
                elif geno == 1:
                    gt = 'het'
                else:
                    if geno != 2:
                        logger.warning(f'{s} genotype is {geno}; not 0, 1, or 2. Assuming this means {s} is heterozygous')
                        gt = 'het'
                    else:
                        if (chrom == 'X' or chrom == 'chrX') and sex[s] == 'Male':
                            gt = 'hemi'
                        else:
                            gt = 'hom'
                        
                # add formatted line to sample dict
                # columns needed for output: gene, genotype, variant, functional effect, hgvsc, hgvsp, rsid, gnomad global AF, gnomad popmax AF, gnomad popmax pop, clinvar clinical significance, clinvar clinical review status (no longer in seqr downloads), number of clinvar stars, clinvar url (also no longer in seqr download), PMID 
                samples[s].append([gene, gt, variant, worst_consequence, hgvsc, hgvsp, rsid, global_af, popmax_af, pop, clinvar_clinsig, '.', clinvar_gold_stars, '.', pmid])

    return samples
        

def write_out(samples: dict, out: str, report: bool):
    """
    Writes out files for each sample for MYOSEQ reports

    :param dict samples: Dictionary of samples (key) and their variants (with frequencies; value)
    :param str out: Name of output file
    :param bool report: Whether output file should be named 'flagged' (REPORT variants') or 'genes' (all rare variants')
    :return: None
    :rtype: None
    """

    logger.info(f'Number of samples found: {len(samples)}')
    for sample in samples:
        if not report:
            outfile = f'{out}/{sample}.genes.txt'
        else:
            outfile = f'{out}/{sample}.flagged.txt'

        # open a file for each sample in specified output directory
        with open(outfile, 'w') as o:

            # write header to output
            o.write('gene_name\tgenotype\tvariant\tfunctional_class\thgvs_c\thgvs_p\trsid\tgnomad_global_af\tgnomad_pop_max_af\tgnomad_pop_max_population\tclinvar_clinsig\tclinvar_clnrevstat\tnumber_of_stars\tclinvar_url\tcomments\n')

            # write variants to output
            for line in samples[sample]:
                o.write('\t'.join(line) + '\n')


def main(args):

    logger.info('Getting sex information for samples')
    sex = get_sex(args.sex)

    logger.info('Getting gnomAD global AF, popmax AF, and popmax pop from hail tsv export')
    variants = get_gnomad_info(args.tsv)

    logger.info('Getting seqr export information')
    samples = get_output_fields(args.seqr, variants, sex)
    
    logger.info('Writing output files')
    write_out(samples, args.out, args.report)


if __name__ == '__main__':
    
    # Get args from command line
    parser = argparse.ArgumentParser(description='Reformats seqr variant tsv export in preparation for downstream pdf generation')
    parser.add_argument('-x', '--sex', help='Input file of samples and their sex information', required=True)
    parser.add_argument('-s', '--seqr', help='Input file of variants exported from seqr', required=True)
    parser.add_argument('-t', '--tsv', help='TSV exported by hail after joining seqr variant loci/alleles to gnomAD hail tables', required=True)
    parser.add_argument('-r', '--report', help='Output file is for variants tagged REPORT', action='store_true')
    parser.add_argument('-o', '--out', help='Output directory', required=True)
    args = parser.parse_args()

    main(args)
