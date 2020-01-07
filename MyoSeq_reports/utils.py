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
    return set(sorted(gene_list, reverse=True))


def get_samples(sampleFile):
    """
    Get requested sample IDs from input file (one sample per line)

    :param str samp: Name of file containing requested samples
    :return: List of sample IDs
    :rtype: list
    """
    samples = []
    with open(sampleFile) as s:
        for line in s:
            samples.append(line.strip())
    return samples


def check_missing_samples(samples, check, fname):
    """
    Checks if any requested samples are missing from a file

    :param list samples: List of strings (requested sample IDs)
    :param list check: List of strings (sample IDs found in file)
    :param str fname: File name
    :return: None
    :rtype: None
    """
    missing = set(samples) - set(check)
    if len(missing) > 0:
        logging.warning('{} not found in {} file'.format(','.join(missing), fname))
