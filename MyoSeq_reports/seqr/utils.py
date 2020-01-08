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
