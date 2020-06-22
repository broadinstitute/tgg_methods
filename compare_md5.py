import argparse
import subprocess
import csv
import logging

from gnomad.utils.file_utils import get_file_stats

from google.cloud import storage

logging.basicConfig(
    format="%(asctime)s): %(message)s", datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger()
logger.setLevel(logging.INFO)


def compare_md5_dict(origin_md5_dict, dest_md5_dict):
    """Compare the md5 origin and destination dictionaries. If md5s do not match,
    document in dictionary. If file is missing, document in set of missing files.
    :param origin_md5_dict: dictionary of paths and md5 values for origin files
    :param dest_md5_dict: dictionary of paths and md5 values for destination files
    :return dictionary of filesets and unequal md5 value, set of missing files
    """
    needs_attention = {}
    files_missing = set([])
    logger.info('In comparison module')
    for key in origin_md5_dict:
        if key in dest_md5_dict:
            if not origin_md5_dict[key]['md5'] == dest_md5_dict[key]['md5']:
                needs_attention[key] = [origin_md5_dict[key]['path'], dest_md5_dict[key]['path']]
        else:
            files_missing.add(origin_md5_dict[key]['path'])
    return needs_attention, files_missing


def create_md5_dict(file_dict)->dict:
    """Create a dictionary for each subfile in a directory with file key 
    being directory + subfile path and values being the full path and md5
    containing md5s.
    :param file_dict: dictionary of directories 
    :return dictionary of files, paths, and md5 value
    """
    md5_dict = {}
    for key in file_dict:
        logger.info(f'Examining the directory: {key}')
        for subfile in file_dict[key]:
            subfile_key = subfile[subfile.find(key.split("/")[-1]):]
            size, int_size, md5 = get_file_stats(subfile)
            md5_dict[subfile_key] = {'path' : subfile, 'md5' : md5}
    return md5_dict


def read_hail_paths(bucket, bucket_file):
    """Read file from google bucket and create sets of origin files 
    and destination files
    :param bucket: top level google bucket with project e.g. "gnomad-public"
    :param bucket_file: Path to file containing columns of origin and destination files
    :return set of origin files, set of destination files
    """
    bucket.get_blob(bucket_file).download_to_filename('paths')
    origin = []
    dest = []
    with open('paths', 'r') as paths:
        csv_input = csv.reader(paths, delimiter='\t', skipinitialspace=True)
        for rows in csv_input:
            origin.append(rows[0])
            dest.append(rows[1])
    return origin, dest


def create_file_dict(directory_list):
    """Generate a dictionary where top level directories are keys and the set of
    all files in each top level directory are the values
    :param directory_list: List of top level directories
    :return file_dict: dictionary of directories and their subfiles
    """
    file_dict = {}
    for directory in directory_list:
        logger.info(f'Building dictionary for {directory}')
       # Remove file string from file name and then convert "/" to "_", no need for new key assignm
        subfiles = (
                subprocess.check_output(["gsutil", "ls", "-r", directory])
                .decode("utf8")
                .strip()
                .split("\n")
            )
        subfiles = [f for f in subfiles if not (f.endswith("/") or f.endswith(":") or f == '')]
        file_dict[directory] = subfiles
    return file_dict


def main(args):
    input_bucket_file = args.input_bucket_file
    output_bucket_file = args.output_bucket_file

    client = storage.Client()
    bucket = client.get_bucket(args.bucket)

    origin, dest = read_hail_paths(bucket, input_bucket_file)
    # Get names of vcf shards
    if len(origin) != len(dest):
        raise ValueError("Number of origins does not equal numbers of destinations")
    else:
        origin_file_dict = create_file_dict(origin)
        dest_file_dict = create_file_dict(dest)
        origin_md5_dict = create_md5_dict(origin_file_dict)
        dest_md5_dict = create_md5_dict(dest_file_dict)
        needs_attention, files_missing = compare_md5_dict(origin_md5_dict, dest_md5_dict)
        
        output = bucket.blob(output_bucket_file)
        results = []
        if len(needs_attention) != 0:
            logger.info(f'md5 mismatch:\n {needs_attention}')
            results.append(f'md5 mismatch: {needs_attention}')
        if len(files_missing) != 0:
            logger.info(f'Files not copied:\n {files_missing}')
            results.append(f'Files not copied: {files_missing}')
        
        output.upload_from_string('\n'.join(results))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-bucket-file", help="File with two columns: origin and destination", required=True)
    parser.add_argument("--output-bucket-file", help="File to output comparison results", required=True)
    parser.add_argument("--bucket", help="Bucket where input and ouptut files live", required=True)

    args = parser.parse_args()
    main(args)