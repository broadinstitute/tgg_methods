import argparse
import sys
from collections import namedtuple
from typing import List, Tuple
from urllib.parse import urlparse

from google.cloud.storage import Client


def _remove_prefix(s: str, prefix: str) -> str:
    if s.startswith(prefix):
        return s[len(prefix) :]
    return s


GCSDiffResult = namedtuple(
    "GCSDiffResult", ["missing_objects", "extra_objects", "changed_objects"]
)


def diff_gcs_directories(
    base_directory_url: str, target_directory_url: str
) -> Tuple[List[str], List[str], List[str]]:
    """
    Compare objects under different GCS prefixes.

    :param base_directory_url: URL for base directory
    :param target_directory_url: URL for target directory

    :returns: Tuple with 3 elements:
        List of objects in base directory that are not present in target directory
        List of objects in target directory that are not present in base directory
        List of objects with different content in base and target directory
    """
    base = urlparse(base_directory_url)
    target = urlparse(target_directory_url)

    if base.scheme != "gs":
        raise ValueError("base_directory_url must be a gs:// URL")

    if target.scheme != "gs":
        raise ValueError("target_directory_url must be a gs:// URL")

    client = Client()

    base_blobs = client.list_blobs(base.hostname, prefix=base.path.strip("/") + "/")
    base_blobs = {
        _remove_prefix(blob.name, base.path.strip("/")): blob for blob in base_blobs
    }

    missing_objects = set(base_blobs.keys())
    extra_objects = []
    changed_objects = []

    target_blobs = client.list_blobs(
        target.hostname, prefix=target.path.strip("/") + "/"
    )

    for blob in target_blobs:
        key = _remove_prefix(blob.name, target.path.strip("/"))

        missing_objects.discard(key)

        try:
            if blob.md5_hash != base_blobs[key].md5_hash:
                changed_objects.append(key)
        except KeyError:
            extra_objects.append(key)

    return GCSDiffResult(list(missing_objects), extra_objects, changed_objects)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("base_directory")
    parser.add_argument("target_directory")
    args = parser.parse_args()

    diff = diff_gcs_directories(args.base_directory, args.target_directory)

    for obj in diff.missing_objects:
        print(f"- {obj}")

    for obj in diff.extra_objects:
        print(f"+ {obj}")

    for obj in diff.changed_objects:
        print(f"* {obj}")

    if diff.missing_objects or diff.extra_objects or diff.changed_objects:
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
