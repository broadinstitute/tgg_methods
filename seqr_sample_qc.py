from gnomad_hail import *
from gnomad_hail.resources.sample_qc import *

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("seqr_sample_qc")
logger.setLevel(logging.INFO)


def main(args):
    hl.init()

    print("HELLO WORLD")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--my-first-argument')
    args = parser.parse_args()

    main(args)
