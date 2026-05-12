#!/usr/bin/env python3
import sys, os, shutil, tempfile, subprocess, argparse, logging

__version__ = "0.6.7"
logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.DEBUG)

class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter): pass

desc = "parallel fastq-dump wrapper, extra args will be passed through"
epi = """DESCRIPTION:
Example: parallel-fastq-dump --sra-id SRR2244401 --threads 4 --outdir out/ --split-files --gzip
"""

parser = argparse.ArgumentParser(description=desc, epilog=epi, formatter_class=CustomFormatter)
parser.add_argument("-s","--sra-id", help="SRA id", action="append")
parser.add_argument("-t","--threads", help="number of threads", default=1, type=int)
parser.add_argument("-O","--outdir", help="output directory", default=".")
parser.add_argument("-T","--tmpdir", help="temporary directory", default=None)
parser.add_argument("-N","--minSpotId", help="Minimum spot id", default=1, type=int)
parser.add_argument("-X","--maxSpotId", help="Maximum spot id", default=None, type=int)
parser.add_argument("-V","--version", help="shows version", action="store_true", default=False)

def pfd(args: argparse.Namespace, srr_id: str, extra_args: list[str]) -> None:
    """Run fastq-dump in parallel by splitting the spot range across threads.

    Each thread writes its chunk to a temporary directory; chunks are then
    concatenated in order into the output directory.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments including threads, outdir, minSpotId,
        maxSpotId, and tmpdir.
    srr_id : str
        SRA run accession to download.
    extra_args : list of str
        Additional flags forwarded verbatim to each fastq-dump invocation.
    """
    tmp_dir = tempfile.TemporaryDirectory(prefix="pfd_", dir=args.tmpdir)
    logging.info(f"tempdir: {tmp_dir.name}")
    n_spots = get_spot_count(srr_id)
    logging.info(f"{srr_id} spots: {n_spots}")
    start = max(args.minSpotId, 1)
    end = min(args.maxSpotId, n_spots) if args.maxSpotId is not None else n_spots
    blocks = split_blocks(start, end, args.threads)
    logging.info(f"blocks: {blocks}")
    ps = []
    for i in range(args.threads):
        d = os.path.join(tmp_dir.name, str(i))
        os.mkdir(d)
        cmd = ["fastq-dump","-N",str(blocks[i][0]),"-X",str(blocks[i][1]),"-O",d]+extra_args+[srr_id]
        logging.info(f"CMD: {' '.join(cmd)}")
        p = subprocess.Popen(cmd)
        ps.append(p)
    wfd = {}
    for i,p in enumerate(ps):
        exit_code = p.wait()
        if exit_code != 0:
            logging.warning(f"fastq-dump error! exit code: {exit_code}")
            sys.exit(1)
        tmp_path = os.path.join(tmp_dir.name, str(i))
        for fo in os.listdir(tmp_path):
            if fo not in wfd: wfd[fo] = open(os.path.join(args.outdir, fo), "wb")
            with open(os.path.join(tmp_path, fo), "rb") as fd:
                shutil.copyfileobj(fd, wfd[fo])
            os.remove(os.path.join(tmp_path, fo))
    for fd in wfd.values(): fd.close()

def split_blocks(start: int, end: int, n_pieces: int) -> list[list[int]]:
    """Divide a spot-ID range into equal-sized blocks for parallel processing.

    The final block absorbs any remainder from integer division.

    Parameters
    ----------
    start : int
        First spot ID (inclusive).
    end : int
        Last spot ID (inclusive).
    n_pieces : int
        Number of blocks to produce.

    Returns
    -------
    list of list of int
        Each inner list is ``[block_start, block_end]`` (both inclusive).
    """
    total = end - start + 1
    avg = total // n_pieces
    out = []
    last = start
    for i in range(n_pieces):
        out.append([last, last + avg - 1])
        last += avg
        if i == n_pieces - 1: out[i][1] += total % n_pieces
    return out

def get_spot_count(sra_id: str) -> int:
    """Return the total spot count for an SRA run using sra-stat.

    Parameters
    ----------
    sra_id : str
        SRA run accession to query.

    Returns
    -------
    int
        Total number of spots.

    Raises
    ------
    IndexError
        If the sra-stat output cannot be parsed.
    """
    cmd = ["sra-stat","--meta","--quick",sra_id]
    logging.info(f"CMD: {' '.join(cmd)}")
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    txt = stdout.decode().rstrip().split("\n")
    total = 0
    try:
        for l in txt: total += int(l.split("|")[2].split(":")[0])
    except IndexError:
        msg = "sra-stat output parsing error!\n--sra-stat STDOUT--\n{}\n--sra-stat STDERR--\n{}"
        raise IndexError(msg.format("\n".join(txt), stderr.decode().rstrip()))
    return total

def partition(f, l: list) -> tuple[list, list]:
    """Partition a list into two groups based on a boolean predicate.

    Parameters
    ----------
    f : callable
        Function that returns True or False for each element.
    l : list
        List to partition.

    Returns
    -------
    tuple
        A 2-tuple ``(matching, not_matching)`` where each element is a list.
    """
    r = ([],[])
    for i in l: r[0 if f(i) else 1].append(i)
    return r

def is_sra_file(path: str) -> bool:
    """Determine whether a path looks like an SRA file or accession.

    Parameters
    ----------
    path : str
        File path or accession string to inspect.

    Returns
    -------
    bool
        True if the path ends with ``.sra`` or contains a recognised accession
        prefix (SRR, ERR, DRR).
    """
    f = os.path.basename(path)
    if f.lower().endswith(".sra"): return True
    if any(x in f.upper() for x in ["SRR","ERR","DRR"]): return True
    return False

def main() -> None:
    """Parse arguments and dispatch parallel fastq-dump for each provided SRA ID."""
    args, extra = parser.parse_known_args()
    if args.version:
        print(f"parallel-fastq-dump : {__version__}")
        subprocess.Popen(["fastq-dump","-V"]).wait()
        sys.exit(0)
    elif args.sra_id:
        extra_srrs, extra_args = partition(is_sra_file, extra)
        args.sra_id.extend(extra_srrs)
        logging.info(f"SRR ids: {args.sra_id}")
        logging.info(f"extra args: {extra_args}")
        if not os.path.isdir(args.outdir) and args.outdir != ".":
            os.makedirs(args.outdir)
        if args.tmpdir and not os.path.isdir(args.tmpdir) and args.tmpdir != ".":
            os.makedirs(args.tmpdir)
        for si in args.sra_id: pfd(args, si, extra_args)
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()
