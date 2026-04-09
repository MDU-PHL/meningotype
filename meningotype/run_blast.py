import subprocess


def seqBLAST(query, db, blast, outfmt, perc_identity=90, evalue='1e-20', num_threads=1, culling_limit=1, dust=None, fasta_data=None, extra_flags=None):
    command = [
        blast,
        "-query", query,
        "-db", db,
        "-evalue", evalue,
        "-outfmt", outfmt,
        "-culling_limit", str(culling_limit),
        "-num_threads", str(num_threads)
    ]
    if blast == 'blastn':
        command += ["-perc_identity", str(perc_identity)]
        if dust is not None:
            command += ["-dust", dust]
    if extra_flags:
        command += extra_flags
    if query == '-':
        result = subprocess.run(
            command,
            input=fasta_data,
            capture_output=True,
            text=True,
            check=True
        )
    else:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=True
        )
    return result.stdout, result.stderr