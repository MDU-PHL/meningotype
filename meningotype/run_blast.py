import subprocess


# NcbiblastnCommandline(query=f, db=allelesdb, task='blastn', perc_identity=90, evalue='1e-20', outfmt='"6 sseqid pident length"', culling_limit='1', num_threads=cpus)

def seqBLAST(query, db, blast, outfmt, perc_identity=90, evalue='1e-20', num_threads=1, culling_limit=1, dust='no', fasta_data=None):
    command = [
			blast,
			"-query", query,
			"-db", db,
			"-perc_identity", str(perc_identity),
			"-evalue", evalue,
			"-outfmt", outfmt,
			"-culling_limit", str(culling_limit),
			"-num_threads", str(num_threads)
		]
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