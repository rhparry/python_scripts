#Rhys Parry r.parry@uq.edu.au Assembly python script for use with UQ HPC. It downloads SRA files to the $TMPDIR generates fastq files and then uses the outputs for MEGAHIT. Requires a accessions file, sra-tools and megahit in the environment. 
import os
import subprocess
import shutil

# Define paths
BASE_OUTPUT_DIR = "/QRISdata/Q5900/Plas_Mary"
FINAL_CONTIGS_DIR = os.path.join(BASE_OUTPUT_DIR, "contigs")
COMPLETED_FILE = os.path.join(BASE_OUTPUT_DIR, "completed.txt")
FAILED_FILE = os.path.join(BASE_OUTPUT_DIR, "failed_accessions.txt")
TMPDIR = os.getenv("TMPDIR", "/tmp")  # Ensure TMPDIR is set

# Create necessary directories
os.makedirs(FINAL_CONTIGS_DIR, exist_ok=True)
os.makedirs(TMPDIR, exist_ok=True)
with open(COMPLETED_FILE, "a") as f:
    pass  # Ensure the file exists
with open(FAILED_FILE, "a") as f:
    pass  # Ensure the file exists

# Define the singularity-based commands
PREFETCH_CMD = (
    "singularity exec -B /sw/local/rocky8.6/noarch/qcif/shpc_modules/"
    "quay.io/biocontainers/sra-tools/3.0.3--h87f3376_0//99-shpc.sh:/.singularity.d/env/99-shpc.sh "
    "/sw/local/rocky8.6/noarch/qcif/containers/quay.io/biocontainers/sra-tools/3.0.3--h87f3376_0/"
    "quay.io-biocontainers-sra-tools-3.0.3--h87f3376_0-sha256:c9f92683e10091c3ef93066b1fcbdeeba89af49242ab778a9c8cc006f6be82a3.sif "
    "/usr/local/bin/prefetch.3"
)
FASTERQ_DUMP_CMD = (
    "singularity exec -B /sw/local/rocky8.6/noarch/qcif/shpc_modules/"
    "quay.io/biocontainers/sra-tools/3.0.3--h87f3376_0//99-shpc.sh:/.singularity.d/env/99-shpc.sh "
    "/sw/local/rocky8.6/noarch/qcif/containers/quay.io/biocontainers/sra-tools/3.0.3--h87f3376_0/"
    "quay.io-biocontainers-sra-tools-3.0.3--h87f3376_0-sha256:c9f92683e10091c3ef93066b1fcbdeeba89af49242ab778a9c8cc006f6be82a3.sif "
    "/usr/local/bin/fasterq-dump.3"
)

# Maximum retries
MAX_RETRIES = 2

# Process each accession in the list
with open("accessions_list.txt", "r") as accessions:
    for accession in accessions:
        accession = accession.strip()
        if not accession:
            continue

        # Check if accession has already been processed
        if os.path.isfile(COMPLETED_FILE):
            with open(COMPLETED_FILE, "r") as completed:
                if accession in completed.read():
                    print(f"Skipping {accession}: already completed.")
                    continue

        retry_count = 0
        while retry_count <= MAX_RETRIES:
            # Download the .sra file
            print(f"Downloading data for {accession} with prefetch.3...")
            if subprocess.run(f"{PREFETCH_CMD} {accession} -O {TMPDIR}", shell=True).returncode != 0:
                print(f"prefetch.3 failed for {accession}. Retrying... ({retry_count})")
                retry_count += 1
                continue

            # Locate the .sra file
            accession_dir = os.path.join(TMPDIR, accession)
            sra_file = os.path.join(accession_dir, f"{accession}.sra")
            if not os.path.exists(sra_file):
                print(f"Error: {sra_file} not found. Retrying... ({retry_count})")
                retry_count += 1
                continue

            # Validate the .sra file
            print(f"Validating {sra_file}...")
            validation_cmd = (
                f"singularity exec -B {TMPDIR}:/mnt "
                "/sw/local/rocky8.6/noarch/qcif/containers/quay.io/biocontainers/sra-tools/3.0.3--h87f3376_0/"
                "quay.io-biocontainers-sra-tools-3.0.3--h87f3376_0-sha256:c9f92683e10091c3ef93066b1fcbdeeba89af49242ab778a9c8cc006f6be82a3.sif "
                f"/usr/local/bin/vdb-validate /mnt/{accession}/{accession}.sra"
            )
            if subprocess.run(validation_cmd, shell=True).returncode != 0:
                print(f"Validation failed for {accession}. Retrying... ({retry_count})")
                retry_count += 1
                continue

            # Run fasterq-dump.3
            print(f"Running fasterq-dump.3 for {accession}...")
            fasterq_cmd = (
                f"{FASTERQ_DUMP_CMD} --split-files --threads 32 --mem 64GB --outdir {accession_dir} {sra_file}"
            )
            if subprocess.run(fasterq_cmd, shell=True).returncode != 0:
                print(f"fasterq-dump.3 failed for {accession}. Retrying... ({retry_count})")
                retry_count += 1
                continue

            # Check for output FASTQ files
            print(f"Checking for FASTQ files for {accession} in {accession_dir}...")
            fastq_1 = os.path.join(accession_dir, f"{accession}_1.fastq")
            fastq_2 = os.path.join(accession_dir, f"{accession}_2.fastq")
            single_fastq = os.path.join(accession_dir, f"{accession}.fastq")

            if not (os.path.isfile(fastq_1) or os.path.isfile(fastq_2) or os.path.isfile(single_fastq)):
                print(f"Error: FASTQ files for {accession} not found. Retrying... ({retry_count})")
                retry_count += 1
                continue

            # Determine single-end or paired-end
            if os.path.isfile(fastq_1) and os.path.isfile(fastq_2):
                print(f"Paired-end reads detected for {accession}.")
                megahit_cmd = f"megahit -1 {fastq_1} -2 {fastq_2} -t 32 -o {accession_dir}/megahit_out"
            elif os.path.isfile(fastq_1):  # Assume single-end if only _1.fastq exists
                print(f"Single-end reads detected for {accession}.")
                megahit_cmd = f"megahit -r {fastq_1} -t 32 -o {accession_dir}/megahit_out"
            else:
                print(f"Error: FASTQ files for {accession} not found. Retrying... ({retry_count})")
                retry_count += 1
                continue

            # Run MEGAHIT
            print(f"Running MEGAHIT for {accession}...")
            if subprocess.run(megahit_cmd, shell=True).returncode == 0:
                contigs_file = os.path.join(accession_dir, "megahit_out", "final.contigs.fa")
                if os.path.isfile(contigs_file):
                    # Modify contig file headers
                    print(f"Modifying headers in contigs for {accession}...")
                    modified_contigs = os.path.join(FINAL_CONTIGS_DIR, f"{accession}_final_contigs.fa")
                    with open(contigs_file, "r") as infile, open(modified_contigs, "w") as outfile:
                        for line in infile:
                            if line.startswith(">"):
                                line = f">{accession}_" + line[1:].replace(" ", "_").replace("=", "_")
                            outfile.write(line)

                    # Log completion
                    with open(COMPLETED_FILE, "a") as completed:
                        completed.write(f"{accession}\n")
                    print(f"Processing completed for {accession}.")
                    break
            else:
                print(f"Assembly failed for {accession}. Retrying... ({retry_count})")
                retry_count += 1

        if retry_count > MAX_RETRIES:
            with open(FAILED_FILE, "a") as failed:
                failed.write(f"{accession}\n")
            print(f"Skipping {accession}: maximum retries exceeded.")

