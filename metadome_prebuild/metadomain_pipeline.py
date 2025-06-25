import os
import sys
import yaml
import argparse
import logging
import subprocess
import json
import tempfile

def parse_arguments():
        parser = argparse.ArgumentParser(description="Complete pipeline to identify protein domains and create annotation files")
        parser.add_argument("--config", required=True, default="config/paths.yaml", help="Path to yaml file containing information about files to download and save.")
        parser.add_argument("--cores", required=True, type=int, default=1, help="Number of available cores")
        parser.add_argument("--working_dir_path",required=True,  help="Path to directory containing the pixi.toml, this script and where analysis will be stored")
        parser.add_argument("--mode",required=True, choices=["variants", "metadome", "both"], default="variants",  help="Calculate file containing variant (variants), input files for MetaDome webserver (MetaDome) or both (both).")
        return parser.parse_args()


def setup_logging():
    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=logging.INFO,
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    

def load_config(path="config/paths.yaml"):
    logging.info(f"Loading config from {path}")
    with open(path) as f:
        return yaml.safe_load(f)


def run_cmd(cmd, log_file=None, shell=False):
    logging.info(f"Running command: {cmd}")
    if log_file:
        with open(log_file, "w") as lf:
            result = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, shell=shell)
    else:
        result = subprocess.run(cmd, shell=shell)
    if result.returncode != 0:
        logging.error(f"Command failed (exit {result.returncode}): {cmd}")
        sys.exit(1)


def step0(cfg, cores, pmf):
    logging.info("Step 0: Downloading raw data...")
    # Links
    prot_gencode_link = cfg['download_links']['gencode']['translations']
    trans_gencode_link = cfg['download_links']['gencode']['transcripts']
    genome_gtf_link = cfg['download_links']['gencode']['gtf']
    genome_link = cfg['download_links']['gencode']['genome']
    refseq_link = cfg['download_links']['gencode']['gencode_refseq']
        
    swissprot_prots_link = cfg['download_links']['uniprot']['swissprot_fasta'] 
    uniprot_isoforms_link = cfg['download_links']['uniprot']['isoforms_fasta']
       
    pfam_hmm_link = cfg['download_links']['pfam']['hmm']
    pfam_dat_link = cfg['download_links']['pfam']['hmm_dat']
    
    # Outputs
    prot_gencode = cfg["raw"]['gencode']['translations']
    trans_gencode = cfg["raw"]['gencode']['transcripts']
    gtf_gencode = cfg["raw"]['gencode']['gtf']
    genome = cfg["raw"]['gencode']['genome']
    refseq = cfg["raw"]['gencode']['gencode_refseq']
    
    swissprot_prots = cfg["raw"]['uniprot']['swissprot_fasta'] 
    uniprot_isoforms = cfg["raw"]['uniprot']['isoforms_fasta']
       
    pfam_hmm_gz = cfg["raw"]['pfam']['hmm_gz']
    pfam_dat_gz = cfg["raw"]['pfam']['hmm_dat_gz']
    
    # Proceed with downloads
    logging.info("Download required files...")

    cmd = f"wget --no-verbose -O {prot_gencode} {prot_gencode_link}"
    run_cmd(cmd, shell=True)

    cmd = f"wget --no-verbose -O {trans_gencode} {trans_gencode_link}"
    run_cmd(cmd, shell=True)

    cmd = f"wget --no-verbose -O {gtf_gencode} {genome_gtf_link}"
    run_cmd(cmd, shell=True)
    
    cmd = f"wget --no-verbose -O {genome} {genome_link}" 
    run_cmd(cmd, shell=True)

    cmd = f"pixi run --manifest-path {pmf} bgzip -d {genome}"
    run_cmd(cmd, shell=True)

    cmd = f"pixi run --manifest-path {pmf} samtools faidx {cfg["raw"]["gencode"]['genome_unc']}"
    run_cmd(cmd, shell=True)

    cmd = f"wget --no-verbose -O {refseq} {refseq_link}"
    run_cmd(cmd, shell=True)

    cmd = f"wget --no-verbose -O {swissprot_prots} {swissprot_prots_link}"
    run_cmd(cmd, shell=True)
    
    cmd = f"wget --no-verbose -O {uniprot_isoforms} {uniprot_isoforms_link}"
    run_cmd(cmd, shell=True)
    
    cmd = f"wget --no-verbose -O {pfam_hmm_gz} {pfam_hmm_link}"
    run_cmd(cmd, shell=True)   
     
    cmd = f"pixi run --manifest-path {pmf} bgzip -d {pfam_hmm_gz}"
    run_cmd(cmd, shell=True)
    
    cmd = f"wget --no-verbose -O {pfam_dat_gz} {pfam_dat_link}"
    run_cmd(cmd, shell=True)
    
    cmd = f"pixi run --manifest-path {pmf} bgzip -d {pfam_dat_gz}"
    run_cmd(cmd, shell=True)
    
    return "Executed step 0"

    
def step1(cfg, cores, pmf):
    logging.info("Step 1: Mapping SwissProt proteins to Gencode translations...")
    raw = cfg["raw"]
    inter = cfg["intermediate"]
    res1 = cfg["results"]["step1_mapping"]
    logs = cfg["logs_dir"]

    logging.info("Combining Uniprot fasta files...")
    cmd = f"zcat {raw['uniprot']['isoforms_fasta']} {raw['uniprot']['swissprot_fasta']} > {inter['uniprot_combined_fasta']}"
    run_cmd(cmd, shell=True)

    logging.info("Extracting human proteins...")
    cmd = f"grep OX=9606 {inter['uniprot_combined_fasta']} | sed 's/^>//;s/ .*//' > {inter['uniprot_ids']}"
    run_cmd(cmd, shell=True)
    cmd = f"pixi run --manifest-path {pmf} seqkit grep -f {inter['uniprot_ids']} {inter['uniprot_combined_fasta']} > {inter['uniprot_human_fasta']}"
    run_cmd(cmd, shell=True)
    
    logging.info("Maing BLAST database with Uniprot human data...")

    cmd = (
        f"pixi run --manifest-path {pmf} makeblastdb -in {inter['uniprot_human_fasta']} "
        f"-parse_seqids -dbtype prot -blastdb_version 5 -out {inter['blast_db_dir']} "
        f"-title uniprot_sprot_canonical_and_varsplic"
    )
    print(cmd)
    run_cmd(cmd, shell=True)

    with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
        tmp_path = tmp_file.name

    run_cmd(f'pixi run --manifest-path {pmf} bgzip -d -c {raw["gencode"]["translations"]} > {tmp_path}', shell=True)

    logging.info("Blasting Gencode sequences to the Uniprot ones...")
    
    cmd = (
        f"pixi run --manifest-path {pmf} blastp -num_threads {cores} "
        f"-db {inter['blast_db_dir']} "
        f"-query {tmp_path} "
        f"-outfmt 10 -evalue 1e-5 "
        f"-out {inter['blast_output']}"
    )
    
        
    print(cmd)
    run_cmd(cmd, shell=True)

    if os.path.exists(tmp_path):
        os.remove(tmp_path)

    logging.info("Extracting perfect matches...")

    cmd = (
        f"awk -F, '$3==100.0 && $7==$9 && $8==$10' "
        f"{inter['blast_output']} > {res1['perfect_matches']}"
    )
    run_cmd(cmd, shell=True)
    
    return "Executed step 1"


def step2(cfg, cores, pmf):
    logging.info("Step 2: Checking concordance and quality of the Uniprot-Gencode matches...")
    cfg2 = cfg["results"]["step2_checks"]
    cfg_raw = cfg["raw"]
    cfg_res1 = cfg["results"]["step1_mapping"]
    scripts = cfg["scripts_dir"]

    cmd = (
        f"pixi run --manifest-path {pmf} python {os.path.join(scripts,'Step2-uniprot-to-gencode-mapping-checks.py')} "
        f"--blast {cfg_res1['perfect_matches']} "
        f"--swissprot {cfg['intermediate']['uniprot_combined_fasta']} "
        f"--gencode_prot {cfg_raw['gencode']['translations']} "
        f"--gencode_trans {cfg_raw['gencode']['transcripts']} "
        f"--n_cores {cores} "
        f"--out_pass {cfg2['pass_csv']} "
        f"--out_fail {cfg2['fail_csv']}"
    )
    run_cmd(cmd, shell=True)
    
    return "Executed step 2"


def step3(cfg, pmf):
    logging.info("Step 3: Extracting MANE IDs and creating annotation file...")
    cfg3 = cfg["results"]["step3_mane_annotation"]
    cfg2 = cfg["results"]["step2_checks"]
    raw = cfg["raw"]
    scripts = cfg["scripts_dir"]

    cmd = (
        f"pixi run --manifest-path {pmf} python {os.path.join(scripts,'Step3-Gencode-Uniprot-matched-MANEmapper-extractsequences.py')} "
        f"--mapper_csv {cfg2['pass_csv']} "
        f"--protein_fasta {raw['gencode']['translations']} "
        f"--transcript_fasta {raw['gencode']['transcripts']} "
        f"--output_csv {cfg3['output_csv']} "
        f"--output_protein_fa {cfg3['matched_protein_fa']} "
    )
    run_cmd(cmd, shell=True)
    
    return "Executed step 3"


def step4(cfg, cores, pmf):
    logging.info("Step 4: Running pfamscan to identify protein domains...")
    raw   = cfg["raw"]["pfam"]
    cfg3  = cfg["results"]["step3_mane_annotation"]
    cfg4  = cfg["results"]["step4_pfam"]
    scripts = cfg["scripts_dir"]
    logs = cfg["logs_dir"]

    pfam_dir = os.path.dirname(raw["hmm"])

    logging.info("Running pfam_scan.pl...")
    if not os.path.exists(cfg4['output_pfamscan_file']): 
        cmd = (
            f"pixi run --manifest-path {pmf} python -u {os.path.join(scripts,'Step4a-pfam_scan.py')} "
            f"--input_pfam {pfam_dir} "
            f"--manifest_path {pmf} "
            f"--n_cores {cores} "
            f"--fasta_to_annotate {cfg3['matched_protein_fa']} "
            f"--output_file {cfg4['output_pfamscan_file']} "
        ) 
        run_cmd(cmd, shell=True)

    logging.info("Parsing pfam_scan.pl results...")

    cmd = (
        f"pixi run --manifest-path {pmf} python -u {os.path.join(scripts,'Step4b-metadomains.py')} "
        f"--pfamscan {cfg4['output_pfamscan_file']} "
        f"--fasta {cfg3['matched_protein_fa']} "
        f"--hmmdatabase {raw['hmm']} "
        f"--n_cores {cores} "
        f"--output_dir {cfg4['annotated_metadomains_dir']}"
    )
    run_cmd(cmd, shell=True)
    
    return "Executed step 4"


def step5(cfg, cores, pmf, mode):
    """Step 5: single-SNV finder, in variant, metadome, or both modes."""
    cfg5 = cfg["results"]["step5_single_snv"]
    raw  = cfg["raw"]["gencode"]
    scripts = cfg["scripts_dir"]

    def run_variant():
        logging.info("Step 5 (variants): Finding genomic coordinates and computing possible SNVs...")
        script = "Step5-single_snv_finder.py"
        cmd = (
            f"pixi run --manifest-path {pmf} python -u {os.path.join(scripts,script)} "
            f"--gtf {raw['gtf']} "
            f"--genome_fa {raw['genome_unc']} "
            f"--n_cores {cores} "
            f"--transcript_csv {cfg['results']['step3_mane_annotation']['output_csv']} "
            f"--output_dir {cfg5['snv_tables_dir']}"
        )
        run_cmd(cmd, shell=True)

    def run_metadome():
        logging.info("Step 5 (MetaDome): Computing genomic, transcript, and protein mappings for MetaDome...")
        script = "MetaDomeStep5-single_snv_finder.py"
        cmd = (
            f"pixi run --manifest-path {pmf} python -u {os.path.join(scripts,script)} "
            f"--gtf {raw['gtf']} "
            f"--genome_fa {raw['genome_unc']} "
            f"--n_cores {cores} "
            f"--transcript_csv {cfg['results']['step3_mane_annotation']['output_csv']} "
            f"--output_dir {cfg5['snv_tables_dir_metadome']}"
        )
        run_cmd(cmd, shell=True)

    if mode in ("variants", "both"):
        run_variant()
    if mode in ("metadome", "both"):
        run_metadome()



def step6(cfg, cores, pmf):
    logging.info("Step 6: Calculating metapositionss...")
    cfg6 = cfg["results"]["step6_metapositions"]
    scripts = cfg["scripts_dir"]

    cmd = (
        f"pixi run --manifest-path {pmf} python {os.path.join(scripts,'Step6-stokholm2tbl.py')} "
        f"--n_cores {cores} "
        f"--input_folder {cfg['results']['step4_pfam']['annotated_metadomains_dir']} "
        f"--output {cfg6['metadomain_positions']}"
    )
    run_cmd(cmd, shell=True)
    
    return "Executed step 6"


def step7(cfg, cores, pmf, mode):
    """Step 7: merge outputs for variants, metadome, or both."""
    scripts = cfg['scripts_dir']
    idmapper = cfg['results']['step3_mane_annotation']['output_csv']
    metaposition = cfg['results']['step6_metapositions']['metadomain_positions']
    gencode_refseq = cfg['raw']['gencode']['gencode_refseq']
    pfam_interpro = cfg['raw']['pfam']['pfam_interpro']
    uniprot_name = cfg['intermediate']['uniprot_ids']
    pfamscan_output= cfg["results"]["step4_pfam"]['output_pfamscan_file']
    metadome_annot = cfg['results']['metadome_data_fields']

    def run_variant():
        script = f"Step7-output_merge_reoder_columns.py"
        cmd = (
            f"pixi run --manifest-path {pmf} python {os.path.join(scripts, script)} "
            f"--idmapper {idmapper} "
            f"--metaposition {metaposition} "
            f"--genomic_folder {cfg['results']['step5_single_snv']['snv_tables_dir']} "
            f"--output {cfg['results']['step7_finaloutput']['final_output']} "
            f"--refseq {gencode_refseq} "
            f"--n_cores {cores}"
        )
        run_cmd(cmd, shell=True)
    
    def run_metadome():
        script = f"MetaDomeStep7-output_merge_reoder_columns.py"
        cmd = (
            f"pixi run --manifest-path {pmf} python {os.path.join(scripts, script)} "
            f"--idmapper {idmapper} "
            f"--metaposition {metaposition} "
            f"--genomic_folder {cfg['results']['step5_single_snv']['snv_tables_dir_metadome']} "
            f"--output {cfg['results']['step7_finaloutput']['final_output_metadome']} "
            f"--uniprot_name {uniprot_name} "
            f"--pfamscan_output {pfamscan_output} "
            f"--n_cores {cores} "
            f"--pfam_interpro {pfam_interpro} "
            f"--genome_build {metadome_annot['genome_build']} "
            f"--source {metadome_annot['source']} "
            f"--GENCODE_version {metadome_annot['GENCODE_version']} "
            f"--PFAM_version {metadome_annot['PFAM_version']} "
            f"--refseq {gencode_refseq}"
        )
        run_cmd(cmd, shell=True)

    if mode in ("variants", "both"):
        logging.info("Step 7 (variants): merging variant SNV outputs...")
        run_variant()
    if mode in ("metadome", "both"):
        logging.info("Step 7 (MetaDome): merging MetaDome outputs...")
        run_metadome()

    return "Executed step 7"


def step8(cfg, pmf):
    logging.info("Step 8: Collecting stats about the metapositions...")
    input_csv = cfg['results']['step7_finaloutput']['final_output']
    scripts = cfg["scripts_dir"]

    cmd = (
        f"pixi run --manifest-path {pmf} python {os.path.join(scripts,'VariantsStep8-stats.py')} "
        f"--input_csv {input_csv}"
    )
    run_cmd(cmd, shell=True)
    
    return "Executed step 8"

def main():

    args = parse_arguments()
    
    # Find the pixi manifest
    pmf = os.path.join(args.working_dir_path, "pixi.toml")
    
    os.chdir(args.working_dir_path)
    
    setup_logging()

    
    # Load config
    cfg = load_config(args.config)
    
    # create all output directories up-front
    for name, path in cfg.get("output_directories", {}).items():
        os.makedirs(path, exist_ok=True)
        logging.debug(f"Created directory for {name}: {path}")

    # Checkpointing setup
    state_file = cfg.get("state_file", ".pipeline_state.json")
    if os.path.exists(state_file):
        with open(state_file) as sf:
            data = json.load(sf)
            completed = set(data.get("completed", []))
    else:
        completed = set()

    # Validate that previously completed steps actually have outputs
    expected_outputs = {
        "step0": [
            cfg["raw"]["gencode"]["translations"],
            cfg["raw"]["gencode"]["transcripts"],
            cfg["raw"]["gencode"]["gtf"],
            cfg["raw"]["uniprot"]["swissprot_fasta"],
            cfg["raw"]["uniprot"]["isoforms_fasta"],
            cfg["raw"]["pfam"]["hmm"],
            cfg["raw"]["pfam"]["hmm_dat"],
        ],
        "step1": [
            cfg["intermediate"]["uniprot_combined_fasta"],
            cfg["intermediate"]["uniprot_human_fasta"],
            cfg["intermediate"]["uniprot_ids"],
            cfg["intermediate"]["blast_output"],
        ],
        "step2": [
            cfg["results"]["step2_checks"]["pass_csv"],
            cfg["results"]["step2_checks"]["fail_csv"],
        ],
        "step3": [
            cfg["results"]["step3_mane_annotation"]["output_csv"],
            cfg["results"]["step3_mane_annotation"]["matched_protein_fa"],
        ],
        "step4": [cfg["results"]["step4_pfam"]["output_pfamscan_file"]],
        "step6": [cfg["results"]["step6_metapositions"]["metadomain_positions"]]
    }

    # Validate outputs from the two modes
    expected_outputs.setdefault("step5", [])
    expected_outputs.setdefault("step7", [])

    if args.mode in ("variants", "both"):
        expected_outputs["step5"].append(
            cfg["results"]["step5_single_snv"]["snv_tables_dir"]
        )
        expected_outputs["step7"].append(
            cfg["results"]["step7_finaloutput"]["final_output"]
        )
    
    if args.mode in ("metadome", "both"):
        expected_outputs["step5"].append(
            cfg["results"]["step5_single_snv"]["snv_tables_dir_metadome"]
        )
        expected_outputs["step7"].append(
            cfg["results"]["step7_finaloutput"]["final_output_metadome"]
        )
        
    for step in list(completed):
        missing = [p for p in expected_outputs.get(step, []) if not os.path.exists(p)]
        if missing:
            logging.warning(f"Outputs for {step} missing ({missing}); will run that step")
            completed.remove(step)

    # Pipeline steps
    steps = [
        ("step0", "Download & prepare inputs",            lambda: step0(cfg, args.cores, pmf)),
        ("step1", "SwissProt id mapping & filtering",     lambda: step1(cfg, args.cores, pmf)),
        ("step2", "Run metadomains script",               lambda: step2(cfg, args.cores, pmf)),
        ("step3", "Annotate with MANE & subset FASTA",    lambda: step3(cfg, pmf)),
        ("step4", "PFAM scan & parse",                    lambda: step4(cfg, args.cores, pmf)),
        ("step5", "Mapper calculator",                    lambda: step5(cfg, args.cores, pmf, args.mode)),
        ("step6", "Convert Stockholm alignments",         lambda: step6(cfg, args.cores, pmf)),
        ("step7", "Merge results",                        lambda: step7(cfg, args.cores, pmf, args.mode)),
    ]
    
    if args.mode in ("variants", "both"):
        steps.append(("step8", "Collect stats",         lambda: step8(cfg, pmf)))

    print(steps)
    
    # Run what needed
    for name, desc, func in steps:
        if name in completed:
            logging.info(f"=== Skipping {desc} ({name}) — already completed ===")
            continue

        logging.info(f"=== Step {name[-1]}: {desc} ===")
        func()

        # record completion
        completed.add(name)
        with open(state_file, "w") as sf:
            json.dump({"completed": list(completed)}, sf)

    logging.info("Finished!")



if __name__ == "__main__":
    main()
