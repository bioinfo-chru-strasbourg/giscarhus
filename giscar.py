from datetime import datetime
from hashlib import md5
from multiprocessing import Pool
from os import makedirs, mkdir, stat, umask
from pathlib import Path
from shutil import copy2, rmtree
from subprocess import run
from time import time


def main(
    max_run_age: int,
    repository: Path,
    tmp_folder: Path,
    singularity_image: Path,
    r_script: Path,
    r_model: Path,
    dataset_validation: Path,
    genome_reference: Path,
    threads: int,
    ignored_samples: list[str],
    manual_purity_bool: bool,
    manual_purity: dict[str, str],
):
    date = datetime.now().strftime("%Y%m%d") + "-" + datetime.now().strftime("%H%M%S")
    if str(dataset_validation) == ".":
        for folder in repository.iterdir():
            if (
                (
                    time() - stat(repository / folder).st_mtime
                    < max_run_age * 24 * 60 * 60
                )
                and (repository / folder / "STARKCopyComplete.txt").is_file()
                and not (repository / folder / "GISCARCcomplete.txt").is_file()
            ):
                folder_path = repository / folder
                launching_process(
                    folder_path,
                    date,
                    tmp_folder,
                    singularity_image,
                    r_script,
                    r_model,
                    genome_reference,
                    threads,
                    ignored_samples,
                    manual_purity_bool,
                    manual_purity,
                )
    else:
        launching_process(
            dataset_validation,
            date,
            tmp_folder,
            singularity_image,
            r_script,
            r_model,
            genome_reference,
            threads,
            ignored_samples,
            manual_purity_bool,
            manual_purity,
        )


def launching_process(
    input: Path,
    date: str,
    tmp_folder: Path,
    singularity_image: Path,
    r_script: Path,
    r_model: Path,
    genome_reference: Path,
    threads: int,
    ignored_samples: list[str],
    manual_purity_bool: bool,
    manual_purity: dict[str, str],
):
    folder = input.name
    repository = input.parent
        

    samples = [
        sample.name
        for sample in (repository / folder).iterdir()
        if sample.is_dir() and sample.name not in ignored_samples
    ]
    samples_fastq = {
        sample: {
            "R1": next((repository / folder).glob(f"{sample}_*_R1_001.fastq.gz")),
            "R2": next((repository / folder).glob(f"{sample}_*_R2_001.fastq.gz")),
        }
        for sample in samples
    }

    if (tmp_folder / "giscar").is_dir():
        rmtree(tmp_folder / "giscar")

    for sample in samples:
        if not (tmp_folder / "giscar" / folder / sample).is_dir():
            makedirs(tmp_folder / "giscar" / folder / sample, mode=0o755)
    tmp_folder = tmp_folder / "giscar" / folder

    samplesheet = find_samplesheet(repository, folder)
    purity_dic = get_purity(samplesheet, manual_purity_bool, manual_purity)

    with Pool(threads) as pool:
        args = [
            (fastq_file, tmp_folder / sample)
            for sample, files in samples_fastq.items()
            for fastq_file in files.values()
        ]
        pool.starmap(worker_multi_md5_copy, args)

    with Pool(threads) as pool:
        args = [
            (
                sample,
                str(next((tmp_folder / sample).glob(f"{sample}_*_R1_001.fastq.gz"))),
                str(next((tmp_folder / sample).glob(f"{sample}_*_R2_001.fastq.gz"))),
                tmp_folder,
                purity_dic,
                singularity_image,
                r_script,
                r_model,
                genome_reference,
            )
            for sample in samples
        ]
        pool.starmap(worker_launch_giscar, args)
    
    hrd_output = list(Path(tmp_folder).glob("*/HRD/R_output/HRD_uniqueData.tsv"))
    merged_output = output_merger(hrd_output, date, folder, tmp_folder)
    copy2(merged_output, repository / folder)

    for sample in samples:
        sample_tmp_folder = tmp_folder / sample
        output_file = sample_tmp_folder / "HRD" / "R_output" / "HRD_uniqueData.tsv"

        makedirs(repository / folder / sample / "GISCAR", exist_ok=True)
        copy2(output_file, repository / folder / sample / "GISCAR" / f"GISCAR_HRD.{sample}.{date}.tsv")

    rmtree(tmp_folder / "giscar")
    with open(repository / folder / "GISCARComplete.txt", "w"):
        pass


def get_purity(samplesheet: Path, manual_purity_bool: bool, manual_purity: dict[str, str]) -> dict:
    if manual_purity_bool is False:
        purity_dic = {}
        with open(samplesheet, "r") as read_file:
            purity_dic = {
                line.split(",")[0]: value
                for line in read_file
                if "TUM#" in line
                for cell in line.strip().split(",")
                if "TUM#" in cell
                for purity in cell.split("!")
                if "TUM#" in purity
                for value in purity.split("#")
                if value.isdigit()
            }
        return purity_dic
    else:
        return manual_purity


def output_merger(hrd_output: list[Path], date: str, folder: str, tmp_folder: Path) -> Path:
    input_files = [str(file) for file in hrd_output]
    output_file = tmp_folder / f"GISCAR_HRD.{folder}.{date}.tsv"
    cmd = ["awk", "FNR==1 && NR!=1 { next; } { print }"] + input_files
    print(" ".join(cmd))
    with open(output_file, "w") as f:
        run(cmd, stdout=f, check=True)
    return output_file


def find_samplesheet(repository: Path, folder: str) -> Path:
    samples = [sample for sample in (repository / folder).iterdir() if sample.is_dir()]
    samplesheet = samples[0] / "STARK" / f"{samples[0].name}.SampleSheet.csv"
    return samplesheet


def worker_multi_md5_copy(source: Path, destination: Path) -> None:
    makedirs(destination, exist_ok=True)
    with open(source, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            md5().update(chunk)
    source_md5 = md5().hexdigest()

    copy2(source, destination)
    copied_file = destination / source.name

    with open(copied_file, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            md5().update(chunk)
    destination_md5 = md5().hexdigest()

    if source_md5 != destination_md5:
        print(
            f"MD5 checksum mismatch for {source.name} after copy to {destination}, please check the file"
        )
        raise ValueError(
            f"MD5 checksum mismatch for {source.name} after copy to {destination}, please check the file"
        )


def worker_launch_giscar(
    sample: str,
    r1_fastq: Path,
    r2_fastq: Path,
    tmp_folder: Path,
    purity_dic: dict[str, str] | dict,
    singularity_image: Path,
    r_script: Path,
    r_model: Path,
    genome_reference: Path,
):
    purity = purity_dic.get(sample)
    if purity is None:
        purity = "60"
    host_tmp = Path("/home1/data/tmp/")
    sample_tmp_folder = tmp_folder / sample

    giscar_cmd = [
        "singularity",
        "run",
        "--bind",
        f"{str(host_tmp)}:/tmp/,{str(tmp_folder)}:{str(tmp_folder)},{str(genome_reference.parent)}:{str(genome_reference.parent)}",
        str(singularity_image),
        "--purity",
        str(float(purity)/100),
        "-o",
        str(sample_tmp_folder),
        "--fastq1",
        str(r1_fastq),
        "--fastq2",
        str(r2_fastq),
        "--genome-reference",
        str(genome_reference),
    ]
    print(" ".join(giscar_cmd))
    run(giscar_cmd, check=True)

    r_input = sample_tmp_folder / "HRD" / "Call"
    r_output = sample_tmp_folder / "HRD" / "R_output"
    makedirs(r_output, mode=0o755)
    r_cmd = [
        "Rscript",
        str(r_script),
        "-C",
        str(r_input),
        "--model",
        str(r_model),
        "-O",
        str(r_output),
    ]
    print(" ".join(r_cmd))
    run(r_cmd, check=True)


if __name__ == "__main__":
    max_run_age: int = 30  # how far back to look for runs
    repository = Path("/home1/L_PROD/HUSDIAGGEN/BRCANESS/")
    tmp_folder = Path("/home1/data/tmp/")
    singularity_image = Path(
        "/home1/L_PROD/NGS/BAS/rauchm/giscar_script/HRD-v4_demo.simg"
    )
    r_script = Path("/home1/L_PROD/NGS/BAS/rauchm/giscar_script/pipelineHRDscorev4.r")
    r_model = Path("/home1/L_PROD/NGS/BAS/rauchm/giscar_script/modFinalv4.RData")
    dataset_validation = Path("/home1/L_PROD/HUSDIAGGEN/BRCANESS/250327_NB551590_0607_AHMGL2AFX7")
    genome_reference = Path("/home1/DB/STARK/genomes/current/hg19.fa")
    original_umask = umask(0o000)
    threads: int = 30
    ignored_samples = ["NTC"]
    manual_purity_bool = False
    manual_purity = {
        "":"",
        "":"",
        "":"",
        "":"",
    }
    main(
        max_run_age,
        repository,
        tmp_folder,
        singularity_image,
        r_script,
        r_model,
        dataset_validation,
        genome_reference,
        threads,
        ignored_samples,
        manual_purity_bool,
        manual_purity,
    )
    umask(original_umask)
