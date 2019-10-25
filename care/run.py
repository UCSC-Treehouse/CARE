import os
import csv
import glob
import json
import shutil
import nbformat
import argparse
import subprocess
from distutils.version import LooseVersion
from nbconvert.preprocessors import ExecutePreprocessor
from nbconvert.preprocessors.execute import CellExecutionError
import warnings
import sys
import traceback





def run_notebook(notebook):
    '''
        Run a Jupyter notebook and save it in-place.
        The notebook is run with a cwd of its parent directory.
        If the notebook errors out, the error output is captured in the saved notebook.
        Takes:
            notebook: relative or absolute path to an .ipynb file that will be overwritten
            with its executed .ipynb version.
        Returns:
            True if the notebook executed fully;
            False if there was a CellExecutionError or a "Kernel died" RuntimeError
    '''
    with open(notebook) as f:
        nb = nbformat.read(f, as_version=4)
    ep = ExecutePreprocessor(timeout=4000, kernel_name='python2')
    ran_ok=False

    try:
        ep.preprocess(nb, {'metadata': {'path': os.path.dirname(notebook)}})
        ran_ok=True
    except CellExecutionError as e:
        # Any deliberately thrown KeyboardInterrupt, etc will be caught here
        print("Error running notebook.")
        parsed_error = str(e).split("(PRINTTHIS)")
        if len(parsed_error) > 0:
            print(parsed_error[1])
        print("Please check the notebook file for more details.")
        print(notebook)
    except RuntimeError as e:
        if str(e) == "Kernel died before replying to kernel_info":
            print( "Kernel died while running notebook! Skipping this sample." ) # TODO retry instead
        else:
            raise
    finally:
        with open(notebook, 'wt') as f:
            nbformat.write(nb, f)
        return ran_ok

def copy_protocol(from_dir, to_dir, copy_beta=False):
    '''
        Copy the appropriate files into the tertiary output dir. This includes all notebooks,
        the treehouse library, annotations.json, and the templates file.
    '''
    shutil.copytree(os.path.join(from_dir, "treehouse"), os.path.join(to_dir, "treehouse"),
                    ignore=shutil.ignore_patterns('*.pyc')) # will create the parent dir as well
    shutil.copy(os.path.join(from_dir, "summary.template"), to_dir)
    shutil.copy(os.path.join(from_dir, "slides.template"), to_dir)
    shutil.copy(os.path.join(from_dir, "annotations.json"), to_dir)
    for notebook in glob.glob(os.path.join(from_dir, "*.ipynb")):
        shutil.copy(notebook, to_dir)
    if copy_beta: # Also copy the beta directory
        shutil.copytree(os.path.join(from_dir, "beta"), os.path.join(to_dir, "beta"),
                        ignore=shutil.ignore_patterns('*.pyc')) # will create the parent dir as well



# Locate the latest rsem_genes.results file for a sample
def locate_expression(samplepath):
    all_rnaseq = sorted(glob.glob(os.path.join(samplepath, "secondary", "ucsc_cgl-rnaseq-cgl-pipeline-*")),
                        key=LooseVersion)
    # Get the file, or most recent if there are more than one.
    if len(all_rnaseq) >= 1:
        result = os.path.join(all_rnaseq[-1], "RSEM", "rsem_genes.results")
    else:
        result = ""
    return result

# Parse a provided manifest and return an array of dicts:
# { sample, disease, rollupfile, alias }
# missing values might be either None or ''
def parse_manifest(path):
    keys = ["sample", "disease", "rollupfile", "alias"]
    try:
        with open(path, "r") as f:
            rows = csv.DictReader(f, fieldnames=keys, dialect='excel-tab')
            result = []
            for row in rows:
                result.append(row)
    except IOError:
        print("  Missing or invalid manifest!")
        result = []
    return result

# Create the sample_info.json for a specific sample
# Should have keys: sample disease alias
def create_auxillary_files(samples_details, sample, output_dir):
    # Get the info for this sample
    created_ok = True
    info = [ x for x in samples_details if x["sample"] == sample][0]
    info["rollup"] = None # Will be populated below if rollupfile exists

    if info["disease"]:
        print("  Found disease '{}' for sample {}.".format(info["disease"], sample))
    if info["alias"]:
        print("  Found compendium alias '{}' for sample {}".format(info["alias"], sample))

    if info["rollupfile"]:
        print("  Loading roll-up cohort file '{}' for sample {}".format(info["rollupfile"], sample))
        try:
            with open(os.path.join("/work/rollup", info["rollupfile"]), "r") as f:
                info["rollup"] = list(map(lambda x: x.strip(), f.readlines()))
        except IOError:
            print("  Couldn't find provided roll-up cohort file! Skipping this sample.")
            created_ok = False
    with open(os.path.join(output_dir, "sample_info.json"), "w") as f:
        json.dump(info, f)
    return created_ok

# If a QC file is present (may be in json or csv format), create it as bam_umend_qc.json in the output dir.
def create_qc_json_file(input_dir, output_dir):
    json_filename = "bam_umend_qc.json"
    found = False 
    # Check for a pre-existing json file 
    new_json_path = os.path.join(input_dir, "secondary", "ucsctreehouse-bam-umend-qc-*",json_filename)
    found_json = sorted(glob.glob(new_json_path), key=LooseVersion)
    if found_json:
        shutil.copy(found_json[-1], output_dir)
        print("  Found secondary QC JSON file.")
        found=True
    else:
        # json not found; look for old format file & convert to json
        old_txt_path= os.path.join(input_dir, "secondary", "ucsc_cgl-rnaseq-cgl-pipeline-*",
                                     "QC", "bamQC", "readDist.txt_*_qc.txt")
        found_txt = sorted(glob.glob(old_txt_path), key=LooseVersion)
        if found_txt:
            try:
                with open(found_txt[-1], "r") as f:
                    print("  Found secondary QC txt file; parsing to JSON.")
                    parsed = next(csv.DictReader(f, dialect='excel-tab'))
                with open(os.path.join(output_dir, json_filename), "w") as f:
                    json.dump(parsed, f)
                found= True 
            except IOError:
                print("  Error trying to parse secondary QC text file; skipping this sample.")
        else:
            print("  No secondary MEND QC file found; skipping this sample.")
    return found
      
def main(raw_args=None):
    
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument("-c", "--clean", action='store_true', default=False,
                        help="Start with a clean copy of the protocol")
    parser.add_argument("-j", "--jupyter", action='store_true', default=False,
                        help="Start jupyter notebook server")
    parser.add_argument("--start", default=False,
                        help="Start the analysis at this notebook number")
    parser.add_argument("--end", "--stop", default=float("inf"),
                        help="After this notebook is run, stop the analysis")
    parser.add_argument("--beta", action="store_true", default=False,
                        help="Also run beta notebooks after the protocol notebooks.")
    dirs_to_link = ["inputs", "outputs", "references", "cohort", "rollup"]
    for linkdir in dirs_to_link:
        parser.add_argument("--{}".format(linkdir), help="Path to {} directory, as mounted".format(linkdir))
    args = parser.parse_args(raw_args)

    print( "Running with arguments: {}".format(args), flush=True)

    if args.jupyter:
        # Set Jovyan's GID to match passed in NB_GID if present
        # (Base jupyter/datascience-notebook:72cca2a7f3ea is too old to support this)
        if os.environ["NB_GID"]:
            gid_bytes = subprocess.check_output(["id", "-g", os.environ["NB_USER"]])
            current_gid = "".join(map(chr, gid_bytes)).strip()
            if current_gid != os.environ["NB_GID"]:
                print("Setting gid to {}".format(os.environ["NB_GID"]))
                current_groupname = "".join(map(chr, subprocess.check_output(["id", "-g", "-n", os.environ["NB_USER"]]))).strip()
                subprocess.check_call(["groupmod", "-g", os.environ["NB_GID"], "-o", current_groupname])
        
        print("Starting Jupyter notebook server...")
        try:
            proc = subprocess.Popen(["start-notebook.sh"], shell=False)
            proc.communicate()
        finally:
            proc.terminate()
            proc.kill()
        return
    
    else:

        # If desired dirs are not present, make symlinks in the work directory
        # for input, output, reference, and cohort. Reason: Docker may not have been able to
        # mount these dirs at the desired spot
        basedir = "/work"
        for linkdir in dirs_to_link:
            try:
                source = getattr(args, linkdir)
                destination = os.path.join(basedir, linkdir)
                file_exists_message = "{} OK (found existing.)".format(destination)
                os.symlink(source, destination)
            except FileExistsError:
                print(file_exists_message)
                # Continue on - its ok to not provide the parameter if the dest already is a file
            except TypeError:
                if os.path.isdir(destination):
                    print(file_exists_message)
                else:
                    print("Please provide either --{} parameter or create the {} dir yourself. Can't find this dir, so exiting.".format(
                        linkdir, destination))
                    exit(1)
            try:
                os.stat(destination)
            except FileNotFoundError:
                print("Can't find {} as the source of the {} dir. Please provide a valid source. Exiting.".format(source, destination))
                exit(1)
        
        # Check for manifest
        manifest_file = "/app/manifest.tsv"
        if os.path.isfile(manifest_file):
            print("Analyzing all samples listed in the manifest...")
            # Load that up
            samples_details = parse_manifest(manifest_file)
            all_samples = map(lambda x: x["sample"], samples_details)
        else:
            print("Couldn't find file manifest.tsv! Please provide a manifest listing your sample IDs.")
            all_samples = []

        for sample_id in all_samples:
            if args.beta:
                print("*** BETA MODE ***\nWill run beta notebooks!")

            input_dir = os.path.join("/work/inputs", sample_id)
            output_dir = os.path.join("/work/outputs", sample_id)

            input_file = locate_expression(input_dir)

            tumormap_alt_input_file = os.path.join(input_dir, "alternate_expression_for_tumormap.tsv")

            if not os.path.isfile(input_file):
                print("{}: No rsem_genes.results file found; skipping.".format(sample_id))
                continue

            # Set up the output dir
            if args.clean:
                if args.start:
                    print("{}: Refusing to start from a later notebook after cleaning the output dir.".format(sample_id))
                    continue
                print("{}: Overwriting existing protocol copy".format(sample_id))
                shutil.rmtree(output_dir, ignore_errors=True)
            if os.path.isdir(output_dir):
                if args.start:
                    print("{}: Running with existing output".format(sample_id))
                else:
                    print("{}: Refusing to overwrite existing output; try --clean.".format(sample_id))
                    continue
            else:
                print("{}: Copying protocol & inputs".format(sample_id))
                copy_protocol("/app", output_dir, args.beta) # if beta is true, copy beta dir also

                # rsem_genes.results
                shutil.copy(input_file, output_dir)

                # Further input files not sourced directly from the protocol dir

                # sample_info.json - roll_up file must be present if requested
                if not create_auxillary_files(samples_details, sample_id, output_dir):
                    continue
                # bam_umend_qc.json - must be present
                if not create_qc_json_file(input_dir, output_dir):
                    continue

                if os.path.isfile(tumormap_alt_input_file):
                    print("  Found alternate_expression_for_tumormap.tsv file to use for Tumormap placement.")
                    shutil.copy(tumormap_alt_input_file, output_dir)
                    
            # Environment variables
            os.environ["TREEHOUSE_SAMPLE_ID"] = sample_id

            # Get the numerical prefix for each notebook and run in that order
            all_notebooks = glob.glob(os.path.join(output_dir, "*.ipynb"))
            notebooks_by_num = {float(os.path.basename(nb).split("_")[0]): nb for nb in all_notebooks}
            print("  Running all notebooks in numerical order")
            all_ran = False
            for notebook_num in sorted(notebooks_by_num.keys()):
                if notebook_num < float(args.start):
                    print("    Notebook {} is before {}; skipping.".format(notebook_num, args.start))
                    continue
                if notebook_num > float(args.end):
                    print("    Notebook {} has completed; halting execution.".format(args.end))
                    break
                print("    {}".format(notebooks_by_num[notebook_num]), flush=True)
                if not run_notebook(notebooks_by_num[notebook_num]):
                    print("    Notebook failed: terminating this sample.")
                    break
                    # TODO - mark this run as tertiary-failed --- somehow
                all_ran = True
            if args.beta and all_ran: # Run beta when regular notebooks have completed succesfully.
                beta_notebooks = glob.glob(os.path.join(output_dir, "beta", "*.ipynb"))
                print("  Running beta notebooks in arbitrary order!")
                for notebook in beta_notebooks:
                    print("    {}".format(notebook))
                    if not run_notebook(notebook):
                        print("    Warning: Beta notebook failed: {}".format(os.path.basename(notebook)))
                        # Run the next one anyhow!
    print("Done!")


if __name__ == '__main__':
    main()
