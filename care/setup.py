import os
import stat
import json
import tarfile
import click
import requests
from contextlib import closing
import run as runjupyter

@click.group()
def cli():
    pass

def download_item(url, local_tgz_filepath):
    '''
    Download TGZ from url and save into local_tgz_filepath
    '''
    # closing required for requests=2.11.2; can omit if we upgrade version
    with closing(requests.get(url, stream=True)) as r:
        r.raise_for_status()
        with open(local_tgz_filepath, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192): 
                if chunk: # filter out keep-alive new chunks
                    f.write(chunk)
    click.echo("Downloaded tgz to {}".format(local_tgz_filepath))

def extract_item(local_tgzfile, final_location):
    '''
    Extract tgz at path local_tgzfile. Must contain exactly one dir named
    the same as the tgz's filename. Extract that dir's contents to the final_location dir.
    '''
    tgz_name = os.path.basename(local_tgzfile)
    expected_tgz_dirname = os.path.splitext(tgz_name)[0] # FILENAME.tgz -> FILENAME

    with tarfile.open(local_tgzfile, "r") as tar:
        # Ensure it contains exactly 1 dir, named as expected
        tar_root_dir = os.path.commonprefix(tar.getnames())
        if tar_root_dir != expected_tgz_dirname:
            click.echo("Error: expected {} to contain exactly one dir named {}".format(url, expected_tgz_dirname))
            raise ValueError
        for name in tar.getnames():
            if "/../" in name:
                click.echo("Error: Found path in tar that refers to parent dir: {}".format(name))
                raise ValueError
        # Move all items out of the root dir (except for the dir itself) and extract directly to dest dir
        members_to_extract = []
        for member in tar.getmembers():
            if member.name != tar_root_dir:
                member.name = member.name.replace(tar_root_dir + "/", "", 1)
                members_to_extract.append(member)
        tar.extractall(path=final_location, members=members_to_extract)
    # And chmod the final_location dir to ugo-r so step0.ipynb doesn't error out
    read_only_mode = (~stat.S_IWUSR) & (~stat.S_IWGRP) & (~stat.S_IWOTH)
    os.chmod(final_location, os.stat(final_location).st_mode & read_only_mode)
    click.echo("Extracted into {}".format(final_location))
        
@cli.command()
def run():
    with open("/app/inputs.json", "r") as f: 
        inputs = json.load(f)
    tgz_dir = inputs["local_tgz_dir"]

    for item, attrs in inputs["inputs"].items():
        click.echo("Checking for {} in {}...".format(item, attrs["dir"]), nl=False)
        tgz_path_if_exists = os.path.join(tgz_dir, attrs["download"])
        if os.access(attrs["dir"], os.R_OK):
            click.echo("Found!")
        elif os.access(tgz_path_if_exists, os.R_OK):
            click.echo("Found TGZ file; extracting! (Found at: {})".format(tgz_path_if_exists))
            extract_item(tgz_path_if_exists, attrs["dir"])
        else:
            url = "/".join([inputs["download_base_url"], attrs["download"]])
            click.echo("Couldn't find {}; downloading from {} and extracting.".format(item, url))
            download_item(url, tgz_path_if_exists)
            extract_item(tgz_path_if_exists, attrs["dir"])

    # Set up output dir
    try:
        os.mkdir(inputs["output_dir"])
        click.echo("Made output dir {}".format(inputs["output_dir"]))
    except OSError as e:
        if e.errno != os.errno.EEXIST:
            raise

    # Then, transfer control to run.py
    runjupyter.main(raw_args=["--rollup", "/work/inputs"])

if __name__ == "__main__":
    cli()
