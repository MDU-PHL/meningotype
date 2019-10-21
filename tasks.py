'''
Automate the build of the container and validation.
This uses invoke

To run:
inv build_container

Find out what tasks are available:
inv --list
'''

from invoke import task
import subprocess, datetime

from mdu_writer.verifications.meningotype_write import WriteMenigotypeVerify

@task
def update_singularity(c):
    date = datetime.datetime.today().strftime("%d_%m_%y")
    c.run("python3 update_meningotype.py")
    c.run("git add *")
    c.run(f"git commit -m updated singulairty {date}")
    c.run("git push")

@task
def build_container(c):
    config = toml.load("config.toml")
    print("Building container")
    c.sudo(f'singularity build salmonella_typing.simg Singularity')
    c.sudo(f"mv {config['container_dir']}/salmonella_typing.simg {config['archive_dir']}") #in preparation for a directory of contianers in the config file
    c.sudo(f"mv salmonella_typing.simg {config['container_dir']}") #in preparation for a directory of contianers in the config file

@task
def run_verification(c):
    config = toml.load("config.toml")
    c.run(f"python3 meningotype/meningotype.py --verify --verification_path {config['verification_path']} --verification_type {config['verification_type']}")

