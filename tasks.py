'''
Automate the build of the container and validation.
This uses invoke

To run:
inv build_container

Find out what tasks are available:
inv --list
'''

from invoke import task
import subprocess, datetime, yaml, toml

# from mdu_writer.verifications.meningotype_write import WriteMenigotypeVerify

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
    date = datetime.datetime.today().strftime("%d_%m_%y")
    archive_dir = f"{config['archive_dir']}/{config['version']}_DB{date}"
    print("Building container")
    c.sudo(f'singularity build salmonella_typing.simg Singularity')
    c.sudo(f"mkdir {archive_dir}") # make archive directory for this image
    c.sudo(f"cp {config['container_dir']}/salmonella_typing.simg {archive_dir}") #copy to the archive directory
    c.sudo(f"mv salmonella_typing.simg {config['container_dir']}") #move to the execution directory

@task
def run_verification(c):
    config = toml.load("config.toml")
    if 'verification_path' in config and 'verification_type' in config:
        c.run(f"python3 meningotype/meningotype.py --verify --verification_path {config['verification_path']} --verification_type {config['verification_type']}")
    else:
        print("Please provide a path for verification data and a verification type.")
@task
def write_verification(c, config_path):
    write = WriteMenigotypeVerify(config_path)
    write.write_doc()

@task
def push_verification(c):
    config = toml.load("config.toml")
    cfg_path = config['verification_path']
    c.run(f"cd {cfg_path}/meningotype.wiki && git add Validation-Report.md && git status && git commit -m 'verification' && git push")

