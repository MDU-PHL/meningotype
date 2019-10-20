'''
Automate the build of the container and validation.
This uses invoke

To run:
inv build_container

Find out what tasks are available:
inv --list
'''

from invoke import task

@task
def clean(c):
    pass

@task(clean)
def build_container(c):
    pass

@task(build_container)
def run_validation(c):
    pass