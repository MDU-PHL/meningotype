'''
Given a new mlst version or DB, update the container
'''

import pathlib
import click
import jinja2
import toml
import pendulum
import subprocess
import shlex


def load_template(name):
    '''
    Return the singularity recipe template as unicode text
    '''
    template = pathlib.Path(name).read_text()
    return template


@click.command()
@click.option("--version", default=None)
@click.option("--mlst_version", default="latest")
@click.option("--author", default=None)
@click.option("-c", "--config", default="config.toml")
def update_meningotype_singularity(version, mlst_version, author, config):
    '''
    Use the config.toml, or override any of the options via the command line
    '''
    # load the params
    config = toml.load(config)
    if version is not None:
        config['version'] = version
    if mlst_version is not None:
        config['mlst_version'] = mlst_version
    if author is not None:
        config['author'] = author
    # load the template
    loader = jinja2.FunctionLoader(load_template)
    env = jinja2.Environment(loader=loader)
    SINGULARITY_RECIPE = env.get_template("_singularity.j2").render(config)
    # create global version
    global_recipe = pathlib.Path("Singularity")
    global_recipe.write_text(SINGULARITY_RECIPE)


if __name__ == "__main__":
    update_meningotype_singularity()
