from os import path
from setuptools import setup
from meningotype import __version__ as version


here = path.abspath(path.dirname(__file__))

with open(path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()


setup(name='meningotype',
	version=version,
	description='In silico serotyping and finetyping (porA and fetA) of Neisseria meningitidis',
	long_description=long_description,
	long_description_content_type="text/markdown",
	classifiers=[
		'Development Status :: 4 - Beta',
		'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
		'Programming Language :: Python :: 3.6',
		'Programming Language :: Python :: 3.7',
		'Programming Language :: Python :: 3.8',
		'Programming Language :: Python :: 3.9',
		'Topic :: Scientific/Engineering :: Bio-Informatics',
		'Topic :: Scientific/Engineering :: Medical Science Apps.',
		'Intended Audience :: Science/Research',
	],
	keywords='microbial genomics Neisseria sequence typing',
	url='https://github.com/MDU-PHL/meningotype',
	author='Jason Kwong',
	author_email='kwongj@gmail.com',
	license='GPLv3',
	packages=['meningotype'],
	install_requires=[
		'argparse',
		'BioPython',
	],
	test_suite='nose.collector',
	tests_require=[],
	entry_points={
		'console_scripts': ['meningotype=meningotype.meningotype:main'],
	},
	include_package_data=True,
	zip_safe=False)
