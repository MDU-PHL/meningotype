from setuptools import setup


def readme():
	with open('README.md') as f:
		return f.read()

def current_version():
	with open('meningotype/VERSION.txt') as f:
		return f.read().rstrip()

setup(name='meningotype',
	version=current_version(),
	description='In silico serotyping and finetyping (porA and fetA) of Neisseria meningitidis',
	long_description=readme(),
	classifiers=[
		'Development Status :: 4 - Beta',
		'License :: OSI Approved :: GPLv3',
		'Programming Language :: Python :: 2.7',
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
	console_scripts=[
		'meningotype/nmen.py',
		'meningotype/ctrA.py',
		'meningotype/menwy.py',
		'meningotype/finetype.py',
		'meningotype/porB.py',
	],
	entry_points={
		'console_scripts': ['meningotype=meningotype.meningotype:main'],
	},
	include_package_data=True,
	zip_safe=False)
