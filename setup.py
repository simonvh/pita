from distutils.core import setup
import setuptools

VERSION = "1.10"
DESCRIPTION = """
pita - pita improves transcript annotation
"""

setup (name = 'pita',
		version = VERSION,
		description = DESCRIPTION,
		author='Simon van Heeringen',
		author_email='s.vanheeringen@ncmls.ru.nl',
		license='MIT',
		packages=[
			'pita'
		],
		scripts=[
			"scripts/pita",
		],
		data_files=[]
)
