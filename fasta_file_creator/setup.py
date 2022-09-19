from setuptools import setup,find_packages

setup(
    name="pw2_fasta_creator",
    version="09.15.2022",
    description="create fasta file for pw2",
    packages=find_packages(),
    author="Katsuya Lex Colon",
    author_email="kcolon@caltech.edu",
    py_modules=['pw2_fasta_creator'],
    entry_points={'console_scripts': ['create_pw2_fasta = pw2_fasta_creator.pw2_fasta_creator:main']}
)