import setuptools

setuptools.setup(
    name="serratus_summarizer",
    version="0.0.1",
    author="The Serratus Project",
    author_email="",
    description="Tools to make summaries of alignments and reference sequences.",
    pacakges=['serratus_summarizer'],
    install_requires=[
        'biopython',
        ],
    classifiers = [
        "Programming Language :: Python :: 2",
    ]
    )
