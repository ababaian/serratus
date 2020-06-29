## serraplace 
Phylogenetic Placement of high-quality serratus assemblies

Building the container:
`docker build --no-cache -t serraplace .`

Though if subsequently you only change the place.sh script, its OK to re-build with cache:
`docker build -t serraplace .`

Then go to a desired output directory on your local machine (ensure that the contig files are present there under contigs/)
mkdir -p workdir && cd workdir

and run the pipeline (this immediately starts place.sh):
`docker run -v $PWD:$PWD -w $PWD serraplace`

