## Serratax

First try at implementing Serratax in a container.

To build the container, copy the Dockerfile to an empty directory and type:  

`docker build --no-cache -t serratax .`

To run the container interactively:  

`docker run -it -v $PWD:$PWD -w $PWD serratax bash`

The rather obscure -v and -w options cause your current directory to be mounted 
inside the container and have the same pathname as it does on the host machine.
At the bash prompt inside the container, you can run:

`serratax input.fasta ./outputdir`

The input.fasta file must be inside the current directory or below. If below,
give the relative path starting with dot slash `./`.
The output will be written to the given output directory, which also must be
in or under your current directory. Exit the container by typing `exit` at
the bash prompt, and you will find outputdir under your current directory
on the host.

No doubt you can do things differently
by adjusting the -v and -w options, but that's beyond my skillset at the time
of writing.   

To run as a stand-alone command without an interactive prompt:  

`docker run -v $PWD:$PWD -w $PWD serratax serratax input.fasta ./outputdir`

Notice that `serratax` appears twice in the above command line. First time, it's the container name.
Second time, it's the command to execute within the container.  

The final prediction is written to `outputdir/tax.final`. It is a tsv file with three fields: 
1. NCBI taxid, 2. NCBI taxonomy name, and 3. Full taxonomy path in my own format (there is no standard for this AFAIK).
Example:  

`694009  Severe acute respiratory syndrome-related coronavirus   family:Coronaviridae,subfamily:Orthocoronavirinae,genus:Betacoronavirus,subgenus:Sarbecovirus,species:Severe acute respiratory syndrome-related coronavirus`

This is a novice attempt at making a container -- it's my very first! The Dockerfile is very short. Dependencies
(python scripts, muscle and usearch binaries) are copied into the container from a tarfile on S3 which I prepared manually.
Better practice might be to get python code from the repo, but for now this is simpler 
and it works. Suggestions for improvements welcome, email me or open an issue in the repo. `robert@drive5.com`
