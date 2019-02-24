# globus
GLObal Biochemical reconstruction Using Sampling

Source files to build Docker container using diamond as the sequence aligner.

GLOBUS can be run in a docker container (https://www.docker.com/get-started)

HOW TO ANNOTATE A GENOME WITH GLOBUS:

After docker is installed and running, pull the GLOBUS image by typing on a terminal window:

$ docker pull platyias/globus:diamond

Once the image is downloaded, create a container by typing:

$ docker run -dit platyias/globus:diamond

You are now ready to annotate your protein sequences with GLOBUS.

First, copy a multi-fasta file (e.g. sequences.fas) to the genomes folder in the docker container:

$ docker cp sequences.fas :/data/genomes/904306.fas

Replace with the corresponding ID of the container you just created; this can be found by typing:

$ docker ps

Second, run GLOBUS by typing

$ docker exec /bin/bash -c "/data/RunGlobus.sh sequences.fas"

Third, after GLOBUS finishes, copy GLOBUS results from the docker container:

$ docker cp :/data/results/sequences.fas/sequences.fas.All_P_sorted+Iden.txt globus_out.txt

Files inside the docker container can also be explored by typing:

$ docker exec -i -t /bin/bash
