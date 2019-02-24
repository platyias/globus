# globus
GLObal Biochemical reconstruction Using Sampling

This repository contains source files to build a Docker image using usearch as the sequence aligner. Due to license restrictions, usearch needs to be manually added to a local docker image by the user.   

If you use GLOBUS, please cite:

Plata, G., Fuhrer, T., Hsiao, T., Sauer, U., Vitkup, D. 2012. Global probabilistic annotation of metabolic networks enables enzyme discovery. Nat. Chem. Biol. 8:848-854

Code in this repository was contributed by Dr. Jeewoen Sin. 

HOW TO ANNOTATE A GENOME WITH GLOBUS:

GLOBUS can be run in a docker container (https://www.docker.com/get-started)

After docker is installed and running, pull the GLOBUS image by typing on a terminal window:

$ docker pull platyias/globus:diamond

Once the image is downloaded, create a container by typing:

$ docker run -dit platyias/globus:diamond

You are now ready to annotate your protein sequences with GLOBUS.

First, copy a multi-fasta file (e.g. sequences.fas) to the genomes folder in the docker container:

$ docker cp sequences.fas <container id>:/data/genomes/904306.fas
  
Replace <container id> with the corresponding ID of the container you just created; this can be found by typing:
  
$ docker ps

Second, run GLOBUS by typing 

$ docker exec <container id> /bin/bash -c "/data/RunGlobus.sh sequences.fas"
  
Third, after GLOBUS finishes, copy GLOBUS results from the docker container:

$ docker cp <container id>:/data/results/sequences.fas/sequences.fas.All_P_sorted+Iden.txt globus_out.txt

Files inside the docker container can also be explored by typing:

$ docker exec -i -t <container id> /bin/bash



