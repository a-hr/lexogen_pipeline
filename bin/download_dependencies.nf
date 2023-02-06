/*
Download all the containers without needing to run the whole pipeline, allowing to run offline afterwards.

TODO:
    * write a parser to detect the "container = {container-name}" pattern
    * create nextflow-like names for each of the containers
    * pull the images in .img format to the "containers/" directory with singularity 
*/