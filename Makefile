clean:
	@rm -rf work/
	@rm -rf output/
	@rm -f slurm*.out
	@rm -f trace*.txt
	@rm -f .nextflow.log*
	@rm -f *.mmd

# singularity and golang must be available in the working environment
pull:
	echo "Pulling containers ..."
	@mkdir -p containers
	@for container in `grep -oP "(?<=container = ').*(?=')" confs/slurm.config`; do \
		echo "Checking if $$container is already downloaded ..."; \
		containerName=`echo $$container | sed 's/\//-/g' | sed 's/:/-/g'`; \
		if [ -f containers/$$containerName.img ]; then \
			echo "$$containerName.img already exists. Skipping download."; \
			echo ""; \
		else \
			echo "Pulling $$container ..."; \
			singularity -s pull --name $$containerName.img --dir containers/ docker://$$container; \
			echo ""; \
		fi; \
	done
	@echo "Done!"