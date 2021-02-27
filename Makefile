
build:
	conda env update --file environment.yml --prune 

test:
	make build
	pytest

precommit:
	conda env export --no-builds > environment.yml
	make test

