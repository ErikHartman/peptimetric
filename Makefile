
build:
	poetry install

test:
	pytest

precommit:
	make test

