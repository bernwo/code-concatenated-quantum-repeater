PYTHON		= python3
PIP_FLAGS	= --extra-index-url=https://${NETSQUIDPYPI_USER}:${NETSQUIDPYPI_PWD}@pypi.netsquid.org --extra-index-url=https://${NETSQUIDPYPI_USER}:${NETSQUIDPYPI_PWD}@pypiqutech.netsquid.org

requirements: _check_variables
	@cat requirements.txt | xargs -n 1 $(PYTHON) -m pip install ${PIP_FLAGS}
_check_variables:
ifndef NETSQUIDPYPI_USER
	$(error Set the environment variable NETSQUIDPYPI_USER before uploading)
endif
ifndef NETSQUIDPYPI_PWD
	$(error Set the environment variable NETSQUIDPYPI_PWD before uploading)
endif

.PHONY: requirements _check_variables
