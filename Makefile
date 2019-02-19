default:
	pip install -e ../ThredgeCorr --no-binary :all:

install:
	pip install ../ThredgeCorr

pypi:
	rm dist/*
	python setup.py sdist
	twine upload dist/*
