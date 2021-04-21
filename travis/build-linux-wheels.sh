#!/bin/bash
set -e -x

# Compile wheels for python 3.7 & 3.8
for PYBIN in /opt/python/*3[78]*/bin; do
    "${PYBIN}/pip" install -r /io/requirements.txt
    "${PYBIN}/pip" wheel /io/ -w dist/
done

# Bundle external shared libraries into the wheels
for whl in dist/*.whl; do
    auditwheel repair "$whl" --plat $PLAT -w /io/dist/
done

# Install packages and test
#for PYBIN in /opt/python/*/bin/; do
#    "${PYBIN}/pip" install ttcrpy --no-index -f /io/wheelhouse
#    (cd "$HOME"; "${PYBIN}/nosetests" pymanylinuxdemo)
#done
