#!/bin/bash
for LANGUAGE in F95; do
  echo "Running ${LANGUAGE} tests"
  make -s test LANGUAGE=${LANGUAGE}
done
