#!/bin/bash
TESTS=$(ls *.x)
NTESTS=$(echo $TESTS | wc -w)
NSUC=0
FAILED=()

echo "Running tests..."
echo "================="
for binary in $TESTS
do
    printf "\n\033[0;36mTest $binary\033[0m\n"
    ./$binary
    STATUS=(1 - $?)
    NSUC=$((NSUC + STATUS))
    if [[ $NSUC == 0 ]]; then
        FAILED+="$binary"
    fi
done

printf "\n$NSUC out of $NTESTS tests successful.\n"
if [[ ${FAILED[@]} -gt 0 ]]; then
    echo "Failed tests: $FAILED"
fi