#!/usr/bin/env bash
set -euo pipefail

_tests="tests/unit"

while [[ $# -gt 0 ]]
do
_key="$1"

case $_key in
    -t|--tests)
        _tests="$2"
        shift 2
        ;;
    *)
        echo "Unknown option: $_key"
        exit 1
        ;;
esac
done

echo "$_tests"

case "$_tests" in
    *unit*)
        _dc_run_opts="--no-deps --rm"
        ;;
    *)
        _dc_run_opts="--rm"
        ;;
esac

_dc_opts="-f docker-compose.yml"
_command="docker-compose $_dc_opts run $_dc_run_opts celery python -m coverage run --source=metadome -m unittest discover -s $_tests"
_report_command="docker-compose $_dc_opts run $_dc_run_opts celery python -m coverage report"

echo "$_command"
eval "$_command"

echo "$_report_command"
eval "$_report_command"
