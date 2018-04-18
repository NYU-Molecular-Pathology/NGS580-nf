#!/bin/bash
# wrapper script for use in the pipeline with Conda installed trimmomatic;
# used for consistency in pipeline command regardless of environment implementation
trimmomatic "$@"
