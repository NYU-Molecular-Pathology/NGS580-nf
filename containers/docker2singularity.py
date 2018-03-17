#!/usr/bin/env python
import os
import sys
import subprocess as sp
"""
script to build Singularity image files
only if
- a known image file does not already exist for that directory
"""
def main():
    # read list of known files
    singularity_images_file = "singularity.images.txt"
    singularity_images = []
    with open(singularity_images_file) as f:
        for line in f.readlines():
            path = line.strip()
            if len(path) < 1:
                continue
            else:
                singularity_images.append({
                'path': path,
                'basename': os.path.basename(path),
                'dirname': os.path.dirname(path)
                })

    dirnames = [x['dirname'] for x in singularity_images]
    basenames = [x['basename'] for x in singularity_images]

    search_dir = sys.argv[1]

    # check if provided dir is already in the list
    if search_dir in dirnames:
        print("Directory {0} already has a known prior image file; Exiting...".format(search_dir))
        return()

    # build Docker command
    print("Directory {0} does not have a known prior image file; creating Singularity image...".format(search_dir))
    command = """
    docker run -v /var/run/docker.sock:/var/run/docker.sock -v $(pwd)/{0}:/output --privileged -t --rm singularityware/docker2singularity stevekm/ngs580-nf:{0}
    """.format(search_dir)
    print("command is:\n{0}".format(command))
    process = sp.Popen(command, stdout = sp.PIPE, stderr = sp.PIPE, shell = True, universal_newlines = True)
    proc_stdout, proc_stderr = process.communicate()
    proc_stdout = proc_stdout.strip()
    proc_stderr = proc_stderr.strip()
    print(proc_stdout)
    print(proc_stderr)

if __name__ == '__main__':
    main()
