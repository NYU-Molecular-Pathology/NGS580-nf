Variant calling container based on [GATK 4.0 Docker container](https://hub.docker.com/r/broadinstitute/gatk3/tags/3.8-0/).

```
$ docker run --rm -ti broadinstitute/gatk:4.0.1.2 bash
root@20c59635538e:/usr# cat /etc/*-release
DISTRIB_ID=Ubuntu
DISTRIB_RELEASE=16.04
DISTRIB_CODENAME=xenial
DISTRIB_DESCRIPTION="Ubuntu 16.04.3 LTS"
NAME="Ubuntu"
VERSION="16.04.3 LTS (Xenial Xerus)"
ID=ubuntu
ID_LIKE=debian
PRETTY_NAME="Ubuntu 16.04.3 LTS"
VERSION_ID="16.04"
HOME_URL="http://www.ubuntu.com/"
SUPPORT_URL="http://help.ubuntu.com/"
BUG_REPORT_URL="http://bugs.launchpad.net/ubuntu/"
VERSION_CODENAME=xenial
UBUNTU_CODENAME=xenial
```

# Software

- GATK 4.0

- Samtools 0.1.19

- Lofreq 2.1.2

- htslib 1.7

- ANNOVAR revision 150617

- vcflib
