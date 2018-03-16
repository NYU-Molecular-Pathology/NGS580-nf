Variant calling container based on [GATK 3.8 Docker container](https://hub.docker.com/r/broadinstitute/gatk3/tags/3.8-0/).

```
$ docker run --rm -ti broadinstitute/gatk3:3.8-0 bash
root@20c59635538e:/usr# cat /etc/*-release
PRETTY_NAME="Debian GNU/Linux 8 (jessie)"
NAME="Debian GNU/Linux"
VERSION_ID="8"
VERSION="8 (jessie)"
ID=debian
HOME_URL="http://www.debian.org/"
SUPPORT_URL="http://www.debian.org/support"
BUG_REPORT_URL="https://bugs.debian.org/"
```
Using `/usr/GenomeAnalysisTK.jar` located in the default working directory

# Software

- GATK 3.8

- Samtools 0.1.19

- Lofreq 2.1.2

- htslib 1.7

- ANNOVAR revision 150617

- vcflib
