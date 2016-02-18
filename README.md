VCFaid installation (using recursive clone)
------------------------------------------

`git clone --recursive https://github.com/tobiasrausch/vcfaid.git`

`cd vcfaid/`

`make all`

VCFaid installation (using package manager)
-------------------------------------------

Install Linux packages:

`apt-get update`

`apt-get install -y build-essential g++ git cmake zlib1g-dev ant libbz2-dev libboost-date-time-dev libboost-program-options-dev libboost-system-dev libboost-filesystem-dev libboost-iostreams-dev`

Install htslib:

`git clone https://github.com/samtools/htslib.git`

`cd htslib && make && make lib-static && cd ..`

Set environment variables for boost libraries and htslib:

`export BOOST_ROOT=/usr`

`export SEQTK_ROOT=<htslib_path>`

Build vcfaid:

`git clone https://github.com/tobiasrausch/vcfaid.git`

`cd vcfaid/ && touch .htslib .boost && make all && cd ..`


Running gq
----------

`./src/gq -g 30 -v output.vcf.gz input.vcf.gz`


Credits
-------
VCFaid takes quite a bit of actual code fragments from [arfer](https://github.com/ekg/arfer).
