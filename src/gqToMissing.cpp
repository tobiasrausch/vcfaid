/*
============================================================================
VCFaid
============================================================================
Copyright (C) 2016 Tobias Rausch

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
============================================================================
Contact: Tobias Rausch (rausch@embl.de)
============================================================================
*/

#define _SECURE_SCL 0
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <fstream>

#define BOOST_DISABLE_ASSERTS
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include "arfer.h"
#include "gq.h"

using namespace vcfaid;

struct Config {
  int32_t gqthreshold;
  boost::filesystem::path outfile;
  boost::filesystem::path vcffile;
};


template<typename TConfig>
inline int32_t 
_setToMissing(TConfig const& c) {

  // Open VCF file
  htsFile* ifile = bcf_open(c.vcffile.string().c_str(), "r");
  bcf_hdr_t* hdr = bcf_hdr_read(ifile);

  // Open output file
  htsFile *fp = hts_open(c.outfile.string().c_str(), "wb");
  bcf_hdr_t *hdr_out = bcf_hdr_dup(hdr);
  bcf_hdr_write(fp, hdr_out);

  bcf1_t* rec = bcf_init();
  while (bcf_read(ifile, hdr, rec) == 0) {
    bcf_unpack(rec, BCF_UN_ALL);
    int ngt = 0;
    int32_t* gt = NULL;
    int ngq = 0;
    int32_t* gq = NULL;
    bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt);
    bcf_get_format_int32(hdr, rec, "GQ", &gq, &ngq);
    for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
      if (gq[i] < c.gqthreshold) {
	gt[i*2] = bcf_gt_missing;
	gt[i*2 + 1] = bcf_gt_missing;
      }
    }
    bcf_update_genotypes(hdr_out, rec, gt, bcf_hdr_nsamples(hdr) * 2);

    // Write record
    bcf_write1(fp, hdr_out, rec);
    
    // Clean-up
    free(gq);
    free(gt);
  }
  bcf_destroy(rec);

  // Close output VCF
  bcf_hdr_destroy(hdr_out);
  hts_close(fp);

  // Build index
  bcf_index_build(c.outfile.string().c_str(), 14);

  // Close VCF
  bcf_hdr_destroy(hdr);
  bcf_close(ifile);
  return 0;
}

int main(int argc, char **argv) {

#ifdef PROFILE
  ProfilerStart("vcfaid.prof");
#endif

  Config c;

  // Parameter
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("gqthreshold,g", boost::program_options::value<int32_t>(&c.gqthreshold)->default_value(20), "GQs below will be GT=./.")
    ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("var.bcf"), "BCF output file")
    ;

  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "input VCF/BCF file")
    ;

  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);

  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);

  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("input-file"))) {
    std::cout << "Usage: " << argv[0] << " [OPTIONS] <input.vcf.gz>" << std::endl;
    std::cout << visible_options << "\n";
    return 1;
  } 
  
  // Check VCF file
  if (!(boost::filesystem::exists(c.vcffile) && boost::filesystem::is_regular_file(c.vcffile) && boost::filesystem::file_size(c.vcffile))) {
    std::cerr << "Input VCF/BCF file is missing: " << c.vcffile.string() << std::endl;
    return 1;
  }

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  int r=_setToMissing(c);

  // End
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;


#ifdef PROFILE
  ProfilerStop();
#endif

  return r;
}
