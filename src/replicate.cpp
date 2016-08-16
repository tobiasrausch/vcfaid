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
#include <htslib/vcf.h>

#include "arfer.h"

using namespace vcfaid;

struct Config {
  typedef std::vector<std::string> TSamples;
  TSamples ctrl;
  TSamples tmrs;

  typedef std::map<std::string, uint32_t> TSampleIdMap;
  TSampleIdMap simap;

  float minBAF;
  int32_t minReplicateSupport;
  boost::filesystem::path outfile;
  boost::filesystem::path vcffile;
  boost::filesystem::path samplefile;
};


template<typename TConfig>
inline int32_t
_checkReplicates(TConfig const& c) {
  // Open VCF file
  htsFile* ifile = bcf_open(c.vcffile.string().c_str(), "r");
  bcf_hdr_t* hdr = bcf_hdr_read(ifile);

  // Open output file
  std::ofstream ofile(c.outfile.string().c_str());
  
  // Parse VCF
  bcf1_t* rec = bcf_init();
  while (bcf_read(ifile, hdr, rec) == 0) {
    bcf_unpack(rec, BCF_UN_ALL);
    int ngt = 0;
    int32_t* gt = NULL;
    int nrv = 0;
    int32_t* rv = NULL;
    int nrr = 0;
    int32_t* rr = NULL;
    int ndv = 0;
    int32_t* dv = NULL;
    int ndr = 0;
    int32_t* dr = NULL;
    bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt);
    bcf_get_format_int32(hdr, rec, "DV", &dv, &ndv);
    bcf_get_format_int32(hdr, rec, "DR", &dr, &ndr);
    bcf_get_format_int32(hdr, rec, "RV", &rv, &nrv);
    bcf_get_format_int32(hdr, rec, "RR", &rr, &nrr);
    bool precise = false;
    if (bcf_get_info_flag(hdr, rec, "PRECISE", 0, 0) > 0) precise = true;
    
    typedef std::vector<double> TBaf;
    typedef std::vector<int32_t> TAltSupport;
    typedef std::vector<bool> TCarrier;
    TBaf ctrlBaf(c.ctrl.size(), -1);
    TAltSupport tmrsAltSupport(c.tmrs.size(), -1);
    TCarrier car(c.ctrl.size(), false);
    
    // Estimate allele frequency
    uint32_t ac[2];
    ac[0] = 0;
    ac[1] = 0;
    for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
      typename TConfig::TSampleIdMap::const_iterator itSI = c.simap.find(hdr->samples[i]);
      if (itSI != c.simap.end()) {
	if ((bcf_gt_allele(gt[i*2]) != -1) && (bcf_gt_allele(gt[i*2 + 1]) != -1)) {
	  int gt_type = bcf_gt_allele(gt[i*2]) + bcf_gt_allele(gt[i*2 + 1]);
	  ++ac[bcf_gt_allele(gt[i*2])];
	  ++ac[bcf_gt_allele(gt[i*2 + 1])];
	  if (c.ctrl[itSI->second] == itSI->first) {
	    if (gt_type != 0) {
	      car[itSI->second] = true;
	      if (precise) ctrlBaf[itSI->second] = (double) rv[i] / (double) (rr[i] + rv[i]);
	      else ctrlBaf[itSI->second] = (double) dv[i] / (double) (dr[i] + dv[i]);
	    }
	  } else {
	    if (precise) tmrsAltSupport[itSI->second] = rv[i];
	    else tmrsAltSupport[itSI->second] = dv[i];
	  }
	}
      }
    }

    // Filter rare SVs using the "tumor" replicate
    double af = (double) ac[1] / (double) (ac[0] + ac[1]);
    if ((af > 0) && (af <= 0.01)) {
      int32_t bestSupport = -1;
      double bestBAF = -1;
      for(uint32_t i = 0; i<car.size(); ++i) {
	if (car[i]) {
	  if (tmrsAltSupport[i] >= bestSupport) {
	    // It's enough to be higher than min BAF
	    if ((ctrlBaf[i] >= c.minBAF) || (ctrlBaf[i] >= bestBAF)) {
	      bestSupport = tmrsAltSupport[i];
	      bestBAF = ctrlBaf[i];
	    }
	  }
	}
      }
      std::string varId = rec->d.id;
      if ((bestBAF >= c.minBAF) && (bestSupport >= c.minReplicateSupport)) ofile << varId << "\t1" << std::endl;
      else ofile << varId << "\t0" << std::endl;
    }


    // Clean-up
    if (dv != NULL) free(dv);
    if (dr != NULL) free(dr);
    if (rv != NULL) free(rv);
    if (rr != NULL) free(rr);
    if (gt != NULL) free(gt);
  }
  bcf_destroy(rec);

  // Close output file
  ofile.close();

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
    ("samples,s", boost::program_options::value<boost::filesystem::path>(&c.samplefile), "sample file")
    ("baf,b", boost::program_options::value<float>(&c.minBAF)->default_value(0.25), "min. B-allele frequency")
    ("support,p", boost::program_options::value<int32_t>(&c.minReplicateSupport)->default_value(2), "min. support in replicate")
    ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.tsv"), "output tsv file")
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
  if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("samples"))) {
    std::cout << "Usage: " << argv[0] << " [OPTIONS] -s <samples.tsv> <input.bcf>" << std::endl;
    std::cout << visible_options << "\n";
    return 1;
  } 
  
  // Check VCF file
  if (!(boost::filesystem::exists(c.vcffile) && boost::filesystem::is_regular_file(c.vcffile) && boost::filesystem::file_size(c.vcffile))) {
    std::cerr << "Input VCF/BCF file is missing: " << c.vcffile.string() << std::endl;
    return 1;
  }

  // Check sample file
  if (!(boost::filesystem::exists(c.samplefile) && boost::filesystem::is_regular_file(c.samplefile) && boost::filesystem::file_size(c.samplefile))) {
    std::cerr << "Sample file is missing " << c.samplefile.string() << std::endl;
    return 1;
  } else {
    // Get samples
    uint32_t pos = 0;
    std::ifstream sampleFile(c.samplefile.string().c_str(), std::ifstream::in);
    if (sampleFile.is_open()) {
      while (sampleFile.good()) {
	std::string sampleFromFile;
	getline(sampleFile, sampleFromFile);
	typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	boost::char_separator<char> sep(",\t ");
	Tokenizer tokens(sampleFromFile, sep);
	Tokenizer::iterator tokIter = tokens.begin();
	if (tokIter != tokens.end()) {
	  std::string control = *tokIter++;
	  if (tokIter != tokens.end()) {
	    std::string tumor = *tokIter;
	    c.ctrl.push_back(control);
	    c.tmrs.push_back(tumor);
	    c.simap.insert(std::make_pair(control, pos));
	    c.simap.insert(std::make_pair(tumor, pos));
	    ++pos;
	  }
	}
      }
      sampleFile.close();
    }
    if ((c.ctrl.empty()) || (c.tmrs.empty())) {
      std::cerr << "No samples specified." << std::endl;
      return 1;
    }
  }
  
  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  int r = _checkReplicates(c);

  // End
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;



#ifdef PROFILE
  ProfilerStop();
#endif

  return r;
}
