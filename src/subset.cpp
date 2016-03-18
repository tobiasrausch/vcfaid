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

#include "gq.h"

using namespace vcfaid;

struct Config {
  boost::filesystem::path bedfile;
  boost::filesystem::path outfile;
  boost::filesystem::path vcffile;
};

template<typename TConfig, typename TScores>
inline void
_parseScores(TConfig const& c, TScores& scores)
{
  typedef typename TScores::mapped_type TMappedType;
  typedef typename TScores::key_type TKeyType;
  std::ifstream bedFile(c.bedfile.string().c_str(), std::ifstream::in);
  if (bedFile.is_open()) {
    while (bedFile.good()) {
      std::string line;
      getline(bedFile, line);
      typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
      boost::char_separator<char> sep(" \t,;");
      Tokenizer tokens(line, sep);
      Tokenizer::iterator tokIter = tokens.begin();
      if (tokIter!=tokens.end()) {
	TKeyType id = *tokIter++;
	TMappedType score = boost::lexical_cast<TMappedType>(*tokIter++);
	scores[id] = score;
      }
    }
  }
}

template<typename TConfig, typename TScores>
inline int32_t 
_processVCF(TConfig const& c, TScores const& scores) {

  // Open VCF file
  htsFile* ifile = bcf_open(c.vcffile.string().c_str(), "r");
  if (ifile == NULL) {
    std::cerr << "VCF file is missing " << c.vcffile.string() << std::endl;
    return 1;
  }
  bcf_hdr_t* hdr = bcf_hdr_read(ifile);

  // Open output file
  htsFile *fp = hts_open(c.outfile.string().c_str(), "wg");
  bcf_hdr_t *hdr_out = bcf_hdr_dup(hdr);
  bcf_hdr_remove(hdr_out, BCF_HL_INFO, "SCORE");
  bcf_hdr_append(hdr_out, "##INFO=<ID=SCORE,Number=1,Type=Float,Description=\"Structural Variant Score.\">");
  bcf_hdr_write(fp, hdr_out);

  bcf1_t* rec = bcf_init();
  while (bcf_read(ifile, hdr, rec) == 0) {
    bcf_unpack(rec, BCF_UN_INFO);
    std::string svid(rec->d.id);
    if (scores.find(svid) != scores.end()) {
      float score = scores.at(svid);
      _remove_info_tag(hdr_out, rec, "SCORE");
      bcf_update_info_float(hdr_out, rec, "SCORE", &score, 1);

      // Write record
      bcf_write1(fp, hdr_out, rec);
    }
  }
  bcf_destroy(rec);

  // Close output VCF
  bcf_hdr_destroy(hdr_out);
  hts_close(fp);

  // Close VCF
  bcf_hdr_destroy(hdr);
  bcf_close(ifile);
  return 0;
}

int main(int argc, char **argv) {

#ifdef PROFILE
  ProfilerStart("subset.prof");
#endif

  Config c;

  // Parameter
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("bed,b", boost::program_options::value<boost::filesystem::path>(&c.bedfile), "bedfile of id & score of variants to keep")
    ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("var.vcf.gz"), "VCF output file")
    ;

  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "input VCF file")
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
    std::cerr << "Input VCF file is missing " << c.vcffile.string() << std::endl;
    return 1;
  }

  // Check VCF file
  if (!(boost::filesystem::exists(c.bedfile) && boost::filesystem::is_regular_file(c.bedfile) && boost::filesystem::file_size(c.bedfile))) {
    std::cerr << "Input BED file is missing " << c.bedfile.string() << std::endl;
    return 1;
  }

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  // Parse selected Ids and Scores
  typedef std::map<std::string, double> TScores;
  TScores scores;
  _parseScores(c, scores);

  // Filter Ids and add scores
  int r=_processVCF(c, scores);

  // End
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;



#ifdef PROFILE
  ProfilerStop();
#endif

  return r;
}
